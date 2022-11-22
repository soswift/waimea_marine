
# This script loads and cleans the pipeline outputs for metabolite and microbe analysis pipelines
# The R package phyloseq is used to organize sample data, taxonomy, and feature counts for both microbes and metabolites
# Phyloseq objects and flat csv tables are passed to subsequent scripts to generate statistical analyses and figures
# Setup environment ----------------------------------------

# Load libraries

library(phyloseq)
library(data.table)
library(eulerr)
library(vegan)
library(pairwiseAdonis)
library(tidyr)
library(biomformat)

#library(speedyseq)

## identify data files

# microbial data files
abundance_file   <- "data/raw/abundance_table_100.shared"

sample_data_file <- "data/raw/new_marine_meta.csv"

taxonomy_file    <- "data/raw/annotations_100.taxonomy"

unifrac_file     <- "data/raw/unifrac_unweighted_100.csv"

tree_file        <- "data/raw/FastTree_100.tre"

# source lots of helper functions for plotting and subsetting
source("src/helper_functions.R")

# source function for reading and cleaning abundance table
source("src/clean_and_query_16S.R")

# source helper function for making phyloseq objects
source("src/make_phyloseq.R")

# source function for writing out phyloseq objects
source("src/physeq_csv_out.R")

# function for reading unifrac dists from flat tables
source("src/read_unifrac.R")

# function for formatting microbial data for mmvec
source("src/format_mmvec.R")


# set color scheme for sample types
sample_type_cols <- c(CCA = "#6469ed",
                      Coral = "#e49c4c",
                      Limu = "#7cc854" )

# 1. Load and clean microbial data ----------------------------
# Load and clean data from the metaflowimcs pipeline.
# Expects standards pipeline output format.
# Data is subset to only include desired marine samples.
# Read counts are subsampled to even depth and relative abundance transformed.

# read in data and sample data, 

marine_data <- clean_16S_tables(abundance_file = abundance_file,
                                taxonomy_file = taxonomy_file,
                                metadata_file = sample_data_file,
                                description = "ST_marine",
                                output_dir = "data/processed",
                                id_column = "sequencing_id",
                                cull = list(min.num = 3,
                                            min.abund = 0.00001,
                                            min.single.abund = 0.001))

marine_phy <- make_phyloseq(abund_file = "data/processed/ST_marine_abundance_table.csv",
                            meta_file  = "data/processed/ST_marine_metadata_table.csv",
                            tax_file   = "data/processed/ST_marine_taxonomy_table.csv",
                            id_column  = "sequencing_id")

rm(marine_data)

# check variation in sequencing depths and subsample to standard value
nrow(otu_table(marine_phy))
length(sample_names(marine_phy))
hist(sample_sums(marine_phy), breaks = 30)

marine_phy_rar <- rarefy_even_depth(marine_phy,
                                    rngseed = 1,
                                    sample.size = 15000,
                                    replace = F)
# check total OTUs
nrow(otu_table(marine_phy_rar))

# check which samples were dropped due to low read counts
marine_data[["metadata"]][!(sequencing_id %in% sample_names(marine_phy_rar))]
length(sample_names(marine_phy_rar))

# convert the OTU counts to relative abundance
final_marine_phy <- transform_sample_counts(marine_phy_rar,
                                            function(x) x / sum(x))

# drop sediment and water samples
final_marine_phy <- subset_samples(final_marine_phy,
                                   sample_type %in% c("Limu","CCA","Coral"))

# Re-classify genera based on our target genera for Limu and Coral
keep_genera <- c("Jania", "Halimeda","Galaxaura", "Porites", "Monitipora", "Pocillopora")

genus_data <- sample_data(final_marine_phy)$genus

sample_data(final_marine_phy)$genus <- ifelse(genus_data %in% keep_genera,
                                              genus_data,
                                              "Other" )

final_marine_phy <- prune_taxa(x = final_marine_phy,
                               taxa = taxa_sums(final_marine_phy) > 0)


saveRDS(final_marine_phy, file = "data/processed/final_marine_phy.rds")

rm(marine_phy_rar)

# write marine microbial sample tables out as flat files
physeq_csv_out(final_marine_phy,
               description = "all_marine_microbe",
               outdir = "data/processed/table_exports/")

# finally, get unifrac distance object for the selected samples
# this helper function completes the distance triangle matrix, matches the samples, and converts the matrix to a dist object
final_unifrac <- read_unifrac(unifrac_file = "data/raw/unifrac_weighted_100.csv", 
                              phyloseq_obj = final_marine_phy)

# write out unifrac subset to marine samples
saveRDS(final_unifrac, "data/processed/final_unifrac.rds")
write.csv(as.matrix(final_unifrac),
          "data/processed/table_exports/all_marine_microbe_unifrac_dist.csv")

# Final_unifrac contains all the relevant marine microbial samples.
# Some of these samples may not have corresponding metabolomics data.
# Incomparable samples will be dropped for joint analyses.


# 2. export data for MMVEC -----------------------------------------

# to match metabolomics data, replace the id used for sequencing with the physical sample name
# from the metadata, pull out the sample barcode for each sequencing id, matching the order in the abundance table

new_ids  <-  marine_data[["metadata"]][ match( colnames( get_taxa( final_marine_phy ) ), sequencing_id ),
                                        sample_barcode]

format_mmvec(phyloseq_obj = final_marine_phy,
             id_column = "Group",
             new_ids = new_ids,
             output_dir = "data/processed",
             description = "marine")

# 3. Load and clean metabolomics data ------------------------------------------------

# Metabolomics was data provided in one big table with peak classification information and peak areas in the same file.
# The sample data file matches the microbial sample data

# metabolomics data files
chem_raw_file    <- "data/raw/Helena_CCA_Coral_Limu_FeatureTable.txt"

chem_sample_file <- "data/raw/new_marine_meta.csv"


# this table combines abundance and metabolite feature metadata
chem_raw <- fread(chem_raw_file)

# read in metabolomics sample data
chem_meta <- fread(chem_sample_file)

# clean up sample genus info
keep_genera <- c("Jania", "Halimeda","Galaxaura", "Porites", "Monitipora", "Pocillopora")

chem_meta[ , genus := ifelse(genus %in% keep_genera, genus, "Other" )]


# format metabolite data using custom cleaning functions
# 3 data versions: 
# 'no filter' includes everything that was run
# 'w_blanks' includes blanks and samples, but drops 'bad samples' that clustered with blanks in the raw data set
# the final version includes only good samples, not bad samples or blanks 

## RAW
# pull out abundance data (i.e. peak areas)
chem_abund_raw <- get_chem_abundance(chem_data = chem_raw, bad_samples = NULL)
# pull out peak metadata (i.e. classifications, networks, etc.)
chem_taxa_raw <- get_chem_peak_data(chem_data = chem_raw)

# arrange metadata
chem_meta_raw <- chem_meta[match(colnames(chem_abund_raw), chem_meta$sample_barcode)]
row.names(chem_meta_raw) <- chem_meta_raw[ , sample_barcode]

chem_phy_raw <- make_chem_phyloseq(chem_abund = chem_abund_raw,
                                        chem_meta = chem_meta_raw,
                                        chem_tax = chem_taxa_raw,
                                        id_column = "sample_barcode")


## WITH BLANKS
# drop bad samples (determined by ordination of raw dataset)
bad_samples <- c(2638,
                 2839,
                 2684,
                 2750,
                 2815,
                 2650,
                 2795,
                 2862,
                 2905)

# pull out abundance data (i.e. peak areas)
chem_abund_w_blanks <- get_chem_abundance(chem_data = chem_raw,
                                          bad_samples = bad_samples)

# pull out peak metadata (i.e. classifications, networks, etc.)
chem_taxa_w_blanks <- get_chem_peak_data(chem_data = chem_raw)

# arrange metadata
chem_meta_w_blanks <- chem_meta[match(colnames(chem_abund_w_blanks), chem_meta$sample_barcode)]
row.names(chem_meta_w_blanks) <- chem_meta_w_blanks[ , sample_barcode]

chem_phy_w_blanks <- make_chem_phyloseq(chem_abund = chem_abund_w_blanks,
                                        chem_meta = chem_meta_w_blanks,
                                        chem_tax = chem_taxa_w_blanks,
                                        id_column = "sample_barcode")

## NO BLANKS
# drop blanks and remove any chemical features not present in the real samples
not_blanks <- !grepl("blank", colnames(chem_abund_w_blanks))

chem_abund <- chem_abund_w_blanks[ , not_blanks]
chem_abund <- chem_abund[rowSums(chem_abund) > 0, ]
write.csv(chem_abund, "data/processed/table_exports/all_marine_metabolite_raw_abundance.csv")

# subset metadata to match cleaned abundance
chem_meta <- chem_meta[match(colnames(chem_abund), chem_meta$sample_barcode)]


# call function to make phyloseq of metabolite data
chem_phy <- make_chem_phyloseq(chem_abund = chem_abund,
                               chem_meta = chem_meta, 
                               chem_tax = chem_taxa_w_blanks,
                               id_column = "sample_barcode")

# convert raw peak areas to relative abundance
chem_phy <- transform_sample_counts(chem_phy,
                                    function(x) x / sum(x))


# write out flat tables for the real samples and the real samples with blanks included
physeq_csv_out(chem_phy_raw,
               description = "raw_marine_metabolite",
               outdir = "data/processed/table_exports/")

physeq_csv_out(chem_phy,
               description = "all_marine_metabolite",
               outdir = "data/processed/table_exports")

physeq_csv_out(chem_phy_w_blanks,
               description = "all_marine_metabolite_w_blanks",
               outdir = "data/processed/table_exports/")

saveRDS(chem_phy_raw, "data/processed/chem_phy_raw.rds")
saveRDS(chem_phy_w_blanks, "data/processed/chem_phy_w_blanks.rds")
saveRDS(chem_phy, "data/processed/chem_phy.rds")

# clean up
rm(chem_abund_w_blanks, 
   chem_taxa_w_blanks,
   chem_meta_w_blanks,
   chem_abund,
   chem_meta)


