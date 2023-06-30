# Analyses related to reviewer responses
library(data.table)
library(vegan)
library(ComplexHeatmap)
library(ggpubr)
library(pairwiseAdonis) 

setwd("~/Documents/Bioinformatics/Projects/SuperTransect/Analysis/github_waimea_marine/")
source("src/colors.R")

theme_set(theme_minimal())
# ----------------------------------------------------------------
#1. Cyanobacteria ----------------------------------------------------
# ----------------------------------------------------------------


# Cyanobacteria are another group of important/dominant primary producers on coral reefs
# any interesting patterns here?


fam_abund = fread("data/processed/microbe_family_abundance_marine_waimea.csv")
sample_dat = fread("data/processed/table_exports/all_marine_microbe_sample_flat_table.csv")

meta.cols = c(
  "family",
  "class",
  "n_samp_obs",
  "sample_type_DA",
  "adj_p_val",
  "FC_LimuVCoral",
  "FC_LimuVCCA",
  "FC_CoralVCCA"
)

cyano_abund = fam_abund[class == "Cyanobacteriia" & n_samp_obs > 30,]

write.csv(cyano_abund,"data/revision/cyanobacteriia_summary.csv")


# make matrix for heatmap
cyano_mat = as.matrix(cyano_abund[ , .SD, .SDcol =!meta.cols], rownames = cyano_abund$family)
cyano_mat = cyano_mat[,colSums(cyano_mat) > 0]

sample_dat = sample_dat[match(colnames(cyano_mat), as.character(sample_dat$sequencing_id)),]

# make dot plots

cyano_dt = as.data.frame(t(cyano_mat))
cyano_dt$sequencing_id = as.integer(row.names(cyano_dt))

dat = merge(sample_dat, cyano_dt, by = "sequencing_id")

violins = lapply(cyano_abund$family, function(a_family){
  
  ggplot(dat, aes(x = sample_type,
                  y = log(.data[[a_family]]),
                  color = sample_type))+
    geom_violin(draw_quantiles = c(0.5))+
    geom_jitter()+
    scale_color_manual(values = sample_type_cols, guide = F)+
    labs(y = "Log Relative Abundance",
         title = a_family,
         subtitle = paste("DA =", cyano_abund[family == a_family, sample_type_DA]))
  
})
names(violins) = cyano_abund$family



cyano_p = ggarrange(plotlist = violins, nrow = 4, ncol = 2, common.legend = T)
cyano_p
ggsave("output/violins/cyano_family_violins.pdf", plot = cyano_p, width = 8, height = 10)

# ----------------------------------------------------------------
# 2. Important Features --------------
# ----------------------------------------------------------------
# These 128 ion features might be listed into a table in the supplemental?
# They seem like important direct information regarding the hypothesis and the interpretation.
rf = read.csv("data/processed/all_metabolite_RF_LM_ord.csv")
metabolite_meta  = read.csv("data/processed/table_exports//all_marine_metabolite_tax_flat_table.csv")


drop.cols = c(
  "X",
  "Y",
  "Z",
  "transp"
)

rf = rf[, !(colnames(rf) %in% drop.cols)]  
rf = rf[rf$RF_selection == "Important",]


metabolite_meta$featureID = metabolite_meta$X

rf = merge(rf, metabolite_meta, all.x = T, all.y = F, by = "featureID")


write.csv(rf, "data/revision/rf_selected_features.csv")

# ----------------------------------------------------------------
# 3. Sample Numbers ----------------------------------------------------------------
# ----------------------------------------------------------------
# How many samples?
meta = fread("data/processed/ST_marine_full_run_metadata_table.csv")

# from each genus
meta [ ,.N, by = genus]

# from each type
meta[ ,.N, by = common_name]

# total
meta[common_name %in% c("Coral","CCA","Algae"), .N]

# summary list of all samples collected
total_samples_collected = meta[common_name %in% c("Coral","CCA","Algae"), .(sample_name, sample_type, genus)]
fwrite(total_samples_collected, "data/revision/total_samples_collected.csv")

# samples that were successfully sequenced
microbe_meta = fread("data/processed/table_exports/all_marine_microbe_sample_flat_table.csv")

total_samples_sequenced = microbe_meta[ , .(sample_name, sample_type, genus)]
fwrite(total_samples_sequenced, "data/revision/total_samples_sequenced.csv")

# ----------------------------------------------------------------
#4. Extraction negatives -------------------------------------
# ----------------------------------------------------------------
# We need to identify extraction negatives from the DNA extraction plates
# The extraction negatives are not included in pipeline outputs?
# Compare total reads and do a quick ordination

meta = fread("data/raw/new_marine_meta.csv")
controls = fread("data/raw/16S_qc/all_sample_negative_control_key.csv")
#abund = fread("data/raw/16S_qc/main/details/all_postprocessing_97.relabund", drop = c("label","numOtus"))
raw_tax = fread("data/raw/16S_qc/Results/raw/details/raw_consensusClassification_97.taxonomy")
abund = fread("data/raw/16S_qc/Results/raw/details/raw_abundanceTable_97.shared", drop = c("label","numOtus"))
sample_dat = fread("data/processed/table_exports/all_marine_microbe_sample_flat_table.csv")


# clean up
setnames(abund, "Group", "sequencing_id")


tax = fread(text = raw_tax$Taxonomy,
            sep = ";")
names(tax) = c(
  "kingdom",
  "phylum",
  "class",
  "order",
  "family",
  "genus",
  "species"
)

# drop the species column because it looks messy
tax[ , species := NULL]

# add the OTU identifier back into the taxonomy table (everything should still be in the same order)
tax[ , OTU := raw_tax$OTU]


# identify PCR and extraction negatives from relevant PCR plates
pcr_plates = unique(meta$pcr_plate)
neg_controls = controls[control_detail %in% c("ExtractionNegative","PCRNegative") &
                          pcr_plate %in% pcr_plates]

# make joined table of real samples and controls
meta$control_detail = "Sample"
meta = meta[!is.na(sequencing_id)]
meta = meta[sequencing_id %in% sample_dat$sequencing_id,]

merge_cols = c("sequencing_id",
               "pcr_plate",
               "control_detail",
               "common_name",
               "sample_id")

blank_meta = rbind(neg_controls[ , ..merge_cols] , meta[ , ..merge_cols]) 

# do some comparisons

## read counts
abund_mat = as.matrix(abund[ , .SD, .SDcols = !"sequencing_id"], rownames = abund$sequencing_id)
read_counts = rowSums(abund_mat)

blank_meta$read_counts = ifelse(as.character(blank_meta$sequencing_id) %in% names(read_counts), 
                                yes = read_counts[as.character(blank_meta$sequencing_id)],
                                no =  NA)


ggplot(blank_meta , aes(x = control_detail, y = read_counts))+
  geom_boxplot()+
  geom_jitter()

blank_meta[ !is.na(read_counts) , mean(read_counts), by = control_detail]

blank_meta[ !is.na(read_counts) , .N, by = control_detail]



## ASV counts
pa_abund = abund[ , lapply(.SD, as.logical), by = "sequencing_id"]
sample_ASV_counts = data.table(counts = rowSums(pa_abund[, .SD,
                                                         .SDcols = !"sequencing_id"]),
                               sequencing_id = abund$sequencing_id)

blank_meta_counts = merge(blank_meta, sample_ASV_counts, all.y = T)
blank_meta_counts[ , mean(counts), by = control_detail]


ggplot(blank_meta_counts, aes(x = control_detail, y = counts, color = control_detail))+
  geom_boxplot()+
  scale_color_manual(values = c("darkgray","black","darkblue"), guide = NULL)+
  labs(x = NULL, y = "Unique ASVs")


## Top Taxa
tax_abund = melt.data.table(abund,
                            id.vars = "sequencing_id",
                            variable.name = "OTU",
                            value.name = "reads")

tax_abund = merge(tax_abund, tax[, .(OTU,genus)], by = "OTU")

tax_abund = tax_abund[ , .(reads = sum(reads)),
                       by = .(genus, sequencing_id)]

tax_abund = merge(tax_abund, blank_meta, by = "sequencing_id")

tax_abund[ , pa := as.integer(as.logical(reads))]

tax_summary = tax_abund[ pa > 0 ,
                         .(OTU_count = sum(reads)),
                         by = c("control_detail","genus")]

control_tax_summary = tax_summary[control_detail %in% c("ExtractionNegative",
                                                        "PCRNegative")][
                                                          order(OTU_count,
                                                                decreasing = T)]



## Ordination
ord_samples = function(samps){
  samps = samps[samps %in% rownames(unifrac)]
  dist_mat = unifrac[samps,samps]
  ord = metaMDS(comm = as.dist(dist_mat))
  out = list(points = as.data.table(ord$points, keep.rownames = "seqID"),
             stress = round(ord$stress, 3))
  return(out) 
}


#----------------------------------------------------------------
# 5. Nitrogen containing compounds ---------------------------------
#----------------------------------------------------------------
# quickly check percentage of N containing compounds in coral vs. cca vs. algae
metabolite_abundance_file = "data/processed/table_exports/all_marine_metabolite_abundance_flat_table.csv"
sample_data_file          = "data/processed/table_exports/all_marine_metabolite_sample_flat_table.csv"
metabolite_metadata_file  = "data/processed/table_exports/all_marine_metabolite_tax_flat_table.csv"

abund = fread(metabolite_abundance_file, header = T)
samp = fread(sample_data_file)
tax = fread(metabolite_metadata_file)

setnames(abund, "V1", "featureID")
setnames(tax, "V1", "featureID")

# pull out compounds in corals vs compounds in algae
types = c("Coral", "Algae", "CCA")

type_feats = list()

for(a_type in types){
  
  type_samples = samp[common_name == a_type, as.character(sample_barcode)]
  
  # pull out features that are >0 in the sample type
  type_feat_names = abund[rowSums(abund[ , .SD, .SDcols=type_samples]) > 0, featureID]

  type_feats[[a_type]] = tax[featureID %in% type_feat_names,]
  
  print(a_type)
  print(nrow(type_feats[[a_type]][grepl("N",Qemistree_csi_smiles)]) / nrow(type_feats[[a_type]]))
  
  }


#----------------------------------------------------------------
  # 6. Metabolite Threshold Filter ---------------------------------
#----------------------------------------------------------------
metabolite_abundance_file = "data/processed/table_exports/all_marine_metabolite_w_blanks_abundance_flat_table.csv"
sample_data_file          = "data/processed/table_exports/all_marine_metabolite_w_blanks_sample_flat_table.csv"
random_forest_and_LM_file  = "data/processed/all_metabolite_RF_LM_ord.csv"
metabolite_taxonomy_file = "data/processed/table_exports/all_marine_metabolite_tax_flat_table.csv"

abund = fread(metabolite_abundance_file, header = T)
samp = fread(sample_data_file)
rf_lm_ord = fread(random_forest_and_LM_file)
tax = fread(metabolite_taxonomy_file)

setnames(abund, "V1", "featureID")
setnames(tax, "V1", "featureID")

# get max in samples and blanks
blank_ids = samp[sample_type == "Blank", sample_barcode]
samp_ids = samp[sample_type != "Blank", sample_barcode]

blank_max = abund[, .(blank_max = max(.SD)),
                  .SDcols = blank_ids,
                  by = featureID]
sample_max = abund[, .(sample_max = max(.SD)), 
                   .SDcols = samp_ids,
                   by = featureID]


# merge outputs
max_dt = merge(blank_max, sample_max, by = "featureID")
max_dt[ , sample_blank_ratio := sample_max/blank_max]

# drop features only found in blanks
max_dt = max_dt[sample_max > 0,]

# identify features where max in samples is less than 3 times higher than max in blanks
bad_feats = max_dt[sample_max/blank_max < 3]
nrow(bad_feats)


# identify statistically significant feats that did not pass this 3:1 threshold
bad_rf_lm_ord = merge(bad_feats,
                      rf_lm_ord,
                      by = "featureID",
                      all.y = F)

bad_rf_lm_ord = merge(bad_rf_lm_ord,
                      tax,
                      by = "featureID",
                      all.y = F)

write.csv(bad_rf_lm_ord, "data/revision/features_failing_threshold_cutoff.csv")

#----------------------------------------------------------------
# 6. Plate effects ---------------------------------
#----------------------------------------------------------------

source("src/helper_functions.R")

# read in
plates = fread("data/revision/Plates_only_MetaboAnalyst_samplesremoved.txt")
rel_abund = fread("data/processed/table_exports/all_marine_metabolite_abundance_flat_table.csv", header = T)
samps = fread("data/processed/table_exports/all_marine_metabolite_sample_flat_table.csv")

# clean
plates = transpose(plates, keep.names = "sample_name", make.names = "row ID")
plates = plates[ , .(sample_name, plate)]
plates[ , sample_name:= gsub("_.*", "", sample_name)]

#waimea_samples = samps[site_name == "WaimeaBay", sample_name]
#plates = plates[!(sample_name %in% waimea_samples),]

rel_abund = transpose(rel_abund, make.names = "V1", keep.names = "sample_name")
rel_mat = as.matrix(rel_abund, rownames = "sample_name")
rel_mat = rel_mat[plates$sample_name,]

samps[,sample_name := as.character(sample_name)]

# Samples --------------------



# dist
d = vegdist(rel_mat, method = "bray")


# ordination
ord = metaMDS(comm = d)
ord_out = list(
  points = as.data.table(ord$points, keep.rownames = "sample_name"),
  stress = round(ord$stress, 3)
)

points = ord_out$points

ord_meta = merge(points, plates, by = "sample_name")
ord_meta = merge(ord_meta, samps, by = "sample_name")


ord_meta[ , .N, by = .(plate, site_name)]


# plot
p = ggplot(ord_meta, aes(x = MDS1,
                     y = MDS2,
                     color = plate,
                     shape = site_name)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("indianred3",
                                "dodgerblue3"))+
  labs(caption = paste("Stress =", ord_out$stress))

ggsave("output/NMDS/metabolite_ordination_no_Waimea_Bay.pdf",
       plot = p,
       scale = 2)
·
adonis2(d ~ plate + sample_type*site_name, data = ord_meta[match(row.names(rel_mat), sample_name)])



# Blanks -----------------------------

# read in
b_abund = fread("data/processed/table_exports/all_marine_metabolite_w_blanks_abundance_flat_table.csv", header = T)
blanks = fread("data/processed/table_exports/all_marine_metabolite_w_blanks_sample_flat_table.csv")
plates = fread("data/revision/Plates_only_MetaboAnalyst.txt")

# clean
plates = transpose(plates, keep.names = "sample_name", make.names = "row ID")
plates = plates[ , .(sample_name, plate)]
plates[ , sample_name:= gsub("_.*", "", sample_name)]

# blanks = blanks[ sample_type == "Blank",]

blanks = merge(blanks, plates, by = "sample_name", all.x = T)

blanks[sample_type == "Blank", 
       plate := gsub("P(.)_.+", "plate \\1",sample_barcode)]

blanks[ ,blank := ifelse(sample_type == "Blank",
                        "Blank",
                        "Sample")]

b_abund = transpose(b_abund,
                    make.names = "V1",
                    keep.names = "sample_name")

b_mat = as.matrix(b_abund, rownames = "sample_name")
b_mat = b_mat[blanks$sample_name,]
b_mat = b_mat[ , colSums(b_mat) > 0]

# convert to relative abundance
#b_rel_mat = apply(b_mat, 2, function(x) x/sum(x))

# dist
b_d = vegdist(b_mat, method = "bray")


# ordination
ord = metaMDS(comm = b_d)
ord_out = list(
  points = as.data.table(ord$points,
                         keep.rownames = "sample_name"),
  stress = round(ord$stress, 3)
)

points = ord_out$points

ord_meta = merge(points,
                 blanks,
                 by = "sample_name")

# plot
set.seed(2021)
b_p = ggplot(ord_meta, aes(x = MDS1,
                     y = MDS2,
                     color = plate,
                     shape = sample_type)) +
  geom_point(size = 3) +
  scale_shape_manual(values = c(19,23,24,10),name = "Sample Type")+
  scale_color_manual(values = c("indianred3",
                                "dodgerblue3"),
                     name = "Extraction Plate")+
  labs(caption = paste("Stress =", ord_out$stress))




# adonis
plates_perm = adonis2(b_mat ~ plate + sample_type * site_name,
        data = ord_meta[match(row.names(b_mat), sample_name)])

write.csv(plates_perm, "data/revision/batch_effects_PERMANOVA.csv")



