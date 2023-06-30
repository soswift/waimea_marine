# Script for generating analyses that rely on entire microbial communities and metabolite profiles
# Analyses include ordinations, permanovas, betadispersion tests, and paired ordination (Procrustes).
# All of these analyses are based on distance matrices generated from relative abundance data.

## Load libraries --------------------------
# General utility
library(phyloseq)
library(data.table)
library(tidyr)
library(biomformat)
# Plotting and analysis
library(eulerr)
library(vegan)

# library(devtools)
# install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis") 
library(pairwiseAdonis) 

library(cowplot)

# source helper functions for plotting and subsetting
source("src/helper_functions.R")

# source function for writing out phyloseq objects
source("src/physeq_csv_out.R")

# function for reading unifrac dists from flat tables
source("src/read_unifrac.R")

# plot colors
# variables =  genus_cols, sample_type_cols, site_cols
source("src/colors.R")


# Load cleaned microbe and metabolite data----------------------------------------------------------------
## microbe
final_marine_phy <- readRDS("data/processed/final_marine_phy.rds")
micro_phy <- final_marine_phy

final_unifrac <- readRDS("data/processed/final_unifrac.rds")

## metabolite
chem_phy <- readRDS("data/processed/chem_phy.rds")
chem_phy_w_blanks <- readRDS("data/processed/chem_phy_w_blanks.rds")
chem_phy_raw <- readRDS("data/processed/chem_phy_raw.rds")

# 4. All metabolite NMDS and calculate permanovas -----------------------------------------------------

# ordination of all samples and blanks to determine 'bad' samples

raw_dist <- vegdist(veganifyOTU(chem_phy_raw),
                    method = "bray")

raw_ord <- ordinate(chem_phy_raw,
                    method = "NMDS",
                    distance = raw_dist)

check_p <- plot_ordination(
  chem_phy_raw,
  raw_ord,
  color = "sample_type",
  shape = "site_name") +
  scale_color_manual(values = sample_type_cols)+
  labs(title = "All Metabolite NMDS, Bray-Curtis",
       subtitle = paste("stress =" ,signif(raw_ord$stress, 2)))+
  stat_ellipse(aes(color = sample_type, shape = NULL),level = 0.5)

ggsave("output/NMDS/raw_metabolites_NMDS.pdf",
       plot = check_p,
       width = 9,
       height = 7)


# try multiple distance algorithms
# initalize a list of distance matrices
chem_dist <- list()

# initialize a list of plots
chem_p <- list()

# NMDS for all samples, by type
chem_dist$canberra <- vegdist(veganifyOTU(chem_phy),
                              method = "canberra")
chem_dist$bray <- vegdist(veganifyOTU(chem_phy),
                          method = "bray")
chem_dist$euclidean <- vegdist(veganifyOTU(chem_phy),
                               method = "euclidean")

# run NMDS on all three distance methods and compare plots
for(a in seq_along(chem_dist) ){
  
  dist_name <- names(chem_dist[a])
  
  chem_ord  <- ordinate(chem_phy,
                        method = "NMDS",
                        distance = chem_dist[[a]])
  
  chem_p[[a]] <- plot_ordination(
    chem_phy,
    chem_ord,
    color = "sample_type",
    shape = "site_name",
    title = paste0("All Metabolite NMDS, ", dist_name)) +
    scale_color_manual(values = sample_type_cols)+
    labs(subtitle = paste("stress =" ,signif(chem_ord$stress, 2)))+
    stat_ellipse(aes(color = sample_type, shape = NULL),level = 0.5)
  
}

# plot bray curtis distance for main figure
ggsave("output/NMDS/all_metabolites_NMDS.pdf",
       plot = chem_p[[2]],
       width = 10,
       height = 7)

g <-arrangeGrob(chem_p[[1]], chem_p[[2]], chem_p[[3]],
                nrow = 1, ncol = 3)

ggsave(paste0("output/NMDS/all_metabolites_NMDS_dists.pdf"),
       plot = g,
       width = 15,
       height = 5)

# run NMDS on metabolites with blanks included
dist_blanks <- vegdist(veganifyOTU(chem_phy_w_blanks),
                       method = "bray")

chem_ord_w_blanks  <- ordinate(chem_phy_w_blanks,
                      method = "NMDS",
                      distance = dist_blanks)

blank_p <- plot_ordination(chem_phy_w_blanks,
                          chem_ord_w_blanks,
                          color = "sample_type",
                          shape = "site_name",
                          title = "All Metabolite NMDS, bray") +
                          scale_color_manual(values = c("darkgray",sample_type_cols))+
                          labs(subtitle = paste("stress =" ,signif(chem_ord_w_blanks$stress, 2)))
                      
ggsave("output/NMDS/all_metabolites_NMDS_w_blanks.pdf",
       plot = blank_p)

# run permanovas on metabolites, modeling by sample type, site, and the interaction
chem_anova_all <-
  do_permanova(
    chem_phy,
    var_name = "sample_type*site_name",
    dist_obj = chem_dist$bray,
    description = "All microbe samples"
  )
write.csv(chem_anova_all,
          "output/permanova/all_metabolite_permanova_by_site_type.csv")


# run type 2 pairwise permanova for all samples just by sample type
chem_anova_pairwise <- pairwise.adonis2(chem_dist$bray ~ sample_type,
                                        data = as(sample_data(chem_phy),
                                                  "data.frame"))
write.csv(chem_anova_pairwise,
          "output/permanova/all_metabolite_pairwise_permanova_by_type.csv")

# calculate beta dispersion by sample type 
chem_sam_dat <- as(sample_data(chem_phy), "data.frame")

chem_groups <- chem_sam_dat[match(labels(chem_dist$bray),
                                  chem_sam_dat$sample_barcode), 
                            "sample_type"]

chem_betadisp <- betadisper(chem_dist$bray, chem_groups)

# write out sample type beta dispersion summary and tukeyHSD test for significance
sink("output/permanova/all_metabolite_betadispersion.txt" )
print(chem_betadisp)
print(TukeyHSD(chem_betadisp))
sink()

# write out distance matrices 
write.csv(as.matrix(chem_dist$bray),
          "data/processed/table_exports/all_marine_metabolite_bray_dist.csv")
write.csv(as.matrix(dist_blanks),
          "data/processed/table_exports/all_marine_metabolite_w_blanks_bray_dist.csv")

# 5. Metabolite by sample type NMDS and permanova --------------------------------------------
# Subset metabolite samples to a single sample type (Limu, Coral, CCA) for within sample type statistical tests and plots.
# Within samples of a given sample type, does the host genus or sampling site correlate with metabolite composition?

sample_types <- unique(chem_phy@sam_data$sample_type)

chem_type_p <- list()

chem_anovas <- list()

for(i in sample_types){
  # subset
  a_phy  <- subset_samples(chem_phy, sample_type == i)
  a_phy  <- prune_taxa(taxa_sums(a_phy) > 0, a_phy)
  
  a_dist <- subset_dist(chem_dist[[1]], a_phy)
  
  # plot by site
  a_ord  <- ordinate(a_phy,
                     method = "NMDS",
                     distance = a_dist)
  
  chem_type_p[[paste0(i,"_site")]] <- plot_ordination(a_phy, a_ord, color = "site_name") +
                                                     stat_ellipse( level = 0.5) +
                                                     theme(aspect.ratio = 1, legend.key.size = unit(0.5, "cm")) +
                                                     # standard site colors
                                                     scale_color_manual(values = site_cols) 
  
  # plot by genus
  chem_type_p[[paste0(i,"_genus")]] <- plot_ordination(a_phy, a_ord, color = "genus") +
                                                    stat_ellipse( level = 0.5) +
                                                    theme(aspect.ratio = 1, legend.key.size = unit(0.5, "cm"))+
                                                    scale_color_manual(values = genus_cols)
  
  # calcluate permanovas for site
  chem_anovas[[paste0(i,"_site")]]<- do_permanova(a_phy, 
                                  var_name = "site_name",
                                  dist_obj = a_dist, 
                                  description =  paste0(i, " site metabolite samples"))
  if(i != "CCA"){
   a_phy <- subset_samples(a_phy, genus != "Other")
   a_dist <- subset_dist(a_dist, a_phy)
   chem_anovas[[paste0(i,"_genus")]]<- do_permanova(a_phy, 
                                  var_name = "genus",
                                  dist_obj = a_dist, 
                                  description =  paste0(i, " genus metabolite samples"))
}
}

# write out permanova results and arranged plots
write.csv(bind_rows(chem_anovas),
          "output/permanova/sample_type_metabolite_permanova_by_site_and_genus.csv")

# convert plot to grob
chem_type_g <-lapply(chem_type_p, ggplotGrob)

# set a standard width for all grobs
std_width <- chem_type_g[[1]]$widths
chem_type_std <- lapply(chem_type_g, function(x) {
                  x$widths <- std_width
                  return(x)})
# arrange ordinations
g <- ggarrange(chem_type_std[[1]], chem_type_std[[3]], chem_type_std[[5]],
                chem_type_std[[2]], chem_type_std[[4]], chem_type_std[[6]],
                nrow = 2, ncol = 3)


ggsave(paste0("output/NMDS/sample_type_metabolites_NMDS.pdf"),
       plot = g, 
       width = 15,
       height = 5)

# 6. All microbe NMDS and permanova ---------------------------------------------------------
# Run all the same analyses (permanova, NMDS) on the microbial data.
# Also test for differences in alpha diversity across sample types and sites.

# unifrac distances
micro_dist <- final_unifrac

# NMDS for all samples, by type
micro_ord  <- ordinate(micro_phy,
                       method = "NMDS",
                       distance = micro_dist)

plot_ordination(micro_phy, 
                micro_ord, 
                color = "sample_type",
                shape = "site_name",
                title = "All Microbe NMDS, Unifrac") +
  scale_color_manual(values = sample_type_cols)+
  stat_ellipse(mapping = aes(shape = NULL))

ggsave("output/NMDS/all_microbes_NMDS.pdf",width = 7, height = 5)


# run permanova on all microbes by sample type, site, and the interaction
micro_anova_all <- do_permanova(micro_phy,
                                  var_name = "sample_type*site_name",
                                  micro_dist,
                                  description = "all microbe samples")
write.csv(micro_anova_all, "output/permanova/all_microbe_permanova_by_site_type.csv")

micro_anova_pairwise <- pairwise.adonis2(micro_dist ~ sample_type,
                                         data = as(sample_data(micro_phy),
                                                   "data.frame"))
write.csv(micro_anova_pairwise, "output/permanova/all_microbe_pairwise_permanova_by_type")

# run betadispersion by sample type on all microbe samples
micro_sam_dat <- as(sample_data(micro_phy), "data.frame")

micro_groups <- micro_sam_dat[match(labels(micro_dist),
                                  micro_sam_dat$sequencing_id), 
                            "sample_type"]

micro_betadisp <- betadisper(micro_dist, micro_groups)

# write out beta dispersion summary and tukey test resulst
sink("output/permanova/all_microbe_betadispersion.txt" )
print(micro_betadisp)
print(TukeyHSD(micro_betadisp))
sink()

# read in alpha diversity indices for the microbial communities
# these were calculated as part of the metaflowmics bioninformatics pipeline
micro_div <-read.table("data/raw/diversity/all_alphaDiversity_100.summary", header = T)

micro_div <- rename(micro_div,   "sequencing_id" ="group" )

micro_div <- merge(micro_div, as.data.frame(as.matrix(sample_data(micro_phy))), by = "sequencing_id")

# model alpha diversity as a function of sample type (limu, coral, CCA), site, and the interaction
chao_lm <- lm(chao ~sample_type*site_name, data = micro_div)
shan_lm <- lm(shannoneven ~ sample_type*site_name, data = micro_div)

chao_anova <- anova(chao_lm)
shan_anova <- anova(shan_lm)

sink("output/alphadiversity/alpha_div_by_type_and_site.txt")
print(chao_anova)
print(shan_anova)
sink()


# 7. Microbe by sample type NMDS and permanovas ---------------------------
# Subset microbial samples to a single sample type (Limu, Coral, CCA) for within sample type statistical tests and plots.
# Within samples of a given sample type, does the host genus or sampling site correlate with microbial composition?

# initialize list of plots and permanovas for microbial samples types
micro_type_p <- list()

micro_anovas <- list()

for(i in sample_types){
  
  # subset
  a_phy  <- subset_samples(micro_phy, sample_type == i)
  a_phy  <- prune_taxa(taxa_sums(a_phy) > 0, a_phy)
  
  # get distances
  a_dist <- subset_dist(micro_dist, a_phy)
  set.seed(2020)
  # plot by site
  a_ord  <- ordinate(a_phy,
                     method = "NMDS",
                     distance = a_dist)
  
  micro_type_p[[paste0(i,"_site")]] <- plot_ordination(a_phy, a_ord, color = "site_name") +
                                                      stat_ellipse( level = 0.5) +
                                                      theme(aspect.ratio = 1, legend.key.size = unit(0.5, "cm")) +
                                                      # standard site colors
                                                      scale_color_manual(values = site_cols) 
                                                    
  # plot by genus
  micro_type_p[[paste0(i,"_genus")]] <- plot_ordination(a_phy, a_ord, color = "genus") +
                                                      stat_ellipse( level = 0.5) +
                                                      theme(aspect.ratio = 1, legend.key.size = unit(0.5, "cm"))+
                                                      scale_color_manual(values = genus_cols)
                                                                                  
  # run permanovas for site
  micro_anovas[[paste0(i,"_site")]] <- do_permanova(
                                        a_phy, 
                                        var_name = "site_name",
                                        dist_obj = a_dist, 
                                        description =  paste0(i, " site microbe samples"))

  # run permanovas for genus (but not for CCA, which lacks genera information)
  if(i != "CCA"){
  a_phy <- subset_samples(a_phy, genus != "Other")
  a_dist <- subset_dist(a_dist, a_phy)
  micro_anovas[[paste0(i,"_genus")]] <- do_permanova(
                                        a_phy, 
                                        var_name = "genus",
                                        dist_obj = a_dist, 
                                        description =  paste0(i, " genus microbe samples"))
  }
  
}

# convert plot to grob
micro_type_g <-lapply(micro_type_p, ggplotGrob)

# set a standard width for all grobs
std_width <- micro_type_g[[1]]$widths
micro_type_std <- lapply(micro_type_g, function(x) {
  x$widths <- std_width
  return(x)})

# arrange plots and write out results
g <-arrangeGrob(micro_type_std[["Limu_site"]], micro_type_std[["Coral_site"]], micro_type_std[["CCA_site"]],
                micro_type_std[["Limu_genus"]], micro_type_std[["Coral_genus"]], micro_type_std[["CCA_genus"]],
                nrow = 2, ncol = 3)

ggsave("output/NMDS/sample_type_microbes_NMDS.pdf",
       plot = g, width = 15, height = 5)

write.csv(bind_rows(micro_anovas),
          "output/permanova/sample_type_microbe_permanova_by_site_and_genus.csv")


# 8. Pair up microbes and metabolites ------------------------------------------
# Run paired analyses on microbes and metabolites, matched up by the physical samples from which data were collected
# These analyses are run on a subset of the total samples where both metabolite and microbial data were available.

# find mismatches between microbial and metabolite samples
micro_sample_bcodes <- micro_phy@sam_data$sample_barcode
print( paste("Total micro samples = ", length(micro_sample_bcodes)))

chem_sample_bcodes  <- chem_phy@sam_data$sample_barcode
print( paste("Total chem samples = ", length(chem_sample_bcodes)))

chem_no_micro <- chem_phy@sam_data[ ! (chem_phy@sam_data$sample_barcode %in% micro_sample_bcodes) ] %>%
  as.matrix() %>%
  as.data.frame()

print( paste("In chem data, but not in micro = ", nrow(chem_no_micro) ))

micro_no_chem <- micro_phy@sam_data[ ! (micro_phy@sam_data$sample_barcode %in% chem_sample_bcodes)]%>%
  as.matrix() %>%
  as.data.frame()

print( paste("In micro data, but not in chem = ", nrow(micro_no_chem)))

# subset both datasets to the samples that overlap
pair_chem_phy <- subset_samples(chem_phy,
                                sample_barcode %in% micro_phy@sam_data$sample_barcode)

pair_micro_phy  <- subset_samples(micro_phy,
                                  sample_barcode %in% chem_phy@sam_data$sample_barcode)

# check number samples in both data sets
nsamples(pair_chem_phy)
nsamples(pair_micro_phy)

# switch microbial barcodes over to metabolite barcodes so that everything has a unified identifier
## pull out otu table and sample data
otus    <- pair_micro_phy@otu_table
samples <- pair_micro_phy@sam_data
taxa    <- pair_micro_phy@tax_table

## change sample ids to metabolite barcodes
samples$sequence_id <- row.names(samples)
row.names(samples)  <- samples$sample_barcode
colnames(otus)      <- samples$sample_barcode

pair_micro_phy  <- phyloseq(samples, otus, taxa)

# update the microbial sample ids on the unifrac distance matrix
pair_micro_dist <- as.matrix(micro_dist)[samples$sequence_id, samples$sequence_id]
colnames(pair_micro_dist) <- samples$sample_barcode
row.names(pair_micro_dist) <- samples$sample_barcode

pair_micro_dist <- as.dist(pair_micro_dist)

rm(samples, otus, taxa)

# re-generate a bray curtis distance matrix for the paired metabolite samples
pair_chem_dist <- vegdist(veganifyOTU(pair_chem_phy), method = "bray")

## check that the microbial and metabolite names match
if (all( sample_names(pair_micro_phy) %in% sample_names(pair_chem_phy) ) ){
  message("chem and micro sample names match")
} else (
  warning("chem and micro sample names are different")
)

# write out paired flat tables
# microbes
physeq_csv_out(pair_micro_phy, 
               description = "paired_marine_microbe",
               outdir = "data/processed/table_exports")

write.csv(as.matrix(pair_micro_dist), 
          "data/processed/table_exports/paired_marine_microbe_unifrac_dist.csv")

saveRDS(pair_micro_phy, "data/processed/table_exports/paired_chem_phyloseq.rds")
write_biom(make_biom(data = otu_table(pair_micro_phy)), "data/processed/paired_microbe.biom")

# metabolites
physeq_csv_out(pair_chem_phy,
               description = "paired_marine_metabolite",
               outdir = "data/processed/table_exports")

write.csv(as.matrix(pair_chem_dist),
          "data/processed/table_exports/paired_marine_metabolite_bray_dist.csv")

saveRDS(pair_chem_phy, "data/processed/table_exports/paired_chem_phyloseq.rds")
write_biom(make_biom(data = otu_table(pair_chem_phy)), "data/processed/paired_metabolite.biom")


# 9. NMDS and Mantel of microbes and metabolites ----------------------------------------------
# Perform paired ordinations, procrustes rotations, and mantel tests.
# Do metabolite profiles and microbial communities tell the same story?

# Initialize plot list
p <- list()

# Ordinate all samples
p$all <- paired_ordination(
  microbe_phy = pair_micro_phy,
  chem_phy = pair_chem_phy,
  description = "All",
  shape = "site_name",
  color = "sample_type"
)

# Separately ordinate samples from each sample type
sample_types <- unique(pair_chem_phy@sam_data$sample_type)

for(i in sample_types){
  # subset
  c_phy  <- subset_samples(pair_chem_phy, sample_type == i)
  c_phy  <- prune_taxa(taxa_sums(c_phy) > 0, c_phy)
  
  m_phy  <- subset_samples(pair_micro_phy, sample_type == i)
  m_phy   <- prune_taxa(taxa_sums(m_phy) > 0, m_phy)
  
  p[[i]] <- paired_ordination(
    microbe_phy = m_phy,
    chem_phy = c_phy,
    description = i,
    shape = "site_name",
    color = "genus"
  )
}

# arrange plots and write out
grid.arrange(p$all$micro      , p$all$chem,     p$all$proc,
             p$Limu$micro     , p$Limu$chem,    p$Limu$proc,
             p$Coral$micro    , p$Coral$chem,   p$Coral$proc,
             p$CCA$micro      , p$CCA$chem,     p$CCA$proc,
             nrow = 4, ncol = 3
)

g <-arrangeGrob(p$all$micro      , p$all$chem,     p$all$proc,
                p$Limu$micro     , p$Limu$chem,    p$Limu$proc,
                p$Coral$micro    , p$Coral$chem,   p$Coral$proc,
                p$CCA$micro      , p$CCA$chem,     p$CCA$proc,
                nrow = 4, ncol = 3
)

ggsave("output/NMDS/combined_microbe_metabolite_procrust.png", g, width = 15, height = 20)

# 10. Identify outlier CCA -----------------------------------------------------
# Some CCA samples are obvious outliers based on ordinations
# Identify these outliers for exclusion from subsequent ordinations.

cca_p <- list()

# Plot metabolite ordination with CCA labelled
chem_ord  <- ordinate(chem_phy,
                      method = "NMDS",
                      distance = chem_dist[[2]])

cca_p[["chem"]] <- plot_ordination(
                                  chem_phy,
                                  chem_ord,
                                  color = "sample_type",
                                  shape = "site_name",
                                  title = paste0("All Metabolite NMDS, ", "bray")
                                ) + 
                                  # standard sample type colors
                                  scale_color_manual(values = sample_type_cols) + 
                                  geom_text(aes(label = ifelse(sample_type == "CCA",sample_barcode,"")))


# Plot microbe ordination with CCA labelled
micro_ord  <- ordinate(micro_phy,
                       method = "NMDS",
                       distance = micro_dist)

cca_p[["micro"]] <- plot_ordination(
  micro_phy,
  micro_ord,
  color = "sample_type",
  shape = "site_name",
  title = "All Microbe NMDS, Unifrac"
) +
  scale_color_manual(values = sample_type_cols) +
  geom_text(aes(label = ifelse(sample_type == "CCA", sample_barcode, "")))


g <-arrangeGrob(grobs = cca_p)

ggsave("output/NMDS/CCA_outlier_comparison.pdf", plot = g)

# plot paired ordination with all CCA removed

no_cca_micro <- subset_samples(pair_micro_phy, sample_type != "CCA")
no_cca_chem <- subset_samples(pair_chem_phy, sample_type != "CCA")

p <- paired_ordination(
  microbe_phy = no_cca_micro,
  chem_phy = no_cca_chem,
  description = "All",
  shape = "site_name",
  color = "sample_type",
  unifrac_dist = pair_micro_dist
)

g <- arrangeGrob(p$micro, p$chem, p$proc, ncol = 3, nrow = 1)

ggsave("output/NMDS/no_cca_procrustes.pdf", plot = g, width = 15)

# 11. All microbe and metabolite heatmaps ------------------------------------------------------
# Make large heatmaps of all microbes and metabolites by sample as an exploratory analysis.

# plot_heatmap(phyloseq_obj = micro_phy, description = "microbe")
# 
# plot_heatmap(phyloseq_obj = chem_phy, description = "metabolite",
#              dist_method = "canberra")

# 12. Metabolite and Microbe Eulerr diagrams -----------------------------------------------
# For corals and algae, generate euler plots to show overlap in composition by sample type.
# Within sample types, generate euler plots to show overlap by genus

# merge all samples by sample type (coral, limu, cca)
all_chem_eul  <- euler_subset(chem_phy, group_by = "sample_type")

all_micro_eul <- euler_subset(micro_phy, group_by = "sample_type")

# plot overlap by sample type

pdf("output/Euler/Metabolite_Sample_Type_Euler.pdf")
print(plot(all_chem_eul, fills = sample_type_cols, quantities = T))
dev.off()

pdf("output/Euler/Microbe_Sample_Type_Euler.pdf")
print(plot(all_micro_eul, fills = sample_type_cols, quantities = T))
dev.off()



# for each sample type (coral/limu) and data type (metabolite, microbe), plot overlap by genus
eul_plots <- list()

# metabolite
eul_plots$chem_coral <- euler_subset(chem_phy, group_by = "genus",
                                     sample_type == "Coral")
eul_plots$chem_limu  <- euler_subset(chem_phy, group_by = "genus",
                                     sample_type == "Limu" & genus != "Other")
# microbe
eul_plots$micro_coral <- euler_subset(micro_phy, group_by = "genus",
                                        sample_type == "Coral")
eul_plots$micro_limu <- euler_subset(micro_phy, group_by = "genus",
                                        sample_type == "Limu" & genus != "Other")
# write out
for(p in names(eul_plots)){
  pdf(paste0("output/Euler/", p, "_genus_euler.pdf"))
  print(plot(eul_plots[[p]],fills = genus_cols, quantities = T))
  dev.off()
}

