## Load libraries --------------------------
# General utility
library(phyloseq)
library(data.table)
library(tidyr)
library(biomformat)
# Plotting and analysis
library(vegan)
library(pairwiseAdonis)

# source helper functions for plotting and subsetting
source("src/helper_functions.R")


# set color scheme for sample types
sample_type_cols <- c(CCA = "#6469ed",
                      Coral = "#e49c4c",
                      Limu = "#7cc854" )

# Load Cleaned Microbe and Metabolite Data----------------------------------------------------------------
final_marine_phy <- readRDS("data/processed/final_marine_phy.rds")
micro_phy        <- final_marine_phy
final_unifrac    <- readRDS("data/processed/final_unifrac.rds")
chem_phy         <- readRDS("data/processed/chem_phy.rds")

# 13. Statistics on Metabolite Networks by Sample Type -----------------------------------------------
# For each metabolite network, run anova and look for fold changes in abundance between sample types
# This information can be used to confirm that networks actually associate with a given sample type.
# Results inform network ordinations, cytoscape outputs, etc. with simple linear statistics. 

# pull out metabolite peak data
peak_data <- as.data.frame(
                      as(tax_table(chem_phy),
                         "matrix") )

# sum relative abundance of features in each network
# phyloseq treats taxonomy as hierarchical, so we need to drop most columns
# we only need 'componentindex', which is the network and 'cluster index' which is a unique id

tax_mat <- as(tax_table(chem_phy),
                  "matrix")
tax_mat <- tax_mat[ , c("componentindex","cluster index")]

net_tax_table <- tax_table(tax_mat)

net_phy <-phyloseq(sample_data(chem_phy),
                   otu_table(chem_phy),
                   net_tax_table)

# sum relative areas within a network
network_merge <- tax_glom(net_phy,
                            "componentindex")

# remove networks with only a single component
network_merge <- subset_taxa(network_merge,
                             componentindex != "  -1")

# the column componentindex identifies the network
unique_networks <- unique(as(tax_table(network_merge),
                             "matrix")[,"componentindex"]) 

# for each network, subset the data and run anova 
net_data <- list()
for(a_network in unique_networks) {

 network_phy <- subset_taxa(network_merge,
                             componentindex == a_network)
  
  # create a map identifying samples by sample types
  sample_map <- setNames(sample_data(network_phy)[["sample_type"]],
                         sample_names(network_phy))
  sample_types <-c("Limu",
                   "CCA",
                   "Coral")
  
  # get mean relative abundance for each sample type present in the network
  mean_RAs <- list() 
  for(a_sample_type in sample_types){
  mean_RAs[[a_sample_type]] <- mean(
                                 otu_table(
                                 subset_samples(network_phy,
                                                sample_type == a_sample_type)))
  }
  
  # prior to modeling, normalize data by transforming relative abundances to arcsin(sqrt(x))
  arcsin_phy <- transform_sample_counts(network_phy,
                                          fun = function(x)
                                          asin(sqrt(x)) )
  # pull out transformed counts of sample_types, clean up
  sample_abunds <- data.frame(t( as(otu_table(arcsin_phy),
                        "matrix")))
  sample_abunds[is.na(sample_abunds)] <- 0
  
  # make table with columns for sample type and abundance
  sample_table <- data.frame(sample_type = sample_map[row.names(sample_abunds)],
                             abundance = sample_abunds[[1]]
                             )

  # run anova on sample type and get relevant outputs
  aov_out <-  aov(abundance ~ sample_type, data = sample_table)
  aov_sum <- summary(aov_out)
  
  p_val  <- aov_sum[[1]][["Pr(>F)"]][[1]]
  sum_sq <- aov_sum[[1]]["Sum Sq"][[1,1]]
  f_val  <- aov_sum[[1]]["F value"][1,1]
  
  # run tukey HSD to test anova resulst
  tuk_out <-TukeyHSD(aov_out)
  
  # assemble the output as a list
  # output list shows for each network, mean RA by sample type, log2 fold changes, anova + tukey results
  net_data[[a_network]] <- list(
           network = a_network,
           limu_RA = mean_RAs$Limu,
           cca_RA = mean_RAs$CCA,
           coral_RA = mean_RAs$Coral,
           FC_LimuVCoral = fc(mean_RAs$Limu, mean_RAs$Coral),
           FC_LimuVCCA = fc(mean_RAs$Limu, mean_RAs$CCA),
           FC_CoralVCCA = fc(mean_RAs$Coral, mean_RAs$CCA),
           tuk_Coral_CCA_diff  = tuk_out[[1]]["Coral-CCA","diff"],
           tuk_Coral_CCA_p     = tuk_out[[1]]["Coral-CCA","p adj"],
           tuk_Limu_CCA_diff   = tuk_out[[1]]["Limu-CCA","diff"],
           tuk_Limu_CCA_p      = tuk_out[[1]]["Limu-CCA","p adj"],
           tuk_Limu_Coral_diff = tuk_out[[1]]["Limu-Coral","diff"],
           tuk_Limu_Coral_p    = tuk_out[[1]]["Limu-Coral","p adj"],
           f_stat = f_val,
           sum_of_sq = sum_sq,
           p_val = p_val
           )

}

# put the fold changes for all networks into a nice data.frame
# note: this does not include networks that only showed up in one sample
network_fold_changes <- do.call(rbind, net_data)
network_fold_changes <-as.data.frame(apply(network_fold_changes, 2, as.numeric))

# adjust p values using 
network_fold_changes$adj_p_val <- p.adjust(network_fold_changes$p_val, method = "BH")

write.csv(network_fold_changes,
          "data/processed/network_anova_and_fold_changes.csv",
          row.names = F)


