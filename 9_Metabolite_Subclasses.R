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
sample_type_cols <- c(CCA   = "#6469ed",
                      Coral = "#e49c4c",
                      Limu  = "#7cc854" )

# Load Cleaned Metabolite Data----------------------------------------------------------------
chem_phy         <- readRDS("data/processed/chem_phy.rds")

# 13. Statistics on Metabolite Subclass by Sample Type -----------------------------------------------
# For each metabolite subclass, run anova and look for fold changes in abundance between sample types
# This information can be used to confirm that subclasss actually associate with a given sample type.
# Results inform subclass ordinations, cytoscape outputs, etc. with simple linear statistics. 

# pull out metabolite peak data
peak_data <- as.data.frame(
  as(tax_table(chem_phy),
     "matrix") )

# sum relative abundance of features in each subclass
# phyloseq treats taxonomy as hierarchical, so we need to drop most columns
# we only need 'Qemistree_sbuclass', which is the subclass
# we can also keep higher level classifications 'class' and 'superclass'
tax_mat <- as(tax_table(chem_phy),
              "matrix")
tax_mat <- tax_mat[ , c("Qemistree_superclass",
                        "Qemistree_class","Qemistree_subclass")]

sub_tax_table <- tax_table(tax_mat)

sub_phy <-phyloseq(sample_data(chem_phy),
                   otu_table(chem_phy),
                   sub_tax_table)

# sum relative areas within a subclass
sub_merge <- tax_glom(sub_phy,
                          "Qemistree_subclass")

# the column 'Qemistree_subclass' identifies the subclass
unique_subclasses <- unique(as(tax_table(sub_merge),
                             "matrix")[,"Qemistree_subclass"]) 

# for each subclass, subset the data and run anova 
subclass_data <- list()
for(a_subclass in unique_subclasses) {
  
  # subset compounds by subclass
  subclass_phy <- subset_taxa(sub_merge,
                             Qemistree_subclass == a_subclass)
  
  # create a map identifying samples by sample types
  sample_map <- setNames(sample_data(sub_phy)[["sample_type"]],
                         sample_names(sub_phy))
  sample_types <-c("Limu",
                   "CCA",
                   "Coral")
  
  # get mean relative abundance for each sample type present in the subclass
  mean_RAs <- list() 
  for(a_sample_type in sample_types){
    mean_RAs[[a_sample_type]] <- mean(
      otu_table(
        subset_samples(subclass_phy,
                       sample_type == a_sample_type)))
  }
  
  # prior to modeling, normalize data by transforming relative abundances to arcsin(sqrt(x))
  arcsin_phy <- transform_sample_counts(subclass_phy,
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
 
  # count how many samples the family was observed in
  n_samp_obs <- nrow(sample_table[sample_table$abundance > 0,])
  
  # run anova on sample type and get relevant outputs
  aov_out <-  aov(abundance ~ sample_type, data = sample_table)
  aov_sum <- summary(aov_out)
  
  p_val  <- aov_sum[[1]][["Pr(>F)"]][[1]]
  sum_sq <- aov_sum[[1]]["Sum Sq"][[1,1]]
  f_val  <- aov_sum[[1]]["F value"][1,1]
  
  # run tukey HSD to test anova results
  tuk_out <-TukeyHSD(aov_out)
  
  # assemble the output as a list
  # output list shows for each subclass, mean RA by sample type, log2 fold changes, anova + tukey results
  subclass_data[[a_subclass]] <- list(
    subclass = a_subclass,
    n_samp_obs = n_samp_obs,
    Limu_MA = mean_RAs$Limu,
    CCA_MA = mean_RAs$CCA,
    Coral_MA = mean_RAs$Coral,
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

# put the fold changes for all subclasss into a nice data.frame
subclass_fold_changes <- do.call(rbind, subclass_data)
subclasses <- row.names(subclass_fold_changes)

subclass_fold_changes <-as.data.frame(apply(subclass_fold_changes, 2, as.numeric))
subclass_fold_changes$subclass <- subclasses


# adjust p values using 
subclass_fold_changes$adj_p_val <- p.adjust(subclass_fold_changes$p_val, method = "BH")

# classify enrichment in a sample type
classify_DA <- function(a_subclass, results_table = subclass_fold_changes){
  results_table <- results_table[results_table$subclass == a_subclass,]
  # vector identifying enrichment in sample types 
  enriched <- c()
  # tests
  if(results_table$tuk_Coral_CCA_p  < 0.05 & results_table$Coral_MA > results_table$CCA_MA |
     results_table$tuk_Limu_Coral_p < 0.05 & results_table$Coral_MA > results_table$Limu_MA){
    enriched <- append(enriched,"Coral")
  }
  if(results_table$tuk_Coral_CCA_p < 0.05 & results_table$CCA_MA > results_table$Coral_MA |
     results_table$tuk_Limu_CCA_p  < 0.05 & results_table$CCA_MA > results_table$Limu_MA){
    enriched <- append(enriched,"CCA")
  }
  if(results_table$tuk_Limu_CCA_p   < 0.05 & results_table$Limu_MA > results_table$CCA_MA |
     results_table$tuk_Limu_Coral_p < 0.05 & results_table$Limu_MA > results_table$Coral_MA){
    enriched <- append(enriched,"Limu")
  }
  if(is.null(enriched)){
    return("NA")
  }else{
    return(paste(enriched, sep = "/", collapse = ""))
  }
}

# get differentially abundant sample types for all features
subclass_fold_changes$sample_type_DA <-
  vapply(subclass_fold_changes$subclass,
         classify_DA,
         FUN.VALUE = character(1))


write.csv(subclass_fold_changes,
          "data/processed/subclass_anova_and_fold_changes.csv",
          row.names = F)


