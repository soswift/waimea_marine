# Previous analyses focused on metabolite features and classes associated withi
# sample type (Limu, CCA, Coral). Repeat these analyses with microbes.
# Sample type association will help us pair up microbes and metabolites in biclusters.

## Load libraries --------------------------
# General utility
library(phyloseq)
library(data.table)

# source helper functions for plotting and subsetting
source("src/helper_functions.R")


# set color scheme for sample types
sample_type_cols <- c(CCA   = "#6469ed",
                      Coral = "#e49c4c",
                      Limu  = "#7cc854" )

# Load Cleaned Microbe and Metabolite Data----------------------------------------------------------------
micro_phy <- readRDS("data/processed/final_marine_phy.rds")

# Summarize percentage of reads in different microbial orders
taxs <- as.data.table(as(tax_table(micro_phy), "matrix"))
sum_taxs<- taxs[ , .(count = .N), "class"][ order(count, decreasing = T)]
sum_taxs[ , perc := count/sum(count)]
sum_taxs

# 13. Statistics on microbe families by Sample Type -----------------------------------------------
# For each microbe class, run anova and look for fold changes in abundance between sample types
# This information can be used to confirm that classes actually associate with a given sample type.
# Results inform ordinations, biclusters, etc. with simple linear statistics. 

# sum relative abundance of features in each class
sub_merge <- tax_glom(micro_phy,
                      "family")

# for each class, subset the data and run anova 
unique_families <- unique(as(tax_table(sub_merge),
                         "matrix")[,"family"]) 
family_data <- list()

for(a_family in unique_families) {
  
  # subset compounds by family
  family_phy <- subset_taxa(sub_merge,
                              family == a_family)
  
  # create a map identifying samples by sample types
  sample_map <- setNames(sample_data(micro_phy)[["sample_type"]],
                         sample_names(micro_phy))
  sample_types <-c("Limu",
                   "CCA",
                   "Coral")
  
  # get mean relative abundance for each sample type present in the family
  mean_RAs <- list() 
  for(a_sample_type in sample_types){
    mean_RAs[[a_sample_type]] <- mean(
      otu_table(
        subset_samples(family_phy,
                       sample_type == a_sample_type)))
  }
  
  # prior to modeling, normalize data by transforming relative abundances to arcsin(sqrt(x))
  arcsin_phy <- transform_sample_counts(family_phy,
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
  
  # run tukey HSD to test anova resulst
  tuk_out <-TukeyHSD(aov_out)
  
  # assemble the output as a list
  # output list shows for each family, mean RA by sample type, log2 fold changes, anova + tukey results
  
  family_data[[a_family]] <- list(
    family = a_family,
    n_samp_obs  = n_samp_obs, 
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

# put the fold changes for all families into a nice data.frame
family_fold_changes <- do.call(rbind, family_data)
families <- row.names(family_fold_changes)

family_fold_changes <-as.data.frame( apply( family_fold_changes, 2, as.numeric))
family_fold_changes$family <- families


# adjust p values using 
family_fold_changes$adj_p_val <- p.adjust(family_fold_changes$p_val, method = "BH")

# classify enrichment in a sample type
classify_DA <- function(a_family, results_table = family_fold_changes){
  results_table <- results_table[results_table$family == a_family,]
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
family_fold_changes$sample_type_DA <-
  vapply(family_fold_changes$family,
         classify_DA,
         FUN.VALUE = character(1))

# stacked barplots
# TODO: clean this up, pull out data frame and use ggplot
plot_bar(sub_merge, fill = "family", x = "sample_type")

# hackily write the lines out then find and replace * with " for a
# copy-paste list :)
family_fold_changes <- as.data.table(family_fold_changes)
top_families <- family_fold_changes[order(n_samp_obs, decreasing = T)][1:20]
writeLines(paste0("*",top_families$family,"*,"))

 fwrite(family_fold_changes,
          "data/processed/microbe_family_anova_and_fold_changes.csv")

# write out table with per-sample abundance for each microbial family ------
fam_export <- family_fold_changes[ ,.(family,
                                     n_samp_obs,
                                     sample_type_DA,
                                     adj_p_val,
                                     FC_LimuVCoral,
                                     FC_LimuVCCA,
                                     FC_CoralVCCA)]
# get class info for each family
match_tax_names <-function(name_vec,
                           from = c("kingdom","phylum","class","order","family","genus","otu"),
                           to   = c("kingdom","phylum","class","order","family","genus","otu"),
                           phy_obj = fam_glom){
  tax_tab  <- phy_obj@tax_table@.Data
  tax_tab  <- cbind(tax_tab, otu = row.names(tax_tab))
  name_key <- tax_tab[, to]
  names(name_key) <- tax_tab[ , from]
  
  return(name_key[name_vec])
  
}
 

fam_export[ , class := match_tax_names(family,
                                       from = "family",
                                       to   = "class",
                                       phy_obj = fam_glom)]

# get abundance table
fam_abund <- as.data.table(as(fam_glom@otu_table, "matrix"), keep.rownames = T)

fam_abund[ , family := match_tax_names(rn,
                                       from = "otu",
                                       to = "family",
                                       phy_obj = fam_glom)]
fam_abund[ ,rn := NULL]

# merge
fam_export <- merge(fam_export, fam_abund)

# get sample table
sample_dat <- as.data.table(as(fam_glom@sam_data, "matrix"), keep.rownames = T)
setnames(sample_dat, "rn", "SampleID")

# write out
fwrite(fam_export, "data/processed/microbe_family_abundance_marine_waimea.csv")
fwrite(sample_dat, "data/processed/microbe_family_samples_marine_waimea.csv")

