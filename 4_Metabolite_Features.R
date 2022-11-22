# 15. LM on Random Forest Filtered Metabolits --------------------------------
# Identify metabolite features that are highly correlated with sample types
# Since there are many features to test, we use several methods to test for sample type association.
# Random forests rapidly assess the importance of many features, but focuses on predicitve value
# Linear models suffer from multiple testing bias and require CLR transformation and zero handling to account for the compositionality of the data
# ANCOM is a more rigorous method of linear modeling of compositiional data. 

# Random Forest:
# The file 'run_RF.R' was run on a computing cluster
# random.Forest was run on the cleaned metabolite data and produced variable importance scores for each metabolite.
# Importance scores indicate that the presence of the metabolite increased the predictive ability of the random forest classification model.


# Linear models on each metabolite:
# As a complementary approach to random forest, linear models were run on all metabolite features.
# To account for compositionality, data were clr transformed. Tukey HSD was used to determine differences in metabolite relative abundance between sample types.

# read in metabolite abundance (relative abundance of peak areas) and sample data (information on samples)
# abundance of metabolites will be used to predict sample type (Limu, Coral, CCA)

library(data.table)
library(ggplot2)
library(compositions)
source("src/helper_functions.R")

metabolite_abundance_file <- "data/processed/table_exports/all_marine_metabolite_abundance_flat_table.csv"
sample_data_file          <- "data/processed/table_exports/all_marine_metabolite_sample_flat_table.csv"
metabolite_metadata_file  <- "data/processed/table_exports//all_marine_metabolite_tax_flat_table.csv"

abund_clean <- read.csv(metabolite_abundance_file, header = T, row.names = 1)
abund_clean <- t(abund_clean)

sam_dat   <- read.csv(sample_data_file, header = T, row.names = 1)
sam_dat$sample_barcode <- paste0("X", sam_dat$sample_barcode)

chem_dat  <- read.csv(metabolite_metadata_file,  header = T, row.names = 1)
chem_dat$featureID <- row.names(chem_dat)


# Read in random forest results (Variable Importance scores)
sample_type_RF <- readRDS("data/processed/sample_type_rf.rds")

# Pull out importance scores for each metabolite
# randomForest scores
rf_imp_scores <- as.data.table(sample_type_RF$importance, keep.rownames = T)
setnames(rf_imp_scores, "rn", "featureID")

# identify 'high scoring' metabolites based on randomForest scores higher than 2 Standard Deviation from mean
rf_imp_scores[ ,RF_selection :=
                 ifelse(MeanDecreaseAccuracy > (mean(MeanDecreaseAccuracy) + 2*sd(MeanDecreaseAccuracy)),
                 "Important", "Not Important"),]

# write out
fwrite(rf_imp_scores, "data/processed/all_metabolite_RF_importance.csv")
write.csv(sample_type_RF$localImportance, "data/processed/all_metabolite_RF_local_importance.csv")

# CLR Transform data prior to modelling -------------------------------------------------------------
# Transform abundances out of the simplex
# add a tiny number to all values in the abundance table
minval <- min(abund_clean[abund_clean > 0])/2
abund_clean <- abund_clean + minval

# perform CLR transformation
norm_abund <-matrix(compositions::clr(abund_clean),
                    nrow = nrow(abund_clean),
                    ncol = ncol(abund_clean),
                    dimnames = dimnames(abund_clean))

write.csv(norm_abund, "data/processed/CLR_transformed_metabolite_abundance.csv")

# norm abund is a matrix with samples as rows and features as columns
norm_abund <- as.data.frame(norm_abund)
norm_abund$sample_type <- as.character(sam_dat$sample_type)
norm_abund$site_name   <- as.character(sam_dat$site_name)


# run_lm() runs linear models on each feature with sample type as independent variable
run_lm <-function(feat_name, lm_abund){
  formula_call <- as.formula(paste0(feat_name, " ~ sample_type")) 
  
  lm_out <- lm(formula_call,
               data = lm_abund)
  return(lm_out)
}

# run linear models on all features
all_feats <- colnames(abund_clean)
lm.all <- lapply( all_feats,
                  run_lm,
                  lm_abund = norm_abund)
names(lm.all) <- all_feats

# summarize_lm() takes a feature name and the output from run_lm() and generates summary stats for that feature
# when run with lapply, generates a list of stats for the provided features

summarize_lm <- function(i, lm.out){
  
  # summarize model output
  lm_sum <- summary(aov(lm.out[[i]]))
  
  # run tukey HSD
  feat_tuk_out <- TukeyHSD(aov(lm.out[[i]]))
  
  summary_out <-
    data.frame(
      featureID = i,
      p_val  = lm_sum[[1]][["Pr(>F)"]][[1]],
      sum_sq = lm_sum[[1]]["Sum Sq"][[1,1]],
      f_val  = lm_sum[[1]]["F value"][1,1],
      tuk_Coral_CCA_diff  = feat_tuk_out[[1]]["Coral-CCA","diff"],
      tuk_Coral_CCA_p     = feat_tuk_out[[1]]["Coral-CCA","p adj"],
      tuk_Limu_CCA_diff   = feat_tuk_out[[1]]["Limu-CCA","diff"],
      tuk_Limu_CCA_p      = feat_tuk_out[[1]]["Limu-CCA","p adj"],
      tuk_Limu_Coral_diff = feat_tuk_out[[1]]["Limu-Coral","diff"],
      tuk_Limu_Coral_p    = feat_tuk_out[[1]]["Limu-Coral","p adj"]
    )
  return(summary_out)
}


# get results for all features
all_lm <- lapply(all_feats,
                 summarize_lm,
                 lm.out = lm.all)

all_results <- do.call("rbind", all_lm)

all_results <- merge(all_results,
                     all_metabolite_scores,
                     all.x = T, by = "featureID")

# adjust p values when considering all features
all_results$adj_p_val <- p.adjust(all_results$p_val, method = "BH")


pdf(file = "output/randomforest/pvalue_vs_RF.pdf")
plot(x= log10(all_results$adj_p_val), y = all_results$VSURF_score,
     xlab =  "LM log10 P-value",
     ylab = "VSURF Importance Score")
plot.new()
plot(x = log10(all_results$adj_p_val), y = all_results$MeanDecreaseAccuracy,
     xlab =  "LM log10 P-value",
     ylab = "RandomForest MDA Score")
dev.off()

# Make log fold change plots -------------------------------------------------
# Visually compare CLR transformed abundances between sample types for the most significantly different features

# Calculate the mean CLR transformed relative abundance by sample type
mean_clr_abund <- as.data.table(norm_abund)[ , lapply(.SD, mean),
                                             by = sample_type,
                                             .SDcols = !"site_name"]
mean_clr_abund <- transpose(mean_clr_abund,
                            keep.names = "featureID",
                            make.names = "sample_type")
setnames(mean_clr_abund,
         old = c("Limu","CCA","Coral"),
         new = c("Limu_MA","CCA_MA","Coral_MA"))

all_results <- merge(all_results,
                     mean_clr_abund,
                     by = "featureID")

# classify each feature for 'enrichment' in a sample type or multiple sample types
# enrichment is based on significance of tukeyHSD test and mean clr(relative abundance)

# classify_DA() takes a single feature and identifies the sample types
# in which that feature is significantly more abundanct

classify_DA <- function(ft, results_table = all_results){
  results_table <- results_table[results_table$featureID == ft,]
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
all_results$sample_type_DA <-
  vapply(as.character(all_results$featureID),
         classify_DA,
         FUN.VALUE = character(1))

# randomly identify some features significantly correlated with sample type
sig_feats <- as.character(sample(all_results$featureID[all_results$adj_p_val < 0.05], 20))
sig_feats <- as.character(all_results[head(order(all_results$adj_p_val), n = 15), "featureID"])

# plot boxplots to visually verify classification
pdf("output/DA/example_DA.pdf")
lapply(sig_feats, function(ft){
  formula_call <- as.formula(paste0(ft, " ~ sample_type")) 
  boxplot(formula_call, data = norm_abund)
  title(sub = paste("Sig. More Abundant in", all_results[all_results$featureID == ft, "sample_type_DA"],
                    "Adj P-val = ", round(all_results[all_results$featureID == ft, "adj_p_val"], 3)))
})
dev.off()

# write out
write.csv(all_results,
          "data/processed/all_metabolites_sample_type_anovas.csv",
          row.names = F)

## Compare statistically (RF or CLR/LM) selected metabolites with MMVEC ordination --------------------------------------

# RF identified metabolites that are important for classifying sample types
# Do these metabolites also significantly separate out in the MMVEC ordination?
source("src/biclust_helper_functions.R")

# set color scheme for sample types
sample_type_cols <- c(CCA       = "#6469ed",
                      Coral     = "#e49c4c",
                      Limu      = "#7cc854",
                      CCALimu   = "#10f1f1",
                      CoralLimu = "#e8f700",
                      CoralCCA  = "#e20ef0")

site_cols <- c(SharksCove       = "#440154FF",
               ThreeTablesBeach = "#3B528BFF",
               TwinRocks        = "#21908CFF", 
               UppersBeach      = "#5DC863FF",
               WaimeaBay        = "#FDE725FF")

genus_cols <- c(# Limu
  Jania =  "darkseagreen",
  Halimeda = "green2",
  Galaxaura = "darkgreen",
  Other = "grey58",
  # Coral
  Pocillopora = "orange2",
  Porites = "gold",
  Monitipora = "khaki3")

# read in RF importance scores and linear model results
lm_results <- fread("data/processed/all_metabolites_sample_type_anovas.csv")
lm_results[sample_type_DA == "NA", sample_type_DA := NA]

rf_results <- fread("data/processed/all_metabolite_RF_importance.csv")

# MMVEC ordination eigenvalues for microbes and metabolites
mic_mmvec_ord <- fread("data/raw/mmvec/microbe_ordination.tsv") 
met_mmvec_ord <- fread("data/raw/mmvec/metabolite_ordination.tsv") 

# clean up microbe and feature names
met_mmvec_ord[ ,featureID:= sub("metabolite",
                                "id_",
                                featureID) ]

# remove outlier
met_mmvec_ord <- met_mmvec_ord[Y < 0.05, ]

# check that features selected by RF are in the ordination
sum(rf_results$featureID %in% met_mmvec_ord$featureID)/length(rf_results$featureID)

# merge ordination scores with RF/linear model information
met_ord <- merge(met_mmvec_ord,
                 lm_results,
                 by = "featureID",
                 all.x = T)
met_ord <- merge(met_ord,
                 rf_results,
                 by = "featureID",
                 all.x = T)

# LM selected metabolites
met_ord$LM_selection <- ifelse(met_ord$adj_p_val < 0.001 &
                                    !is.na(met_ord$adj_p_val) &
                                    !is.na(met_ord$sample_type_DA) &
                                    any(met_ord[, c("Coral_MA", "CCA_MA", "Limu_MA")] > 1),
                                yes = "Significant",
                                no = "Not Significant")

met_ord[ , transp:= ifelse(LM_selection == "Not Significant", 0.2 , 0.8)]

# plot the ordination with high RF scored metabolites highlighted
p1 <- ggplot(data = met_ord,
           aes(x = X, y = Y, color = RF_selection))+
  geom_point(alpha = met_ord$transp) +
  theme_minimal()

p1

ggsave("output/randomforest/RF_variables_mmvec_biplot.pdf", 
       plot = p, width = 9, height = 8)

# plot the ordination with LM significant metabolites colored by sample type class
# RF selected metabolites are filled
p2 <- ggplot(data = met_ord,
             aes(x = X, y = Y, shape = RF_selection, col = sample_type_DA))+
  geom_point(alpha = met_ord$transp) +
  scale_shape_manual(values = c(19,1))+
  scale_color_manual(values = sample_type_cols, na.value = "gray")+
  theme_minimal()

p2

ggsave("output/randomforest/LM_variables_mmvec_biplot.pdf", 
       plot = p2, width = 9, height = 8)

# plot the ordination with only LM significant metabolites shown
p3 <- ggplot(data = met_ord[met_ord$LM_selection ==  "Significant",],
  aes(x = X, y = Y, shape = LM_selection, col = sample_type_DA))+
  geom_point() +
  theme_minimal()
p3
ggsave("output/randomforest/LM_variables_mmvec_biplot_only_significant.pdf", 
       plot = p3, width = 9, height = 8)

p4 <-  ggplot(data = met_ord[met_ord$LM_selection ==  "Significant",],
                   aes(x = X, y = Y, shape = LM_selection, col = sample_type_DA))+
  geom_point() +
  theme_minimal()

# Test correlations of X, Y, and Z of ordination with sample type classification
# For X axis
x_met_lm_mod <- lm(X  ~ sample_type_DA, data = met_ord, na.action = na.omit)
x_met_lm_aov <-anova(x_met_lm_mod)
x_met_lm_aov

# For Y axis
y_met_lm_mod <- lm(Y  ~ sample_type_DA, data = met_ord, na.action = na.omit)
y_met_lm_aov <-anova(y_met_lm_mod)
y_met_lm_aov

# For Z axis
z_met_lm_mod <- lm(Z  ~ sample_type_DA, data = met_ord, na.action = na.omit)
z_met_lm_aov <-anova(z_met_lm_mod)
z_met_lm_aov


# Check what classes of compounds were selected as important
met_ord[RF_selection == "Important", .N, by = "Qemistree_class"]

met_ord[RF_selection == "Important" & Qemistree_class == "Glycerophospholipids",featureID]

fwrite(met_ord, "data/processed/all_metabolite_RF_LM_ord.csv")

