# Plot mmvec ordination with specific microbial families highlighted

library(data.table)
library(ggplot2)
library(ggrepel)
library(ggpubr)


# plot colors
# variables =  genus_cols, sample_type_cols, site_cols
source("src/colors.R")
sample_type_cols = sample_type_cols[c("CCA","Coral","Limu")]

##
# --- Input files ---- #
##

# lm: linear model results indicating the enrichment of metabolite features in a sample type
# rf: random forest results indicating importance of metabolite features in predicting sample type
# mmvec: PCoAs produced with mmvec based on metabolite/microbe co-occurence scores
# micro_family:linear model results indicating enrichment of microbial families in a sample type
# micro_tax: taxonomic information for OTUs
lm_file = "data/processed/all_metabolites_sample_type_anovas.csv"
rf_file = "data/processed/all_metabolite_RF_importance.csv"

mmvec_met_file = "data/raw/mmvec/metabolite_ordination.tsv"
met_tax_file = "data/processed/table_exports/all_marine_metabolite_tax_flat_table.csv"

mmvec_micro_file = "data/raw/mmvec/microbe_ordination.tsv"
micro_family_file = "data/processed/microbe_family_anova_and_fold_changes.csv"
micro_tax_file = "data/processed/table_exports/all_marine_microbe_tax_flat_table.csv"

##
# --- Import and organize --- # 
##

lm_results = fread(lm_file)
rf_results = fread(rf_file)

met_mmvec_ord = fread(mmvec_met_file ) 
met_tax = fread(met_tax_file)

mic_mmvec_ord = fread(mmvec_micro_file, header = T)
micro_tax = fread(micro_tax_file, header = T)


# clean
lm_results[sample_type_DA == "NA", sample_type_DA := NA]
met_mmvec_ord[ ,featureID:= sub("metabolite",
                                "id_",
                                featureID) ]

setnames(micro_tax, "V1", "OTU_ID")
setnames(met_tax, "V1", "featureID")

mic_mmvec_ord[ , X := as.numeric(X)]
mic_mmvec_ord[ , Y := as.numeric(Y)]


# remove outlier
met_mmvec_ord = met_mmvec_ord[Y < 0.05, ]

# check the precentage of features selected by RF that are in the mvvec ordination
sum(rf_results$featureID %in% met_mmvec_ord$featureID)/length(rf_results$featureID)

# merge ordination scores with RF/linear model information
met_ord = merge(met_mmvec_ord,
                 lm_results,
                 by = "featureID",
                 all.x = T)
met_ord = merge(met_ord,
                 rf_results,
                 by = "featureID",
                 all.x = T)

met_ord = merge(met_ord,
                met_tax,
                by = "featureID",
                all.x = T)
met_ord[, LibraryID := gsub("Spectral Match to (.+) from NIST14","\\1",LibraryID)]


mic_ord = merge(mic_mmvec_ord,
                micro_tax,
                by = "OTU_ID")



##
# --- Exploratory Plots ----- #
##
# The mmvec ordinations include microbe-metabolite associations that are not
# particularly strong orinteresting. We want to highlight microbe-metabolite interactions that
# are specific to particular sample types. We can indicate the importance of these interactions
# in three ways:
# 1) random forest 'importance' scores identifying metabolite features that can predict sample type
# 2) linear models for each metabolite feature identifying enrichment in a sample type
# 3) linear models for microbial families identifying microbial families that are associated with sample types

# organize data for plotting
met_ord$LM_selection = ifelse(met_ord$adj_p_val < 0.001 &
                                !is.na(met_ord$adj_p_val) &
                                !is.na(met_ord$sample_type_DA) &
                                any(met_ord[, c("Coral_MA", "CCA_MA", "Limu_MA")] > 1),
                              yes = "Significant",
                              no = "Not Significant")

met_ord[ , lm_transp:= ifelse(LM_selection == "Not Significant" | ! (sample_type_DA %in% c("CCA", "Coral", "Limu")), 0.2 , 0.8)]
met_ord[ , rf_transp:= ifelse(RF_selection != "Important", 0.2 , 0.8)]



## Basic plots to orient to dataset
# Plot 1: RF selected metabolites highlight
p1 = ggplot(data = met_ord,
            aes(x = X, y = Y, color = RF_selection))+
            geom_point(alpha = met_ord$rf_transp) +
            theme_minimal()

p1

ggsave("output/randomforest/RF_variables_mmvec_biplot.pdf", 
       plot = p1, width = 9, height = 8)



# Plot 2: LM selected metabolites highlight
p2 = ggplot(data = met_ord,
             aes(x = X, y = Y, shape = RF_selection, col = sample_type_DA))+
            geom_point(alpha = met_ord$lm_transp) +
            scale_shape_manual(values = c(19,1))+
            scale_color_manual(values = sample_type_cols, na.value = "gray")+
            theme_minimal()

p2

ggsave("output/randomforest/LM_variables_mmvec_biplot.pdf", 
       plot = p2, width = 9, height = 8)



# Plot 3: LM selected metabolites only
p3 = ggplot(data = met_ord[met_ord$LM_selection ==  "Significant",],
             aes(x = X, y = Y, shape = LM_selection, col = sample_type_DA))+
  geom_point() +
  theme_minimal()
p3
ggsave("output/randomforest/LM_variables_mmvec_biplot_only_significant.pdf", 
       plot = p3, width = 9, height = 8)


##
# --- Final Plot ----- #
##

# points of interest
# color by our two interesting groups of compounds
# label library matches that were lm significant 

metabolite_highlight = c(
  "Lyso-PAF C-18",
  "1-Stearoyl-2-hydroxy-sn-glycero-3-phosphocholine",
  "1-(9Z-Octadecenoyl)-sn-glycero-3-phosphocholine",
  "17(18)-EpETE",
  "20-Hydroxy-(5Z,8Z,11Z,14Z)-eicosatetraenoic acid",
  "Fucoxanthin")

met_colors = met_ord[Qemistree_subclass %in% c("Glycerophosphocholines") |
                       Qemistree_direct_parent %in% c("Long-chain fatty acids")|
                       LibraryID == "Fucoxanthin"]

met_labels = met_colors[!(LibraryID %in% c("missing", "N/A"))
                        & LM_selection == "Significant"
                        ]

met_labels[, LibraryID := gsub("Spectral Match to (.+) from NIST14","\\1",LibraryID)]

# organize microbes for plotting
microbe_highlight = c(
  "Kiloniellaceae" = "white",
  "Saprospiraceae" = "black"
)




# Color metabolites by sample type association
fp = ggplot(data = met_ord)+

    stat_ellipse(data = met_ord[sample_type_DA %in% c("CCA","Coral","Limu")],
                 aes(x = X,
                     y = Y,
                     fill = sample_type_DA),
                 alpha = 0.6,
                 geom = "polygon")+

    geom_point(aes(x = X,
                   y = Y),
               alpha = 0.2,
               size = 3,
               pch = 20,
               color = "#1A1A1A"
    )+
    geom_segment(data = mic_ord[family %in% names(microbe_highlight)],
               aes(x = 0,
                   y = 0,
                   xend = 0.0002*X,
                   yend = 0.0002*Y,
                   color = family),
               arrow = arrow(length = unit(3, "mm"),
                             type = "closed"),
               alpha = 0.9,
               size = 1)+
    geom_label_repel(data = met_labels[LibraryID %in% metabolite_highlight],
                     aes(x = X, y = Y, label = LibraryID),
                      nudge_y = -0.01)+
     scale_fill_manual(values = sample_type_cols)+
     scale_color_manual(values = microbe_highlight)+
     theme_minimal()

fp

ggsave(plot = fp, "output/mmvec_ordination/final_ord.pdf", width = 15, height = 10)


##
# --- Supplemental Plots ----- #
##

# Make versions with additional families for supplemental figure


alt_families = c(
  # algae
  "Rhodobacteraceae",
  "Flavobacteriaceae",
  "Vibrionaceae",
  # CCA
  "Tenderiaceae",
  # coral
  "Rhizobiales_unclassified",
  "Nitrosopumilaceae",
  "Nitrosococcaceae",
  "Nitrospiraceae",
  "Burkholderiaceae"



)

p_list = list()

for( fam in alt_families){
  
p_list[[fam]] =  ggplot(data = met_ord)+
    
                    stat_ellipse(data = met_ord[sample_type_DA %in% c("CCA","Coral","Limu")],
                                 aes(x = X,
                                     y = Y,
                                     fill = sample_type_DA),
                                 alpha = 0.7,
                                 geom = "polygon")+
                    
                    # geom_point(aes(x = X,
                    #                y = Y),
                    #            alpha = 0.1,
                    #            size = 3,
                    #            pch = 20,
                    #            color = "#1A1A1A"
                    # )+
                    geom_segment(data = mic_ord[family %in% fam],
                                 aes(x = 0,
                                     y = 0,
                                     xend = 0.0002*X,
                                     yend = 0.0002*Y),
                                 arrow = arrow(length = unit(3, "mm"),
                                               type = "closed"),
                                 alpha = 0.9,
                                 size = .3,
                                 color = "black")+
                    # geom_label_repel(data = met_labels[LibraryID %in% metabolite_highlight],
                    #                  aes(x = X, y = Y, label = LibraryID),
                    #                  nudge_y = -0.01)+
                    scale_fill_manual(values = sample_type_cols, name = "Metabolite Sample Type Enrichment")+
                    theme_minimal()+
                    labs(title = fam)

}

ggpubr::ggarrange(plotlist = p_list, common.legend = T,ncol = 3, nrow = 3)
ggsave("output/mmvec_ordination/supplemental_ords.pdf", width = 15, height = 12)
