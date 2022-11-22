# A few subclasses of metabolites are highlighted for their biological relevance
# This script generates basic barplots showing relative abundance of each subclass in the sample types
library(data.table)
library(ggpubr)

sample_type_cols <- c(CCA       = "#6469ed",
                      Coral     = "#e49c4c",
                      Limu      = "#7cc854")

# Metabolites -----------------------------------
# Read in data
metabolite_abundance_file <- "data/processed/table_exports/all_marine_metabolite_abundance_flat_table.csv"
sample_data_file          <- "data/processed/table_exports/all_marine_metabolite_sample_flat_table.csv"
metabolite_metadata_file  <- "data/processed/table_exports//all_marine_metabolite_tax_flat_table.csv"
subclass_fold_change_file <- "data/processed/subclass_anova_and_fold_changes.csv"


interesting_subclasses <- c(
  "Amino acids, peptides, and analogues",
  "Fatty acids and conjugates",
  "Eicosanoids",
  "Lineolic acids and derivatives",
  "Glycerophosphocholines",
  "Terpene lactones",
  "Diterpenoids",
  "Estrane steroids",
  "Sesquiterpenoids",
  "Purines and purine derivatives",
  "Monoterpenoids",
  "Glycerophosphoethanolamines",
  "Triterpenoids",
  "Steroid esters"
)

new_interesting_subclasses <- c(
  "Carbonyl compounds",
  "Eicosanoids",
  "Estrane steroids",
  "Fatty acid esters",
  "Fatty acids and conjugates",
  "Fatty alcohols",
  "Glycerophosphocholines",
  "Glycerophosphoethanolamines",
  "Lineolic acids and derivatives",
  "Monoradylglycerols",
  "Monoterpenoids",
  "Purines and purine derivatives",
  "Sesquiterpenoids",
  "Terpene lactones"
)



# read in metabolite data
abund_clean <- fread(metabolite_abundance_file,
                        header = T)

sam_dat   <- fread(sample_data_file,
                      header = T,
                      colClasses = "character",
                      select = c("sample_barcode","sample_type"))

chem_dat  <- fread(metabolite_metadata_file,
                   header = T,
                   select = c("V1", "Qemistree_subclass"))

subclass_fc <- fread(subclass_fold_change_file,
                     header = T)

setnames(chem_dat, "V1","id")
setnames(abund_clean, "V1", "id")


# Set up table for violins
# For each sample, get the mean relative abundance of metabolites in each subclass

abund_clean[1:5,1:5]
chem_dat[1:5,1:2]

mean_abund <- merge(abund_clean, chem_dat)[ , lapply(.SD, mean),
                                            by = "Qemistree_subclass",
                                            .SDcols = !"id"]
mean_abund <- transpose(mean_abund,
                        make.names = "Qemistree_subclass",
                        keep.names = "sample_barcode")

mean_abund <- merge(mean_abund, sam_dat)

# check that subclass names are accurate
new_interesting_subclasses %in% names(mean_abund)

# Make box plots
bx_plts <- lapply(new_interesting_subclasses, function(a_subclass){

  ggviolin(data = mean_abund,
           x = "sample_type",
           y = a_subclass,
           color = "sample_type",
           palette = sample_type_cols,
           add = c("jitter","mean"),
           add.params = list(size =0.5, alpha = 0.5),
           xlab = "Sample type",
           ylab = "Mean relative abundance",
           title = a_subclass)
})

# Convert to grobs
stack_grobs <-lapply(bx_plts, ggplotGrob)

# standardize widths
std_width <- stack_grobs[[1]]$widths
stack_std <- lapply(stack_grobs, function(x) {
  x$widths <- std_width
  return(x)})

# arrange
g <- ggarrange(plotlist=stack_std,
               nrow = 5,
               ncol = 3,
               common.legend = F)
g

# save composite plot
ggsave(filename = "output/violins/subclass_violins.pdf",
       width = 6,
       height = 8,
       scale = 2)

# write out stats summary for the selected subclasses
subclass_fc <- subclass_fc[subclass %in% interesting_subclasses]
fwrite(subclass_fc, "output/violins/subclass_violin_summary_stats.csv")


# Microbes ----------------------------------------
# remove metabolite data
rm(list = ls())

sample_type_cols <- c(CCA       = "#6469ed",
                      Coral     = "#e49c4c",
                      Limu      = "#7cc854")

# identify files 
microbe_abundance_file    <- "data/processed/table_exports/all_marine_microbe_abundance_flat_table.csv"
microbe_abundance_file    <- "data/processed/table_exports/all_marine_microbe_abundance_flat_table.csv"
micro_sample_data_file    <- "data/processed/table_exports/all_marine_microbe_sample_flat_table.csv"
microbe_metadata_file     <- "data/processed/table_exports/all_marine_microbe_tax_flat_table.csv"
class_fold_change_file    <- "data/processed/microbe_class_anova_and_fold_changes.csv"

# read in microbe data
abund_clean <- fread(microbe_abundance_file,
                     header = T)

sam_dat   <- fread(micro_sample_data_file,
                   header = T,
                   colClasses = "character",
                   select = c("sequencing_id","sample_type"))

micro_dat  <- fread(microbe_metadata_file,
                   header = T,
                   select = c("V1", "class"))

class_fc <- fread(class_fold_change_file,
                     header = T)

setnames(micro_dat, "V1","otu")
setnames(abund_clean, "V1", "otu")


# Set up table for violins
# For each sample, get the mean relative abundance of microbes in each subclass

abund_clean[1:5,1:5]
micro_dat[1:5,1:2]

mean_abund <- merge(abund_clean, micro_dat)[ , lapply(.SD, mean),
                                            by = "class",
                                            .SDcols = !"otu"]
mean_abund <- transpose(mean_abund,
                        make.names = "class",
                        keep.names = "sequencing_id")

mean_abund <- merge(mean_abund, sam_dat, by = "sequencing_id")


# Identify interesting microbes based on sample type association and fold change

top_classes <- class_fc[(abs(FC_LimuVCoral) > 1 |
                                   abs(FC_LimuVCCA)   > 1 |
                                   abs(FC_CoralVCCA)  > 1) &
                                  sample_type_DA != "NA"]
top_classes <- top_classes[order(n_samp_obs, decreasing = T)][1:20]

# check that classnames are accurate
interesting_classes <- top_classes$class
interesting_classes %in% names(mean_abund)

# Make box plots
bx_plts <- lapply(interesting_classes, function(a_class){
  
  ggviolin(data = mean_abund,
            x = "sample_type",
            y = a_class,
            color = "sample_type",
            palette = sample_type_cols,
            add = "jitter",
            add.params = list(size =0.5, alpha = 0.5),
            xlab = "Sample type",
            ylab = "Mean relative abundance",
            title = a_class)
  
})

# Convert to grobs
stack_grobs <-lapply(bx_plts, ggplotGrob)

# standardize widths
std_width <- stack_grobs[[1]]$widths
stack_std <- lapply(stack_grobs, function(x) {
  x$widths <- std_width
  return(x)})

# arrange
g <- ggarrange(plotlist=stack_std,
               nrow = 5,
               ncol = 4,
               common.legend = F)
g

# save composite plot
ggsave(filename = "output/violins/microbe_class_violins.pdf",
       width = 8,
       height = 8,
       scale = 2)

# write out stats summary for the selected subclasses
class_fc <- class_fc[class %in% interesting_classes]
fwrite(class_fc, "output/violins/microbe_class_violin_summary_stats.csv")

