# Run Random Forest on Metabolites Data to Find Good Predictors of Sample Type
library(randomForest)
# set seed for parallel
set.seed(2020, "L'Ecuyer-CMRG")

# read in metabolite abundance (relative abundance of peak areas) and sample data (information on samples)
# abundance of metabolites will be used to predict sample type (Limu, Coral, CCA)

metabolite_abundance_file <- "data/processed/table_exports/all_marine_metabolite_abundance_flat_table.csv"
sample_data_file <- "data/processed/table_exports/all_marine_metabolite_sample_flat_table.csv"

abund_raw <- read.csv(metabolite_abundance_file,
                  header = T,
                  row.names = 1)
sam_dat <- read.csv(sample_data_file,
                    header = T,
                    row.names = 1)

# clean and arrange abundance data for random forest 
abund <- as.data.frame(t(abund_raw))

row.names(sam_dat) <- paste0("X",row.names(sam_dat))

abund_clean <- abund[ row.names(sam_dat), ]

all(row.names(abund_clean) == row.names(sam_dat))


# run randomForest to get local importance variables (a.k.a. casewise importance)
rf_dat <- abund_clean
rf_dat <- sam_dat$sample_type

print("running RF")
tm <- proc.time()
rf.out <- randomForest(sample_type ~ .,data = rf_dat, ntree = 500, localImp = T, do.trace = 5)
proc.time() - tm
saveRDS(rf.out, "out/sample_type_rf.rds")


