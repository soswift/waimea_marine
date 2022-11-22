# Run a bajillion correlations in parallel
library(parallel)
library(compositions)

metabolite_file <- "data/processed/table_exports/paired_marine_metabolite_abundance_flat_table.csv"
microbe_file    <- "data/processed/table_exports/paired_marine_microbe_abundance_flat_table.csv"
cutoff <- 20

n_cores <- 4

# read in
chem_abund  <- read.csv(metabolite_file, header = T, row.names = 1)
micro_abund <- read.csv(microbe_file, header = T, row.names = 1)

# order by samples
chem_abund  <- chem_abund[ , order(colnames(chem_abund))]
micro_abund <- micro_abund[ , order(colnames(micro_abund))]

# check samples match
all(colnames(micro_abund) == colnames(chem_abund))


# change each microbe into a list
micro_list <- data.frame(t(micro_abund))

# get_hardcor() runs correlations for a microbe on each row of a metabolite matrix
get_hardcor <- function(a_microbe, min_n) {
  apply(X = chem_abund, MARGIN = 1, function(x) {
    # determine total number of data points where either microbe or metabolite is present
    b <- (x >0) | (a_microbe > 0)
    # if total data points is too few, return NA
    if(sum(b) > min_n){
    # run correlations only on samples where the microbe or metabolite occured
    return(cor(x[b], a_microbe[b], method = "spearman"))
    }else{
      return(NA)
    }
  })
}

# run on all microbes
all_cors <- mclapply(X = micro_list,
                      FUN = get_hardcor, 
                      min_n = cutoff,
                      mc.cores = n_cores)
saveRDS(all_cors, "data/processed/all_cors_cutoff_20.rds")

