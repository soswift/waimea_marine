library(data.table)

# Reading in Data ------------------------------------------
chem_peak_file <- "Metabolites_FeatureTable.csv"
chem_meta_file <- "MS_samples_waimea_supertransect - Full Metadata Waimea_Supertransect.csv"

# read metabolite table
chem_raw <- read.csv(chem_peak_file,
                     header = F,
                     stringsAsFactors = F)
# check
chem_raw[1:5,1:5]

# assign sample names as column names
colnames(chem_raw) <- chem_raw[2, ]

# drop sampe type assignments (we'lladd them back in later)
chem_raw <- chem_raw[-c(1:2), ]

# read metabolite metadata
chem_meta <- read.csv(chem_meta_file,
                      stringsAsFactors = F)
# check
chem_meta[1:5,1:5]

# Cleaning Data ------------------------------------------

# function for cleaning metabolite peak data
get_chem_abundance <- function(chem_data = chem_raw){
  # pull out the abundances
  abund_cols <-
    colnames(chem_data)[grepl("Peak area", colnames(chem_data))]
  
  chem_abund <- chem_data[, colnames(chem_data) %in% abund_cols]
  
  # clean up abundance column names
  sample_names <-
    gsub("(....)_.._.+", "\\1", colnames(chem_abund))
  
  sample_names <-
    gsub(".mz.+","",sample_names)
  
  # clean up abundance row names:
  # change row names to cluster index, paste "id_" to prevent indexing confusion
  id_names<- paste0("id_", chem_raw$`cluster index`)

  # convert abundance data to matrix
  chem_abund <- as.matrix(chem_abund)
 
  # make sure all values are numeric (not character strings)
  chem_abund <- t(apply(chem_abund, 1, as.numeric))
  colnames(chem_abund) <- sample_names
  rownames(chem_abund) <- id_names
  return(chem_abund)
}

# pull out/clean peak data from raw format
chem_abund <- get_chem_abundance(chem_data = chem_raw)

# check
chem_abund[1:5,1:5]

# store sample names as a vector
chem_sample_names <- colnames(chem_abund)

# Relative abundance and Summed Relative Abundance -------------------------------

# for each column, divide each value in the column by the column sum
# transpose back to features as rows, samples as columns
chem_RA <- apply(chem_abund, 2, function(column_val) column_val/sum(column_val))
colnames(chem_RA) <- chem_sample_names

# check
chem_RA[1:5,1:5]

# function get_sums() will sum values across samples based on a category in the metadata (e.g. sample_type)

get_sums <- function(abundance, group_col, id_col, meta, new_name){
  abund_dt <- as.data.table(t(abundance), keep.rownames = T)
  abund_dt$group_vec <- meta[ match(abund_dt$rn, meta[[id_col]]), group_col ]
  abund_dt[, rn:=NULL]
  abund_dt <- melt.data.table(abund_dt,
                              id.vars = "group_vec", 
                              variable.name = "ID",
                              value.name = "abund")
  abund_dt[ , sum(abund), by = .(group_vec, ID)]
  abund_dt <- dcast(abund_dt, ID ~ group_vec,
                    value.var = "abund",
                    fun.aggregate = sum)
  setnames(abund_dt, "ID", new_name)
  return(abund_dt)
}

# calculate sample_type sums of relative abundance
chem_sums <- get_sums(abundance = chem_RA,
                       group_col = "sample_type",
                       id_col = "sample_barcode",
                       meta = chem_meta,
                       new_name = "featureID")

# check
chem_sums[1:5,1:5]

# Blanks are not present in the metadata, so the column is called "NA"
# let's change that to "Blank"
setnames(chem_sums, "NA", "Blank")

# check
chem_sums[1:5,1:5]

# Write Out Files As CSV -------------------------------------------------

write.csv(chem_abund,
          "Qtree_clean_metabolite_peak_data.csv",
          row.names = F)
write.csv(chem_RA,
          "Qtree_metabolite_relative_abundance.csv",
          row.names = F)
write.csv(chem_sums,
          "Qtree_metabolite_sample_type_sums.csv",
          row.names = F)

