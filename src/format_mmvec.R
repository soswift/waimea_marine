# convert cmaiki pipeline abundance table to MMVEC format

format_mmvec <- function(phyloseq_obj, id_column = "Group", new_ids, output_dir, description){
  
  # Format abundance table ----------
  abundance_table <- as(otu_table(phyloseq_obj), "matrix")

  # assign new sample ids as column names
  colnames(abundance_table) <- as.character(new_ids)

  # write out as csv file
  write.csv(abundance_table,
         file.path(output_dir,
                   paste0(description,"_abund_mmvec.csv")),
                   row.names = T)
  
  
  print(abundance_table[1:5,1:5])
  
  # Format taxonomy table ----------
  # assumes first two columns are "OTU" and "Size", which are not needed
  
  taxonomy_table <- as(tax_table(phyloseq_obj), "matrix")
  OTU_names <- row.names(taxonomy_table)
  
  taxonomy_table <- as.data.table(taxonomy_table)
  
  # columns to concatenate
  columns_to_smoosh <- colnames(taxonomy_table)
  columns_to_smoosh
  
  # prefixes to add to each column
  prefixes <- c("k__","p__", "c__","o__","f__","g__")
  
  # paste prefixes onto each column
  for (i  in seq_along(prefixes)) {
    taxonomy_table[ , (columns_to_smoosh[i]) := lapply(.SD, function(x, pref)
      {
      paste(pref, x, sep = "")
      },
      pref=prefixes[i]),
      .SDcols = columns_to_smoosh[i]]
  }
  # smoosh columns 
  taxonomy_table[, Taxon:= do.call(paste, c(.SD, sep = ";")), .SDcols = columns_to_smoosh]
  
  # select OTU and new taxonomy column
  taxonomy_table[ , OTU:= OTU_names]
  taxonomy_clean <- taxonomy_table[, .(OTU, Taxon)]
  taxonomy_clean$Taxon <- paste(taxonomy_clean$Taxon, "s__",sep = ";")
  
  # rename OTU to Feature ID
  setnames(taxonomy_clean, "OTU", "Feature ID")
  
  # add prefixes to each taxonomic rank
  
  fwrite(taxonomy_clean,
         file.path(output_dir,
                   paste0(description,"_taxonomy_mmvec.tsv")),
                   sep = "\t",
                   row.names = F)
  
  taxonomy_clean[1:5]
}