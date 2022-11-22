# physeq_flat_out() writes all the tables in a phyloseq out as flat .csv files
# provide a phyloseq object, a descriptive string for file neams, file separator, and output directory

physeq_csv_out <-
  function(phyloseq_obj,
           description = NULL,
           outdir = ".") {
    
    mats_out <- list()
    
    mats_out$abundance    <- as(otu_table(phyloseq_obj), "matrix")
    
    if (!is.null(phyloseq_obj@sam_data)) {
      mats_out$sample <- as(sample_data(phyloseq_obj), "matrix")
    }
    
    if (!is.null(phyloseq_obj@tax_table)) {
      mats_out$tax    <- as(tax_table(phyloseq_obj), "matrix")
    }
    
    for(i in names(mats_out)){
      
      write.csv(mats_out[[i]],
                  file.path(outdir, paste(description, i, "flat_table.csv", sep = "_")))
    }
      
  }
