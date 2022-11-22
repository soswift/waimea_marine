# read in unifrac distances and match to samples in processed phyloseq object
# read_unifrac() reads in a flat unifrac .csv matrix file.
# it then subsets the matrix to match sample names in a phyloseq object before converting to a dist object

read_unifrac <- function(unifrac_file, phyloseq_obj){
  
  unifrac <- read.csv(unifrac_file, header = T, row.names = 1)
  
  colnames(unifrac) <- sub("X", "", colnames(unifrac))
  
  unifrac <-
    unifrac[sample_names(phyloseq_obj), sample_names(phyloseq_obj)]
  
  unifrac <- as.dist(as.matrix(unifrac))
  
  print(unifrac)
  
}