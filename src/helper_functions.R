# set seed
set.seed( 100 )
# set ggplot themes
  # load libraries
  ## analysis
  library(phyloseq)
  library(vegan)
  ## plotting
  library(gridExtra)
  library(ggplot2)
  library(viridisLite)
  library(pheatmap)
  library(ggpubr)
  
theme_set(theme_minimal())
  ## utility
  library(dplyr)
  library(data.table)
  library(usedist)
  library(dendsort)

# define functions -------------------------------------------------------------------------

# fc() calculates fold change while handling zeroes
fc <- function(x,y){
  if(y==0 && x == 0){
    return(0)
  }else if(y == 0){
    return(Inf)
  }else{
    return(log2(x/y))
  }
}

# function get_chem_abundance() separates feature abundance from feature metadata
# specifically for table provided by Helena
get_chem_abundance <- function(chem_data = chem_raw, bad_samples = bad_samples){

  # pull out the abundances
  abund_cols <-
    colnames(chem_data)[grepl("Peak area", colnames(chem_data))]
  chem_abund <- chem_data[, ..abund_cols]
  
  ## clean up abundance column names
  colnames(chem_abund) <-
    gsub("(....)_.._.+", "\\1", colnames(chem_abund))
  
  colnames(chem_abund) <-
    gsub(".mz.+","",colnames(chem_abund))
  
  good_cols <-
    colnames(chem_abund)[!(colnames(chem_abund) %in% bad_samples)]
  chem_abund <- chem_abund[, ..good_cols]
  chem_abund <- as.matrix(chem_abund)
  
  # change row names to cluster index, paste "id_" to prevent indexing confusion
  row.names(chem_abund)<- paste0("id_", chem_raw$`cluster index`)
  
  # drop features that aren't present in the dataset anymore
  chem_abund <- chem_abund[rowSums(chem_abund) > 0,]
  return(chem_abund)
}

# get_chem_peak_data()

get_chem_peak_data <- function(chem_data = chem_raw, feature_names = row.names(chem_abund)){
  peak_data_cols <- colnames(chem_data)[!grepl("Peak area", colnames(chem_data))]
  chem_peaks <- chem_raw[ , ..peak_data_cols]
  chem_peaks <- as.matrix(chem_peaks)
  row.names(chem_peaks) <- paste0("id_", chem_raw$`cluster index`)
  return(chem_peaks)
}


# function make_chem_phyloseq() makes a phyloseq object from metabolomics data

make_chem_phyloseq <-
  function(chem_abund = chem_abund, chem_meta = chem_meta, chem_tax = chem_tax, id_column = "sample_barcode") {
    # make 'otu table' out of chemical abundance
    abund_phy <- otu_table(chem_abund, taxa_are_rows = T)
    
    sample_names <- chem_meta[[id_column]]
    # make sample data
    sample_phy <- sample_data(chem_meta)
    sample_names(sample_phy) <-sample_names
    # make 'taxonomy' table using data about peaks (e.g. networks)
    taxa_phy <- tax_table(chem_tax)
    
    phy <- phyloseq(abund_phy, sample_phy, taxa_phy)
    
    return(phy)
  }


## from phyloseq source code, import veganifyOTU function
veganifyOTU <- function(physeq){
  if(taxa_are_rows(physeq)){physeq <- t(physeq)}
  return(as(otu_table(physeq), "matrix"))
}


# function do_permanova() pulls info from phyloseq and returns permanova results
do_permanova <- function(phyloseq_obj, var_name, dist_obj, description){
  set.seed(1)
  a_meta <- as(sample_data(phyloseq_obj), "data.frame")
  
  a_formula<- as.formula(paste0("dist_obj ~ ", var_name ))
  aov_results <- as.data.frame(adonis2(a_formula, data = a_meta))

  aov_results$description <- description
  return(aov_results)
}


# function plot_NMDS() will ordinate and plot a phyloseq object using a standard set of parameters
plot_NMDS <- function(a_phy,
                      desc, 
                      color_var, 
                      shape_var=NULL,
                      dist_obj,
                      dist_method = "bray"){
  
  set.seed(1)
  a_ord  <- ordinate(a_phy,
                     method = "NMDS",
                     distance = dist_obj)
  
  p <- plot_ordination(a_phy, a_ord,
                       color = color_var,
                       shape = shape_var) +
                       stat_ellipse( level = 0.5)+
                       theme(aspect.ratio = 1,
                             legend.key.size = unit(0.5, "cm"))
  
  ggsave(paste0("output/NMDS/",
                i,"_",desc,"_NMDS_",dist_method,"_by_", color_var,
                ".pdf"),
         plot = p,
         width = 7,
         height = 5)
  return(p)
}

# function subset_dist() transforms distance object to matrix, subsets, and changes back

subset_dist <- function(dist_obj, phyloseq_obj){
  sample_subset<- sample_names(phyloseq_obj)
  
  # get distances
  a_mat <- as.matrix(dist_obj)[sample_subset, sample_subset]
  a_dist <- as.dist(a_mat)
  return(a_dist)
  
}

# function plot_heatmap()... plots a heatmap! Currently, using pheatmap. Will update to complexHeatmap

plot_heatmap <-
  function(phyloseq_obj,
           description=NULL,
           dist_method = "bray"){
    
    meta <- sample_data(phyloseq_obj)
    otus <- t(as(otu_table(phyloseq_obj), "matrix"))
    
    # update names to include sample types
    new_names <-
      paste(meta$sequencing_id, meta$sample_type, sep = "_")
    row.names(otus) <- new_names
    
    # select annotation columns from meta
    annotation_row = data.frame(
      Type = meta$sample_type,
      Site =  meta$site_name,
      Genus = meta$genus
    )
    
    annotation_row <- droplevels(annotation_row)
    row.names(annotation_row) = new_names
    
    # assign colors
    
    ## define a function for grabbing rcolors
    get_colors <- function(table, column, pal) {
      color_pal <- viridis(nlevels(table[[column]]), option = pal)
      names(color_pal) <- levels(table[[column]])
      return(color_pal)
    }
    
    type_cols <- get_colors(annotation_row, "Type", "D")
    
    site_cols <- get_colors(annotation_row, "Site", "C")
    
    genus_cols <- get_colors(annotation_row, "Genus", "E")
    
    # generate bray distance matrix, cluster, sort
    message(paste0("generating ", dist_method, " distance matrix for ", description, "rows"))
    message("Column clusters based on euclidean distance")
    clust_rows <-
      hclust(vegdist(otus, method = dist_method), method = "average")
    
    # sort dendrogram so outgroups are more apparent
    clust_rows <- as.hclust(dendsort(as.dendrogram(clust_rows)))
    
    
    message("plotting")
    # generate heatmap
    pdf(
      file  = paste0("output/heatmap/", description, "_sample_heatmap.pdf"),
      width = 15,
      height = 10
    )
    
    pheatmap(
      otus,
      #clustering_distance_rows = drows,
      annotation_row = annotation_row,
      show_colnames = F,
      annotation_names_row = T,
      treeheight_row = 100,
      fontsize_row = 4,
      scale = "column",
      cluster_rows = clust_rows,
      clustering_method = "average",
      color = viridis(n = 50, option = "A"),
      annotation_colors = list(Type = type_cols,
                               Site = site_cols,
                               Genus = genus_cols)
    )
    dev.off()
    
  }

# function paired_ordination will plot nmds for microbial and metabolite data and return a list of plots

paired_ordination <- function(microbe_phy,
                              chem_phy,
                              description = NULL,
                              method = "NMDS", 
                              chem_dist_method = "bray", 
                              microbe_dist_method = "Unifrac",
                              color = "genus",
                              shape = "site_name",
                              unifrac_dist = pair_micro_dist)
{
  # Ordinations #
  
  # Microbe 
  micro_dat <- sample_data(microbe_phy)
  
  ## distance
  if(microbe_dist_method == "Unifrac"){
    
    ST.dist<- subset_dist(unifrac_dist, microbe_phy)
    
  } else(
    ST.dist <- vegdist(veganifyOTU(microbe_phy), method = microbe_dist_method)
  )
  
  ## ordinate
  message("Ordinating microbes")
  ST.ord <- ordinate(microbe_phy, method = method, ST.dist)
  
  # plot microbe ordination
  p_micro = plot_ordination(
    microbe_phy,
    ST.ord,
    type = "samples",
    color = color,
    shape = shape,
    title = paste(method, "Microbe", microbe_dist_method, description, "Samples", sep = ", ")
  )
  
  # scale point shapes
  p_micro <- p_micro + scale_shape_discrete(solid = T) 
  
  # if color is continuous, set uniform gradient color
  if(is.numeric(micro_dat[[color]]) | is.integer(micro_dat[[color]])){
    p_micro <- p_micro + uniform_gradient_color()
  }
  
  print(p_micro)
  
  
  # Metabolite
  chem_dat <- sample_data(chem_phy) 
  
  ## distance
  ST_chem.dist <- vegdist(veganifyOTU(chem_phy), method = chem_dist_method)
  
  ## ordinate
  message("Ordinating metabolites")
  
  if(method == "NMDS"){
    # implements distance calculations internally
    ST_chem.ord <- ordinate(chem_phy, method = method, chem_dist_method)
  } else{
    # provide the distance object
    ST_chem.ord <- ordinate(chem_phy, method = method, ST_chem.dist)
  }
  
  # plot metabolite ordination
  p_chem = plot_ordination(
    chem_phy,
    ST_chem.ord,
    type = "samples",
    color = color,
    shape = shape,
    title = paste(method, "Metabolite", chem_dist_method, description, "Samples", sep = ", ")
  )
  
  ## scale point shapes
  p_chem <- p_chem + scale_shape_discrete(solid = T) 
  ## if color is continuous, set uniform gradient color
  
  if(is.numeric(chem_dat[[color]]) | is.integer(chem_dat[[color]])){
    p_chem <- p_chem + uniform_gradient_color()
  }
  
  print(p_chem)
  
  
  # Combined Analysis #
  
  # Mantel
  message("Calculating Mantel")
  
  ## arrange MS distance matrix to match the microbe sample order
  ids <- colnames(ST.dist)
  chem_dist_mat <- ST_chem.dist %>% as.matrix()
  chem_dist_mat <- chem_dist_mat[ids, ids]
  chem_dist <- as.dist(chem_dist_mat)
  
  ## calculate mantel stat
  mant_score <- mantel(ST.dist, ST_chem.dist)
  
  
  
  # Procustes
  message("Procrustes rotation") 
  ## pull out ordination point coordinates and join to metadata
  
  pull_vec <- function(ord, metadata) {
    ord_df        <- ord$points %>% as.data.frame()
    colnames(ord_df) <-c("points_x","points_y")
    ord_df$safe_id   <- row.names(ord$points)
    
    metadata$safe_id <- row.names(metadata)
    ord_df <-left_join(ord_df, metadata, by = "safe_id")
    ord_dt <- as.data.table(ord_df)
    setkey(ord_dt, safe_id)
    setorder(ord_dt, safe_id)
    return(ord_dt)
  }
  
  # data table of microbe ordination
  ord_micro.dt <- pull_vec(ST.ord, micro_dat)
  # data table of chem ordination
  ord_chem.dt <- pull_vec(ST_chem.ord, chem_dat)
  
  ## calculate procustes
  proc <- ade4::procuste(dfX = ord_micro.dt[, list(points_x, points_y)],
                         dfY = ord_chem.dt[, list(points_x, points_y)])
  
  
  # add new procustes coords to data tables
  ord_micro.dt[, c("proc.x","proc.y") := proc$tabX]
  ord_chem.dt[, c("proc.x","proc.y") := proc$tabY]
  
  ord_micro.dt[, data_type := "Microbes"]
  ord_chem.dt[, data_type := "Metabolites"]
  
  cols_to_keep <- c("proc.x",
                    "proc.y",
                    "safe_id",
                    "site_name",
                    "sample_type",
                    "genus",
                    "data_type")
  
  ord_comb.dt <- rbindlist(list(ord_micro.dt[, cols_to_keep, with = F],
                                ord_chem.dt[, cols_to_keep, with = F]))
  ord_comb.dt$plot_color <- ord_comb.dt[[color]]
  
  # generate ggplot
  p_proc <- ggplot(ord_comb.dt, aes(proc.x, proc.y,
                                    color = plot_color,
                                    shape = data_type)) +
    
    geom_point(size = 1) +
    scale_shape_manual(values = c(19,1))+ 
    geom_line(aes(group = safe_id), color = "gray", alpha = 0.3)+
    labs(title = paste(description, ": Procruste Rotation of Ordinations"),
         caption = paste("Mantel Statistic = ", round(mant_score$statistic, 3),
                         "Mantel Significance = ", round(mant_score$signif, 3)),
         x = "",
         y = ""
    )
  if(is.numeric(ord_comb.dt$plot_color) | is.integer(ord_comb.dt$plot_color)){
    p_proc <- p_proc + uniform_gradient_color()
  }
  
  print(p_proc)
  # return plots as list
  
  plot_list <- list(p_micro,p_chem, p_proc)
  names(plot_list) <- c( "micro", "chem", "proc")
  return(plot_list)
  
}

# merge_by_group() merges a phyloseq object and outputs the otu table
merge_by_group <- function(phyloseq_obj, group){
  # merge samples by a group, summing read acounts
  group_phy <- merge_samples(phyloseq_obj, group, fun = sum)
  
  # pull out matrix
  group_otus <- as(otu_table(group_phy), "matrix")
  
  return(group_otus)
}

# euler_subset() plots a eulerr plot for a phyloseq object
# plot circles defined using 'group_by', data subset using expression passed with "..."

euler_subset <- function(physeq, group_by = "genus", ...) {
  type_merge <-
    suppressWarnings(
      merge_by_group(subset_samples(physeq, ...), group = group_by))
  
  type_merge[type_merge > 0] <- 1
  type_merge <- t(type_merge)
  type_eul <- euler(type_merge)
  return(type_eul)
}


