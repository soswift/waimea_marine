# functions that help with the correlation bicluster
library(data.table)
library(ComplexHeatmap)
library(dendsort)

# get_means() aggregates a matrix and generates sums by a metadata category
get_means <- function(abundance, group_col, id_col, meta, new_name){
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
                    fun.aggregate = mean)
  setnames(abund_dt, "ID", new_name)
  return(abund_dt)
}

# top_levels() makes updates a vector so that it only shows the top n categories
# useful for limiting colors in a plot
top_levels <- function(Vec, N = 7, exclude = NA){
  top_levels <-names( sort( table( Vec[!(grepl(exclude, Vec))]), decreasing = T)[1:N])
  Vec[!(Vec %in% top_levels) | Vec %in% exclude] <- "other"
  return(Vec)
}

# scale_dat() normalizes mmvec data by centering to mean for each metabolite
scale_dat <- function(mat_dat){
  mat_dat <- t(scale(t(mat_dat)))
  mat_dat[is.na(mat_dat)] <- 0
  return(mat_dat)
}

# get_top_type() returns the name of the highest sum for each row
# assumes the first column is an ID name, subsequent columns are sums by type
get_top_type <- function(type_means){
  type_names <- colnames(type_means[ , -1])
  
  # determine if any of the types is greater than the other two combined
  top_types <-t(apply( type_means[ , -1], 1, function(x) sapply(x, function(y) y > sum(x) - y)))
  top_types <-apply(top_types, 1, function(x) ifelse(any(x),type_names[x],NA))
  names(top_types) <- type_means[[1]]
  return(top_types)
}

# type_match() takes the assigned types of microbes/metabolites and a matrix showing correlations between the two
# if types are the same, returns type name, if different returns "NS"
# results are filtered based on a correlation cutoff point (e.g. r = 0.4)
type_match <- function(mi_type, me_type, cor_mat, cutoff = 0.4){
  
  # full matrices of microbe x metabolite indicating types for each
  mim <- t(matrix(rep(mi_type, times =  length(me_type)),
                  ncol = length(me_type),
                  dimnames = list(names(mi_type), names(me_type))))
  
  mem <- matrix(rep(me_type, length(mi_type)),
                ncol = length(mi_type),
                dimnames = list(names(me_type), names(mi_type)))
  # check equal
  all(colnames(mim) == colnames(mem))
  
  # if types in the matrices don't match, replace with "NS" for 'Not Same'
  mimem <- mem
  mimem[mimem != mim] <- "NS"
  
  # match to correlation matrix
  mimem <- mimem[row.names(cor_mat), colnames(cor_mat)]
  
  # replace value with "NA" if below cutoff
  mimem[cor_mat < cutoff] <- NA
  return(mimem)
}

# quick_brew_pal() generates a palette based on sorted unique values of x.
# x is a character vector, pal is the brewer palette.
# brewer.pal() has a minimum n of 3, so subset if n < 3
# Returns a named vector palette
quick_brew_pal <- function(x, pal){
  uni_x <- sort(unique(x))
  a_pal    <- brewer.pal(length(uni_x), pal)[1:length(uni_x)] 
  names(a_pal) <- uni_x
  a_pal["other"] <- "lightgray"
  return(a_pal)
}


# generate_phymap() makes a complex Heatmap that 'clusters' using phylogenetic trees (or equivalent qemistree)
# takes matrix where metabolites are rows, otus are columns, trees for both as dendrograms.
# If either tree is null, defaults to hclust.
# 'meta' files are used for annotation.
# returns a list where the first item is the heatmap plot, the second is the underlying matrix arranged by the dendrograms
generate_phymap <- function(
  # dendrograms
  micro_tree = fastree_raw,
  chem_tree = qemistree_raw,
  # correlation matrix
  correlation_mat = all_cors_mat,
  # annotation data
  micro_meta = asv_table,
  chem_meta = tip_data,
  # number of discrete annotation colors to display
  n_col = 7,
  # optionally specify microbial class and family
  # provide as list of vectors with names "class" and "family
  # otherwise, set to NA
  taxa_highlight = list(family = fams, class = top_3_classes),
  # barchart sums
  micro_b = micro_means,
  chem_b = chem_means,
  # add boxplots
  box_plot = F,
  # excluded annotations grep
  excluded = "unclass|uncult|-1",
  # main heatmap colors
  type_cols = DA_cols
){
  

  ## Microbe Tree
  # if tree is provided
  # subset microbe tree to ASVs that are present in correlation matrix and vice versa
  if(!is.null(micro_tree)){
    mat_asvs         <- colnames(correlation_mat)
    micro_tree_asvs  <- mat_asvs[ mat_asvs %in% micro_tree$tip.label]
    micro_tree_clean <- keep.tip(micro_tree, micro_tree_asvs)
    correlation_mat  <- correlation_mat[ , micro_tree_asvs]
    
    # microbe tree as dendrogram
    micro_dendro <- as.dendrogram.phylo(micro_tree_clean)
    
    # if tree is null, then use hclust to generate dendrogram
  } else {
    dist_mat <- t(correlation_mat)
    dist_mat[is.na(dist_mat)] <- 0
    micro_dendro <- as.dendrogram( hclust ( dist( dist_mat ) ) )
  }
  
  ## Metabolite Tree
  # if tree is provided
  # subset metabolite tree to features in the correlation matrix and vice versa
  if(!is.null(chem_tree)){
    mat_feats       <- row.names(correlation_mat)
    chem_tree_feats <-
      mat_feats[mat_feats %in% chem_tree$tip.label]
    
    chem_tree_clean <- keep.tip(chem_tree, chem_tree_feats)
    correlation_mat <- correlation_mat[chem_tree_feats , ]
    
    # metabolite tree as dendrogram
    chem_dendro  <- as.dendrogram.phylo(chem_tree_clean)
    
    # if tree is null, use hclust to generate dendrogram
  } else {
    dist_mat <- correlation_mat
    dist_mat[is.na(dist_mat)] <- 0
    chem_dendro <- as.dendrogram( hclust ( dist( dist_mat ) ) )
  }
  
  # IMPORTANT: manually order the correlation matrix to match the dendrograms
  correlation_mat <- correlation_mat[ labels(chem_dendro),
                                      labels(micro_dendro)]
  
  ## Metadata
  # arrange metadata and barplots to match the correlation matrix and dendrograms
  chem_meta <- chem_meta[ match( labels(chem_dendro),
                                 chem_meta$featureID) , ]
  chem_b <- as.matrix(chem_b[match(labels(chem_dendro),
                                   chem_b$featureID) , .(CCA,Coral,Limu)])
  
  micro_meta <- micro_meta[ match(labels(micro_dendro),
                                  micro_meta$OTU_ID) , ]
  micro_b <- as.matrix( micro_b[ match(labels(micro_dendro),
                                       micro_b$OTU_ID) , .(CCA,Coral,Limu) ])
  
  ## Colors
  # set color scheme for heatmap data and annotations
  # type_cols is used for sample type (stacked barplots)
  # other colors (e.g. 'direct_prt_cols') used for discrete row/col annotations
  # note: need to be named vectors, no NA names allowed
  
  # the heatmap color scheme depends on the correlation matrix data type
  if(class(correlation_mat[1]) == "numeric"){
    col_fun <- circlize::colorRamp2( c(0, 1),
                                     c("white","black"))
  }
  if(class(correlation_mat[1]) == "character"){
    col_fun <- type_cols
  }
  
  # metabolite annotation top levels and color palettes
  dir_par_levels <- top_levels(chem_meta$direct_parent,
                               N = n_col,
                               exclude = excluded)
  network_levels <- top_levels(chem_meta$componentindex,
                               N = n_col,
                               exclude = excluded)
  
  dir_par_cols <- quick_brew_pal(dir_par_levels, "Dark2")
  network_cols <- quick_brew_pal(network_levels, "Pastel2")
  
  # microbe annotation top levels and palettes
  if(all(is.na(taxa_highlight))){
  family_levels <- top_levels(micro_meta$family,
                             N = n_col,
                             exclude = excluded)
  class_levels <- top_levels(micro_meta$class,
                             N = n_col,
                             exclude = excluded)
  }else{
    # keep only specified microbial taxanomic levels
    family_levels <- micro_meta$family
    class_levels <-  micro_meta$class
    
    family_levels[!(family_levels %in% taxa_highlight$family)] <- "other"
    class_levels[!(class_levels %in% taxa_highlight$class)] <- "other"
  }
  
  family_cols   <- quick_brew_pal(family_levels, "Set1")
  class_cols    <- quick_brew_pal(class_levels, "Pastel1") 
  
  ## Barplots
  # relativize values for stacked barplots
  chem_b  <- t(apply(chem_b, 1,
                     FUN =  function(x) x/sum(x)))
  micro_b <- t(apply(micro_b, 1,
                     FUN =  function(x) x/sum(x)))
  
  ## Boxplots
  # generate boxplot values
  if(isTRUE(box_plot)){
    # use abund_list() to generate boxplot data that matches heatmap matrix
    box_abunds <- abund_by_group(cor_mat = correlation_mat)
    
    # set relative annotation width and names for plot
    # 
    ann_width = c(1,1,5,10, 10)
    ann_names = c(T,T,F,F,F)
    
    # identify separate wrapper functions of box_plots() for metabolites and microbes
    # these will be called in columnAnnotation/rowAnnotation
    chem_boxplot <- function(index){
      box_plots(index = index, 
                type = "chem",
                abund_list = box_abunds,
                box_direction = "horizontal")
    }
    micro_boxplot <- function(index){
      box_plots(index = rev(index),
                type = "micro",
                abund_list = box_abunds,
                box_direction = "vertical")
    }
  }else{
    # if box_plot == F, don't generate box plots as row/column annotations, instead set to NULL
      micro_boxplot = NULL
      chem_boxplot = NULL
      
      # relative annotation width in plot
      ann_width = c(1,1,4,10)
      ann_names = c(T,T,F,F)
    }
  
  
  ## Metabolite Annotations
  # library IDs
  library_ids <- chem_meta$Qemistree_ms2_library_match
  library_ids[library_ids == "missing"] <- ""
  
  LM_selection <- chem_meta$LM_selection
  LM_selection <- ifelse(LM_selection == "Significant", "*LM","")
  
  sample_type_DA <- chem_meta$sample_type_DA
  sample_type_DA[is.na(sample_type_DA)] <- ""
  
  metabolite_text <- paste(sample_type_DA, LM_selection, library_ids)
  
  # metabolite row annotation
  ha_row <- rowAnnotation(
    # discrete chemical taxonomic groups
    DirectParent = dir_par_levels,
    Network      = network_levels,
    # stacked barplots
    Sample_Type = anno_barplot(
      apply(chem_b, 2, as.numeric),
      gp = gpar(fill = type_cols,
                col = type_cols,
                option = "A"
                ),
      border = F,
      bar_width = 0.7
    ),
    # boxplot
    boxplot = chem_boxplot,
    
    
    # annotations
    metabolite_text = anno_text(metabolite_text, gp = gpar(fontsize = 20)),
    
    # row annotation parameters
    annotation_width = ann_width,
    width = unit(4, "in"),
    show_annotation_name = ann_names,
    annotation_name_gp = gpar(fontsize = 20),
    col = list(DirectParent = dir_par_cols,
               Network = network_cols)
  )
  
  # microbe column annotation
  ha_col <- columnAnnotation(
    # boxplots
    boxplot = micro_boxplot,
 
     # stacked barplots
    Sample_type = anno_barplot(
      micro_b,
      gp = gpar(
        fill = type_cols,
        col = type_cols,
        option = "A"),
      border = F,
      bar_width = 0.7
    ),
  
    # discrete microbial taxonomic groups
    Class = class_levels,
    Family = family_levels,
   
    # column annotation parameters
    annotation_height = rev(ann_width[-length(ann_width)]),
    height = unit(4, "in"),
    show_annotation_name = rev(ann_names[-length(ann_width)]),
    annotation_name_gp = gpar(fontsize = 20),
    col = list(Class = class_cols,
                Family = family_cols)
 
  )
  
  # generate heatmap
  ht <- Heatmap(
    correlation_mat,
    cluster_rows = chem_dendro,
    cluster_columns = micro_dendro,
    
    # row parameters (metabolite)
    row_dend_width = unit(3, "in"),
    row_names_gp = gpar(fontsize = 6),
    show_row_names = F,
    
    # column parameters (microbe)
    column_dend_height = unit(3, "in"),
    column_dend_side = "bottom",
    column_names_gp = gpar(fontsize = 6),
    show_column_names = F,
    
    # colors
    na_col = "white",
    col = col_fun,
    rect_gp = gpar(col = NA),
    
    # annotations
    top_annotation =  ha_col,
    right_annotation = ha_row,
    
    # size
    width = unit(20, "in")
  )
  
  ht_draw <- draw(ht,
                  heatmap_legend_side = "bottom",
                  annotation_legend_side = "bottom")
  
  
  print(paste0("Total tips in Qemistree:","  ",length(labels(chem_dendro))))
  print(paste0("Total tips in 16S phylogeny:", "  ", length(labels(micro_dendro))))
  
  # return heatmap plot, data matrix, microbe info, and metabolite info
  out <- list(heatmap = ht_draw,
              matrix = correlation_mat,
              microbes = micro_meta,
              metabolites = chem_meta)
  
  return(out)
}


# save_heatmap() saves a heatmap as png with specified filename
save_heatmap <- function(heatmap,
                         filename,
                         outdir = "output/heatmap/",
                         csv = F,
                         file_type = c("png","pdf")){
  # save png image
  if("png" %in% file_type){
  png(filename = paste0(outdir,filename,".png"),
      width = 35, height = 25, res = 300, units = "in" )
  print(heatmap[[1]])
  dev.off()
  }
  if("pdf" %in% file_type){
    pdf(file = paste0(outdir,filename,".pdf"),
        width = 35, height = 25)
    print(heatmap[[1]])
    dev.off()
  }
  if(isTRUE(csv)){
    # write data
    write.csv(heatmap[2], paste0(outdir, filename, ".csv"))
  }
}

# sub_cor_mat() subsets corrlation matrix based on metadata groups
sub_cor_mat <- function(cor_mat = z_cat_mat,
                        chem_col = "class",
                        micro_col = "order",
                        chem_val = NULL,
                        micro_val = NULL) {
  
  if(!is.null(chem_val)){
  # subset by specified metabolite group
  cor_mat <- cor_mat[row.names(cor_mat) %in%
                     tip_data$featureID[tip_data[[chem_col]] %in% chem_val], ]
  }
  
  if(!is.null(micro_val)){
  # subset by specified microbe group
  cor_mat <- cor_mat[ , colnames(cor_mat) %in%
                        asv_table$OTU_ID[ asv_table[[micro_col]] %in% micro_val] ]
  }
  # drop empty rows/columns
  cor_mat <- cor_mat[apply(cor_mat, 1, function(x) !(all(is.na(x)))),
                      apply(cor_mat, 2, function(x) !(all(is.na(x))))]
  
  if (dim(cor_mat)[1] == 0) {
    stop("no rows, something is up with your metabolite subset")
  }
  if (dim(cor_mat)[2] == 0) {
    stop("no columns, something is up with your microbe subset")
  }
  return(cor_mat)
}

# abund_by_group() pulls out abundance matrices for metabolites/microbes for each sample group
# takes:
# 1) correlation matrix where microbes are columns and metabolites are rows
# 2) abundance tables for microbes/metabolites where samples are columns
# 3) metadata table with informaiton on sample groups
# id col matches 
abund_by_group <- function(cor_mat = z_cat_mat,
                           m_abund = micro_abund,
                           c_abund = chem_abund,
                           meta = sample_dat,
                           group = "sample_type",
                           sample_col = "sample_barcode"
) {
  # match abundance matrices to the correlation matrix
  m_abund <- m_abund[match(colnames(cor_mat), row.names(m_abund)), ]
  c_abund <- c_abund[match(row.names(cor_mat), row.names(c_abund)), ]
  
  # transform by log10 
  m_abund <- log10(1e-06 + m_abund)
  c_abund <- log10(1e-06 + c_abund)
  # vector of grouping ids
  group_vec <- meta[[group]]
  # for each group, get sample names
  group_samples <- lapply(unique(group_vec), function(x) {
    meta[[sample_col]][group_vec == x]
  })
  names(group_samples) <- unique(group_vec)
  # pull out microbe and metabolite abundances for each group of samples
  micro_list  <- lapply(group_samples, function(grp) {
      micro = m_abund[ , colnames(m_abund) %in% grp]
  }
  )
  chem_list <-  lapply(group_samples, function(grp){  
    chem = c_abund[ , colnames(c_abund) %in% grp]
  }
  )
  abund_list <- list(chem = chem_list,
                     micro = micro_list)
  
  # return list of microbe and metabolite abundance matrices
  return(abund_list)
}

# box_plots() generates metabolite and microbe box plots for annotating heatmap
# gets index passed from complexHeatmap (i.e. which row/column is being plotted)
# specify range using 'rg' so plots are legible
# specify type ("chem" or "micro")
# box_direction can be "vertical" (for columns) or "horizontal (for rows)
box_plots = function(index, type, abund_list = abund_list, box_direction = "horizontal") {
  # determine number of rows/columns
  nr = length(index)

  # subset data to type (chem or micro) and determine range of values
  type_abund <- abund_list[[type]]
  rg <<- range(abund_list)
  
  # define plot area for boxplots
  # note: if plotting for columns, the index order is reversed
  if(box_direction == "horizontal"){
  pushViewport(viewport(xscale = rg,
                        yscale = c(0.5, nr + 0.5)))
    index_vec <- seq_along(index)
    grid.xaxis()
  }
  if(box_direction == "vertical"){
  pushViewport(viewport(xscale = c(0.5, nr + 0.5),
                        yscale = rg))
    index_vec <- seq_along(index)
    grid.yaxis()
  }
  
  # plot the three boxplots, then shift to next row/colunn
  for(i in index_vec) {
    # CCA
    grid.boxplot(type_abund$CCA[index[i], ],
                 pos = nr-i+1 - 0.2,
                 box_width = 0.2, 
                 gp = gpar(fill = "blue"),
                 direction = box_direction,
                 pch = NA,
                 outline = F)
    # Coral
    grid.boxplot(type_abund$Coral[index[i], ],
                 pos = nr-i+1 + 0,
                 box_width = 0.2, 
                 gp = gpar(fill = "orange"),
                 direction = box_direction,
                 pch = NA,
                 outline = F)
    # Limu
    grid.boxplot(type_abund$Limu[index[i], ],
                 pos = nr-i+1 + 0.2,
                 box_width = 0.2, 
                 gp = gpar(fill = "green"),
                 direction = box_direction,
                 pch = NA,
                 outline = F)
  }
  popViewport()
}


# Graph of microbe abundance vs. metabolite abundance should look pretty
# graph_pair() takes a microbe/metabolite pair and plots x/y corellation
# pairs_dat = data table with columns "OTU_ID" and "featureID" 
# chem_dat = abundance matrix for metabolites, columns are samples
# micro_dat = abundance matrix for microbes, columns are samples
# type_cols = named vector of colors, named by categories to be plotted
graph_pair <-
  function(pairs_dat,
           chem_dat = chem_abund,
           micro_dat = micro_abund,
           metadata = sample_dat, 
           get_mean = F,
           type_cols = DA_cols) {
    
    # pull out microbe/metabolite pair
    OTU = pairs_dat$OTU_ID
    Feature = sub("id_","",pairs_dat$featureID)
    
    # get sample metadata
    samples = colnames(micro_dat)
    row.names(metadata) <- metadata$sample_barcode
    sample_types <-metadata[samples, "sample_type"]
    
    # If no mean requested, graph the pair based on raw abundance
    if(get_mean == F){
      
      pair <- data.frame(
        chem = chem_dat[row.names(chem_dat) %in% Feature, samples],
        micro = micro_dat[row.names(micro_dat) %in% OTU, samples],
        sample_type = sample_types)
    } 
    # If mean abundance is requested, 
    # calculate the mean abundance of across all microbe/metabolite pairs
    if(get_mean == T) {
      Feature = unique(Feature)
      # get mean metabolite and microbe abundance if can
      if(sum(row.names(chem_dat) %in% Feature) > 1){
        chem_mean  <- apply(chem_dat[row.names(chem_dat) %in% Feature, samples], 2, mean)
        Feature = "Mean Rel. Abund"
        print("Calculating feature mean")
      }else{ 
        warning("Only 1 feature in the data, not calculating mean")
        chem_mean <- chem_dat[row.names(chem_dat) %in% Feature, samples]
        }
     
       if(sum(row.names(micro_dat) %in% OTU) > 1){
         micro_mean <- apply(micro_dat[row.names(micro_dat) %in% OTU, samples], 2, mean)
         OTU = "Mean Microbe Rel. Abund"
         print("Calculating microbe mean")
       }else{ 
         warning("Only 1 microbe in the data, not calculating mean")
         micro_mean <- micro_dat[row.names(micro_dat) %in% OTU, samples]
         }
      
      pair <-  data.frame(
                          chem = chem_mean,
                          micro = micro_mean,
                          sample_type = sample_types
                        )
      
    }
    
   
    # calculate correlation
    spear_cor <- round( cor(x = pair$chem,
                           y = pair$micro,
                           method = "spearman"),
                       2)

    # plot   
   p <- ggplot( data = pair,
            aes(x = log(chem + 1e-6),
                y = log(micro + 1e-6)))+
      geom_point(aes(col = sample_types), size = 3)+
      labs(
            subtitle = paste0(
              "Microbe: ",  pairs_dat$family, "\n",
              "Metabolite: ", pairs_dat$Qemistree_ms2_library_match,
              "  ", pairs_dat$featureID, "\n",
              "Mean scaled MMVEC: ", round(mean(pairs_dat$z_mmvec), 2), 
              "   RF: ", pairs_dat$RF_selection
              ),
            x = paste(pairs_dat$Qemistree_ms2_library_match, "Log RA"),
            y = paste(pairs_dat$family, "Log RA"))+
     scale_color_manual(values = type_cols)+
     coord_equal()+
     theme_minimal()+
     theme(axis.title =element_text(size=22),
           panel.border = element_rect(colour = "darkgray", size = 1, fill = NA))
   print(p)
  }


# ht_wrapper() subsets, makes a heatmap, and saves it
# we want to plot lots of zoomed in heatmaps and this wrapper function
# makes it easy to plot a bunch of heatmaps for subsets of the metabolite data
ht_wrapper <-function(chem_val = "Organic compounds",
                      chem_col = "kingdom",
                      micro_col = "kingdom",
                      micro_val = "Bacteria",
                      cat_matrix = z_cat_mat,
                      ...){
  print(chem_val)
  sub_mat <- sub_cor_mat(cor_mat = cat_matrix,
                         chem_col = chem_col,
                         micro_col = micro_col,
                         chem_val = chem_val,
                         micro_val = micro_val)
  # generate heatmap
  ht_sub <- generate_phymap( correlation_mat = sub_mat, ...)
  save_heatmap(ht_sub, paste(chem_val,micro_val,"bicluster", sep = "_"))
  return(ht_sub)
}