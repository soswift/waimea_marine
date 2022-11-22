## Generate a two-way heatmap of metabolites and microbes
## Microbes are arrange by fasttree phylogeny, metabolites are arranged by qemistree phylogeny
## Match up all the tip names with heatmap data (correlations)
## Pass phylogenies to complexHeatmap as dendrograms
  
  library(ape)
  library(data.table)
  library(dendextend)
  library(ComplexHeatmap)
  library(phylogram)
  library(circlize)
  library(ggplot2)
  library(RColorBrewer)
  
  # source functions
  source("src/biclust_helper_functions.R")
  
# Identify data files -------------------------------------------------
  # metabolite tree tip key, feature classification, relative abundance data, and RandomForest/LM/PCoA results
  tip_data_file   <- "data/raw/qemistree/Fingerprints to features.tsv"
  peak_data_file  <- "data/processed/table_exports/all_marine_metabolite_tax_flat_table.csv"
  chem_abund_file <- "data/processed/table_exports/paired_marine_metabolite_abundance_flat_table.csv"
  rf_lm_ord_file  <- "data/processed/all_metabolite_RF_LM_ord.csv"
  
  # microbial sample data, taxonomy, relative abundance data
  sample_data_file <- "data/processed/table_exports/paired_marine_microbe_sample_flat_table.csv"
  asv_tax_file     <- "data/processed/table_exports/paired_marine_microbe_tax_flat_table.csv"
  micro_abund_file <- "data/processed/table_exports/paired_marine_microbe_abundance_flat_table.csv"
  
  # microbial and metabolite phylogenies
  micro_tree_file <- "data/raw/FastTree_100.tre"
  chem_tree_file  <- "data/raw/qemistree/qemistree.tree"
  
  # mmvec results 
  mmvec_file      <- "data/raw/mmvec/Ranks_result.tsv"
  
  # set parameters
  # cutoff for median mmvec score to include a given metabolite feature
  mmvec_cutoff = 2
  
  DA_cols <-c(CCA = "#6469ed",
              Coral = "#e49c4c",
              Limu = "#7cc854",
              NS = "#808080")
  
  
 # Read in flat files -----------------------------
  
  # key to matching qemistree tip labels to other metabolomics data(mmvec etc.)
  tip_to_feature <- fread(tip_data_file, key = "id")
  
  # additional feature information (networks, classifications, etc.) metabolites
  peak_data <- fread(peak_data_file)
  setnames(peak_data,"V1", "featureID")
  peak_data[ , featureID:= gsub("id_","",featureID)]

  peak_data[ , Qemistree_ms2_library_match := gsub("Spectral Match to|from NIST14","",Qemistree_ms2_library_match)]
  
  # Random Forest, CLR-LM classification, and ordination values
  rf_lm_ord <-fread(rf_lm_ord_file)
  rf_lm_ord[ ,featureID := sub("id_","",featureID)]
  
  # merge the qemistree tips with metabolite peak information
  # this allows feature information to be mapped onto the tree
  tip_data <- merge(tip_to_feature, peak_data, by = "featureID")
  tip_data <- merge(tip_data, rf_lm_ord[ , .(featureID,
                                             adj_p_val,
                                             sample_type_DA,
                                             X, Y, Z,
                                             MeanDecreaseAccuracy,
                                             RF_selection,
                                             LM_selection)],
                    by = "featureID", all = T)
  
  tip_data$componentindex <- sub(" +", "", tip_data$componentindex)
  qem_id_map <- setNames(tip_data$featureID, tip_data$id)
  
  
  
# Read and Organize Metabolite Qemistree Tree -----------------------------
  # read in qemistree tree file 
  qemistree_raw <- read.tree(file = chem_tree_file)
  
  # update tip labels to 'featureID' to match the rest of our metabolite data
  qemistree_raw$tip.label <- unname(qem_id_map[qemistree_raw$tip.label])
  qem_tip_data <- tip_data[featureID %in% qemistree_raw$tip.label]
  qemistree_raw <- drop.tip(qemistree_raw, tip_data[class == "unclassified", featureID])
  
# Read Microbe FastTree Phylogenetic Tree ------------------------
  # read in fastree file
  fastree_raw <- read.tree(micro_tree_file)
  
# Read in ASV and metabolite Abundance Data-----------------------
  # microbe data
  asv_table <- fread(asv_tax_file, header = T)
  setnames(asv_table, "V1", "OTU_ID")
  
  micro_abund <- read.csv(micro_abund_file, row.names = 1)
  micro_abund <- as.matrix(micro_abund)
  colnames(micro_abund) <- sub("X","",colnames(micro_abund))
  
  sample_dat <- read.csv(sample_data_file)
  
  # metabolite data
  chem_abund <- read.csv(chem_abund_file, row.names = 1)
  chem_abund  <- as.matrix(chem_abund)
  colnames(chem_abund) <- sub("X","",colnames(chem_abund))
  row.names(chem_abund) <- sub("id_","",row.names(chem_abund))
  
  # match sample order and check
  micro_abund <- micro_abund[ , colnames(chem_abund)]
  all(colnames(micro_abund) == colnames(chem_abund))
  
  # get sample type sums for microbes and metabolites (sourced from biclust_helpers.R)
  # use these for barplots on the heatmap and to identify realtionsips betwen sample types
  micro_means <- get_means(micro_abund,
                          group_col = "sample_type",
                          id_col = "sample_barcode",
                          meta = sample_dat,
                          new_name = "OTU_ID")
  
  chem_means <- get_means(chem_abund,
                          group_col = "sample_type",
                          id_col = "sample_barcode",
                          meta = sample_dat,
                         new_name = "featureID")
  
# Read MMVEC and correlation tables -------------------------------
  # MMVEC data (matrix showing pairwise mmvec scores)
  mmvec_table <- fread(mmvec_file,
                       key = "featureid")
  
  mmvec_table[ , featureid := sub("metabolite", "", featureid)]
  setnames(mmvec_table, "featureid", "featureID")
  
# Filter Matrix For Heatmap -------------------------
  ## Subset correlation matrices to match features in qemistree
  # and filter out uninteresting metabolites
  
  mmvec_mat <- as( mmvec_table[ , .SD, .SDcols = !"featureID"], "matrix")
  row.names(mmvec_mat) <- mmvec_table$featureID

  mmvec_mat <- mmvec_mat[row.names(mmvec_mat) %in% qemistree_raw$tip.label ,
                         colnames(mmvec_mat) %in% asv_table$OTU_ID]
  
  # Subset mmvec data to rows that had a high median raw MMVEC score
  mmvec_mat <- mmvec_mat[apply(mmvec_mat, 1, median) >= mmvec_cutoff , ]
  
  # scale the mmvec data (z-score across metabolites)
  # mmvec scores tend to be uniformly high or low for a given metabolite.
  # scaling allows for better interpretation across metabolites on a heatmap
  z_mmvec_mat <- scale_dat(mmvec_mat)
  
# Generate Categorical Matrix -------------------------------------
  # use get_top_type() to figure out the top type for each microbe and metabolite (based on mean relative abundance)
  # type_match() takes this information and generates a categorical matrix
  # each cell in the matrix gets a value based on two criteria:
  #1) is the co-occurance score above a specified cutoff (e.g. 1, 0.4, etc.)
  #2) were the microbe and metabolite most abundant in the same sample type or different sample types?
  
  # calculate sample type with highest mean relative abundance
  micro_top_type <- get_top_type(micro_means)
  chem_top_type  <- get_top_type(chem_means)
  
  # find matches
  # z_cat_mat is a character matrix indicating enrichment in sample types
  # for each microbe and metabolite pair with a co-occurance score above the cutoff
  z_cat_mat <- type_match(mi_type = micro_top_type,
                          me_type = chem_top_type,
                          cor_mat = z_mmvec_mat,
                          cutoff = 1)
  # alternatively
  # use the LM sample type association instead of highest mean rel abund
  # each row will be colored by the results of the LM analysis
  chem_sample_DAs <- tip_data[match(row.names(z_mmvec_mat),featureID), sample_type_DA]
  chem_sample_DAs[chem_sample_DAs == ""] <- "AllTypes"

  alt_cat_mat <- matrix(data = rep(chem_sample_DAs, ncol(z_mmvec_mat)),
                        nrow = nrow(z_mmvec_mat),
                        ncol = ncol(z_mmvec_mat))
  alt_cat_mat[z_mmvec_mat < 1] <- ""
  row.names(alt_cat_mat) <- row.names(z_mmvec_mat)
  colnames(alt_cat_mat) <- colnames(z_mmvec_mat)
  alt_cat_mat[, apply(alt_cat_mat, 2, function(x) ifelse(all(x == ""), F,T))]
  
  
  
  
# Summarize data that made it through MMVEC filters -------------------------
  # write out summary of mmvec filtered metabolites
  # sum mmvec features by subclass and indicate if there is a library match
  mmvec_summary <- tip_data[featureID %in% row.names(z_cat_mat)]
  mmvec_summary <- mmvec_summary[ ,
                                  .(N = .N, 
                                    lib_match = ifelse(
                                      any(ms2_library_match != "missing"),
                                      yes = "Yes",
                                      no = "No")),
                                  by = c("kingdom",
                                         "superclass",
                                         "class",
                                         "subclass")]
  mmvec_summary <- mmvec_summary[order(N, decreasing = T)]
  fwrite(mmvec_summary, "data/processed/mmvec_feature_summary.csv")
  
  
  # Pull out networks that are actually in the correlation matrix
  z_mmvec_nets <- tip_data[tip_data$featureID %in% row.names(z_mmvec_mat),
                           .N, by =  componentindex][order(N, decreasing = T)]
  
  nets_we_have <- z_mmvec_nets[N > 3 & componentindex != "-1", componentindex]
  
# Big Bicluster -----------------------------------
  
  # Identify microbial taxonomic groups that should be displayed on the bicluster
  top_3_classes <-
    unique(top_levels(asv_table[OTU_ID %in% colnames(z_mmvec_mat), class],
                      N = 3, exclude = "unclass|uncult|-1"))
  
  unique(top_levels(asv_table[OTU_ID %in% colnames(z_mmvec_mat), family],
                    N = 20, exclude = "unclass|uncult|-1"))
  fams <- c("Rhodobacteraceae", # Macroalgae - Alphaproteobacteria
            "Flavobacteriaceae", # Macroalgae - Bacterioidia
            "Saprospiraceae", #Macroalgae - Bacterioidia
            "Kiloniellaceae", # Coral/CCA - Alphaproteobacteria
            "Cyclobacteriaceae", # # Macroalgae - Bacterioidia
            "Woeseiaceae", # no association
            "Burkholderiaceae" # Coral - Gammaproteobacteria
            )
  
  # generate full bicluster
  # see generate_phymap() in biclust_helper_functions.R for more information
  # many variables (e.g. trees, sample data) are passed via the default args
   # ht_cat_mmvec <- generate_phymap(correlation_mat = z_cat_mat)
   # save_heatmap(ht_cat_mmvec, "sample_type_heat   map_zmmvec_qemistree")
  # 
  
# Lots of Zoomed in Biclusters -------------------------
  # focus biclusters on specific groups of microbes and metabolites
  
  # subclasses with the most RF 'Important' features in them
  tip_data[RF_selection == "Important", .N, by = subclass][order(N, decreasing = T)]
  
  # subclasses of chemicals showing interesting sample type associations
  cool_subclasses<- c(
    "Amino acids, peptides, and analogues",
    "Fatty acids and conjugates",
    "Eicosanoids",
    "Lineolic acids and derivatives",
    "Glycerophosphocholines",
    "Fatty acid esters"
    )
  
  
  # list of microbes that were differentially abundant between sample types
  diff_microbes <- c(
    "Alphaproteobacteria",
    "Bacteroidia",
    "Gammaproteobacteria"
  )
  
  
  # make biclusters for each subclass of metabolite, only including differentially abundant microbial classes
  # output is a list including
  # 1. the heatmap plot
  # 2. the underlying data matrix
  # 3. microbe information
  # 4. metabolite information

  cool_subclass_hts <- lapply(cool_subclasses, 
                               ht_wrapper,
                               chem_col = "subclass", 
                               micro_val = fams,
                               micro_col = "family",
                               cat_matrix = z_cat_mat)
  

  
  names(cool_subclass_hts) <- cool_subclasses
  
  
  # Long chain fatty acids
  lc_fatty_acid_ht <- ht_wrapper("Long-chain fatty acids",
             chem_col = "direct_parent", 
             micro_val = fams,
             micro_col = "family",
             cat_matrix = z_cat_mat)
  
  
# Identifying Noteworthy Metabolites ------------------------------------------
  ## merge all of our data for microbe metabolite relationships together
  # 1) coerce all data to long format and standardize column names
  # 2) merge all data together and filter pairwise associations by various criteria
  # 3) Validate putative associates with linear correlation plots.
  
  # z-scored mmvec table 
  z_mmvec_table <- data.table(z_mmvec_mat, keep.rownames = T)
  z_mmvec_table_long <- melt.data.table(z_mmvec_table, id.vars = "rn")
  
  setnames(z_mmvec_table_long,
           c("rn","variable","value"),
           c("featureID","OTU_ID","z_mmvec"))
  
  # raw mmvec table
  mmvec_table_long <- melt.data.table(mmvec_table, id.vars = "featureID")
  setnames(mmvec_table_long,
           c("variable","value"),
           c("OTU_ID", "raw_mmvec"))
  
  # put it all together so we can filter microbe/metabolite associations by various criteria
  pair_dat <- merge(z_mmvec_table_long, mmvec_table_long,
                    all.x = T,
                    all.y = F,
                    by = c("featureID","OTU_ID"))
  
  pair_dat <- merge(pair_dat, asv_table,
                    all.x = T,
                    all.y = F,
                    by = "OTU_ID")
  
  pair_dat <- merge(pair_dat, tip_data,
                    all.x = T,
                    all.y = F,
                    by = "featureID")
  
  # graphs of mean relative abundance of a feature across mean relative abundance of microbial families
  
  ## Filtering criteria:
  # plot features and microbes that are in the glycerophosphocholines bilcutser
  glyc_mat <- cool_subclass_hts$Glycerophosphocholines$matrix
  glyc_pairs <- pair_dat[featureID %in% row.names(glyc_mat) &
                          OTU_ID %in% colnames(glyc_mat)]
  
  glyc_pairs <- glyc_pairs[!is.na(z_mmvec)]
  
      pdf("output/Correlations/all_glycero_feature_family_cors.pdf")
      for(feat in unique(glyc_pairs$featureID)){
        for(fam in unique(glyc_pairs$family)){
          feat_fam <- glyc_pairs[family == fam & featureID == feat]
          if(nrow(feat_fam) > 0){
          graph_pair(feat_fam, get_mean = T)
          }
        }
      }
      dev.off()
    
  # plot features and microbes that are in the fatty acids bicluster
  fat_mat <- lc_fatty_acid_ht$matrix
  fat_pairs <- pair_dat[featureID %in% row.names(fat_mat) &
                          OTU_ID %in% colnames(glyc_mat)]  
  fat_pairs <- fat_pairs[!is.na(z_mmvec)]
    
  pdf("output/Correlations/all_fatty_feature_family_cors.pdf")
    for(feat in unique(fat_pairs$featureID)){
      for(fam in unique(fat_pairs$family)){
        feat_fam <- fat_pairs[family == fam & featureID == feat]
        if(nrow(feat_fam) > 0){
          graph_pair(feat_fam, get_mean = T)
        }
      }
    }
    dev.off()
    
    
    
    # Burkholderiaceae, 1−Hexadecyl−sn−glycero−3−phosphocholine (coral, almost same as Lyso-PAF)
    pdf("output/Correlations/Burkholderiaceae_lyso-PAF.pdf")
    graph_pair(glyc_pairs[family=="Burkholderiaceae" & featureID =="25696"], get_mean = T)
    dev.off()
    
    # Kiloniellaceae, 1−Stearoyl−2−hydroxy−sn−glycero−3−phosphocholine
    pdf("output/Correlations/Kiloniellaceae_1−Stearoyl−2−hydroxy−sn−glycero−3−phosphocholine.pdf")
    graph_pair(glyc_pairs[family=="Kiloniellaceae" & featureID =="7266"], get_mean = T)
    dev.off()
  
    # Flavobacteriaceae, 17(18)−EpETE 5263
    pdf("output/Correlations/Flavobacteriaceae_17(18)-EpETE.pdf")
    graph_pair(fat_pairs[family=="Flavobacteriaceae" & featureID =="5263"], get_mean = T)
    dev.off()
    
    # Rhodobacteraceae, 20−Hydroxy−(5Z,8Z,11Z,14Z)−eicosatetraenoic acid
    pdf("output/Correlations/Rhodobacteriaceae-hydroxy-eicosatetranoic_acid.pdf")
    graph_pair(fat_pairs[family=="Rhodobacteraceae" & featureID =="70"], get_mean = T)
    dev.off()
    
    # Flavobacteriaceae, 8-HEPE acid
    pdf("output/Correlations/Flavobacteriaceae 8-HEPE.pdf")
    graph_pair(fat_pairs[family=="Flavobacteriaceae" & featureID =="252"], get_mean = T)
    dev.off()
    
    
  # Determine feature ID name of metabolite from linoleic acid plot
  coral_feature <- tip_data[grepl("Lyso-PAF C-18",LibraryID) , featureID]
  coral_pair <- pair_dat[ family %in% "Rhizobiaceae" & featureID %in% coral_feature]
  
  pdf("output/Correlations/Lyso-paf_rhizobiaceae_coral.pdf")
  graph_pair(coral_pair, get_mean = T)
  dev.off()

  
  
  
