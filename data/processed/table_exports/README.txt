This folder contains flat table exports used for data analysis for the Waimea
marine microbe/metabolomics manuscript.

-------------------------------
all_marine_microbe =
	All available marine microbial samples (i.e. 16S). Includes some
samples that we do not have metabolomics data for. Used for microbe specific
analyses (NMDS, permanova). Blanks removed. Tables with blanks included available on request. 
To link sample names (rownames of abundance table) to sample data, use the 'sequencing_id' column.

-------------------------------

all_marine_metabolite = 
	All available marine metabolomics samples. Includes some samples that
we do not have 16S data for. Used for metabolomic specic analyses (NMDS, PCoA,
Networking, etc.). Blanks removed. Note(!) link sample names (rownames of abundance table) using the 'sample_barcode', which is different from the microbe data. There's a table labled 'tax'. It's not taxonomy, per se, but it contains all of the other data we have about the chemical features (e.g. network, classifications, etc.). 

-------------------------------
paired_marine_microbe/paired_mariine_metabolite =
	Microbial samples subset to only include samples that overlap between
the 16S data and the metabolomics data. Used for paired analyses (procrustes,
mantel, mmvec). Blanks removed. Since these data are paired, sample names (rownames of abundance tables) are both set to 'sample_barcode' so everything matches up. 

