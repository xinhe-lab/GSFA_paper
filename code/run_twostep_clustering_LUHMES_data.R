#!/usr/bin/env Rscript

## Run two-step clustering analysis for LUHMES data

library(optparse)
suppressPackageStartupMessages(library(data.table))
dyn.load('/software/geos-3.7.0-el7-x86_64/lib64/libgeos_c.so') # attach the geos lib for Seurat
suppressPackageStartupMessages(library(Seurat))
require(reshape2)
require(dplyr)
# library(doParallel)
library(foreach)

# Process the command-line arguments --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
parser <- OptionParser()
parser <- add_option(parser,c("--outdir","-o"),type="character",default=NULL)
out    <- parse_args(parser)
outdir <- out$outdir

cat("outdir=", outdir, "\n")

if(!dir.exists(outdir)){
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
}

# 0. Load input data ####
cat("Load input data ... \n")
data_dir <- "/project2/xinhe/kevinluo/GSFA/data/"
combined_obj <- readRDS(file.path(data_dir,"LUHMES_cropseq_data_seurat.rds"))
feature.names <- data.frame(fread(file.path(data_dir, "LUHMES_GSM4219575_Run1_genes.tsv.gz"),
                                  header = FALSE), stringsAsFactors = FALSE)
metadata <- combined_obj@meta.data
table(metadata$orig.ident)

combined_obj

# 1. Data preprocessing ####
cat("Data preprocessing ... \n")

cat("QC ... \n")
# The number of unique genes detected in each cell.
range(combined_obj$nFeature_RNA)

# The total number of molecules detected within a cell
range(combined_obj$nCount_RNA)

# The percentage of reads that map to the mitochondrial genome
range(combined_obj$percent_mt)

cat("Data filtering ... \n")
# filter cells with a library size > 20000 or more than 10% of total read counts from mitochondria genes, as in GSFA.
combined_obj <- subset(combined_obj,
                       subset = percent_mt < 10 & nCount_RNA < 20000)

paste0("Genes: ", dim(combined_obj)[1])
paste0("Cells: ", dim(combined_obj)[2])

# Normalizing the data
combined_obj <- NormalizeData(combined_obj, normalization.method = "LogNormalize", scale.factor = 10000)

# Identification of highly variable features (feature selection)
combined_obj <- FindVariableFeatures(combined_obj, selection.method = "vst", nfeatures = 1000)

# Select the 6k genes used for GSFA
scaled_gene_matrix_in_gsfa <- combined_obj@assays$RNA@scale.data
dim(scaled_gene_matrix_in_gsfa)
selected_gene_ids <- rownames(scaled_gene_matrix_in_gsfa)
cat(length(selected_gene_ids), "genes selected. \n")

all_gene_ids <- rownames(combined_obj)
# Regress out regressed out the differences in experimental batch, unique UMI count, library size, and percentage of mitochondrial gene expression.
# scaled to obtain gene level z-scores.
covariates <- c('orig.ident', 'nCount_RNA', 'nFeature_RNA', 'percent_mt')
combined_obj <- ScaleData(combined_obj, vars.to.regress = covariates, features = all_gene_ids)
saveRDS(combined_obj, file = file.path(outdir, "LUHMES_seurat_processed_data.rds"))

# 2. Perform linear dimensional reduction
combined_obj <- RunPCA(combined_obj, features = VariableFeatures(object = combined_obj))
pdf(file.path(outdir, "LUHMES_seurat_pca.pdf"), width = 5, height = 5)
ElbowPlot(combined_obj, ndims = 50)
dev.off()

# 3. Run non-linear dimensional reduction (UMAP/tSNE) ####
combined_obj <- RunUMAP(combined_obj, dims = 1:30)

# 4. Cluster the cells ####
combined_obj <- FindNeighbors(combined_obj, dims = 1:30)
combined_obj <- FindClusters(combined_obj)
cluster_labels <- Idents(combined_obj)
levels(combined_obj)

# Look at cluster IDs of the first 5 cells
head(cluster_labels, 5)
saveRDS(combined_obj, file = file.path(outdir, "LUHMES_seurat_clustered.rds"))

# 5. Finding differentially expressed features (cluster biomarkers) ####
combined_obj <- readRDS(file.path(outdir, "LUHMES_seurat_clustered.rds"))

cat("Run DE test using MAST...\n")
cat(length(levels(combined_obj)), "clusters.\n")
de.markers <- foreach(i=levels(combined_obj), .packages="Seurat") %do% {
  FindMarkers(combined_obj, ident.1 = i, test.use = "MAST", features = selected_gene_ids)
}
saveRDS(de.markers, file = file.path(outdir, "LUHMES_seurat_MAST_DEGs.rds"))

# session info
sessionInfo()
