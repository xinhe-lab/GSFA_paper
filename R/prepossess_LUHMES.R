library(data.table)
library(tidyverse)
library(Matrix)
library(Seurat)
library(GSFA)
library(ggplot2)
theme_set(theme_bw() + theme(plot.title = element_text(size = 14, hjust = 0.5),
                             axis.title = element_text(size = 14),
                             axis.text = element_text(size = 12),
                             legend.title = element_text(size = 13),
                             legend.text = element_text(size = 12),
                             panel.grid.minor = element_blank())
)
library(gridExtra)
data_dir <- "/project2/xinhe/yifan/Factor_analysis/LUHMES/GSE142078_raw/"

filename_tb <- data.frame(experiment = c("Run1", "Run2", "Run3"),
                          prefix = c("GSM4219575_Run1", "GSM4219576_Run2", "GSM4219577_Run3"),
                          stringsAsFactors = F)
seurat_lst <- list()
guide_lst <- list()
for (run in 1:3){
  experiment <- filename_tb$experiment[run]
  prefix <- filename_tb$prefix[run]
  print(prefix)
  
  feature.names <- data.frame(fread(paste0(data_dir, prefix, "_genes.tsv.gz"),
                                    header = FALSE), stringsAsFactors = FALSE)
  barcode.names <- data.frame(fread(paste0(data_dir, prefix, "_barcodes.tsv.gz"),
                                    header = FALSE), stringsAsFactors = FALSE)
  barcode.names$V2 <- sapply(strsplit(barcode.names$V1, split = "-"),
                             function(x){x[1]})
  # Load the gene count matrix (gene x cell):
  dm <- readMM(file = paste0(data_dir, prefix, "_matrix.mtx.gz"))
  rownames(dm) <- feature.names$V1
  colnames(dm) <- barcode.names$V2
  print(dim(dm))
  # Load the meta data of cells:
  metadata <- data.frame(fread(paste0(data_dir, prefix, "_Cell_Guide_Lookup.csv.gz"),
                               header = T, sep = ','), check.names = F)
  metadata$target <- sapply(strsplit(metadata$sgRNA, split = "_"),
                            function(x){x[1]})
  metadata <- metadata %>% filter(CellBarcode %in% barcode.names$V2)
  targets <- unique(metadata$target)
  targets <- targets[order(targets)]
  # Make a cell by perturbation matrix:
  guide_mat <- data.frame(matrix(nrow = nrow(metadata),
                                 ncol = length(targets)))
  rownames(guide_mat) <- metadata$CellBarcode
  colnames(guide_mat) <- targets
  for (i in targets){
    guide_mat[[i]] <- (metadata$target == i) * 1
  }
  guide_lst[[run]] <- guide_mat
  
  dm.cells_w_gRNA <- dm[, metadata$CellBarcode]
  real_gene_bool <- startsWith(rownames(dm.cells_w_gRNA), "ENSG")
  dm.cells_w_gRNA <- dm.cells_w_gRNA[real_gene_bool, ]
  print("Dimensions of final gene expression matrix:")
  print(dim(dm.cells_w_gRNA))
  
  dm.seurat <- CreateSeuratObject(dm.cells_w_gRNA, project = paste0("LUHMES_", experiment))
  dm.seurat <- AddMetaData(dm.seurat, metadata = guide_mat)
  seurat_lst[[run]] <- dm.seurat
}

combined_obj <- merge(seurat_lst[[1]], c(seurat_lst[[2]], seurat_lst[[3]]),
                      add.cell.ids = c("Run1", "Run2", "Run3"),
                      project = "LUHMES")
print("Dimensions of merged gene expression matrix:")
print(dim(combined_obj))
# Number of cells with each perturbed gene:
colSums(combined_obj@meta.data[, -c(1:3)])

# Seurat QC ####
MT_genes <- feature.names %>% filter(startsWith(V2, "MT-")) %>% pull(V1)
combined_obj[['percent_mt']] <- PercentageFeatureSet(combined_obj, features = MT_genes)
VlnPlot(combined_obj, features = c('nFeature_RNA','nCount_RNA','percent_mt'), pt.size = 0.3)
combined_obj <- subset(combined_obj, subset = percent_mt < 10 & nCount_RNA < 20000)
VlnPlot(combined_obj, features = c('nFeature_RNA','nCount_RNA','percent_mt'), pt.size = 0.3)

print("Dimensions of merged gene expression matrix after QC:")
print(dim(combined_obj))

# Deviance residual transformation ####
dev_res <- deviance_res_transform(t(as.matrix(combined_obj@assays$RNA@counts)))
top_gene_index <- select_top_devres_genes(dev_res, num_top_genes = 6000)
dev_res_filtered <- dev_res[, top_gene_index]

covariate_df <- data.frame(batch = factor(combined_obj$orig.ident),
                           lib_size = combined_obj$nCount_RNA,
                           umi_count = combined_obj$nFeature_RNA,
                           percent_mt = combined_obj$percent_mt)
dev_res_corrected <- covariate_removal(dev_res_filtered, covariate_df)
# Cell x gene normalized expression matrix for GSFA input (Y):
scaled.gene_exp <- scale(dev_res_corrected)
# Cell x perturbation matrix for GSFA input (G):
G_mat <- combined_obj@meta.data[, 4:18]
G_mat <- as.matrix(G_mat)

# Seurat visualization (just to explore the data, not necessary for GSFA) ####
gene_detection <- rowMeans(combined_obj@assays$RNA@counts > 0)
combined_obj <- NormalizeData(combined_obj, normalization.method = "LogNormalize", scale.factor = 1e4)
combined_obj <- ScaleData(combined_obj,
                          features = names(gene_detection)[gene_detection > 0.1],
                          vars.to.regress = c("orig.ident", "nFeature_RNA", "nCount_RNA", "percent_mt"))
# # Save data:
# saveRDS(combined_obj, "processed_data/seurat_obj.merged_scaled_detect_01.corrected_new.rds")
# saveRDS(combined_obj@meta.data, "processed_data/merged_metadata_new.rds")

## PCA + UMAP
combined_obj <- RunPCA(combined_obj, features = kept_features, verbose = F)
ElbowPlot(combined_obj, ndims = 20)
combined_obj <- FindNeighbors(combined_obj, reduction = "pca")
combined_obj <- RunUMAP(combined_obj, reduction = "pca", dims = 1:8)

## Batch effect
DimPlot(combined_obj, reduction = "umap", group.by = "orig.ident")

## Cell cycle effects
CC_genes_df <- data.frame(fread("/project2/xinhe/shared_data/cell_cycle_genes.csv"))
norm.rna <- combined_obj@assays$RNA@data
for (i in 1:ncol(CC_genes_df)){
  gene_ids <- feature.names %>% filter(V2 %in% na.omit(CC_genes_df[, i])) %>% pull(V1)
  phase <- names(CC_genes_df)[i]
  phase.exp <- Matrix::colMeans(norm.rna[gene_ids, ])
  combined_obj[[phase]] <- phase.exp
}
FeaturePlot(combined_obj, features = names(CC_genes_df))

## Marker gene expression
marker_genes <- c("TOP2A", "NEUROD1" , "STMN2", "MAP2")
for (m in marker_genes){
  gene_id <- feature.names$V1[feature.names$V2 == m]
  combined_obj[[m]] <- norm.rna[gene_id, ]
}
FeaturePlot(combined_obj, features = marker_genes,
            min.cutoff = 1,
            cols = c("lightgrey", "pink1", "blue"))

## gRNA distribution
library(grid)
targets <- names(combined_obj@meta.data)[4:18]
plot_data <- DimPlot(combined_obj, reduction = 'umap', group.by = 'Nontargeting')
plot.lst <- list()
for (i in targets){
  plot_data$data[[i]] <- factor(combined_obj[[i]][, 1])
  p <- ggplot(plot_data$data, aes_string(x = "UMAP_1", y = "UMAP_2", color = i)) +
    geom_point(alpha = 0.6, size = 0.3) +
    scale_color_manual(values = c("lightgrey", "red")) +
    labs(title = i) +
    theme(legend.position = "None",
          axis.title = element_blank(),
          panel.grid = element_blank())
  plot.lst[[i]] <- p
}
args <- c(plot.lst, list(nrow = 3,
                         left = textGrob("UMAP 1", rot = 90, vjust = 0.5,
                                         gp = gpar(fontsize = 14, fontface = 'bold')),
                         bottom = textGrob("UMAP 2", gp = gpar(fontsize = 14, fontface = 'bold'))))
do.call(grid.arrange, args)
