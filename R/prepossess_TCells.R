library(data.table)
library(tidyverse)
library(Matrix)
library(Seurat)
library(GSFA)
## Download "GSE119450_RAW.tar" from GEO: GSE119450 and decompress it
## Change the following data directory to where that folder is:
data_dir <- "/project2/xinhe/yifan/Factor_analysis/Stimulated_T_Cells/GSE119450_RAW/"

filename_tb <- 
  data.frame(experiment = c("D1S", "D2S", "D1N", "D2N"),
             prefix = c("GSM3375488_D1S", "GSM3375490_D2S", "GSM3375487_D1N", "GSM3375489_D2N"),
             stringsAsFactors = F)
seurat_lst <- list()
guide_lst <- list()

for (i in 1:4){
  experiment <- filename_tb$experiment[i]
  prefix <- filename_tb$prefix[i]
  print(paste0("Loading data of ", experiment))
  
  feature.names <- data.frame(fread(paste0(data_dir, experiment, "/genes.tsv"),
                                    header = FALSE), stringsAsFactors = FALSE)
  barcode.names <- data.frame(fread(paste0(data_dir, experiment, "/barcodes.tsv"),
                                    header = FALSE), stringsAsFactors = FALSE)
  barcode.names$V2 <- sapply(strsplit(barcode.names$V1, split = "-"),
                             function(x){x[1]})
  # Load the gene count matrix (gene x cell) and annotate the dimension names:
  dm <- readMM(file = paste0(data_dir, experiment, "/matrix.mtx"))
  rownames(dm) <- feature.names$V1
  colnames(dm) <- barcode.names$V2
  print(dim(dm))
  
  # Load the meta data of cells:
  metadata <- data.frame(fread(paste0(data_dir, experiment, "/",
                                      prefix, "_CellBC_sgRNA.csv.gz"),
                               header = T, sep = ','), check.names = F)
  # Example of "metadata" content:
  #            Cell.BC                   gRNA.ID sample.name UMI.count
  # 1 AAAACAAACTGCAGCC ES.sg42.NonTarget.CTRL125     D1_Stim         1
  # 2 AAACCTGAGAGCTGGT              ES.sg21.LAG3     D1_Stim         1
  # 3 AAACCTGAGATACACA           ES.sg5.C10orf54     D1_Stim         1
  metadata$gene_target <- sapply(strsplit(metadata$gRNA.ID, split = "[.]"),
                                 function(x){x[3]})
  metadata$guide <- sapply(strsplit(metadata$gRNA.ID, split = "[.]"),
                           function(x){paste0(x[2], ".", x[3])})
  print(nrow(metadata))
  metadata <- metadata %>% filter(Cell.BC %in% barcode.names$V2)
  targets <- unique(metadata$gene_target)
  targets <- targets[order(targets)]
  # Make a cell by perturbation matrix:
  guide_mat <- data.frame(matrix(nrow = nrow(metadata),
                                 ncol = length(targets)))
  rownames(guide_mat) <- metadata$Cell.BC
  colnames(guide_mat) <- targets
  for (m in targets){
    guide_mat[[m]] <- (metadata$gene_target == m) * 1
  }
  guide_lst[[experiment]] <- guide_mat
  
  # Only keep cells with gRNA info:
  dm.cells_w_gRNA <- dm[, metadata$Cell.BC]
  print("Dimensions of final gene expression matrix:")
  print(dim(dm.cells_w_gRNA))
  
  dm.seurat <- CreateSeuratObject(dm.cells_w_gRNA, project = paste0("TCells_", experiment))
  dm.seurat <- AddMetaData(dm.seurat, metadata = guide_mat)
  seurat_lst[[experiment]] <- dm.seurat
}

combined_obj <- merge(seurat_lst[[1]], 
                      c(seurat_lst[[2]], seurat_lst[[3]], seurat_lst[[4]]),
                      add.cell.ids = filename_tb$experiment,
                      project = "T_cells_all_merged")
print("Dimensions of merged gene expression matrix:")
dim(combined_obj)

# Seurat QC ####
MT_genes <- feature.names %>% filter(startsWith(V2, "MT-")) %>% pull(V1)
combined_obj[['percent_mt']] <- PercentageFeatureSet(combined_obj, features = MT_genes)
VlnPlot(combined_obj, features = c('nFeature_RNA','nCount_RNA','percent_mt'), pt.size = 0.3)
combined_obj <- subset(combined_obj, subset = percent_mt < 10 & nFeature_RNA > 500)

print("Dimensions of merged gene expression matrix after QC:")
dim(combined_obj)
# [1] 33694 24955

# Deviance residual transformation ####
dev_res <- deviance_residual_transform(t(as.matrix(combined_obj@assays$RNA@counts)))
top_gene_index <- select_top_devres_genes(dev_res, num_top_genes = 6000)
dev_res_filtered <- dev_res[, top_gene_index] # 24955 x 6000 matrix

## Regress out library size, unique UMI count, and MT-gene percentage:
covariate_df <- data.frame(lib_size = combined_obj$nCount_RNA,
                           umi_count = combined_obj$nFeature_RNA,
                           percent_mt = combined_obj$percent_mt)
dev_res_corrected <- covariate_removal(dev_res_filtered, covariate_df)

## Cell x gene normalized expression matrix for GSFA input (Y):
scaled.gene_exp <- scale(dev_res_corrected)
sample_names <- colnames(combined_obj@assays$RNA@counts)
gene_names <- rownames(combined_obj@assays$RNA@counts)
rownames(scaled.gene_exp) <- sample_names
colnames(scaled.gene_exp) <- gene_names[top_gene_index]

## Cell x perturbation matrix for GSFA input (G):
G_mat <- combined_obj@meta.data[, 4:24]
G_mat <- as.matrix(G_mat)
# saveRDS(G_mat, "../data/TCells/perturbation_matrix.rds")

## Number of cells with each perturbed gene:
colSums(G_mat)

# Cell group info:
group <- combined_obj$orig.ident
# saveRDS(group, "../data/TCells/sample_groups.rds")
group <- (group %in% c("TCells_D1S", "TCells_D2S")) * 1
## For two group GSFA,
## Group 0: unstimulated cells, Group 1: stimulated cells


# Seurat visualization (just to explore the data, not necessary for GSFA) ####
combined_obj <- NormalizeData(combined_obj, normalization.method = "LogNormalize", scale.factor = 1e4)
gene_detection <- rowMeans(combined_obj@assays$RNA@counts > 0)
kept_features <- names(gene_detection)[gene_detection > 0.1]
combined_obj <- ScaleData(combined_obj,
                          features = kept_features,
                          vars.to.regress = c("nFeature_RNA", "nCount_RNA", "percent_mt"))

## PCA + UMAP:
combined_obj <- RunPCA(combined_obj, features = kept_features, verbose = F)
ElbowPlot(combined_obj, ndims = 20)
combined_obj <- FindNeighbors(combined_obj, reduction = "pca")
combined_obj <- RunUMAP(combined_obj, reduction = "pca", dims = 1:10)

## Donor and stimulation effects:
DimPlot(combined_obj, reduction = "umap", group.by = "orig.ident",
        cols = c("red3", "steelblue4", "forestgreen", "orchid4"))

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
marker_genes <- c("IL7R", "CCR7" , "MKI67", "GZMB")
for (m in marker_genes){
  gene_id <- feature.names$V1[feature.names$V2 == m]
  combined_obj[[m]] <- norm.rna[gene_id, ]
}
FeaturePlot(combined_obj, features = marker_genes,
            min.cutoff = 1,
            cols = c("lightgrey", "pink1", "blue"))

## gRNA distribution
library(ggplot2)
theme_set(theme_bw() + theme(plot.title = element_text(size = 14, hjust = 0.5),
                             axis.title = element_text(size = 14),
                             axis.text = element_text(size = 12),
                             legend.title = element_text(size = 13),
                             legend.text = element_text(size = 12),
                             panel.grid.minor = element_blank())
)
library(gridExtra)
library(grid)

targets <- names(combined_obj@meta.data)[4:24]
plot_data <- DimPlot(combined_obj, reduction = 'umap', group.by = 'NonTarget')
plot.lst <- list()
for (i in names(guide_mat)){
  plot_data$data[[i]] <- factor(combined_obj[[i]][, 1])
  p <- ggplot(plot_data$data, aes_string(x = "UMAP_1", y = "UMAP_2", color = i)) +
    geom_point(alpha = 0.6, size = 0.3) +
    scale_color_manual(values = c("lightgrey", "red")) +
    labs(title = i) +
    theme(legend.position = "None",
          axis.title = element_blank())
  plot.lst[[i]] <- p
}
args <- c(plot.lst,
          list(left = textGrob("UMAP 1", rot = 90, vjust = 0.5,
                               gp = gpar(fontsize = 14, fontface = 'bold')),
               bottom = textGrob("UMAP 2", gp = gpar(fontsize = 14, fontface = 'bold'))))
do.call(grid.arrange, args)
