library(Seurat)
library(Matrix)
library(parallel)
library(svd)
library(optparse)
library(GSFA)

SCALER <- 1294

option_list <- list(
  make_option(
    "--out_folder", action = "store", default = NA, type = 'character',
    help = "Folder where the outputs are going to be saved [required]"
  ),
  make_option(
    "--rep", action = "store", default = 1, type = 'integer',
    help = "The index of the rep of random dataset to simulate, default is %default [optional]"
  )
)
opt <- parse_args(OptionParser(option_list = option_list))
rep <- opt$rep

setwd(opt$out_folder)

generate.data <- function(seed) {
  pbmc2 <- readRDS("../pbmc68k_Bcell_seuratData.rds")
  
  mat <- as(pbmc2@assays$RNA@data, "dgTMatrix")
  set.seed(seed)
  sel.genes <- list()
  sel.cells <- list()
  gRNA.effect <- list() # Sign of gRNA effects {-1,1}
  #gRNA.effect <- c(0.1,0.3,0.5,0.7,0.9)
  high.expr.genes <- rownames(mat)[rowSums(mat>0)>700]
  gRNA.cells <- rmultinom(ncol(mat),1,c(0.1,0.1,0.1,0.7))
  gRNA.genes <- sample(high.expr.genes,900,replace=F)
  sigma <- c(0.3, 0.4, 0.5)
  #gRNA.effect <- rep(0.4,5)
  #Ncells <- 1:5 * 500
  for (i in 1:3) {
    sel.genes[[i]] <- gRNA.genes[((i-1)*300+1):(i*300)]
    sel.cells[[i]] <- colnames(mat)[gRNA.cells[i,]==1]
    #sel.dir[[i]] <- rbinom(length(sel.genes[[i]]), 1, .5) * 2 - 1
    gRNA.effect[[i]] <- rnorm(length(sel.genes[[i]]), mean=0, sd=sigma)
    mat[sel.genes[[i]],sel.cells[[i]]] <- mat[sel.genes[[i]],sel.cells[[i]]] + gRNA.effect[[i]]
  }
  mat1 <- t(exp(mat)-1) * pbmc2$lib.size / SCALER
  mat2 <- round(mat1)
  mat2[mat2<0] <- 0
  
  dev_mat2 <- deviance_residual_transform(as.matrix(mat2))
  rownames(dev_mat2) <- rownames(mat2)
  colnames(dev_mat2) <- colnames(mat2)
  covariate_df <- data.frame(lib_size = rowSums(mat2),
                             umi_count = rowSums(mat2>0),
                             percent_mt = pbmc2$percent.mt)
  dev_res_corrected <- covariate_removal(dev_mat2, covariate_df)
  scaled.gene_exp <- scale(dev_mat2)
  
  saveRDS(list(original_dat=as(pbmc2@assays$RNA@counts, "dgTMatrix"), modified_dat=t(mat2), gRNA.sd=sigma, #gRNA.sign=sel.dir,
               GSFA_Y=scaled.gene_exp, sel.genes=sel.genes, gRNA.cells=sel.cells, effects=gRNA.effect),
          paste0("pbmc68k_600C300G_disjoint_scale1294_randomNorm_sigma0.3to0.5_rep",seed,".rds"))
}

generate.data(rep)

