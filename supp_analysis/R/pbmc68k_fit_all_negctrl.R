library(Seurat)
library(Matrix)
library(parallel)
library(svd)
library(optparse)
library(GSFA)

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
rep <- opt$rep # from 1 to 300

setwd(opt$out_folder)

pbmc <- readRDS(paste0("../batch_simulation/pbmc68k_600C300G_disjoint_scale1294_randomNorm_sigma0.3to0.5_rep",rep,".rds"))
dat <- as.matrix(pbmc$modified_dat)[colnames(pbmc$GSFA_Y),]
oracle <- pbmc
sim <- pbmc
G <- matrix(0, nrow=ncol(dat), ncol=4)
G[colnames(dat) %in% pbmc$gRNA.cells[[1]],1] <- 1
G[colnames(dat) %in% pbmc$gRNA.cells[[2]],2] <- 1
G[colnames(dat) %in% pbmc$gRNA.cells[[3]],3] <- 1

### GSFA fit to the data
set.seed(rep+100)
nogRNA <- colnames(dat)[!(colnames(dat) %in% unlist(pbmc$gRNA.cells))]
negctrl <- nogRNA[rbinom(length(nogRNA), 1, 0.2/0.7)==1]
G[colnames(dat) %in% negctrl,4] <- 1

dat1 <- sim$GSFA_Y
#rownames(dat) <- colnames(sim$modified_dat)
#colnames(dat) <- rownames(sim$modified_dat)

fit_count <- fit_gsfa_multivar(Y = dat1, G = G, K = 20,
                               prior_type = "mixture_normal", init.method = "svd",
                               niter = 3000, used_niter = 1000,
                               verbose = T, return_samples = F)
saveRDS(fit_count, paste0("GSFA_fit_negctrl_rep",rep,".rds"))

deg <- list()
for (i in 1:3){ 
  deg[[i]] <- rownames(fit_count$lfsr)[fit_count$lfsr[,i]<.05]
}

FNR.gsfa <- c(1-sum(deg[[1]] %in% oracle$sel.genes[[1]])/length(oracle$sel.genes[[1]]),
              1-sum(deg[[2]] %in% oracle$sel.genes[[2]])/length(oracle$sel.genes[[2]]),
              1-sum(deg[[3]] %in% oracle$sel.genes[[3]])/length(oracle$sel.genes[[3]]))
FDP.gsfa <- c(1-sum(deg[[1]] %in% oracle$sel.genes[[1]])/length(deg[[1]]),
              1-sum(deg[[2]] %in% oracle$sel.genes[[2]])/length(deg[[2]]),
              1-sum(deg[[3]] %in% oracle$sel.genes[[3]])/length(deg[[3]]))

### MAST and Wilcoxon in Seurat fit to the data
expr <- CreateSeuratObject(dat, min.features = 0, min.cells = 0)
expr$g1 <- 0
expr$g1[G[,1]==1] <- 1
expr$g1[G[,4]==1] <- 2
expr$g2 <- 0
expr$g2[G[,2]==1] <- 1 
expr$g2[G[,4]==1] <- 2
expr$g3 <- 0
expr$g3[G[,3]==1] <- 1
expr$g3[G[,4]==1] <- 2
expr <- NormalizeData(expr, normalization.method = "LogNormalize", scale.factor = 10000)

pred.markers.MAST <- list()
deg.MAST.bonf <- list()
deg.MAST.bh <- list()
fdp.MAST.bonf <- numeric(3)
fdp.MAST.bh <- numeric(3)
fnr.MAST.bonf <- numeric(3)
fnr.MAST.bh <- numeric(3)
for (i in 1:3) {
  pred.markers.MAST[[i]] <- FindMarkers(expr, ident.1 = 1, ident.2 = 2, group.by = paste0("g",i), min.pct=0, test.use = "MAST", logfc.threshold = 0)
  pred.markers.MAST[[i]]$q_BH <- p.adjust(pred.markers.MAST[[i]]$p_val, method = "BH")
  deg.MAST.bonf[[i]] <- rownames(pred.markers.MAST[[i]])[pred.markers.MAST[[i]]$p_val_adj<.05]
  deg.MAST.bh[[i]] <- rownames(pred.markers.MAST[[i]])[pred.markers.MAST[[i]]$q_BH<.05]
  fdp.MAST.bonf[i] <- 1-sum(deg.MAST.bonf[[i]] %in% pbmc$sel.genes[[i]])/length(deg.MAST.bonf[[i]])
  fnr.MAST.bonf[i] <- 1-sum(deg.MAST.bonf[[i]] %in% pbmc$sel.genes[[i]])/length(pbmc$sel.genes[[i]])
  fdp.MAST.bh[i] <- 1-sum(deg.MAST.bh[[i]] %in% pbmc$sel.genes[[i]])/length(deg.MAST.bh[[i]])
  fnr.MAST.bh[i] <- 1-sum(deg.MAST.bh[[i]] %in% pbmc$sel.genes[[i]])/length(pbmc$sel.genes[[i]])
}

pred.markers.wilcox <- list()
deg.wilcox.bonf <- list()
deg.wilcox.bh <- list()
fdp.wilcox.bonf <- numeric(3)
fdp.wilcox.bh <- numeric(3)
fnr.wilcox.bonf <- numeric(3)
fnr.wilcox.bh <- numeric(3)
for (i in 1:3) {
  pred.markers.wilcox[[i]] <- FindMarkers(expr, ident.1 = 1, ident.2 =2, group.by = paste0("g",i), min.pct=0, test.use = "wilcox", logfc.threshold = 0)
  pred.markers.wilcox[[i]]$q_BH <- p.adjust(pred.markers.wilcox[[i]]$p_val, method = "BH")
  deg.wilcox.bonf[[i]] <- rownames(pred.markers.wilcox[[i]])[pred.markers.wilcox[[i]]$p_val_adj<.05]
  deg.wilcox.bh[[i]] <- rownames(pred.markers.wilcox[[i]])[pred.markers.wilcox[[i]]$q_BH<.05]
  fdp.wilcox.bonf[i] <- 1-sum(deg.wilcox.bonf[[i]] %in% pbmc$sel.genes[[i]])/length(deg.wilcox.bonf[[i]])
  fnr.wilcox.bonf[i] <- 1-sum(deg.wilcox.bonf[[i]] %in% pbmc$sel.genes[[i]])/length(pbmc$sel.genes[[i]])
  fdp.wilcox.bh[i] <- 1-sum(deg.wilcox.bh[[i]] %in% pbmc$sel.genes[[i]])/length(deg.wilcox.bh[[i]])
  fnr.wilcox.bh[i] <- 1-sum(deg.wilcox.bh[[i]] %in% pbmc$sel.genes[[i]])/length(pbmc$sel.genes[[i]])
}



saveRDS(list(dea.MAST=pred.markers.MAST, dea.wilcox=pred.markers.wilcox,
             fdp.wilcox.bh=fdp.wilcox.bh, fnr.wilcox.bh=fnr.wilcox.bh,
             fdp.wilcox.bonf=fdp.wilcox.bonf, fnr.wilcox.bonf=fnr.wilcox.bonf,
             fdp.MAST.bh=fdp.MAST.bh, fnr.MAST.bh=fnr.MAST.bh,
             fdp.MAST.bonf=fdp.MAST.bonf, fnr.MAST.bonf=fnr.MAST.bonf,
             deg.gsfa=deg, fdp.gsfa=FDP.gsfa, fnr.gsfa=FNR.gsfa,
             true.deg=pbmc$sel.genes, gRNA.effect=oracle$effects),
        paste0("All_fit_pbmc68k_600C300G_disjoint_scale1294_randomNorm_sigma0.3to0.5_negctrl_",rep,".rds"))

