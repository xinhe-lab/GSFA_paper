library(foreach)
library(GSFA)
library(Seurat)
library(parallel)

#sim <- readRDS("pi_0.05_simulated_data/simulated_dataset.rep_1.rds")
setwd("~/Documents/Yifan_GSFA/seurat_fdp_power_pi20/")
seeds <- readRDS("../1K_random_seeds.rds")

param_pi <- 0.2 # 0.05, 0.1 or 0.2

params <- list(N = 4000, P = 6000, K = 10, M = 6, G_prob = 0.05,
               sigma_w2_true = rep(0.5, 10), psi_true = 1)
params$beta_true <- matrix(0, nrow = params$M + 1, ncol = params$K)
params$beta_true[params$M + 1, ] <- 0.5 # Intercept
params$beta_true[1, 1] <- 0.1
params$beta_true[2, 2] <- 0.2
params$beta_true[3, 3] <- 0.3
params$beta_true[4, 4] <- 0.4
params$beta_true[5, 5] <- 0.5
params$beta_true[6, 6] <- 0.6

params$pi_true <- rep(param_pi, params$K)

# cores=detectCores()
# cl <- makeCluster(cores[1]-1) #not to overload your computer
# registerDoParallel(cl)

seurat.eval <- function(expr, sim, asso, thres=.05) {
  # Associate clusters to guides with P < thres
  # If thres<=0, pick the smallest p values.
  proxies <- list()
  
  if (thres>0) {
    proxies <- foreach(i=1:6) %do% which(asso[,i]<thres)
  } else {
    proxies <- NULL
  }
  # Differential gene testing for each proxy clusters.
  # pred.markers <- foreach(i=3:6, .packages="Seurat") %dopar% { # Only run FindMarkers for guide3-guide6
  #   if(length(proxies[[i]]) > 0) {
  #     FindMarkers(expr, ident.1 = proxies[[i]]-1, min.pct=0, test.use = "MAST", logfc.threshold = 0)
  #   }
  # }
  pred.markers <- mclapply(proxies, function(x) {if(length(x) > 0) {FindMarkers(expr, ident.1 = x-1, min.pct=0, test.use = "MAST", logfc.threshold = 0)}},
                           mc.cores=getOption("mc.cores", 14L))
  
  fdp.seurat <- numeric(6) # Compute false positive rates
  power.seurat <- numeric(6)
  
  for (i in 1:6) {
    if (length(pred.markers[[i]]) > 0) {
      gene.num <- as.numeric(substr(rownames(pred.markers[[i]]), 5, nchar(rownames(pred.markers[[i]]))))
      fdp.seurat[i] <- 1 - sum(gene.num[pred.markers[[i]]$p_val_adj<0.05] %in% which(sim$sim_data$F[,i]==0)) / length(gene.num[pred.markers[[i]]$p_val_adj<0.05])
      power.seurat[i] <- sum(gene.num[pred.markers[[i]]$p_val_adj<0.05] %in% which(sim$sim_data$F[,i]==1)) / sum(sim$sim_data$F[,i]==1)
    } else {
      fdp.seurat[i] <- NaN
      power.seurat[i] <- NaN
    }
  }
  list(perf=cbind(fdp.seurat, power.seurat), proxy=proxies, markers=pred.markers)
}

for(rep in 1:300) {
  cat("Running simulation data rep ", rep, "\n")
  set.seed(seeds[rep])
  
  sim <- list()
  sim$params <- params
  sim_data <- normal_data_sim(N = params$N, P = params$P, K = params$K, M = params$M,
                              beta_true = params$beta_true,
                              pi_true = params$pi_true,
                              sigma_w2_true = params$sigma_w2_true,
                              psi_true = params$psi_true, G_prob = params$G_prob, offset = T)
  sim$sim_data <- poisson_count_sim_var_scale(sim_data = sim_data)
  rownames(sim$sim_data$count) <- paste0("cell", 1:nrow(sim$sim_data$count))
  colnames(sim$sim_data$count) <- paste0("gene", 1:ncol(sim$sim_data$count))
  
  expr <- CreateSeuratObject(t(sim$sim_data$count), min.features = 0, min.cells = 0)
  expr <- NormalizeData(expr, normalization.method = "LogNormalize", scale.factor = 10000)
  all.genes <- rownames(expr)
  expr <- FindVariableFeatures(expr, selection.method = "vst", nfeatures = 2000)
  expr <- ScaleData(expr, features = all.genes)
  expr <- RunPCA(expr, features = VariableFeatures(object = expr))
  #ElbowPlot(expr) # Different random seeds look the same
  expr <- FindNeighbors(expr, dims = 1:12) # More PCA to capture the factor
  expr <- FindClusters(expr)
  clust.label <- expr@active.ident
  
  clust.memb <- matrix(0, nrow=nlevels(clust.label), ncol=length(clust.label))
  clust.memb[as.numeric(clust.label) + (1:length(clust.label)-1)*nlevels(clust.label)] <- 1
  assoc <- as.data.frame(clust.memb %*% sim$sim_data$G)
  assoc$other <- rowSums(clust.memb) - rowSums(assoc)
  colnames(assoc) <- c(paste0("guide",1:6),"noguides")
  rownames(assoc) <- paste0("cluster", 1:nlevels(clust.label))
  
  # chi.p <- assoc
  # clust.total <- rowSums(assoc)
  # guide.total <- colSums(assoc)
  # cell.total <- sum(assoc)
  # for (i in 1:nrow(assoc)) {
  #   for (j in 1:ncol(assoc)) {
  #     b <- clust.total[i] - assoc[i,j]
  #     c <- guide.total[j] - assoc[i,j]
  #     d <- cell.total - assoc[i,j] - b - c
  #     tab <- matrix(c(assoc[i,j], b, c, d), nrow=2)
  #     chi.p[i,j] <- chisq.test(tab)$p.value
  #   }
  # }
  
  corr.p <- matrix(0, nrow=nrow(clust.memb), ncol=6)
  for(i in 1:nrow(corr.p)) {
    for (j in 1:ncol(corr.p)) {
      corr.p[i,j] <- cor.test(clust.memb[i,], sim$sim_data$G[,j])$p.value
    }
  }
  corr.p <- matrix(p.adjust(corr.p, method="BH"), ncol=ncol(corr.p), nrow=nrow(corr.p))
  perf.cor <- seurat.eval(expr, sim, corr.p, .05)
  
  saveRDS(list(corr_assoc=perf.cor), paste0("cor_seurat_fdp_power_pi20_rep_",rep,".rds"))
}


