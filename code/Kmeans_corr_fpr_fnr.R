library(foreach)
library(GSFA)
library(Seurat)
#library(optparse)
# option_list <- list(
#   make_option(
#     "--out_folder", action = "store", default = NA, type = 'character',
#     help = "Folder where the outputs are going to be saved [required]"
#   ),
#   make_option(
#     "--pi", action = "store", default = 0.2, type = 'double',
#     help = "The density of factors to simulated, must be between 0 and 1, default is %default [optional]"
#   ),
#   make_option(
#     "--rep", action = "store", default = 1, type = 'integer',
#     help = "The index of the rep of random dataset to simulate, default is %default [optional]"
#   )
# )
# opt <- parse_args(OptionParser(option_list = option_list))

#out_folder <- opt$out_folder # "unguided.normal_mixture.pi_0.1/single_datasets/"
#param_pi <- 0.2 # 0.05, 0.1 or 0.2
#rep <- opt$rep # from 1 to 300
#reps <- c(rep, rep+100, rep+200)
reps <- 1:300

#sim <- readRDS("pi_0.05_simulated_data/simulated_dataset.rep_1.rds")
setwd("~/Documents/Yifan_GSFA/kmeans_fdp_power_pi10/")
seeds <- readRDS("../1K_random_seeds.rds")

param_pi <- 0.1 # 0.05, 0.1 or 0.2

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

for(r in reps) {
  cat("Running simulation data rep ", r, "\n")
  out_file <- paste0("kmeans_corr_fpr_fnr_rep_", r, ".rds")
  
  set.seed(seeds[r])
  
  sim_data <- normal_data_sim(N = params$N, P = params$P, K = params$K, M = params$M,
                              beta_true = params$beta_true,
                              pi_true = params$pi_true,
                              sigma_w2_true = params$sigma_w2_true,
                              psi_true = params$psi_true, G_prob = params$G_prob, offset = T)
  #sim$sim_data <- poisson_count_sim_var_scale(sim_data = sim_data)
  #rownames(sim$sim_data$count) <- paste0("cell", 1:nrow(sim$sim_data$count))
  #colnames(sim$sim_data$count) <- paste0("gene", 1:ncol(sim$sim_data$count))
  
  kres <- kmeans(sim_data$Y, iter.max = 30, nstart=3L, centers = 10)
  clust.label <- kres$cluster
  
  clust.memb <- matrix(0, nrow=max(clust.label), ncol=length(clust.label))
  clust.memb[as.numeric(clust.label) + (1:length(clust.label)-1)*max(clust.label)] <- 1
  # assoc <- as.data.frame(clust.memb %*% sim$sim_data$G)
  # assoc$other <- rowSums(clust.memb) - rowSums(assoc)
  # colnames(assoc) <- c(paste0("guide",1:6),"noguides")
  # rownames(assoc) <- paste0("cluster", 1:max(clust.label))
  
  corr.p <- matrix(0, nrow=nrow(clust.memb), ncol=6)
  for(i in 1:nrow(corr.p)) {
    for (j in 1:ncol(corr.p)) {
      corr.p[i,j] <- cor.test(clust.memb[i,], sim_data$G[,j])$p.value
    }
  }
  corr.p <- matrix(p.adjust(corr.p, method="BH"), ncol=ncol(corr.p), nrow=nrow(corr.p))
  
  colnames(sim_data$Y) <- paste0("gene", 1:ncol(sim_data$Y))
  rownames(sim_data$Y) <- paste0("cell", 1:nrow(sim_data$Y))
  
  expr <- CreateSeuratObject(t(sim_data$Y), min.cells = 0, min.features = 0)
  Idents(expr) <- clust.label
  
  guide.deg <- list()
  proxy <- list()
  for (i in 1:ncol(corr.p)) {
    proxy[[i]] <- which(corr.p[,i]<0.05)
    if (length(proxy[[i]]>0)) {
      guide.deg[[i]] <- FindMarkers(expr, ident.1 = proxy[[i]], test.use = "t")
    }
  }
  theta.true <- (sim_data$F * sim_data$U) %*% t(as.matrix(params$beta_true[1:6,]))
  fdp.kmeans <- numeric(6)
  power.kmeans <- numeric(6)
  sel <- sapply(guide.deg, length)
  markers.num <- list()
  for (i in seq_along(guide.deg)) {
    if (length(guide.deg[[i]]) > 0) {
      pred.markers <- rownames(guide.deg[[i]])[guide.deg[[i]]$p_val_adj<.05]
      gene.num <- as.numeric(substr(pred.markers, 5, nchar(pred.markers)))
      markers.num[[i]] <- gene.num
      fdp.kmeans[i] <- sum(gene.num %in% which(theta.true[,i]==0)) / length(gene.num)
      power.kmeans[i] <- sum(gene.num %in% which(abs(theta.true[,i])>0)) / sum(abs(theta.true[,i])>0)
    }
  }
  saveRDS(list(fdp=fdp.kmeans, power=power.kmeans, Nassocs=sel, markers=markers.num), out_file)
}


### Messed up in the last section. Need to recalculate the fdp and power.
setwd("~/Documents/Yifan_GSFA/kmeans_fdp_power_pi20/")
perf <- list.files(pattern = "rds")
reps <- as.numeric(substr(perf, 25, nchar(perf)-4))
#reps <- as.numeric(substr(perf, 31, nchar(perf)-4))
corr.fdp <- matrix(0, ncol=6, nrow=300)
corr.power <- matrix(0, ncol=6, nrow=300)
sel <- matrix(0, nrow=300, ncol=6)
#seeds <- readRDS("../1K_random_seeds.rds")

# Summarize Kmeans performance
for (i in seq_along(perf)) {
  temp <- readRDS(perf[i])
  corr.fdp[reps[i],] <- temp$fdp
  corr.power[reps[i],] <- temp$power
  if (length(temp$Nassocs>0)) {
    for (j in seq_along(temp$Nassocs)){
      sel[reps[i],j] <- temp$Nassocs[j]
    }
  }
}
saveRDS(list(fdp=corr.fdp, power=corr.power, sel=sel), "../kmeans_fdp_power_pi20.rds")

# Summarize Seurat performance
for (i in seq_along(perf)) {
  temp <- readRDS(perf[i])
  corr.fdp[reps[i],] <- temp$corr_assoc$perf[,1]
  corr.power[reps[i],] <- temp$corr_assoc$perf[,2]
  sel[reps[i], ] <- sapply(temp$corr_assoc$proxy, length)
}
saveRDS(list(fdp=corr.fdp, power=corr.power, sel=sel), "../seurat_fdp_power_pi20.rds")



# #This is for recalculation
# for (i in seq_along(perf)) {
#   temp <- readRDS(perf[i])
#   cat("running", i, "data")
#   set.seed(seeds[reps[i]])
#   sim_data <- normal_data_sim(N = params$N, P = params$P, K = params$K, M = params$M,
#                               beta_true = params$beta_true,
#                               pi_true = params$pi_true,
#                               sigma_w2_true = params$sigma_w2_true,
#                               psi_true = params$psi_true, G_prob = params$G_prob, offset = T)
#   theta.true <- (sim_data$F * sim_data$U) %*% t(as.matrix(params$beta_true[1:6,]))
#   temp.sel <- sapply(temp$markers, length)
#   for (j in which(temp.sel>0)) {
#     sel[reps[i],j] <- temp.sel[j]
#     corr.fdp[reps[i],j] <- 1 - sum(temp$markers[[j]] %in% which(theta.true[,j]!=0)) / length(temp$markers[[j]])
#     corr.power[reps[i],j] <- sum(temp$markers[[j]] %in% which(theta.true[,j]!=0)) / sum(theta.true[,j]!=0)
#   }
# }
saveRDS(list(fdp=corr.fdp, Nassoc=sel, power=corr.power), "kmeans_fdp_power.rds")

sel.kmeans <- (fnr==1) & (fpr==0)
sel = (fpr.gsfa==0) & (fnr.gsfa==1)
load("../GSFA_pi0.05_count_fpr_fnr.Rdata")

perf <- data.frame(method=c(rep("GSFA",sum(!sel)), rep("Kmeans", sum(!sel.kmeans))), FPR=c(fpr.gsfa[!sel], fpr[!sel.kmeans]), power=1-c(fnr.gsfa[!sel], fnr[!sel.kmeans]),
                   guide=c(paste0("guide",c(rep(1:6,colSums(!sel)), rep(1:6,colSums(!sel.kmeans))))))
p1 <- ggplot(perf, aes(guide, FPR, fill=method)) + geom_boxplot() + ylab("False positive rates") + xlab("")
p2 <- ggplot(perf, aes(guide, power, fill=method)) + geom_boxplot() + ylab("Power")
grid.arrange(p1,p2,nrow=2)



### Ploting GSFA vs two-step clustering in two scenarios
setwd("~/Documents/Yifan_GSFA/")
GSFA.stats <- list(count=list(), normal=list())

load("GSFA_fdp_power/GSFA_pi0.05_count_fdp_power.Rdata")
GSFA.stats$count$pi5 <- list(fdp=fdp.gsfa.count, power=power.gsfa.count)
load("GSFA_fdp_power/GSFA_pi0.1_count_fdp_power.Rdata")
GSFA.stats$count$pi10 <- list(fdp=fdp.gsfa.count, power=power.gsfa.count)
load("GSFA_fdp_power/GSFA_pi0.2_count_fdp_power.Rdata")
GSFA.stats$count$pi20 <- list(fdp=fdp.gsfa.count, power=power.gsfa.count)

load("GSFA_fdp_power/GSFA_pi0.05_normal_fdp_power.Rdata")
GSFA.stats$normal$pi5 <- list(fdp=fdp.gsfa.normal, power=power.gsfa.normal)
load("GSFA_fdp_power/GSFA_pi0.1_count_fdp_power.Rdata")
GSFA.stats$normal$pi10 <- list(fdp=fdp.gsfa.normal, power=power.gsfa.normal)
load("GSFA_fdp_power/GSFA_pi0.2_count_fdp_power.Rdata")
GSFA.stats$normal$pi20 <- list(fdp=fdp.gsfa.normal, power=power.gsfa.normal)

seurat.pi5 <- readRDS("two_step_clustering_fdp_power/Seurat_corr_fdp_power_nassoc.RDS")
seurat.pi10 <- readRDS("two_step_clustering_fdp_power/seurat_fdp_power_pi10.rds")
seurat.pi20 <- readRDS("two_step_clustering_fdp_power/seurat_fdp_power_pi20.rds")
kmeans.pi5 <- readRDS("two_step_clustering_fdp_power/kmeans_fdp_power.rds")
kmeans.pi10 <- readRDS("two_step_clustering_fdp_power/kmeans_fdp_power_pi10.rds")
kmeans.pi20 <- readRDS("two_step_clustering_fdp_power/kmeans_fdp_power_pi20.rds")

guides <- factor(c(rep(0.1,300), rep(0.2,300), rep(0.3,300), rep(0.4,300), rep(0.5,300), rep(0.6,300)), levels=1:6*0.1)
dat.normal <- data.frame(Performance=c(c(fdp.gsfa.normal[!is.na(fdp.gsfa.normal)]), c(corr.fdp[sel>0]), c(power.gsfa.normal[!is.na(fdp.gsfa.normal)]), c(corr.power[sel>0])),
                        metric=c(rep("Observed FDP", sum(!is.na(fdp.gsfa.normal))+sum(sel>0)), rep("Power", sum(!is.na(fdp.gsfa.normal))+sum(sel>0))),
                        guides=c(guides[!is.na(fdp.gsfa.normal)], guides[sel>0], guides[!is.na(fdp.gsfa.normal)], guides[sel>0]),
                        method=c(rep("GSFA",sum(!is.na(fdp.gsfa.normal))), rep("clustering",sum(sel>0)),rep("GSFA",sum(!is.na(fdp.gsfa.normal))), rep("clustering",sum(sel>0))),
                        scenario="Normal scenario")
dat.count <- data.frame(Performance=c(c(fdp.gsfa.count[!is.na(fdp.gsfa.count)]), c(seurat$fdp[seurat$Nassoc>0]), c(power.gsfa.count[!is.na(fdp.gsfa.count)]), c(seurat$power[seurat$Nassoc>0])),
                         metric=c(rep("Observed FDP", sum(!is.na(fdp.gsfa.count))+sum(seurat$Nassoc>0)), rep("Power", sum(!is.na(fdp.gsfa.count))+sum(seurat$Nassoc>0))),
                         guides=c(guides[!is.na(fdp.gsfa.count)], guides[seurat$Nassoc>0], guides[!is.na(fdp.gsfa.count)], guides[seurat$Nassoc>0]),
                         method=c(rep("GSFA",sum(!is.na(fdp.gsfa.count))), rep("clustering",sum(seurat$Nassoc>0)),rep("GSFA",sum(!is.na(fdp.gsfa.count))), rep("clustering",sum(seurat$Nassoc>0))),
                         scenario="Count scenario")
ggplot(rbind(dat.count, dat.normal), aes(y=Performance, x=guides, col=method)) + 
  geom_boxplot() + facet_wrap(~metric) +
  facet_grid(scenario~metric) + 
  xlab("True Effect Size of Perturbation")
