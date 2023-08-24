library(ggplot2)
library(gridExtra)
library(ggplot2)
library(GSFA)
theme_set(
  theme_bw() +
    theme(plot.title = element_text(size = 14, hjust = 0.5),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 12),
          legend.title = element_text(size = 13),
          legend.text = element_text(size = 12)
    )
)
setwd("~/Documents/Yifan_GSFA/corr_chisq_fpr_fnr/")
#setwd("~/Documents/Yifan_GSFA/kmeans/")
perfs <- list.files()
perfs <- perfs[endsWith(perfs, "rds")]
perfs <- perfs[grepl("FDR", perfs)]
reps <- as.numeric(substr(perfs, 21, nchar(perfs)-4))
corr.fdp <- matrix(0, nrow=300, ncol=6)
corr.power <- matrix(0, nrow=300, ncol=6)
Nassoc <- matrix(0, nrow=300, ncol=6)
sel <- matrix(0, nrow=300, ncol=6)
seeds <- readRDS("../1K_random_seeds.rds")
# Need to recalculate FPR and FNR
# reps <- substr(perfs, 1, nchar(perfs)-4)
# reps <- sapply(strsplit(reps,"_"), function(x) {as.numeric(x[4])})
# 
# for (i in seq_along(perfs)) {
#   temp <- readRDS(perfs[i])
#   Nassoc[reps[i],] <- colSums(temp$corr_fdr<.05)
# }

param_pi <- 0.05
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

for (i in seq_along(perfs)) {
  temp <- readRDS(perfs[i])
  cat("running", i, "data")
  set.seed(seeds[reps[i]])
  sim_data <- normal_data_sim(N = params$N, P = params$P, K = params$K, M = params$M,
                              beta_true = params$beta_true,
                              pi_true = params$pi_true,
                              sigma_w2_true = params$sigma_w2_true,
                              psi_true = params$psi_true, G_prob = params$G_prob, offset = T)
  theta.true <- (sim_data$F * sim_data$U) %*% t(as.matrix(params$beta_true[1:6,]))
  sel[reps[i],] <- sapply(temp$corr_assoc$proxy, length)
  for (j in which(sel[reps[i],]>0)) {
    genes <- rownames(temp$corr_assoc$markers[[j]])[temp$corr_assoc$markers[[j]]$p_val_adj<.05]
    genes <- as.numeric(substr(genes, 5, nchar(genes)))
    corr.fdp[reps[i],j] <- 1 - sum(genes %in% which(theta.true[,j]!=0)) / length(genes)
  }
  
  corr.power[reps[i],] <- 1 - temp$corr_assoc$perf[,"fnr.seurat"]
  Nassoc[reps[i],] <- sapply(temp$corr_assoc$proxy, length)
}
saveRDS(list(FNR=corr.fnr, FPR=corr.fnr, Nassoc=Nassoc), "../Seurat_corr_fnr_fpr_nassoc.RDS")
saveRDS(list(power=corr.power, fdp=corr.fdp, Nassoc=Nassoc), "../Seurat_corr_fdp_power_nassoc.RDS")
temp <- readRDS("../Seurat_corr_fnr_fpr_nassoc.RDS")





load("../GSFA_pi0.05_count_fpr_fnr.Rdata")
gsfa.assocperguide <- 300-colSums((fpr.gsfa==0) & (fnr.gsfa==1))
seurat.assocperguide <- colSums(Nassoc>0)
#kmeans.assocperguide <- colSums(Nassoc>0)
#Nassocperguide <- data.frame(method=c(rep("GSFA",6), rep("Seurat",6)), guide=rep(paste0("guide",1:6),2), Nassoc=c(gsfa.assocperguide, seurat.assocperguide))
Nassocperguide <- data.frame(method=c(rep("GSFA",6), rep("Kmeans",6)), guide=rep(paste0("guide",1:6),2), Nassoc=c(gsfa.assocperguide, kmeans.assocperguide))
Nassocperguide$percent <- Nassocperguide$Nassoc / 3
ggplot(Nassocperguide, aes(x=guide, y=percent, fill=method)) + geom_bar(stat="identity", position="dodge") + 
  ylab("Percentage of guide associations detected") + 
  geom_text(aes(y=percent+1.5, label=sprintf("%1.1f%%", percent)),position=position_dodge(width=1))

sel = (fpr.gsfa==0) & (fnr.gsfa==1)
fpr = corr.fpr
fnr = corr.fnr
perf <- data.frame(method=c(rep("GSFA",sum(!sel)), rep("Seurat", sum(Nassoc>0))), FPR=c(fpr.gsfa[!sel], fpr[Nassoc>0]), power=1-c(fnr.gsfa[!sel], fnr[Nassoc>0]),
                   guide=c(paste0("guide",c(rep(1:6,colSums(!sel)), rep(1:6,colSums(Nassoc>0))))))
p1 <- ggplot(perf, aes(guide, FPR, fill=method)) + geom_boxplot() + ylab("False positive rates") + xlab("")
p2 <- ggplot(perf, aes(guide, power, fill=method)) + geom_boxplot() + ylab("Power")
grid.arrange(p1,p2,nrow=2)

boxplot(c(fpr.gsfa[!sel]), c(corr.fpr[Nassoc>0]), names=c("GSFA","Seurat_correlation"), ylab="False positive rate")
boxplot(c(fnr.gsfa[!sel]), c(corr.fnr[Nassoc>0]), names=c("GSFA","Seurat_correlation"), ylab="False negative rate")
load("../")

temp <- readRDS(perfs[10])
pred.markers <- temp$corr_assoc$markers
gene.num <- as.numeric(substr(rownames(pred.markers[[2]]), 5, nchar(rownames(pred.markers[[2]]))))
sum(gene.num[pred.markers[[2]]$p_val_adj>0.05] %in% which(sim$sim_data$F[,i]==1)) / sum(sim$sim_data$F[,2]==1)
