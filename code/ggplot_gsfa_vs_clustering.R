library(GSFA)
library(ggplot2)


theme_set(
  theme_bw() +
    theme(plot.title = element_text(size = 14, hjust = 0.5),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 12),
          legend.title = element_text(size = 13),
          legend.text = element_text(size = 12)
    )
)
setwd("~/Documents/Yifan_GSFA/")
GSFA.count <- list()
GSFA.normal <- list()

load("GSFA_fdp_power/GSFA_pi0.05_count_fdp_power.Rdata")
GSFA.count$pi5 <- list(fdp=fdp.gsfa.count, power=power.gsfa.count)
load("GSFA_fdp_power/GSFA_pi0.1_count_fdp_power.Rdata")
GSFA.count$pi10 <- list(fdp=fdp.gsfa.count, power=power.gsfa.count)
load("GSFA_fdp_power/GSFA_pi0.2_count_fdp_power.Rdata")
GSFA.count$pi20 <- list(fdp=fdp.gsfa.count, power=power.gsfa.count)

load("GSFA_fdp_power/GSFA_pi0.05_normal_fdp_power.Rdata")
GSFA.normal$pi5 <- list(fdp=fdp.gsfa.normal, power=power.gsfa.normal)
load("GSFA_fdp_power/GSFA_pi0.1_count_fdp_power.Rdata")
GSFA.normal$pi10 <- list(fdp=fdp.gsfa.normal, power=power.gsfa.normal)
load("GSFA_fdp_power/GSFA_pi0.2_count_fdp_power.Rdata")
GSFA.normal$pi20 <- list(fdp=fdp.gsfa.normal, power=power.gsfa.normal)

seurat.pi5 <- readRDS("two_step_clustering_fdp_power/Seurat_corr_fdp_power_nassoc.RDS")
seurat.pi10 <- readRDS("two_step_clustering_fdp_power/seurat_fdp_power_pi10.rds")
seurat.pi20 <- readRDS("two_step_clustering_fdp_power/seurat_fdp_power_pi20.rds")
kmeans.pi5 <- readRDS("two_step_clustering_fdp_power/kmeans_fdp_power.rds")
kmeans.pi10 <- readRDS("two_step_clustering_fdp_power/kmeans_fdp_power_pi10.rds")
kmeans.pi20 <- readRDS("two_step_clustering_fdp_power/kmeans_fdp_power_pi20.rds")

dat.kmeans.FDP <- data.frame(Performance=c(kmeans.pi5$fdp[kmeans.pi5$Nassoc>0], kmeans.pi10$fdp[kmeans.pi10$sel>0], kmeans.pi20$fdp[kmeans.pi20$sel>0]),
                         metric="Observed FDP", method="clustering", scenario="Normal scenario", 
                         pi=factor(c(rep("Factor Density=0.05", sum(kmeans.pi5$Nassoc>0)),rep("Factor Density=0.1", sum(kmeans.pi10$sel>0)), rep("Factor Density=0.2", sum(kmeans.pi20$sel>0))),
                                   levels=c("Factor Density=0.05","Factor Density=0.1","Factor Density=0.2")))
dat.kmeans.power <- data.frame(Performance=c(kmeans.pi5$power[kmeans.pi5$Nassoc>0], kmeans.pi10$power[kmeans.pi10$sel>0], kmeans.pi20$power[kmeans.pi20$sel>0]),
                             metric="Power", method="clustering", scenario="Normal scenario", 
                             pi=factor(c(rep("Factor Density=0.05", sum(kmeans.pi5$Nassoc>0)),rep("Factor Density=0.1", sum(kmeans.pi10$sel>0)), rep("Factor Density=0.2", sum(kmeans.pi20$sel>0))),
                                       levels=c("Factor Density=0.05","Factor Density=0.1","Factor Density=0.2")))

dat.seurat.FDP <- data.frame(Performance=c(seurat.pi5$fdp[seurat.pi5$Nassoc>0], seurat.pi10$fdp[seurat.pi10$sel>0], seurat.pi20$fdp[seurat.pi20$sel>0]),
                             metric="Observed FDP", method="clustering", scenario="Count scenario", 
                             pi=factor(c(rep("Factor Density=0.05", sum(seurat.pi5$Nassoc>0)),rep("Factor Density=0.1", sum(seurat.pi10$sel>0)), rep("Factor Density=0.2", sum(seurat.pi20$sel>0))),
                                       levels=c("Factor Density=0.05","Factor Density=0.1","Factor Density=0.2")))
dat.seurat.power <- data.frame(Performance=c(seurat.pi5$power[seurat.pi5$Nassoc>0], seurat.pi10$power[seurat.pi10$sel>0], seurat.pi20$power[seurat.pi20$sel>0]),
                               metric="Power", method="clustering", scenario="Count scenario", 
                               pi=factor(c(rep("Factor Density=0.05", sum(seurat.pi5$Nassoc>0)),rep("Factor Density=0.1", sum(seurat.pi10$sel>0)), rep("Factor Density=0.2", sum(seurat.pi20$sel>0))),
                                         levels=c("Factor Density=0.05","Factor Density=0.1","Factor Density=0.2")))

dat.GSFA.count.FDP <- data.frame(Performance=c(GSFA.count$pi5$fdp[!is.na(GSFA.count$pi5$fdp)], GSFA.count$pi10$fdp[!is.na(GSFA.count$pi10$fdp)], GSFA.count$pi20$fdp[!is.na(GSFA.count$pi20$fdp)]),
                             metric="Observed FDP", method="GSFA", scenario="Count scenario", 
                             pi=factor(c(rep("Factor Density=0.05", sum(!is.na(GSFA.count$pi5$fdp))),rep("Factor Density=0.1", sum(!is.na(GSFA.count$pi10$fdp))), rep("Factor Density=0.2", sum(!is.na(GSFA.count$pi20$fdp)))),
                                       levels=c("Factor Density=0.05","Factor Density=0.1","Factor Density=0.2")))
dat.GSFA.count.power <- data.frame(Performance=c(GSFA.count$pi5$power[!is.na(GSFA.count$pi5$fdp)], GSFA.count$pi10$power[!is.na(GSFA.count$pi10$fdp)], GSFA.count$pi20$power[!is.na(GSFA.count$pi20$fdp)]),
                           metric="Power", method="GSFA", scenario="Count scenario", 
                           pi=factor(c(rep("Factor Density=0.05", sum(!is.na(GSFA.count$pi5$fdp))),rep("Factor Density=0.1", sum(!is.na(GSFA.count$pi10$fdp))), rep("Factor Density=0.2", sum(!is.na(GSFA.count$pi20$fdp)))),
                                     levels=c("Factor Density=0.05","Factor Density=0.1","Factor Density=0.2")))

dat.GSFA.normal.FDP <- data.frame(Performance=c(GSFA.normal$pi5$fdp[!is.na(GSFA.normal$pi5$fdp)], GSFA.normal$pi10$fdp[!is.na(GSFA.normal$pi10$fdp)], GSFA.normal$pi20$fdp[!is.na(GSFA.normal$pi20$fdp)]),
                                 metric="Observed FDP", method="GSFA", scenario="Normal scenario", 
                                 pi=factor(c(rep("Factor Density=0.05", sum(!is.na(GSFA.normal$pi5$fdp))),rep("Factor Density=0.1", sum(!is.na(GSFA.normal$pi10$fdp))), rep("Factor Density=0.2", sum(!is.na(GSFA.normal$pi20$fdp)))),
                                           levels=c("Factor Density=0.05","Factor Density=0.1","Factor Density=0.2")))
dat.GSFA.normal.power <- data.frame(Performance=c(GSFA.normal$pi5$power[!is.na(GSFA.normal$pi5$fdp)], GSFA.normal$pi10$power[!is.na(GSFA.normal$pi10$fdp)], GSFA.normal$pi20$power[!is.na(GSFA.normal$pi20$fdp)]),
                                   metric="Power", method="GSFA", scenario="Normal scenario", 
                                   pi=factor(c(rep("Factor Density=0.05", sum(!is.na(GSFA.normal$pi5$fdp))),rep("Factor Density=0.1", sum(!is.na(GSFA.normal$pi10$fdp))), rep("Factor Density=0.2", sum(!is.na(GSFA.normal$pi20$fdp)))),
                                             levels=c("Factor Density=0.05","Factor Density=0.1","Factor Density=0.2")))


# plot Count scenario
ggplot(rbind(dat.seurat.FDP, dat.seurat.power, dat.GSFA.count.FDP, dat.GSFA.count.power), 
       aes(y=Performance, x=method, col=metric)) + 
  geom_boxplot() + facet_grid(metric~pi) +
  xlab("") + theme(legend.position = "none") +
  scale_y_continuous(breaks=0:5*0.2)

# plot Normal scenario
ggplot(rbind(dat.kmeans.FDP, dat.kmeans.power, dat.GSFA.normal.FDP, dat.GSFA.normal.power), 
       aes(y=Performance, x=method, col=metric)) + 
  geom_boxplot() + facet_grid(metric~pi) +
  xlab("") + theme(legend.position = "none") + 
  scale_y_continuous(breaks=0:5*0.2)
