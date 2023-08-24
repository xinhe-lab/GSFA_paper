library(optparse)
option_list <- list(
  make_option(
    "--out_folder", action = "store", default = NA, type = 'character',
    help = "Folder where the outputs are going to be saved [required]"
  ),
  make_option(
    "--pi", action = "store", default = 0.2, type = 'double',
    help = "The density of factors to simulated, must be between 0 and 1, default is %default [optional]"
  ),
  make_option(
    "--rep", action = "store", default = 1, type = 'integer',
    help = "The index of the rep of random dataset to simulate, default is %default [optional]"
  )
)
opt <- parse_args(OptionParser(option_list = option_list))

out_folder <- opt$out_folder # "unguided.normal_mixture.pi_0.1/single_datasets/"
param_pi <- opt$pi # 0.05, 0.1 or 0.2
rep <- opt$rep # from 1 to 300
reps <- c(rep, rep+100, rep+200)

#library(tidyverse)
library(foreach)
library(GSFA)
#Rcpp::sourceCpp('/project2/xinhe/yifan/Factor_analysis/rcpp/unguided_GFSA_mixture_normal_prior.cpp')
Rcpp::sourceCpp('unguided_GFSA_mixture_normal_prior.cpp')
#setwd("/project2/xinhe/yifan/Factor_analysis/simulations/")
seeds <- readRDS("1K_random_seeds.rds")

#out_folder <- paste0("simulation_results/multi_factors_n_markers/realistic_count_beta_0.1_to_0.6/",out_folder)
#dir.create(out_folder, recursive = TRUE)

params <- list(N = 4000, P = 6000, K = 10, M = 6, G_prob = 0.05,
               sigma_w2_true = rep(0.5, 10), psi_true = 1)
params$beta_true <- matrix(0, nrow = params$M + 1, ncol = params$K)
params$beta_true[1, 1] <- 0.4
params$beta_true[1, 2] <- 0.4
params$beta_true[1, 3] <- 0.4
params$beta_true[2, 2] <- 0.4
params$beta_true[2, 3] <- 0.4
params$beta_true[2, 4] <- 0.4
params$beta_true[3, 3] <- 0.4
params$beta_true[3, 4] <- 0.4
params$beta_true[3, 5] <- 0.4
params$beta_true[4, 4] <- 0.4
params$beta_true[4, 5] <- 0.4
params$beta_true[4, 6] <- 0.4
params$beta_true[5, 5] <- 0.4
params$beta_true[5, 6] <- 0.4
params$beta_true[5, 1] <- 0.4
params$beta_true[6, 6] <- 0.4
params$beta_true[6, 1] <- 0.4
params$beta_true[6, 2] <- 0.4

params$pi_true <- rep(param_pi, params$K)

for (rep in reps) {
set.seed(seeds[rep])

out_file <- paste0(out_folder, "fit_result.rep_", rep, ".rds")
print(paste0("Simulate the ", rep, "th rep of random count-based datasets and perform unguided GSFA."))
print(paste0("Results will be saved at ", out_file))
  
sim_data <- normal_data_sim(N = params$N, P = params$P, K = params$K, M = params$M,
                            beta_true = params$beta_true,
                            pi_true = params$pi_true,
                            sigma_w2_true = params$sigma_w2_true,
                            psi_true = params$psi_true, G_prob = params$G_prob, offset = T)
sim_data <- poisson_count_sim_var_scale(sim_data = sim_data)
#scaled_data <- deviance_residual_transform(sim_data$count)
#sim_data$scaled <- scaled_data

fit_unguided <- gsfa_gibbs_cpp(Y = sim_data$Y, K = 10,
                               prior_type = "mixture_normal", initialize = "svd",
                               prior_s=50, prior_r=0.2,
                               prior_sb=5, prior_rb=0.2,
                               niter = 3000, ave_niter = 1000,
                               verbose = T, return_samples = F)
# Identy DEGs based on Pr(F=1)>0.95
deg.factor <- foreach(i=1:10) %do% {which(fit_unguided$posterior_means$F_pm[,i]>0.95)}
# Associate factors to guide with marginal linear regression
temp <- as.data.frame(cbind(sim_data$G, fit_unguided$posterior_means$Z_pm))
colnames(temp) <- c(paste0("guide", 1:6), paste0("factor", 1:10))
beta.pval <- matrix(0, nrow=6, ncol=10)
for (f in 1:6) {
  beta.pval[f,] <- foreach(i=1:10, .combine = cbind) %do% {
    form.cur <- as.formula(paste0("factor",i, "~guide",f))
    res <- lm(form.cur, data=temp)
    summary(res)$coefficients[2,4]
  }
}
# aggregate DEGs from factors to guides.
beta.qval <- matrix(p.adjust(c(beta.pval),"BH"), nrow=6)
proxies <- foreach(i=1:6) %do% which(beta.qval[i,]<0.05)
deg.guide <- foreach(i=1:6) %do% unique(unlist(deg.factor[proxies[[i]]]))

# Compute ground truth theta and compute FPR and FNR
theta.true <- (sim_data$F * sim_data$U) %*% t(as.matrix(params$beta_true[1:6,]))
fpr.fa <- foreach(i=1:6, .combine=c) %do% {sum(deg.guide[[i]] %in% which(theta.true[,i]==0)) / sum(theta.true[,i]==0)}
fnr.fa <- 1 - foreach(i=1:6, .combine=c) %do% {sum(deg.guide[[i]] %in% which(abs(theta.true[,i])>0)) / sum(abs(theta.true[,i])>0)}

### Use GSFA for comparison
cat("Running guided GSFA\n")
fit_normal <- fit_gsfa_multivar(Y = sim_data$Y, G = sim_data$G, K = 10,
                               prior_type = "mixture_normal", init.method = "svd",
                               prior_w_s = 50, prior_w_r = 0.2,
                               prior_beta_s = 5, prior_beta_r = 0.2,
                               niter = 3000, used_niter = 1000,
                               verbose = T, return_samples = F)

#theta.true <- (sim_data$F * sim_data$U) %*% t(as.matrix(params$beta_true[1:6,]))
fpr.gsfa <- foreach(i=1:6, .combine=c) %do% {sum(which(fit_normal$lfsr[,i]<0.05) %in% which(theta.true[,i]==0)) / sum(theta.true[,i]==0)}
fnr.gsfa <- foreach(i=1:6, .combine=c) %do% {sum(which(fit_normal$lfsr[,i]>0.05) %in% which(abs(theta.true[,i])>0)) / sum(abs(theta.true[,i])>0)}

res_list <- list(G0 = fit_unguided,
                 G1 = fit_normal,
                 seed = seeds[rep],
                 params = params,
                 fpr.gsfa = fpr.gsfa,
                 fnr.fa = fnr.fa,
                 fnr.gsfa = fnr.gsfa,
                 fpr.fa = fpr.fa)
saveRDS(res_list, file=out_file)
}
