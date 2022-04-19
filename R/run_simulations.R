## This script performs simulation under one factor density setting
## for one rep of randomly simulated data (for both normal-based and count-based
## scenarios).
## Jobs can be submitted in parallel by varying the input arguments to this script.
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

out_folder <- opt$out_folder # "count_pi_0.05"
param_pi <- opt$pi # 0.05, 0.1 or 0.2
rep <- opt$rep # from 1 to 500

library(tidyverse)
library(GSFA)
seeds <- readRDS("../data/simulations/random_seeds.rds")
dir.create(out_folder)
if (!endsWith(out_folder, "/")){
  out_folder <- paste0(out_folder, "/")
}
out_file <- paste0(out_folder, "simulation_result.pi_", param_pi, ".rep_", rep, ".rds")
print(paste0("Simulate the ", rep, "th rep of random count-based datasets and perform GSFA."))
print(paste0("Results will be saved at ", out_file))

set.seed(seeds[rep])

params <- list(N = 4000, P = 6000, K = 10, M = 6, G_prob = 0.05,
               sigma_w2_true = rep(0.5, 10), psi_true = 1)
params$beta_true <- matrix(0, nrow = params$M + 1, ncol = params$K)
params$beta_true[1, 1] <- 0.1
params$beta_true[2, 2] <- 0.2
params$beta_true[3, 3] <- 0.3
params$beta_true[4, 4] <- 0.4
params$beta_true[5, 5] <- 0.5
params$beta_true[6, 6] <- 0.6
params$beta_true[params$M + 1, ] <- 0.5 # Intercept
params$pi_true <- rep(param_pi, params$K)

sim_data <- normal_data_sim(N = params$N, P = params$P, K = params$K, M = params$M,
                            beta_true = params$beta_true,
                            pi_true = params$pi_true, 
                            sigma_w2_true = params$sigma_w2_true,
                            psi_true = params$psi_true, G_prob = params$G_prob, offset = T)
sim_data <- poisson_count_sim_var_scale(sim_data = sim_data)
scaled_data <- deviance_residual_transform(sim_data$count)
sim_data$scaled <- scaled_data

fit_normal <- fit_gsfa_multivar(Y = sim_data$Y, G = sim_data$G, K = 10,
                                prior_type = "mixture_normal", init.method = "svd",
                                prior_w_s = 50, prior_w_r = 0.2,
                                prior_beta_s = 5, prior_beta_r = 0.2,
                                niter = 3000, used_niter = 1000,
                                verbose = T, return_samples = F)

fit_count <- fit_gsfa_multivar(Y = sim_data$scaled, G = sim_data$G, K = 10,
                               prior_type = "mixture_normal", init.method = "svd",
                               prior_w_s = 50, prior_w_r = 0.2,
                               prior_beta_s = 5, prior_beta_r = 0.2,
                               niter = 3000, used_niter = 1000,
                               verbose = T, return_samples = F)

res_list <- list(gsfa_normal = fit_normal,
                 gsfa_count = fit_count,
                 sim_data = sim_data,
                 seed = seeds[rep],
                 params = params)
saveRDS(res_list, out_file)

