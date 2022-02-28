## Recommend to run this script on a computing cluster that allows for
## ~ 50GB memory and ~ 5.5 hours of no interruption.
## GSFA results were generated under default parameters with the output of
## "scaled.gene_exp", "G_mat" and "group" from "preprocess_TCells.R" as GSFA input.
## Gibbs sampling was performed for 2000 iterations in a first run,
## and then carried on for another 2000 iterations in a second run
## using the "restart" option.
## If limited run time and memory can be allocated, set "niter" to a smaller
## number and run more segments.
## Four output files:
## A full Gibbs sampling object; posterior mean estimates of parameters; 
## two LFSR matrices (one for each cell group).
library(GSFA)
library(optparse)
option_list <- list(
  make_option(
    "--expression_file", action = "store", default = NA, type = 'character',
    help = "RDS file of processed sample by gene expression matrix for GSFA input (Y) [required]"
  ),
  make_option(
    "--perturbation_file", action = "store", default = NA, type = 'character',
    help = "RDS file of sample by perturbation matrix for GSFA input (G) [required]"
  ),
  make_option(
    "--sample_group_file", action = "store", default = NA, type = 'character',
    help = "RDS file of sample group information for GSFA by group [required]"
  ),
  make_option(
    "--out_folder", action = "store", default = NA, type = 'character',
    help = "Directory of output folder (must end with a \'/\') [required]"
  ),
  make_option(
    "--restart", action = "store", default = FALSE,
    help = "Flag to resume GSFA on previous result, default is %default"
  ),
  make_option(
    "--previous_res", action = "store", default = NA, type = 'character',
    help = "RDS file storing previous GSFA result to resume GSFA on [required when restart=TRUE]"
  ),
  make_option(
    "--permute", action = "store", default = FALSE,
    help = "Flag to perform permutation on the cells before GSFA, default is %default"
  ),
  make_option(
    "--perm_num",  action = "store", default = 1, type = 'integer',
    help = "Permutation index from 1 to 10, default is %default [required when permute==TRUE]"
  ),
  make_option(
    "--init_method", action = "store", default = "svd", type = 'character',
    help = "Type of initialization method to use, can be \"svd\" or \"random\", default is %default"
  ),
  make_option(
    "--out_suffix", action = "store", default = "svd", type = 'character',
    help = "Suffix of the output file, reflecting the initializtion method, default is %default"
  ),
  make_option(
    "--K",  action = "store", default = 20, type = 'integer',
    help = "Number of factors to use, default is %default [auto-detected when restart=TRUE]"
  ),
  make_option(
    "--niter", action = "store", default = 2000, type = 'integer',
    help = "Number of iterations to sample, default is %default"
  ),
  make_option(
    "--average_niter", action = "store", default = 1000, type = 'integer',
    help = "Number of iterations to average over to obtain the posterior means, default is %default"
  ),
  make_option(
    "--random_seed", action = "store", default = 92629, type = 'integer',
    help = "Set a random seed for Gibbs sampling, default is %default [required for all types of initializations]"
  )
)
opt <- parse_args(OptionParser(option_list = option_list))

expression_file <- opt$expression_file
perturbation_file <- opt$perturbation_file
sample_group_file <- opt$sample_group_file
out_folder <- opt$out_folder
init_method <- opt$init_method # "svd"
out_suffix <- opt$out_suffix
K <- opt$K # 20
niter <- opt$niter
average_niter <- opt$average_niter # 1000
random_seed <- opt$random_seed # 92629

stopifnot(file.exists(expression_file) & file.exists(perturbation_file) &
            file.exists(sample_group_file) & dir.exists(out_folder))
print("Loading input data ...")
scaled.gene_exp <- readRDS(expression_file)
scaled.gene_exp <- t(scaled.gene_exp)
G_mat <- readRDS(perturbation_file)
sample_group <- readRDS(sample_group_file)
stopifnot(rownames(scaled.gene_exp) == rownames(G_mat) & 
            rownames(G_mat) == names(sample_group))

sample_group <- (sample_group %in% c("TCells_D1S", "TCells_D2S")) * 1
## 0: unstimulated, 1: stimulated

if (opt$restart){
  stopifnot(file.exists(opt$previous_res))
  prev_fit <- readRDS(opt$previous_res)
}

if (!opt$permute){
  set.seed(random_seed)
  out_dir <- paste0(out_folder, "gibbs_obj_k", K, ".", out_suffix, ".rds")
} else {
  stopifnot(opt$perm_num %in% 1:10)
  seeds <- c(49553, 72704, 11932, 56826, 49707, 33357, 93747, 95392, 96675, 38186)
  random_seed <- seeds[opt$perm_num]
  print(paste("Permuting cell orders under seed", random_seed, "..."))
  set.seed(random_seed)
  ## Permutation within each cell group:
  new_cell_order <- rep(NA, ncol(scaled.gene_exp))
  new_cell_order[which(sample_group==0)] <- gtools::permute(which(sample_group==0))
  new_cell_order[which(sample_group==1)] <- gtools::permute(which(sample_group==1))
  scaled.gene_exp <- scaled.gene_exp[, new_cell_order]
  out_dir <- paste0(out_folder, "gibbs_obj_k", K, ".", paste0("perm_", opt$perm_num), 
                    ".", out_suffix, ".rds")
}

if (file.exists(out_dir)){
  warnings("Output file exists and will be overwritten!")
}

if (opt$restart == F){
  print(paste0("Performing GSFA with ", K, " factors and ", niter, " iterations."))
  print(paste0("Results will be saved at ", out_dir))
  fit <- fit_gsfa_multivar_2groups(Y = scaled.gene_exp, G = G_mat, 
                                   group = sample_group, K = K,
                                   prior_type = "mixture_normal", init.method = init_method,
                                   prior_w_s = 50, prior_w_r = 0.2,
                                   prior_beta_s = 20, prior_beta_r = 0.2,
                                   niter = niter, used_niter = average_niter,
                                   verbose = T, return_samples = T)
} else {
  print(paste0("Resuming GSFA on top of a previous result saved at ", opt$previous_res,
               " for ", niter, " iterations."))
  print(paste0("Results will be saved at ", out_dir))
  fit <- fit_gsfa_multivar_2groups(Y = scaled.gene_exp, G = G_mat, 
                                   group = sample_group, fit0 = prev_fit,
                                   prior_type = "mixture_normal", init.method = init_method,
                                   prior_w_s = 50, prior_w_r = 0.2,
                                   prior_beta_s = 20, prior_beta_r = 0.2,
                                   niter = niter, used_niter = average_niter,
                                   verbose = T, return_samples = T)
}
print("Finished! Saving the results...")
saveRDS(fit, file = out_dir)

saveRDS(fit$posterior_means, file = sub(pattern = ".rds", replacement = ".PM.rds", x = out_dir))

lfsr0_mat <- fit$lfsr0
lfsr0_mat <- lfsr0_mat[, -ncol(lfsr0_mat)]
print("Group 0, # of genes that pass LFSR < 0.05:")
print(colSums(lfsr0_mat < 0.05))
saveRDS(lfsr0_mat, file = sub(pattern = ".rds", replacement = ".lfsr_mat.group0.rds", x = out_dir))
  
lfsr1_mat <- fit$lfsr1
lfsr1_mat <- lfsr1_mat[, -ncol(lfsr1_mat)]
print("Group 1, # of genes that pass LFSR < 0.05 (computed from all factors):")
print(colSums(lfsr1_mat < 0.05))
saveRDS(lfsr1_mat, file = sub(pattern = ".rds", replacement = ".lfsr_mat.group1.rds", x = out_dir))
