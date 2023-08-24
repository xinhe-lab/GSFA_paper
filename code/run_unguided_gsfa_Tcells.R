library(GSFA)
library(data.table)
library(tidyverse)
library(optparse)

option_list <- list(
  make_option(
    "--out_folder", action = "store", default = NA, type = 'character',
    help = "Directory of output folder (must end with a \'/\') [required]"
  ),
  make_option(
    "--condition", action = "store", default = 'stimulated', type = 'character',
    help = "Treatment condition: stimulated or unstimulated. %default [optional]"
  ),
  make_option(
    "--permute", action = "store", default = FALSE,
    help = "Flag to perform permutation on the cells before GSFA, default is %default [optional]"
  ),
  make_option(
    "--perm_num",  action = "store", default = 1, type = 'integer',
    help = "Permutation index, default is %default [required when permute==TRUE]"
  ),
  make_option(
    "--init_method", action = "store", default = "svd", type = 'character',
    help = "Type of initialization method to use, default is %default [not required when restart=TRUE]"
  ),
  make_option(
    "--out_suffix", action = "store", default = "svd", type = 'character',
    help = "Suffix of the output file, reflecting the initializtion method, default is %default [not required when restart=TRUE]"
  ),
  make_option(
    "--K",  action = "store", default = 20, type = 'integer',
    help = "Number of factors to use, default is %default [auto-detected when restart=TRUE]"
  ),
  make_option(
    "--niter", action = "store", default = 1000, type = 'integer',
    help = "Number of iterations to sample, default is %default [optional]"
  ),
  make_option(
    "--average_niter", action = "store", default = 500, type = 'integer',
    help = "Number of iterations to average over to obtain the posterior means, default is %default [optional]"
  ),
  make_option(
    "--random_seed", action = "store", default = 92629, type = 'integer',
    help = "Set a random seed for Gibbs sampling, default is %default [required for all types of initializations]"
  ),
  make_option(
    "--store_all_samples", action = "store", default = TRUE,
    help = "Flag to store samples throughout Gibbs iterations, can be turned off if storage is limited, default is %default"
  )
)
opt <- parse_args(OptionParser(option_list = option_list))

out_folder <- opt$out_folder # "gsfa_output_detect_01/dev_res_corrected/"
condition <- opt$condition # "stimulated", "unstimulated"
init_method <- opt$init_method # "svd", "random" or "given"
out_suffix <- opt$out_suffix # "svd", "rand_01"
K <- opt$K
niter <- opt$niter
average_niter <- opt$average_niter
random_seed <- opt$random_seed
return_samples <- opt$store_all_samples

Rcpp::sourceCpp('/home/kaixuan/projects/GSFA_analysis/code/unguided_GFSA_mixture_normal_prior.cpp')

setwd("/project2/xinhe/kevinluo/GSFA/unguided_GSFA/Stimulated_T_Cells/")

print("Loading inputs for GSFA...")
# Processed gene expression matrix:
scaled.gene_exp <- readRDS("processed_data/deviance_residual.all_T_cells_merged_top_6k.batch_uncorrected.rds")
# Cell meta data:
metadata <- readRDS("processed_data/metadata.all_T_cells_merged.rds")
# Select cell group:
sample_group <- endsWith(metadata$orig.ident, "S") * 1 # 0: unstimulated, 1: stimulated

if(condition == "stimulated"){
  scaled.gene_exp <- scaled.gene_exp[, which(sample_group == 1)]
  metadata <- metadata[which(sample_group == 1), ]
}else if(condition == "unstimulated"){
  scaled.gene_exp <- scaled.gene_exp[, which(sample_group == 0)]
  metadata <- metadata[which(sample_group == 0), ]
}else{
  stop("condition needs to be stimulated or unstimulated!")
}

# Perturbation info:
G_mat <- metadata[, 4:24]
G_mat <- as.matrix(G_mat)
KO_names <- colnames(G_mat)
negctrl_index <- which(KO_names == "NonTarget")

if(!dir.exists(out_folder)){dir.create(out_folder)}

if (!opt$permute){
  stopifnot(rownames(G_mat) == colnames(scaled.gene_exp))
  set.seed(random_seed)
  out_dir <- paste0(out_folder, "All.gibbs_obj_k", K, ".", out_suffix, ".rds")
  if (file.exists(out_dir)){
    warnings("Output file exists and will be overwritten!")
  }
} else {
  stopifnot(opt$perm_num %in% 1:10)
  seeds <- c(49553, 72704, 11932, 56826, 49707, 33357, 93747, 95392, 96675, 38186)
  print(paste("Permuting cell orders under seed", seeds[opt$perm_num], "..."))
  set.seed(seeds[opt$perm_num])
  new_cell_order <- sample(ncol(scaled.gene_exp))
  # new_cell_order <- rep(NA, ncol(scaled.gene_exp))
  # new_cell_order[which(sample_group==0)] <- gtools::permute(which(sample_group==0))
  # new_cell_order[which(sample_group==1)] <- gtools::permute(which(sample_group==1))
  scaled.gene_exp <- scaled.gene_exp[, new_cell_order]
  if (opt$restart == F){
    out_dir <- paste0(out_folder, "All.gibbs_obj_k", K, ".", paste0("perm_", opt$perm_num), ".rds")
  } else {
    out_dir <- paste0(out_folder, "All.gibbs_obj_k", K, ".", paste0("perm_", opt$perm_num), ".restart.rds")
  }
}

print(paste0("Performing unguided GSFA with ", K, " factors and ", niter, " iterations."))
print(paste0("Results will be saved at ", out_dir))

Y <- t(scaled.gene_exp)
dim(Y)

fit <- gsfa_gibbs_cpp(Y = Y, K = K,
                      prior_type = "mixture_normal", initialize = init_method,
                      prior_s=50, prior_r=0.2,
                      prior_sb=20, prior_rb=0.2,
                      niter = niter, ave_niter = average_niter,
                      verbose = T, return_samples = return_samples)

factor_names <- paste0("Factor_", 1:ncol(fit$posterior_means$Z_pm))
if (is.null(rownames(Y))){
  sample_names <- 1:nrow(Y)
} else {
  sample_names <- rownames(Y)
}
if (is.null(colnames(Y))){
  gene_names <- 1:ncol(Y)
} else {
  gene_names <- colnames(Y)
}

rownames(fit$posterior_means$Z_pm) <- sample_names
colnames(fit$posterior_means$Z_pm) <- factor_names
rownames(fit$posterior_means$F_pm) <- gene_names
colnames(fit$posterior_means$F_pm) <- factor_names
rownames(fit$posterior_means$W_pm) <- gene_names
colnames(fit$posterior_means$W_pm) <- factor_names
class(fit) <- c("gsfa_fit", "list")

print("Finished! Saving the results...")
if (return_samples){
  saveRDS(fit, file = out_dir)
}
names(fit)

big_items <- names(fit)[endsWith(names(fit), "samples")]
for (i in big_items){
  fit[[i]] <- NULL
}
saveRDS(fit, file = sub(pattern = ".rds", replacement = ".light.rds", x = out_dir))

