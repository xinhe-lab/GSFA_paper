#!/usr/bin/env Rscript

## Run SCEPTRE analysis

## ---- load-packages------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
dyn.load('/software/geos-3.7.0-el7-x86_64/lib64/libgeos_c.so') # attach the geos lib for Seurat
library(optparse)
suppressPackageStartupMessages(library(tidyverse))
library(cowplot)
library(Matrix)
library(Seurat)
library(sceptre)

# Process the command-line arguments --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
parser <- OptionParser()
parser <- add_option(parser,c("--datadir","-d"),type="character",default=NULL)
parser <- add_option(parser,c("--outdir","-o"),type="character",default=NULL)
parser <- add_option(parser,c("--regularization","-r"),type="double",default=0.1)
parser <- add_option(parser,c("--B", "-b"),type="integer",default=1000)
parser <- add_option(parser,c("--seed","-s"),type="integer",default=4)
parser <- add_option(parser,c("--npairs","-n"),type="integer",default=NULL)
parser <- add_option(parser,"--permute", action="store_true",type="logical",default=FALSE)
parser <- add_option(parser,"--perm_num", type="integer",default=1)
out    <- parse_args(parser)
datadir           <- out$datadir
outdir            <- out$outdir
regularization    <- out$regularization
B                 <- out$B
seed              <- out$seed
npairs            <- out$npairs
permute           <- out$permute
perm_num          <- out$perm_num

cat("datadir=", datadir, "\n")
cat("outdir=", outdir, "\n")
cat("regularization=", regularization, "\n")
cat("B=", B, "\n")
cat("seed=", seed, "\n")
cat("permute=", permute, "\n")
cat("perm_num=", perm_num, "\n")

## ----load-data-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat('Load data matrices...\n')
load(file.path(datadir, "data.matrices.RData"))
gene_gRNA_group_pairs <- readRDS(file.path(datadir, "gene.gRNA.group.pairs.rds"))

cat(sprintf('Dimenstion of gene expression matrix: %d rows %d columns.\n', nrow(gene_matrix), ncol(gene_matrix)))
cat(sprintf('Dimenstion of combined perturbation matrix: %d rows %d columns.\n', nrow(combined_perturbation_matrix), ncol(combined_perturbation_matrix)))
cat(sprintf('Dimenstion of covariate matrix: %d rows %d columns.\n', nrow(covariate_matrix), ncol(covariate_matrix)))
cat(sprintf('Dimenstion of gene gRNA-group pairs: %d rows %d columns.\n', nrow(gene_gRNA_group_pairs), ncol(gene_gRNA_group_pairs)))

table(gene_gRNA_group_pairs$pair_type)

# permutation
if (!permute){
  stopifnot(colnames(combined_perturbation_matrix) == colnames(gene_matrix))
  set.seed(seed)
  if(!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
}else{
  # Perform one of ten permutation tests
  stopifnot(perm_num %in% 1:10)
  seeds <- c(49553, 72704, 11932, 56826, 49707, 33357, 93747, 95392, 96675, 38186)
  random_seed <- seeds[perm_num]
  print(paste("Permuting cell orders under seed", random_seed, "..."))
  set.seed(random_seed)
  original_cell_ids <- colnames(combined_perturbation_matrix)
  new_cell_order <- sample(ncol(combined_perturbation_matrix))
  combined_perturbation_matrix <- combined_perturbation_matrix[, new_cell_order]
  colnames(combined_perturbation_matrix) <- original_cell_ids
  outdir <- paste0(outdir, "/permutation/", "perm_", perm_num)
  if(!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
}

save(list = c("gene_matrix", "combined_perturbation_matrix", "covariate_matrix"),
     file = file.path(outdir, 'data.matrices.RData'))

if(!is.null(npairs)){
  cat("Test run on", npairs, "gene-gRNA pairs ...\n")
  gene_gRNA_group_pairs <- gene_gRNA_group_pairs %>% sample_n(npairs)
}

saveRDS(gene_gRNA_group_pairs, file.path(outdir, "gene.gRNA.group.pairs.rds"))

## ----run-sceptre------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

if(!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

cat('Running SCEPTRE... \n')
cat(sprintf('Testing %d gene gRNA-group pairs. \n', nrow(gene_gRNA_group_pairs)))

timing <- system.time(
  result <- run_sceptre_high_moi(gene_matrix = gene_matrix,
                                 combined_perturbation_matrix = combined_perturbation_matrix,
                                 covariate_matrix = covariate_matrix,
                                 gene_gRNA_group_pairs = gene_gRNA_group_pairs,
                                 side = 'both',
                                 storage_dir = outdir,
                                 regularization_amount = regularization,
                                 B = B,
                                 seed = seed,
                                 full_output = FALSE))

saveRDS(result, file.path(outdir, 'sceptre.result.rds'))
cat(sprintf("Computation took %0.2f seconds. \n",timing["elapsed"]))

## ----neg-control------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
neg_control_p_vals <- result %>% filter(pair_type == 'negative_control') %>% pull(p_value)
pdf(file.path(outdir, 'sceptre.neg.control.qqplot.pdf'))
qq_plot <- make_qq_plot(neg_control_p_vals)
plot(qq_plot)
dev.off()

## ----FDR--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat('Compute FDR... \n')
candidate_pair_results <- result %>% filter(pair_type == 'candidate')
candidate_pair_results_p_adj <- candidate_pair_results %>%
  mutate(p_val_adj = p.adjust(p_value, method = 'BH'))

saveRDS(candidate_pair_results_p_adj, file.path(outdir, 'sceptre.candidate.pair.results.rds'))

## ----select-sig-pairs-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat('Discovery set... \n')
discovery_set <- candidate_pair_results_p_adj %>% filter(p_val_adj <= 0.1)
saveRDS(discovery_set, file.path(outdir, 'sceptre.discovery.set.results.rds'))

cat('Results saved to', outdir, '\n')

## -----Session Info -------------------------------------------------------------------------------------------------------------------------------
sessionInfo()
