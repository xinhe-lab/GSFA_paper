#!/usr/bin/env Rscript

## Run SCEPTRE analysis for LUHMES data

## ---- load-packages------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
dyn.load('/software/geos-3.7.0-el7-x86_64/lib64/libgeos_c.so') # attach the geos lib for Seurat
library(optparse)
suppressPackageStartupMessages(library(tidyverse))
library(cowplot)
library(Matrix)
library(sceptre)
library(Seurat)

# Process the command-line arguments --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
parser <- OptionParser()
parser <- add_option(parser,c("--datadir","-d"),type="character",default=NULL)
parser <- add_option(parser,c("--outdir","-o"),type="character",default=NULL)
parser <- add_option(parser,c("--regularization","-r"),type="double",default=0.1)
parser <- add_option(parser,c("--B", "-b"),type="integer",default=1000)
parser <- add_option(parser,c("--seed","-s"),type="integer",default=4)
out    <- parse_args(parser)
datadir           <- out$datadir
outdir            <- out$outdir
regularization    <- out$regularization
B                 <- out$B
seed              <- out$seed

cat("datadir=", datadir, "\n")
cat("outdir=", outdir, "\n")
cat("regularization=", regularization, "\n")
cat("B=", B, "\n")
cat("seed=", seed, "\n")

## ----Prepare input data for SCEPTRE-----------------------------------------
ptm <- proc.time()

if(!dir.exists(datadir)) dir.create(datadir, recursive = TRUE)
if(!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

LUHMES_data <- readRDS('/project2/xinhe/yifan/Factor_analysis/shared_data/LUHMES_cropseq_data_seurat.rds')

# gene expression data
gene_matrix <- LUHMES_data@assays$RNA@counts
# gene-by-cell expression matrix
gene_matrix[1:10, 1:3]
dim(gene_matrix)

# cell metadata
metadata <- LUHMES_data@meta.data
covariate_matrix <- metadata[,c('orig.ident', 'nCount_RNA', 'nFeature_RNA', 'percent_mt')]
dim(covariate_matrix)
head(covariate_matrix)

# Perturbation matrix (a binary matrix of perturbations,  rows are gRNA groups and columns are cell barcodes)
combined_perturbation_matrix <- t(metadata[,4:18])
dim(combined_perturbation_matrix)
combined_perturbation_matrix[1:10,1:3]
range(combined_perturbation_matrix)

# Specify the gene-gRNA group pairs to test for association
# We include the 6k genes used for GSFA in this analysis
# Normalized and scaled data used for GSFA, the rownames of which are the 6k genes used for GSFA
scaled_gene_matrix <- LUHMES_data@assays$RNA@scale.data
dim(scaled_gene_matrix)
selected_gene_id <- rownames(scaled_gene_matrix)
all(selected_gene_id %in% rownames(gene_matrix))

gRNA_group <- rownames(combined_perturbation_matrix)
pairs <- expand.grid(selected_gene_id, gRNA_group)
gene_gRNA_group_pairs <- data.frame(gene_id = pairs$Var1, gRNA_group = pairs$Var2, pair_type = "candidate")
gene_gRNA_group_pairs[gene_gRNA_group_pairs$gRNA_group == "Nontargeting", "pair_type"] <- "negative_control"
table(gene_gRNA_group_pairs$pair_type)
table(gene_gRNA_group_pairs$gRNA_group)
dim(gene_gRNA_group_pairs)

save(list = c("gene_matrix", "combined_perturbation_matrix", "covariate_matrix"),
     file = file.path(datadir, 'data.matrices.RData'))
saveRDS(gene_gRNA_group_pairs, file.path(datadir, "gene.gRNA.group.pairs.rds"))

run.time1 <- proc.time() - ptm
cat("Preprocessing time: \n")
print(run.time1)

rm(LUHMES_data)
rm(scaled_gene_matrix)

## ----load-data-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat('Load data matrices...\n')
load(file.path(datadir, "data.matrices.RData"))
gene_gRNA_group_pairs <- readRDS(file.path(datadir, "gene.gRNA.group.pairs.rds"))

cat(sprintf('Dimenstion of gene expression matrix: %d rows %d columns.\n', nrow(gene_matrix), ncol(gene_matrix)))
cat(sprintf('Dimenstion of combined perturbation matrix: %d rows %d columns.\n', nrow(combined_perturbation_matrix), ncol(combined_perturbation_matrix)))
cat(sprintf('Dimenstion of covariate matrix: %d rows %d columns.\n', nrow(covariate_matrix), ncol(covariate_matrix)))
cat(sprintf('Dimenstion of gene gRNA-group pairs: %d rows %d columns.\n', nrow(gene_gRNA_group_pairs), ncol(gene_gRNA_group_pairs)))

table(gene_gRNA_group_pairs$pair_type)

## ----run-SCEPTRE------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

cat('Running SCEPTRE... \n')
cat(sprintf('Testing %d gene gRNA-group pairs. \n', nrow(gene_gRNA_group_pairs)))
ptm <- proc.time()

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

head(result, 10)

## ----neg-control------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# neg_control_p_vals <- result %>% filter(pair_type == 'negative_control') %>% pull(p_value)
# qq_plot <- make_qq_plot(neg_control_p_vals)
# plot(qq_plot)


## ----FDR--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat('Compute FDR... \n')
candidate_pair_results <- result %>% filter(pair_type == 'candidate')
candidate_pair_results_p_adj <- candidate_pair_results %>%
  mutate(p_val_adj = p.adjust(p_value, method = 'BH'))

head(candidate_pair_results_p_adj)

saveRDS(candidate_pair_results_p_adj, file.path(outdir, 'sceptre.candidate.pair.results.rds'))

## ----select-sig-pairs-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
discovery_set <- candidate_pair_results_p_adj %>% filter(p_val_adj <= 0.1)

head(discovery_set)
saveRDS(discovery_set, file.path(outdir, 'sceptre.discovery.set.results.rds'))

run.time2 <- proc.time() - ptm
cat("Run time for SCEPTRE: \n")
print(run.time2)

cat('SCEPTRE results saved to', outdir, '\n')

## ----- Session Info -----
sessionInfo()
