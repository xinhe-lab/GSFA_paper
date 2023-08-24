#!/usr/bin/env Rscript

## Run MUSIC analysis for LUHMES data

## ---- load-packages------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(optparse)
library(data.table)
library(Matrix)
library(MUSIC)
dyn.load('/software/geos-3.7.0-el7-x86_64/lib64/libgeos_c.so') # attach the geos lib for Seurat
library(Seurat)

# Process the command-line arguments --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
parser <- OptionParser()
parser <- add_option(parser,c("--outdir","-o"),type="character",default=NULL)
parser <- add_option(parser,c("--ncpus","-c"),type="integer",default=1)
out    <- parse_args(parser)
outdir <- out$outdir
ncpus  <- out$ncpus

cat("outdir=", outdir, "\n")
cat("ncpus=", ncpus, "\n")

## ----Prepare input data for MUSIC -----------------------------------------
cat("Run MUSIC analysis for LUHMES data... \n")

ptm <- proc.time()

if(!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
setwd(outdir)
dir.create("./music_output", recursive = TRUE, showWarnings = FALSE)

# 0. Load input data ####
cat("Load input data ... \n")
feature.names <- data.frame(fread(paste0("/project2/xinhe/yifan/Factor_analysis/LUHMES/GSE142078_raw/GSM4219575_Run1_genes.tsv.gz"),
                                  header = FALSE), stringsAsFactors = FALSE)

# combined_obj <- readRDS("processed_data/seurat_obj.merged_scaled_detect_01.corrected_new.rds")
combined_obj <- readRDS("/project2/xinhe/yifan/Factor_analysis/shared_data/LUHMES_cropseq_data_seurat.rds")
expression_profile <- combined_obj@assays$RNA@counts
rownames(expression_profile) <- feature.names$V2[match(rownames(expression_profile),
                                                       feature.names$V1)]
targets <- names(combined_obj@meta.data)[4:18]
targets[11] <- "CTRL"
perturb_information <- apply(combined_obj@meta.data[4:18], 1,
                             function(x){ targets[which(x > 0)] })

# 1. Data preprocessing ####
cat("Data preprocessing ... \n")

cat("QC ... \n")
crop_seq_list <- Input_preprocess(expression_profile, perturb_information)
crop_seq_qc <- Cell_qc(crop_seq_list$expression_profile,
                       crop_seq_list$perturb_information,
                       species = "Hs", plot = F)

cat("Data imputation ... \n")
crop_seq_imputation <- Data_imputation(crop_seq_qc$expression_profile,
                                       crop_seq_qc$perturb_information,
                                       cpu_num = ncpus)
# saveRDS(crop_seq_imputation, "music_output/music_imputation.merged.rds")

cat("Data filtering ... \n")
crop_seq_filtered <- Cell_filtering(crop_seq_imputation$expression_profile,
                                    crop_seq_imputation$perturb_information,
                                    cpu_num = ncpus)
saveRDS(crop_seq_filtered, "music_output/music_filtered.merged.rds")

run.time1 <- proc.time() - ptm
cat("Preprocessing time: \n")
print(run.time1)

# 2. Model building -----------------------------------------

# Obtain highly dispersion differentially expressed genes
cat("Model building ... \n")
ptm <- proc.time()

crop_seq_vargene <- Get_high_varGenes(crop_seq_filtered$expression_profile,
                                      crop_seq_filtered$perturb_information, plot = T)
# saveRDS(crop_seq_vargene, "music_output/music_vargene.merged.rds")

# Fit topic model
# default 4:6 topics
cat("Fitting model with 4:6 topics ... \n")
system.time(
  topic_model_list <- Get_topics(crop_seq_vargene$expression_profile,crop_seq_vargene$perturb_information,topic_number=c(4:6)))

# set the number to 20 just to be comparable to GSFA.
cat("Fitting model with 20 topics ... \n")
system.time(
  topic_20 <- Get_topics(crop_seq_vargene$expression_profile,
                         crop_seq_vargene$perturb_information,
                         topic_number = 20))
saveRDS(topic_20, "music_output/music_merged_20_topics.rds")

# Summarize the results

# Summarize the results under 20 topics to be comparable to GSFA

# Gene ontology annotations for top topics
topic_res <- readRDS("music_output/music_merged_20_topics.rds")
topic_func <- Topic_func_anno(topic_res$models[[1]], species = "Hs")
saveRDS(topic_func, "music_output/topic_func.rds")

# Perturbation effect prioritizing
# calculate topic distribution for each cell.
distri_diff <- Diff_topic_distri(topic_res$models[[1]],
                                 topic_res$perturb_information,
                                 plot = T)

t_D_diff_matrix <- dcast(distri_diff %>% dplyr::select(knockout, variable, t_D_diff),
                         knockout ~ variable)
rownames(t_D_diff_matrix) <- t_D_diff_matrix$knockout
t_D_diff_matrix$knockout <- NULL

# Overall perturbation effect ranking list.
distri_diff <- readRDS(file.path(res_dir, "music_output/distri_diff.rds"))

rank_overall_result <- Rank_overall(distri_diff)
print(rank_overall_result)

# topic-specific ranking list.
rank_topic_specific_result <- Rank_specific(distri_diff)
head(rank_topic_specific_result, 10)

# perturbation correlation.
perturb_cor <- Correlation_perturbation(distri_diff,
                                        cutoff = 0.5, gene = "all", plot = T,
                                        plot_path = file.path(res_dir, "music_output/correlation_network_20_topics.pdf"))

head(perturb_cor, 10)

run.time2 <- proc.time() - ptm
cat("Run time for MUSIC: \n")
print(run.time2)

## ----- Session Info -----
sessionInfo()
