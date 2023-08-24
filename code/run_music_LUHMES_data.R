library(data.table)
library(Seurat)
library(MUSIC)
library(ggplot2)
theme_set(theme_bw() + theme(plot.title = element_text(size = 14, hjust = 0.5),
                             axis.title = element_text(size = 14),
                             axis.text = element_text(size = 13),
                             legend.title = element_text(size = 13),
                             legend.text = element_text(size = 12),
                             panel.grid.minor = element_blank())
)
library(ComplexHeatmap)

setwd("/project2/xinhe/kevinluo/GSFA/music_analysis/LUHMES")
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
                                       cpu_num = 10)
saveRDS(crop_seq_imputation, "music_output/music_imputation.merged.rds")

cat("Data filtering ... \n")
crop_seq_filtered <- Cell_filtering(crop_seq_imputation$expression_profile,
                                    crop_seq_imputation$perturb_information,
                                    cpu_num = 10)
saveRDS(crop_seq_filtered, "music_output/music_filtered.merged.rds")

# 2. Model building ####
cat("Model building ... \n")

crop_seq_vargene <- Get_high_varGenes(crop_seq_filtered$expression_profile,
                                      crop_seq_filtered$perturb_information, plot = T)
saveRDS(crop_seq_vargene, "music_output/music_vargene.merged.rds")

## Get_topics() can take up to a few hours to finish,
## depending on the size of data
cat("Fitting model with 5 topics ... \n")

system.time(
  topic_1 <- Get_topics(crop_seq_vargene$expression_profile,
                        crop_seq_vargene$perturb_information,
                        topic_number = 5))
saveRDS(topic_1, "music_output/music_merged_5_topics.rds")

cat("Fitting model with 10 topics ... \n")
system.time(
  topic_2 <- Get_topics(crop_seq_vargene$expression_profile,
                        crop_seq_vargene$perturb_information,
                        topic_number = 10))
saveRDS(topic_2, "music_output/music_merged_10_topics.rds")

cat("Fitting model with 15 topics ... \n")
system.time(
  topic_3 <- Get_topics(crop_seq_vargene$expression_profile,
                        crop_seq_vargene$perturb_information,
                        topic_number = 15))
saveRDS(topic_3, "music_output/music_merged_15_topics.rds")

cat("Fitting model with 20 topics ... \n")
system.time(
  topic_4 <- Get_topics(crop_seq_vargene$expression_profile,
                        crop_seq_vargene$perturb_information,
                        topic_number = 20))
saveRDS(topic_4, "music_output/music_merged_20_topics.rds")

# 3. Pick optimal number of topics ####
cat("Select optimal number of topics... \n")

topic_model_list <- list()
topic_model_list$models <- list()
topic_model_list$perturb_information <- topic_1$perturb_information
topic_model_list$models[[1]] <- topic_1$models[[1]]
topic_model_list$models[[2]] <- topic_2$models[[1]]
topic_model_list$models[[3]] <- topic_3$models[[1]]
topic_model_list$models[[4]] <- topic_4$models[[1]]
optimalModel <- Select_topic_number(topic_model_list$models,
                                    plot = T,
                                    plot_path = "music_output/select_topic_number_5_to_20.pdf")
# The paper said "the larger the score, the better the selected topic number".
# But we probably are going to set the number to 20 just to be comparable to GSFA.
