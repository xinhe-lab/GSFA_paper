factor_matrix_regression <- function(factors, G_mat){
  # factors: nxk, G_mat: nxm
  pval_tb <- matrix(nrow = ncol(factors), ncol = ncol(G_mat))
  beta_reg_tb <- pval_tb
  if (is.null(colnames(G_mat))){
    colnames(G_mat) <- 1:ncol(G_mat)
  }
  colnames(pval_tb) <- paste0("pval-", colnames(G_mat))
  colnames(beta_reg_tb) <- paste0("beta_reg-", colnames(G_mat))
  for (i in 1:ncol(factors)){
    for (j in 1:ncol(G_mat)){
      regress_mat <- data.frame(G = G_mat[, j], Z = factors[, i])
      lm.summary <- summary(lm(Z ~ G, data = regress_mat))
      pval_tb[i, j] <- lm.summary$coefficients[2, 4]
      beta_reg_tb[i, j] <- lm.summary$coefficients[2, 1]
    }
  }
  return(list(pval = pval_tb, beta = beta_reg_tb))
}

make_gibbs_res_tb <- function(gibbs_obj, G, compute_pve = FALSE, cell_indx = NULL){
  names(gibbs_obj) <- sub(pattern = "[.]", replacement = "_", x = names(gibbs_obj))
  if (!is.null(cell_indx)){
    gibbs_obj$Z_pm <- gibbs_obj$Z_pm[cell_indx, ]
    G <- G[cell_indx, ]
  }
  K <- ncol(gibbs_obj$Z_pm)
  res_tb <- data.frame(index = 1:K, pi = gibbs_obj$pi_pm[, 1])
  if (compute_pve){
    sum_var <- rep(0, K)
    for (i in 1:K) {
      mat_tmp <- outer(gibbs_obj$Z_pm[, i], gibbs_obj$W_pm[, i])
      sum_var[i] <- sum(apply(mat_tmp, 2, var))
    }
    res_tb$sum_var <- sum_var
  }
  if (is.null(ncol(G))){
    if (length(G) != nrow(gibbs_obj$Z_pm)){
      stop("Number of samples in G and in the gibbs result do not match!")
    }
    res_tb$beta_reg <- rep(NA, K)
    res_tb$pval <- rep(NA, K)
    for (i in 1:K){
      lm.summary <- summary(lm(gibbs_obj$Z_pm[, i] ~ G))
      res_tb$beta_reg[i] <- lm.summary$coefficients[2, 1]
      res_tb$pval[i] <- lm.summary$coefficients[2, 4]
    }
  } else {
    if (nrow(G) != nrow(gibbs_obj$Z_pm)){
      stop("Number of samples in genotype and in the gibbs object do not match!")
    }
    res_regress <- factor_matrix_regression(gibbs_obj$Z_pm, G)
    res_tb <- cbind(res_tb, res_regress$beta, res_regress$pval)
  }
  return(res_tb)
}

dotplot_beta_PIP <- function(beta_pip_matrix, beta_pm_matrix,
                             marker_names, reorder_markers = marker_names,
                             reorder_factors = NULL,
                             exclude_offset = TRUE,
                             inverse_factors = TRUE,
                             return_dataframe = FALSE,
                             color_lgd_title = "Estimated effect size"){
  # Both 'beta_pip_matrix' and 'beta_pm_matrix' should be factor by guide/marker matrices,
  # Dots will be colored by effect size and sized by PIP value.
  if (exclude_offset){
    beta_pip_matrix <- beta_pip_matrix[, -ncol(beta_pip_matrix)]
    beta_pm_matrix <- beta_pm_matrix[, -ncol(beta_pm_matrix)]
  }
  rownames(beta_pip_matrix) <- 1:nrow(beta_pip_matrix)
  colnames(beta_pip_matrix) <- marker_names
  beta_pip_matrix <- beta_pip_matrix[, reorder_markers]
  beta_pip_df <- as.data.frame(beta_pip_matrix)
  beta_pip_df$Factor <- paste0("Factor ", 1:nrow(beta_pip_df))
  if (!is.null(reorder_factors)){
    beta_pip_df <- beta_pip_df[reorder_factors, ]
  }
  beta_pip_plot_df <- reshape2::melt(beta_pip_df, value.name = "PIP")
  
  rownames(beta_pm_matrix) <- 1:nrow(beta_pm_matrix)
  colnames(beta_pm_matrix) <- marker_names
  beta_pm_matrix <- beta_pm_matrix[, reorder_markers]
  beta_pm_df <- as.data.frame(beta_pm_matrix)
  beta_pm_df$Factor <- paste0("Factor ", 1:nrow(beta_pm_df))
  if (!is.null(reorder_factors)){
    beta_pm_df <- beta_pm_df[reorder_factors, ]
  }
  beta_pm_plot_df <- reshape2::melt(beta_pm_df, id.var = "Factor",
                                    variable.name = "Perturbation",
                                    value.name = "Estimated effect size")
  # beta_pm_plot_df$PIP <- beta_pip_plot_df$PIP
  beta_pm_plot_df <- beta_pm_plot_df %>%
    mutate(PIP = beta_pip_plot_df$PIP,
           Perturbation = factor(Perturbation, levels = reorder_markers))
  if (inverse_factors){
    beta_pm_plot_df <- beta_pm_plot_df %>%
      mutate(Factor = factor(Factor, levels = beta_pip_df$Factor[nrow(beta_pip_df):1]))
  } else {
    beta_pm_plot_df <- beta_pm_plot_df %>%
      mutate(Factor = factor(Factor, levels = beta_pip_df$Factor[1:nrow(beta_pip_df)]))
  }
  beta_pm_plot_df$`Estimated effect size`[beta_pm_plot_df$`Estimated effect size` > 0.6] <- 0.6
  beta_pm_plot_df$`Estimated effect size`[beta_pm_plot_df$`Estimated effect size` < -0.6] <- -0.6
  plot_out <- ggplot(beta_pm_plot_df) +
    geom_point(aes(x = Perturbation, y = Factor,
                   size = PIP, color = `Estimated effect size`)) +
    scale_color_gradientn(limits = c(-0.6, 0.6),
                          colours = c("purple3", "purple2", "grey90", "darkorange", "darkorange1"),
                          breaks = seq(-0.6, 0.6, 0.3)) +
    # scale_color_gradient2(low = "purple3", mid = "grey90", high = "darkorange1") +
    # scale_color_gradientn(colors = c("purple", "grey90", "darkorange1"),
    #                       values = scales::rescale(c(-0.6, 0, 0.6))) +
    theme_void() +
    theme(axis.text.x = element_text(size = 13, angle = 90, hjust = 1),
          axis.text.y = element_text(size = 13, hjust = 1),
          legend.title = element_text(size = 13),
          legend.text = element_text(size = 12),
          plot.margin = margin(10, 10, 10, 12)) +
    labs(color = color_lgd_title)
  return(plot_out)
}

dotplot_effectsize <- function(effect_matrix, lfsr_matrix,
                               reorder_markers = colnames(effect_matrix),
                               color_lgd_title = "Estimated effect size",
                               size_lgd_title = "LFSR",
                               max_score = 0.2,
                               min_score = -0.2,
                               by_score = 0.1){
  # Both inputs should be gene by guide/marker matrices,
  # Dots will be colored by effect size and sized by 1-LFSR value.
  lfsr_binning <- function(lfsr){
    if (lfsr <= 0.05){
      return("0 - 0.05")
    } else if (lfsr <= 0.25){
      return("0.05 - 0.25")
    } else {
      return("> 0.25")
    }
  }
  
  lfsr_matrix <- lfsr_matrix[, reorder_markers]
  effect_matrix <- effect_matrix[, reorder_markers]
  
  pip_df <- as.data.frame(lfsr_matrix)
  pip_df$gene <- rownames(lfsr_matrix)
  pip_plot_df <- reshape2::melt(pip_df, variable.name = "Perturbation",
                                value.name = "LFSR")
  
  effect_df <- as.data.frame(effect_matrix)
  effect_df$gene <- rownames(effect_matrix)
  effect_plot_df <- reshape2::melt(effect_df, variable.name = "Perturbation",
                                   value.name = "Effect_size")
  
  combined_plot_df <- effect_plot_df %>%
    mutate(LFSR = pip_plot_df$LFSR,
           gene = factor(gene, levels = rownames(effect_matrix)),
           Perturbation = factor(Perturbation, levels = reorder_markers))
  combined_plot_df$Effect_size[combined_plot_df$Effect_size> max_score] <- max_score
  combined_plot_df$Effect_size[combined_plot_df$Effect_size< min_score] <- min_score
  combined_plot_df <- combined_plot_df %>%
    rowwise() %>%
    mutate(LFSR_bin = lfsr_binning(LFSR)) %>%
    mutate(LFSR_bin = factor(LFSR_bin, levels = c("> 0.25", "0.05 - 0.25", "0 - 0.05")))
  plot_out <- ggplot(combined_plot_df) +
    geom_point(aes(x = Perturbation, y = gene,
                   size = LFSR_bin, color = Effect_size)) +
    scale_color_gradientn(limits = c(min_score, max_score),
                          colours = c("blue3", "blue", "grey90", "red", "red3"),
                          breaks = seq(min_score, max_score, by_score)) +
    theme_void() +
    theme(axis.text.x = element_text(size = 13, angle = 90, hjust = 1),
          axis.text.y = element_text(size = 13),
          legend.title = element_text(size = 13),
          legend.text = element_text(size = 12)) +
    labs(color = color_lgd_title, size = size_lgd_title)
  return(plot_out)
}

complexplot_perturbation_factor <- function(gamma_pm, beta_pm, 
                                            marker_names,
                                            reorder_markers = marker_names,
                                            reorder_factors = NULL){
  rownames(gamma_pm) <- marker_names
  colnames(gamma_pm) <- paste0("Factor ", 1:ncol(gamma_pm))
  rownames(beta_pm) <- marker_names
  colnames(beta_pm) <- paste0("Factor ", 1:ncol(beta_pm))
  
  if (!is.null(reorder_factors)){
    pip_mat <- gamma_pm[reorder_markers, reorder_factors]
    effect_size_mat <- beta_pm[reorder_markers, reorder_factors]
  }
 
  lgd_list <- list()
  col_fun <- circlize::colorRamp2(breaks = seq(-0.3, 0.3, 0.3),
                                  colors = c("purple3", "grey90", "darkorange1"))
  lgd_list[["effectsize"]] <- Legend(title = "Association\neffect size",
                                     title_gp = gpar(fontsize = 13, fontface = "bold"),
                                     at = seq(-0.3, 0.3, 0.3),
                                     col_fun = col_fun,
                                     labels_gp = gpar(fontsize = 12))
  
  pip_tic_values <- seq(0.25, 1, 0.25)
  pip_tic_labels <- c("0.25", "0.50", "0.75", "1.00")
  
  lgd_list[["pip"]] <- 
    Legend(title = "PIP",
           title_gp = gpar(fontsize = 13, fontface = "bold"),
           labels = pip_tic_labels,
           labels_gp = gpar(fontsize = 12),
           grid_height = unit(6.5, "mm"),
           grid_width = unit(6, "mm"),
           graphics = list(
             function(x, y, w, h) grid.circle(x, y, r = (pip_tic_values[1] + 0.2) * unit(2, "mm"),
                                              gp = gpar(fill = "black")),
             function(x, y, w, h) grid.circle(x, y, r = (pip_tic_values[2] + 0.2) * unit(2, "mm"),
                                              gp = gpar(fill = "black")),
             function(x, y, w, h) grid.circle(x, y, r = (pip_tic_values[3] + 0.2) * unit(2, "mm"),
                                              gp = gpar(fill = "black")),
             function(x, y, w, h) grid.circle(x, y, r = (pip_tic_values[4] + 0.2) * unit(2, "mm"),
                                              gp = gpar(fill = "black"))
           ))
  
  map1 <- Heatmap(effect_size_mat,
                  name = "Association Effect Size",
                  col = col_fun,
                  rect_gp = gpar(type = "none"),
                  cell_fun = function(j, i, x, y, width, height, fill) {
                    grid.rect(x = x, y = y, width = width, height = height, 
                              gp = gpar(col = NA, fill = NA))
                    grid.circle(x = x, y = y,
                                r = (pip_mat[i, j] + 0.2) * unit(2, "mm"),
                                gp = gpar(fill = col_fun(effect_size_mat[i, j]), col = NA))
                  },
                  border_gp = gpar(col = "black"),
                  row_title = "Perturbations",
                  row_title_gp = gpar(fontsize = 16),
                  column_title = "Factors",
                  column_title_gp = gpar(fontsize = 16),
                  cluster_rows = F, cluster_columns = F,
                  show_heatmap_legend = F,
                  row_names_gp = gpar(fontsize = 12),
                  column_names_gp = gpar(fontsize = 12),
                  column_names_rot = 45,
                  column_names_side = "top",
                  column_title_side = "bottom")
  draw(map1, annotation_legend_list = lgd_list)
}

complexplot_gene_factor <- function(genes_df, interest_df,
                                    F_pm, W_pm, reorder_factors = NULL){
  interest_df <- interest_df[interest_df$gene_name %in% genes_df$Name, ]
  interest_df$type <- factor(interest_df$type, levels = unique(interest_df$type))
  rownames(F_pm) <- genes_df$Name
  colnames(F_pm) <- paste0("Factor ", 1:ncol(F_pm))
  rownames(W_pm) <- genes_df$Name
  colnames(W_pm) <- paste0("Factor ", 1:ncol(W_pm))
  
  if (!is.null(reorder_factors)){
    F_pm <- F_pm[, reorder_factors]
    W_pm <- W_pm[, reorder_factors]
  }
  effect_size_mat <- W_pm[interest_df$gene_name, ]
  pip_mat <- F_pm[interest_df$gene_name, ]
  
  lgd_list <- list()
  col_fun <- circlize::colorRamp2(breaks = seq(-0.6, 0.6, 0.3),
                                  colors = c("purple3", "purple2", "grey90", "darkorange", "darkorange1"))
  lgd_list[["effectsize"]] <- Legend(title = "Gene loading",
                                     title_gp = gpar(fontsize = 13, fontface = "bold"),
                                     at = seq(-0.6, 0.6, 0.3),
                                     col_fun = col_fun,
                                     labels_gp = gpar(fontsize = 12))
  
  pip_tic_values <- seq(0.25, 1, 0.25)
  pip_tic_labels <- c("0.25", "0.50", "0.75", "1.00")
  
  lgd_list[["pip"]] <- 
    Legend(title = "PIP",
           title_gp = gpar(fontsize = 13, fontface = "bold"),
           labels = pip_tic_labels,
           labels_gp = gpar(fontsize = 12),
           grid_height = unit(6.5, "mm"),
           grid_width = unit(6, "mm"),
           graphics = list(
             function(x, y, w, h) grid.circle(x, y, r = (pip_tic_values[1] + 0.2) * unit(2, "mm"),
                                              gp = gpar(fill = "black")),
             function(x, y, w, h) grid.circle(x, y, r = (pip_tic_values[2] + 0.2) * unit(2, "mm"),
                                              gp = gpar(fill = "black")),
             function(x, y, w, h) grid.circle(x, y, r = (pip_tic_values[3] + 0.2) * unit(2, "mm"),
                                              gp = gpar(fill = "black")),
             function(x, y, w, h) grid.circle(x, y, r = (pip_tic_values[4] + 0.2) * unit(2, "mm"),
                                              gp = gpar(fill = "black"))
           ))
  
  marker_colormap <- structure(RColorBrewer::brewer.pal(length(levels(interest_df$type)), "Set3"),
                               names = levels(interest_df$type))
  lgd_list[["Marker"]] <-  Legend(title = "Marker annotation",
                                  title_gp = gpar(fontsize = 13, fontface = "bold"),
                                  labels = levels(interest_df$type),
                                  labels_gp = gpar(fontsize = 12),
                                  at = levels(interest_df$type),
                                  legend_gp = gpar(fill = marker_colormap))
  right_annot <- rowAnnotation(Marker = interest_df$type,
                               col = list(Marker = marker_colormap),
                               annotation_legend_param = list(
                                 Marker = list(
                                   title = "Marker annotation",
                                   at = levels(interest_df$type),
                                   labels = levels(interest_df$type)
                                 )
                               ),
                               show_annotation_name = F,
                               show_legend = F)
  
  map1 <- Heatmap(effect_size_mat,
                  name = "Gene loading",
                  col = col_fun,
                  rect_gp = gpar(type = "none"),
                  cell_fun = function(j, i, x, y, width, height, fill) {
                    grid.rect(x = x, y = y, width = width, height = height, 
                              gp = gpar(col = NA, fill = NA))
                    grid.circle(x = x, y = y,
                                r = (pip_mat[i, j] + 0.2) * unit(2, "mm"),
                                gp = gpar(fill = col_fun(effect_size_mat[i, j]), col = NA))
                  },
                  border_gp = gpar(col = "black"),
                  row_title = "Marker genes",
                  row_title_gp = gpar(fontsize = 16),
                  column_title = "Factors",
                  column_title_gp = gpar(fontsize = 16),
                  cluster_rows = F, cluster_columns = F,
                  right_annotation = right_annot,
                  show_heatmap_legend = F,
                  row_names_gp = gpar(fontsize = 12),
                  column_names_gp = gpar(fontsize = 12),
                  column_names_rot = 45,
                  column_names_side = "top",
                  column_title_side = "bottom")
  draw(map1, annotation_legend_list = lgd_list)
}

complexplot_gene_perturbation <- function(genes_df, interest_df,
                                          targets = NULL,
                                          lfsr_mat, lfsr_name = "LFSR",
                                          lfsr_cutoff = 0.05,
                                          effect_mat, effect_name = "GSFA\nsummarized effect",
                                          score_break = seq(-0.2, 0.2, 0.1),
                                          color_break = c("blue3", "blue", "grey90", "red", "red3")){
  interest_df <- interest_df[interest_df$gene_name %in% genes_df$Name, ]
  interest_df$type <- factor(interest_df$type, levels = unique(interest_df$type))
  if (all(startsWith(rownames(lfsr_mat), "ENSG"))){
    rownames(lfsr_mat) <- genes_df$Name[match(rownames(lfsr_mat), genes_df$ID)]
  }
  if (is.null(rownames(effect_mat))){
    rownames(effect_mat) <- genes_df$Name
  } else if (all(startsWith(rownames(effect_mat), "ENSG"))){
    rownames(effect_mat) <- genes_df$Name[match(rownames(effect_mat), genes_df$ID)]
  }
  if (is.null(colnames(effect_mat))){
    colnames(effect_mat) <- colnames(lfsr_mat)
  }
  if (is.null(targets)){
    num_signif_genes <- colSums(lfsr_mat < lfsr_cutoff)
    targets <- names(num_signif_genes)[which(num_signif_genes > 0)]
  }
  selected_effect_mat <- effect_mat[interest_df$gene_name, targets]
  
  selected_lfsr_mat <- lfsr_mat[interest_df$gene_name, targets]
  binned_size_mat <- matrix(0.6, 
                            nrow = nrow(selected_lfsr_mat),
                            ncol = ncol(selected_lfsr_mat))
  rownames(binned_size_mat) <- rownames(selected_lfsr_mat)
  colnames(binned_size_mat) <- colnames(selected_lfsr_mat)
  binned_size_mat[selected_lfsr_mat <= 0.05] <- 1
  binned_size_mat[selected_lfsr_mat > 0.25] <- 0.2
  
  lgd_list <- list()
  col_fun <- circlize::colorRamp2(breaks = score_break,
                                  colors = color_break)
  lgd_list[["effectsize"]] <- Legend(title = effect_name,
                                     at = score_break,
                                     col_fun = col_fun)
  
  lfsr_tic_values <- c(0.2, 0.6, 1)
  lfsr_tic_labels <- c("> 0.25", "0.05 - 0.25", "0 - 0.05")
  
  lgd_list[["LFSR"]] <- 
    Legend(title = lfsr_name,
           labels = lfsr_tic_labels,
           grid_height = unit(6, "mm"),
           grid_width = unit(6, "mm"),
           graphics = list(
             function(x, y, w, h) grid.circle(x, y, r = (lfsr_tic_values[1] + 0.2) * unit(2, "mm"),
                                              gp = gpar(fill = "black")),
             function(x, y, w, h) grid.circle(x, y, r = (lfsr_tic_values[2] + 0.2) * unit(2, "mm"),
                                              gp = gpar(fill = "black")),
             function(x, y, w, h) grid.circle(x, y, r = (lfsr_tic_values[3] + 0.2) * unit(2, "mm"),
                                              gp = gpar(fill = "black"))
           ))
  
  marker_colormap <- structure(RColorBrewer::brewer.pal(length(levels(interest_df$type)), "Set3"),
                               names = levels(interest_df$type))
  lgd_list[["Marker"]] <-  Legend(title = "Marker annotation",
                                  labels = levels(interest_df$type),
                                  at = levels(interest_df$type),
                                  legend_gp = gpar(fill = marker_colormap))
  right_annot <- rowAnnotation(Marker = interest_df$type,
                               col = list(Marker = marker_colormap),
                               annotation_legend_param = list(
                                 Marker = list(
                                   title = "Marker annotation",
                                   at = levels(interest_df$type),
                                   labels = levels(interest_df$type)
                                 )
                               ),
                               show_annotation_name = F,
                               show_legend = F)
  
  map1 <- Heatmap(selected_effect_mat,
                  name = effect_name,
                  col = col_fun,
                  rect_gp = gpar(type = "none"),
                  cell_fun = function(j, i, x, y, width, height, fill) {
                    grid.rect(x = x, y = y, width = width, height = height, 
                              gp = gpar(col = NA, fill = NA))
                    grid.circle(x = x, y = y,
                                r = (binned_size_mat[i, j] + 0.2) * unit(2, "mm"),
                                gp = gpar(fill = col_fun(selected_effect_mat[i, j]), col = NA))
                  },
                  border_gp = gpar(col = "black"),
                  row_title = "Marker genes",
                  column_title = "Perturbations",
                  cluster_rows = F, cluster_columns = F,
                  right_annotation = right_annot,
                  show_heatmap_legend = F,
                  row_names_gp = gpar(fontsize = 10.5),
                  column_names_rot = 45,
                  column_names_side = "top",
                  column_title_side = "bottom")
  draw(map1, annotation_legend_list = lgd_list)
}

plot_pairwise.corr_heatmap <- function(input_mat_1, input_mat_2 = NULL,
                                       name_1 = NULL, name_2 = NULL,
                                       corr_type = c("pearson", "jaccard", "prop_overlap"),
                                       return_corr = FALSE,
                                       label_size = 8,
                                       color_vec = NULL){
  # Please store samples in the columns of 'input_mat'.
  if (is.null(input_mat_2)){
    input_mat_2 <- input_mat_1
  }
  if (is.null(name_2)){
    name_2 <- name_1
  }
  stopifnot(nrow(input_mat_1) == nrow(input_mat_2))
  corr_mat <- matrix(nrow = ncol(input_mat_1), ncol = ncol(input_mat_2))
  
  for (i in 1:ncol(input_mat_1)){
    for (j in 1:ncol(input_mat_2)){
      vec_1 <- input_mat_1[, i]
      vec_2 <- input_mat_2[, j]
      if (corr_type == "pearson"){
        corr_mat[i, j] <- cor(vec_1, vec_2, method = "pearson")
      } else if (corr_type == "jaccard"){
        corr_mat[i, j] <- sum(vec_1 * vec_2) / sum((vec_1 + vec_2) > 0)
      } else {
        corr_mat[i, j] <- sum(vec_1 * vec_2) / sum(vec_1)
      }
    }
  }
  
  if (is.null(colnames(input_mat_1))){
    rownames(corr_mat) <- 1:ncol(input_mat_1)
  } else {
    rownames(corr_mat) <- colnames(input_mat_1)
  }
  if (is.null(colnames(input_mat_2))){
    colnames(corr_mat) <- 1:ncol(input_mat_2)
  } else {
    colnames(corr_mat) <- colnames(input_mat_2)
  }
  
  if (corr_type == "pearson"){
    legend_name <- "Pearson Correlation"
    if (is.null(color_vec)){
      colormap <- circlize::colorRamp2(breaks = c(-1, 0, 1),
                                       colors = c("blue", "white", "red"))
    } else {
      colormap <- circlize::colorRamp2(breaks = c(-1, 0, 1),
                                       colors = color_vec)
    }
    
  }
  if (corr_type == "jaccard" | corr_type == "prop_overlap"){
    legend_name <- ifelse(corr_type == "jaccard",
                          "Jaccard Index", "% of Shared Non-Zero Genes")
    if (is.null(color_vec)){
      colormap <- circlize::colorRamp2(breaks = c(0, 0.3, 1),
                                       colors = c("black", "purple", "gold"))
    } else {
      colormap <- circlize::colorRamp2(breaks = c(0, 0.5, 1),
                                       colors = color_vec)
    }
  }
  
  ht <- Heatmap(corr_mat,
                col = colormap,
                name = legend_name,
                row_title = name_1, column_title = name_2,
                cluster_rows = F, cluster_columns = F,
                row_names_gp = gpar(fontsize = label_size),
                column_names_gp = gpar(fontsize = label_size))
  draw(ht)
  if (return_corr){
    return(corr_mat)
  }
}

source("/project2/xinhe/yifan/GTEx/scripts/qqplot_uniform.R")
summ_pvalues <- function(pvalues, title_text = NULL){
  requireNamespace("gridExtra", quietly = TRUE)
  # distribution histogram
  plot1 <- histogram(pvalues, col = 'grey', type = "count",
                     xlim = c(0, 1), breaks = 50,
                     main = "p-value distribution", xlab = "p-value", ylab = "Count")
  # uniform qq-plot
  plot2 <- qqunif.plot(pvalues, main = "p-value Q-Q plot")
  gridExtra::grid.arrange(plot1, plot2, ncol = 2, top = title_text)
}

barplot_top_enrich_terms <- function(enrich_df, fdr_cutoff = 0.05, FC_cutoff = 2,
                                     size_cutoff = 20,
                                     terms_of_interest = NULL,
                                     top_num = 5,
                                     str_wrap_length = 25,
                                     FC_max = 8,
                                     pval_max = 10){
  enrich_df <- enrich_df %>%
    filter(FDR < fdr_cutoff, enrichmentRatio >= FC_cutoff) %>%
    filter(size >= size_cutoff) %>%
    arrange(-enrichmentRatio)
  if (is.null(terms_of_interest)){
    top_terms_df <- enrich_df %>% slice(1:top_num)
  } else {
    top_terms_df <- enrich_df %>% filter(description %in% terms_of_interest)
  }
  top_terms_df$description_short <- str_wrap(top_terms_df$description, width = str_wrap_length)
  top_terms_df$description_short <- factor(top_terms_df$description_short, 
                                     levels = top_terms_df$description_short)
  top_terms_df$pValue[top_terms_df$pValue < 10^(-pval_max)] <- 10^(-pval_max)
  top_terms_df$enrichmentRatio[top_terms_df$enrichmentRatio > FC_max] <- FC_max
  plot_out <- ggplot(top_terms_df) +
    geom_bar(aes(x = description_short, y = enrichmentRatio, fill = -log10(pValue)),
             stat="identity") +
    coord_flip() +
    scale_y_continuous(limits = c(0, FC_max)) +
    scale_fill_gradientn(limits = c(2, pval_max),
                         colors = c("blue", "red"),
                         breaks = seq(2, pval_max, 2)) +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_text(size = 13))
  return(plot_out)
}
