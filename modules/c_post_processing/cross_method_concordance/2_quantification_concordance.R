#!/usr/bin/env Rscript

# ===============================================
# CROSS-METHOD CONCORDANCE - QUANTIFICATION CONCORDANCE
# ===============================================
# Computes pairwise Spearman and Pearson correlations between methods,
# identifies discordant genes, and generates concordance heatmaps.
#
# Outputs:
#   - tables/pairwise_spearman_per_sample.csv
#   - tables/pairwise_pearson_per_sample.csv
#   - tables/median_correlation_matrix_spearman.csv
#   - tables/median_correlation_matrix_pearson.csv
#   - tables/discordant_genes.csv
#   - figures/method_concordance_heatmap_spearman.png
#   - figures/method_concordance_heatmap_pearson.png
#   - figures/per_sample_correlation_boxplot.png
#   - figures/discordant_genes_heatmap.png

source(file.path(Sys.getenv("CONCORDANCE_SCRIPT_DIR", "."), "0_concordance_config.R"))

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

cat("\n=== STEP 2: Quantification Concordance ===\n\n")

# Load harmonized data
if (!file.exists(HARMONIZED_RDS)) {
  stop("Harmonized data not found. Run 1_load_matrices.R first.")
}
data <- readRDS(HARMONIZED_RDS)
tpm_matrices <- data$tpm_matrices
common_genes <- data$common_genes
common_samples <- data$common_samples
methods <- names(tpm_matrices)
n_methods <- length(methods)

cat("Methods:", paste(sapply(methods, get_short_name), collapse = ", "), "\n")
cat("Genes:", length(common_genes), "| Samples:", length(common_samples), "\n\n")

# -----------------------------------------------
# 2.1 Pairwise correlations (per sample)
# -----------------------------------------------
# For each sample, correlate the gene expression vectors between each pair of methods.
# Use log2(TPM+1) to reduce skewness from highly expressed genes.

cat("--- Computing pairwise correlations per sample ---\n")

method_pairs <- combn(methods, 2, simplify = FALSE)
pair_names <- sapply(method_pairs, function(p) paste(get_short_name(p[1]), "vs", get_short_name(p[2])))

spearman_per_sample <- matrix(NA, nrow = length(common_samples), ncol = length(method_pairs),
                               dimnames = list(common_samples, pair_names))
pearson_per_sample <- matrix(NA, nrow = length(common_samples), ncol = length(method_pairs),
                              dimnames = list(common_samples, pair_names))

# Pre-compute log2(TPM+1) once per method (avoids redundant log2 per pair × sample)
log2_matrices <- lapply(tpm_matrices, function(mat) log2(mat + 1))

for (i in seq_along(method_pairs)) {
  m1 <- method_pairs[[i]][1]
  m2 <- method_pairs[[i]][2]
  mat1 <- log2_matrices[[m1]]
  mat2 <- log2_matrices[[m2]]

  # Vectorized: identify genes with nonzero expression in at least one method
  nonzero_mask <- (mat1 > 0) | (mat2 > 0)  # genes × samples logical matrix
  nonzero_per_sample <- colSums(nonzero_mask)

  for (s in common_samples) {
    if (nonzero_per_sample[s] < 10) next
    nz <- nonzero_mask[, s]
    spearman_per_sample[s, i] <- cor(mat1[nz, s], mat2[nz, s], method = "spearman",
                                      use = "pairwise.complete.obs")
    pearson_per_sample[s, i] <- cor(mat1[nz, s], mat2[nz, s], method = "pearson",
                                     use = "pairwise.complete.obs")
  }
}

# Save per-sample correlations
write.csv(data.frame(Sample = rownames(spearman_per_sample), spearman_per_sample, check.names = FALSE),
          file.path(TABLES_DIR, "pairwise_spearman_per_sample.csv"), row.names = FALSE)
write.csv(data.frame(Sample = rownames(pearson_per_sample), pearson_per_sample, check.names = FALSE),
          file.path(TABLES_DIR, "pairwise_pearson_per_sample.csv"), row.names = FALSE)

cat("  Saved per-sample correlation tables\n")

# -----------------------------------------------
# 2.2 Median correlation matrices (methods x methods)
# -----------------------------------------------

cat("--- Computing median correlation matrices ---\n")

short_names <- sapply(methods, get_short_name)

build_median_cor_matrix <- function(per_sample_mat) {
  cor_mat <- matrix(1, nrow = n_methods, ncol = n_methods,
                     dimnames = list(short_names, short_names))
  for (i in seq_along(method_pairs)) {
    m1_idx <- which(methods == method_pairs[[i]][1])
    m2_idx <- which(methods == method_pairs[[i]][2])
    med_cor <- median(per_sample_mat[, i], na.rm = TRUE)
    cor_mat[m1_idx, m2_idx] <- med_cor
    cor_mat[m2_idx, m1_idx] <- med_cor
  }
  return(cor_mat)
}

median_spearman <- build_median_cor_matrix(spearman_per_sample)
median_pearson <- build_median_cor_matrix(pearson_per_sample)

write.csv(data.frame(Method = rownames(median_spearman), median_spearman, check.names = FALSE),
          file.path(TABLES_DIR, "median_correlation_matrix_spearman.csv"), row.names = FALSE)
write.csv(data.frame(Method = rownames(median_pearson), median_pearson, check.names = FALSE),
          file.path(TABLES_DIR, "median_correlation_matrix_pearson.csv"), row.names = FALSE)

cat("  Spearman median correlations:\n")
print(round(median_spearman, 3))
cat("\n  Pearson median correlations:\n")
print(round(median_pearson, 3))

# -----------------------------------------------
# 2.3 Method concordance heatmaps
# -----------------------------------------------

cat("\n--- Generating concordance heatmaps ---\n")

plot_concordance_heatmap <- function(cor_mat, cor_type, output_file) {
  col_fun <- colorRamp2(
    seq(min(cor_mat, na.rm = TRUE), 1, length.out = 100),
    colorRampPalette(c("#FFFFFF", "#D1E5F0", "#67A9CF", "#2166AC"))(100)
  )

  # Annotation with cell values — use 2 decimal places to avoid overlap
  cell_fun <- function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.2f", cor_mat[i, j]), x, y, gp = gpar(fontsize = 12, fontface = "bold"))
  }

  ht <- Heatmap(cor_mat,
    name = paste("Median", cor_type),
    col = col_fun,
    cell_fun = cell_fun,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_dend = TRUE,
    show_column_dend = TRUE,
    row_names_gp = gpar(fontsize = 13),
    column_names_gp = gpar(fontsize = 13),
    column_names_rot = 45,
    column_title = paste0("Cross-Method Concordance (Median ", cor_type, " Correlation)"),
    column_title_gp = gpar(fontsize = 16, fontface = "bold"),
    heatmap_legend_param = list(
      title = paste("Median\n", cor_type),
      legend_height = unit(5, "cm")
    ),
    width = unit(14, "cm"),
    height = unit(14, "cm")
  )

  png(output_file, width = 1600, height = 1300, res = 150)
  draw(ht, padding = unit(c(30, 30, 25, 40), "mm"))
  dev.off()
  cat("  Saved:", output_file, "\n")
}

plot_concordance_heatmap(median_spearman, "Spearman",
                         file.path(FIGURES_DIR, "method_concordance_heatmap_spearman.png"))
plot_concordance_heatmap(median_pearson, "Pearson",
                         file.path(FIGURES_DIR, "method_concordance_heatmap_pearson.png"))

# -----------------------------------------------
# 2.4 Per-sample correlation boxplot
# -----------------------------------------------

cat("--- Generating per-sample correlation boxplot ---\n")

png(file.path(FIGURES_DIR, "per_sample_correlation_boxplot.png"),
    width = max(1800, 220 * ncol(spearman_per_sample)), height = 1200, res = 150)

par(mar = c(18, 7, 6, 4))
boxplot(spearman_per_sample,
        las = 2, cex.axis = 1.0,
        main = "Spearman Correlation per Sample (All Method Pairs)",
        cex.main = 1.4,
        ylab = "Spearman rho", cex.lab = 1.2,
        col = colorRampPalette(c("#CE93D8", "#4A148C"))(ncol(spearman_per_sample)),
        outline = TRUE, notch = FALSE)
abline(h = 0.9, lty = 2, col = "red")
abline(h = 0.95, lty = 3, col = "darkgreen")
legend("bottomleft", legend = c("rho=0.90", "rho=0.95"), lty = c(2, 3),
       col = c("red", "darkgreen"), cex = 1.0, bg = "white")

dev.off()
cat("  Saved: per_sample_correlation_boxplot.png\n")

# -----------------------------------------------
# 2.5 Identify discordant genes
# -----------------------------------------------
# A gene is "discordant" if its expression varies widely across methods.
# Metric: for each gene, compute the coefficient of variation (CV) of its
# mean expression (across samples) estimated by each method.

cat("\n--- Identifying discordant genes ---\n")

# Compute mean log2(TPM+1) per gene per method — reuse pre-computed log2 matrices
mean_expr_per_method <- sapply(log2_matrices, rowMeans)
colnames(mean_expr_per_method) <- sapply(methods, get_short_name)

# CV across methods for each gene (vectorized SD using rowMeans/rowSums)
gene_means <- rowMeans(mean_expr_per_method)
n_m <- ncol(mean_expr_per_method)
gene_sds <- sqrt(rowSums((mean_expr_per_method - gene_means)^2) / (n_m - 1))
gene_cv <- ifelse(gene_means > 0, gene_sds / gene_means, 0)

# Flag genes with high CV
discordant_mask <- gene_cv > DISCORDANCE_CV_THRESHOLD & gene_means > log2(CONCORDANCE_MIN_EXPR + 1)
n_discordant <- sum(discordant_mask)
cat("  Discordant genes (CV >", DISCORDANCE_CV_THRESHOLD, "):", n_discordant, "\n")

# Build discordant gene table
discordant_df <- data.frame(
  Gene_ID = names(gene_cv)[discordant_mask],
  CV_across_methods = round(gene_cv[discordant_mask], 4),
  Mean_log2TPM = round(gene_means[discordant_mask], 4),
  mean_expr_per_method[discordant_mask, , drop = FALSE],
  check.names = FALSE,
  stringsAsFactors = FALSE
)
discordant_df <- discordant_df[order(-discordant_df$CV_across_methods), ]

write.csv(discordant_df, file.path(TABLES_DIR, "discordant_genes.csv"), row.names = FALSE)
cat("  Saved: discordant_genes.csv (", nrow(discordant_df), "genes )\n")

# Also save concordant summary (all genes with their CV)
all_genes_concordance <- data.frame(
  Gene_ID = names(gene_cv),
  CV_across_methods = round(gene_cv, 4),
  Mean_log2TPM = round(gene_means, 4),
  mean_expr_per_method,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
all_genes_concordance <- all_genes_concordance[order(-all_genes_concordance$CV_across_methods), ]
write.csv(all_genes_concordance, file.path(TABLES_DIR, "all_genes_concordance.csv"), row.names = FALSE)

# -----------------------------------------------
# 2.6 Discordant genes heatmap (top 50)
# -----------------------------------------------

if (n_discordant > 0) {
  cat("--- Generating discordant genes heatmap ---\n")

  top_n <- min(50, n_discordant)
  top_discordant <- head(discordant_df$Gene_ID, top_n)

  disc_mat <- mean_expr_per_method[top_discordant, , drop = FALSE]

  col_fun <- colorRamp2(
    seq(0, max(disc_mat, na.rm = TRUE), length.out = 100),
    colorRampPalette(c("#F7F7F7", "#CE93D8", "#7B1FA2", "#4A148C"))(100)
  )

  ht <- Heatmap(disc_mat,
    name = "Mean\nlog2(TPM+1)",
    col = col_fun,
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    show_row_names = top_n <= 50,
    row_names_gp = gpar(fontsize = if (top_n > 30) 8 else 10),
    row_names_max_width = unit(8, "cm"),
    column_names_gp = gpar(fontsize = 13),
    column_names_rot = 45,
    column_title = paste0("Top ", top_n, " Discordant Genes Across Methods"),
    column_title_gp = gpar(fontsize = 15, fontface = "bold"),
    right_annotation = rowAnnotation(
      CV = anno_barplot(gene_cv[top_discordant],
                        gp = gpar(fill = "#EF8A62"),
                        width = unit(3, "cm")),
      annotation_name_gp = gpar(fontsize = 10)
    ),
    heatmap_legend_param = list(
      title = "Mean\nlog2(TPM+1)",
      legend_height = unit(5, "cm")
    )
  )

  fig_height <- max(1000, 200 + top_n * 28)
  png(file.path(FIGURES_DIR, "discordant_genes_heatmap.png"),
      width = 1600, height = fig_height, res = 150)
  draw(ht, padding = unit(c(30, 25, 25, 40), "mm"))
  dev.off()
  cat("  Saved: discordant_genes_heatmap.png\n")

  # -----------------------------------------------
  # 2.6b Z-score scaled (0–10) heatmap for discordant genes
  # -----------------------------------------------
  # Normalizes each gene across methods so all genes share a common 0–10 scale,
  # making it easy to compare method agreement regardless of absolute expression.

  cat("--- Generating Z-score scaled heatmap for discordant genes ---\n")

  zscore_mat <- t(scale(t(disc_mat)))
  zscore_mat[is.nan(zscore_mat)] <- 0  # zero-SD genes (identical across methods)
  row_mins <- apply(zscore_mat, 1, min, na.rm = TRUE)
  row_maxs <- apply(zscore_mat, 1, max, na.rm = TRUE)
  row_range <- row_maxs - row_mins
  zscore_scaled <- (zscore_mat - row_mins) / ifelse(row_range == 0, 1, row_range) * 10
  # Zero-SD genes get flat 5.0 (midpoint)
  zscore_scaled[row_range == 0, ] <- 5.0

  col_fun_z <- colorRamp2(
    c(0, 5, 10),
    c("#2166AC", "#F7F7F7", "#B2182B")
  )

  cell_fun_z <- function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.1f", zscore_scaled[i, j]), x, y,
              gp = gpar(fontsize = if (top_n > 30) 7 else 9, fontface = "bold"))
  }

  ht_z <- Heatmap(zscore_scaled,
    name = "Z-Score\n(0-10)",
    col = col_fun_z,
    cell_fun = cell_fun_z,
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    show_row_names = top_n <= 50,
    row_names_gp = gpar(fontsize = if (top_n > 30) 8 else 10),
    row_names_max_width = unit(8, "cm"),
    column_names_gp = gpar(fontsize = 13),
    column_names_rot = 45,
    column_title = paste0("Top ", top_n, " Discordant Genes \u2014 Z-Score Scaled (0\u201310)"),
    column_title_gp = gpar(fontsize = 15, fontface = "bold"),
    right_annotation = rowAnnotation(
      CV = anno_barplot(gene_cv[top_discordant],
                        gp = gpar(fill = "#EF8A62"),
                        width = unit(3, "cm")),
      annotation_name_gp = gpar(fontsize = 10)
    ),
    heatmap_legend_param = list(
      title = "Z-Score\n(0-10)",
      at = c(0, 2.5, 5, 7.5, 10),
      labels = c("0 (low)", "2.5", "5 (avg)", "7.5", "10 (high)"),
      legend_height = unit(5, "cm")
    )
  )

  png(file.path(FIGURES_DIR, "discordant_genes_zscore_heatmap.png"),
      width = 1600, height = fig_height, res = 150)
  draw(ht_z, padding = unit(c(30, 25, 25, 40), "mm"))
  dev.off()
  cat("  Saved: discordant_genes_zscore_heatmap.png\n")

  # Save Z-score table
  zscore_df <- data.frame(
    Gene_ID = rownames(zscore_scaled),
    CV_across_methods = round(gene_cv[top_discordant], 4),
    zscore_scaled,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  write.csv(zscore_df, file.path(TABLES_DIR, "discordant_genes_zscore.csv"), row.names = FALSE)
  cat("  Saved: discordant_genes_zscore.csv\n")
}

# -----------------------------------------------
# 2.7 Scatter plots: pairwise method comparison (top expressed genes)
# -----------------------------------------------

cat("--- Generating pairwise scatter plots ---\n")

# Create a multi-panel scatter for all method pairs
n_pairs <- length(method_pairs)
ncols_plot <- min(n_pairs, 5)
nrows_plot <- ceiling(n_pairs / ncols_plot)

png(file.path(FIGURES_DIR, "pairwise_scatter_plots.png"),
    width = 600 * ncols_plot, height = 600 * nrows_plot, res = 150)
par(mfrow = c(nrows_plot, ncols_plot), mar = c(6, 6, 5, 3), oma = c(2, 2, 2, 2), cex = 0.9)

for (i in seq_along(method_pairs)) {
  m1 <- method_pairs[[i]][1]
  m2 <- method_pairs[[i]][2]

  x <- rowMeans(log2_matrices[[m1]])
  y <- rowMeans(log2_matrices[[m2]])

  # Color by expression level
  cols <- ifelse(x + y > 0, adjustcolor("#7B1FA2", alpha.f = 0.15), "grey80")

  plot(x, y,
       pch = 16, cex = 0.4, col = cols,
       xlab = get_short_name(m1),
       ylab = get_short_name(m2),
       cex.lab = 1.1, cex.main = 1.1,
       main = paste0("r=", round(cor(x, y, method = "pearson"), 3),
                      " | rho=", round(cor(x, y, method = "spearman"), 3)))
  abline(0, 1, col = "red", lty = 2)
}

dev.off()
cat("  Saved: pairwise_scatter_plots.png\n")

# -----------------------------------------------
# Save concordance results for report generation
# -----------------------------------------------

concordance_results <- list(
  median_spearman = median_spearman,
  median_pearson = median_pearson,
  spearman_per_sample = spearman_per_sample,
  pearson_per_sample = pearson_per_sample,
  discordant_df = discordant_df,
  gene_cv = gene_cv,
  mean_expr_per_method = mean_expr_per_method,
  n_discordant = n_discordant
)

saveRDS(concordance_results, file.path(OUTPUT_DIR, "concordance_results.rds"))
cat("\n[DONE] Concordance results saved\n")
