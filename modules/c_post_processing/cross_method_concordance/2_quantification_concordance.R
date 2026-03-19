#!/usr/bin/env Rscript

# ===============================================
# CROSS-METHOD CONCORDANCE - QUANTIFICATION CONCORDANCE
# ===============================================
# Computes pairwise Spearman correlations between methods and generates
# a median Spearman concordance heatmap.
#
# Outputs:
#   - tables/pairwise_spearman_per_sample.csv
#   - tables/median_correlation_matrix_spearman.csv
#   - figures/method_concordance_heatmap_spearman.png

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

if (n_methods < 2) {
  stop("Need at least 2 methods for pairwise concordance. Found: ", n_methods)
}

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

  # Vectorized skip detection — guard against missing names
  nonzero_counts <- nonzero_per_sample[common_samples]
  nonzero_counts[is.na(nonzero_counts)] <- 0L
  skipped_samples <- common_samples[nonzero_counts < CORRELATION_MIN_GENES]
  valid_samples <- setdiff(common_samples, skipped_samples)

  # Vectorized Spearman: rank transform per column, then compute Pearson on ranks.
  # Spearman(x,y) = Pearson(rank(x), rank(y)). This avoids N individual cor() calls.
  if (length(valid_samples) > 0) {
    # Zero out non-expressed genes per sample so they don't affect ranking
    # (set to NA, then rank with na.last="keep" to exclude them)
    m1_valid <- mat1[, valid_samples, drop = FALSE]
    m2_valid <- mat2[, valid_samples, drop = FALSE]
    nz_valid <- nonzero_mask[, valid_samples, drop = FALSE]
    m1_valid[!nz_valid] <- NA
    m2_valid[!nz_valid] <- NA

    # Column-wise rank transform (each sample ranked independently)
    r1 <- apply(m1_valid, 2, rank, na.last = "keep")
    r2 <- apply(m2_valid, 2, rank, na.last = "keep")

    # Pearson correlation on ranks = Spearman (vectorized per column)
    # Center each column, compute dot-product correlation
    r1_centered <- sweep(r1, 2, colMeans(r1, na.rm = TRUE))
    r2_centered <- sweep(r2, 2, colMeans(r2, na.rm = TRUE))
    r1_centered[is.na(r1_centered)] <- 0
    r2_centered[is.na(r2_centered)] <- 0

    num <- colSums(r1_centered * r2_centered)
    den <- sqrt(colSums(r1_centered^2) * colSums(r2_centered^2))
    den[den == 0] <- 1  # guard against zero-variance
    spearman_per_sample[valid_samples, i] <- num / den
  }

  if (length(skipped_samples) > 0) {
    cat("  Warning: Skipped", length(skipped_samples), "samples for pair",
        get_short_name(m1), "vs", get_short_name(m2),
        "(fewer than", CORRELATION_MIN_GENES, "expressed genes):",
        paste(head(skipped_samples, 5), collapse = ", "),
        if (length(skipped_samples) > 5) "..." else "", "\n")
  }
}

# Save per-sample correlations
write.csv(data.frame(Sample = rownames(spearman_per_sample), spearman_per_sample, check.names = FALSE),
          file.path(TABLES_DIR, "pairwise_spearman_per_sample.csv"), row.names = FALSE)

cat("  Saved per-sample Spearman correlation table\n")

# -----------------------------------------------
# 2.2 Median correlation matrices (methods x methods)
# -----------------------------------------------

cat("--- Computing median Spearman correlation matrix ---\n")

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

write.csv(data.frame(Method = rownames(median_spearman), median_spearman, check.names = FALSE),
          file.path(TABLES_DIR, "median_correlation_matrix_spearman.csv"), row.names = FALSE)

cat("  Spearman median correlations:\n")
print(round(median_spearman, 3))

# -----------------------------------------------
# 2.3 Median Spearman concordance heatmap
# -----------------------------------------------

cat("\n--- Generating Median Spearman concordance heatmap ---\n")

# Replace any NaN/Inf values in the correlation matrix with NA before visualization
median_spearman[!is.finite(median_spearman) & row(median_spearman) != col(median_spearman)] <- NA
min_cor <- min(median_spearman, na.rm = TRUE)
if (!is.finite(min_cor) || min_cor >= 1) min_cor <- 0.9
if (all(is.na(median_spearman[row(median_spearman) != col(median_spearman)]))) {
  cat("  WARNING: All pairwise correlations are NA — check input data quality\n")
}
col_fun <- colorRamp2(
  seq(min_cor, 1, length.out = 100),
  colorRampPalette(c("#FEE090", "#E0F3F8", "#91BFDB", "#4575B4"))(100)
)

cell_fun <- function(j, i, x, y, width, height, fill) {
  val <- median_spearman[i, j]
  label <- if (is.finite(val)) sprintf("%.2f", val) else "N/A"
  grid.text(label, x, y, gp = gpar(fontsize = 12, fontface = "bold"))
}

ht <- Heatmap(median_spearman,
  name = "Median Spearman",
  col = col_fun,
  cell_fun = cell_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_dend = TRUE,
  show_column_dend = TRUE,
  row_names_gp = gpar(fontsize = 13),
  column_names_gp = gpar(fontsize = 13),
  column_names_rot = 45,
  column_title = "Cross-Method Concordance (Median Spearman Correlation)",
  column_title_gp = gpar(fontsize = 16, fontface = "bold"),
  heatmap_legend_param = list(
    title = "Median\n Spearman",
    legend_height = unit(5, "cm")
  ),
  width = unit(14, "cm"),
  height = unit(14, "cm")
)

.dev_open <- FALSE
tryCatch({
  png(file.path(FIGURES_DIR, "method_concordance_heatmap_spearman.png"),
      width = 1600, height = 1300, res = 150)
  .dev_open <- TRUE
  draw(ht, padding = unit(c(30, 30, 25, 40), "mm"))
  dev.off()
  .dev_open <- FALSE
  cat("  Saved: method_concordance_heatmap_spearman.png\n")
}, error = function(e) {
  if (.dev_open) try(dev.off(), silent = TRUE)
  cat("  Error generating concordance heatmap:", e$message, "\n")
})

# -----------------------------------------------
# Save concordance results for report generation
# -----------------------------------------------

concordance_results <- list(
  median_spearman = median_spearman,
  spearman_per_sample = spearman_per_sample
)

saveRDS(concordance_results, file.path(OUTPUT_DIR, "concordance_results.rds"))
cat("\n[DONE] Concordance results saved\n")
