#!/usr/bin/env Rscript

# ===============================================
# CROSS-METHOD CONCORDANCE - QUANTIFICATION CONCORDANCE
# ===============================================
# Computes pairwise Spearman correlations between methods and generates
# a median Spearman concordance heatmap.
#
# Big O: O(M^2 × S × G × log(G)) where M=methods, S=samples, G=genes.
#   - M^2: pairwise method combinations
#   - S × G × log(G): per-sample column-wise rank transform
#   Pre-computed log2 and nonzero masks reduce constant factor by ~2x.
#
# Outputs:
#   - tables/pairwise_spearman_per_sample.csv
#   - tables/median_correlation_matrix_spearman.csv
#   - figures/method_concordance_heatmap_spearman.png

# Skip re-sourcing when running under concordance_batch_dispatcher.R (already loaded)
if (!exists(".CONC_BATCH_CONFIG_LOADED") || !isTRUE(.CONC_BATCH_CONFIG_LOADED)) {
  .conc_dir <- Sys.getenv("CONCORDANCE_SCRIPT_DIR", {
    if (nzchar(Sys.getenv("WF_MANAGED_ENV", "")))
      stop("[QUANT_CONCORDANCE] CONCORDANCE_SCRIPT_DIR is required under workflow manager (WF_MANAGED_ENV is set).")
    "."
  })
  source(file.path(.conc_dir, "0_concordance_config.R"))

  suppressPackageStartupMessages({
    library(ComplexHeatmap)
    library(circlize)
    library(grid)
  })
}

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

# Pre-compute short names once — reused below in pair_names, build_median_cor_matrix, etc.
# O(M) vapply; avoids redundant get_short_name() calls on lines 60, 131, 151.
short_names <- vapply(methods, get_short_name, character(1))
cat("Methods:", paste(short_names, collapse = ", "), "\n")
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
pair_names <- vapply(method_pairs, function(p) paste(short_names[p[1]], "vs", short_names[p[2]]), character(1))

spearman_per_sample <- matrix(NA, nrow = length(common_samples), ncol = length(method_pairs),
                               dimnames = list(common_samples, pair_names))

# Pre-compute log2(TPM+1) once per method (avoids redundant log2 per pair × sample)
# O(methods × genes × samples) total — done once, reused across all pairs
log2_matrices <- lapply(tpm_matrices, function(mat) log2(mat + 1))

# Nonzero masks computed on-demand per pair below — avoids storing M full G×S
# logical matrices simultaneously (only 2 needed at a time).

for (i in seq_along(method_pairs)) {
  m1 <- method_pairs[[i]][1]
  m2 <- method_pairs[[i]][2]
  mat1 <- log2_matrices[[m1]]
  mat2 <- log2_matrices[[m2]]

  # On-demand nonzero mask for this pair — O(G×S) bitwise OR, no persistent storage
  nonzero_mask <- (mat1 > 0) | (mat2 > 0)
  nonzero_per_sample <- colSums(nonzero_mask)

  # Vectorized skip detection — O(S) logical mask replaces O(S log S) setdiff
  nonzero_counts <- nonzero_per_sample[common_samples]
  nonzero_counts[is.na(nonzero_counts)] <- 0L
  valid_mask_s <- nonzero_counts >= CORRELATION_MIN_GENES
  valid_samples <- common_samples[valid_mask_s]
  skipped_samples <- common_samples[!valid_mask_s]

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
    # O(genes × samples × log(genes)) per matrix.
    # IMPORTANT: must use base R rank(na.last="keep") here — colRanks includes NA
    # positions in ranking, inflating non-NA rank values. Post-hoc NA masking would
    # remove NA ranks but non-NA ranks remain wrong (ranked 1..N instead of 1..K
    # where K = non-NA count). This changes Spearman correlation results.
    r1 <- apply(m1_valid, 2, rank, na.last = "keep")
    r2 <- apply(m2_valid, 2, rank, na.last = "keep")

    # Pearson correlation on ranks = Spearman (vectorized per column)
    # Center each column — t(t(r)-cm) broadcasts without allocating a full G×S
    # intermediate vector (rep() created one). O(G×S) with better memory reuse.
    cm1 <- colMeans(r1, na.rm = TRUE)
    cm2 <- colMeans(r2, na.rm = TRUE)
    # Guard: if a column is all-NA, colMeans returns NaN → replace with 0 to avoid
    # NaN propagation through the entire column during centering.
    cm1[!is.finite(cm1)] <- 0
    cm2[!is.finite(cm2)] <- 0
    # sweep(MARGIN=2) avoids two intermediate G×S transpose allocations per variable
    r1_centered <- sweep(r1, 2, cm1, "-")
    r2_centered <- sweep(r2, 2, cm2, "-")
    # NOTE: Zeroing NAs biases correlation vs pairwise-complete Spearman when
    # many genes are unexpressed in one method. Acceptable for cross-method
    # concordance ranking (relative, not absolute) but not for formal inference.
    r1_centered[is.na(r1_centered)] <- 0
    r2_centered[is.na(r2_centered)] <- 0

    num <- colSums(r1_centered * r2_centered)
    # Explicit multiply avoids ^ S3 dispatch overhead on matrix operand
    den <- sqrt(colSums(r1_centered * r1_centered) * colSums(r2_centered * r2_centered))
    den[den == 0] <- 1  # guard against zero-variance
    spearman_per_sample[valid_samples, i] <- num / den
  }

  if (length(skipped_samples) > 0) {
    cat("  Warning: Skipped", length(skipped_samples), "samples for pair",
        short_names[m1], "vs", short_names[m2],
        "(fewer than", CORRELATION_MIN_GENES, "expressed genes):",
        paste(head(skipped_samples, 5), collapse = ", "),
        if (length(skipped_samples) > 5) "..." else "", "\n")
  }
}

# Save per-sample correlations
.spearman_out <- data.frame(Sample = rownames(spearman_per_sample), spearman_per_sample, check.names = FALSE)
if (.conc_use_dt) {
  data.table::fwrite(.spearman_out, file.path(TABLES_DIR, "pairwise_spearman_per_sample.csv"))
} else {
  write.csv(.spearman_out, file.path(TABLES_DIR, "pairwise_spearman_per_sample.csv"), row.names = FALSE)
}
rm(.spearman_out)

cat("  Saved per-sample Spearman correlation table\n")

# -----------------------------------------------
# 2.2 Median correlation matrices (methods x methods)
# -----------------------------------------------

cat("--- Computing median Spearman correlation matrix ---\n")

build_median_cor_matrix <- function(per_sample_mat) {
  cor_mat <- matrix(1, nrow = n_methods, ncol = n_methods,
                     dimnames = list(short_names, short_names))
  # O(P) with O(1) index lookup via precomputed named vector (was O(P × M) with which())
  method_indices <- setNames(seq_along(methods), methods)
  for (i in seq_along(method_pairs)) {
    m1_idx <- method_indices[method_pairs[[i]][1]]
    m2_idx <- method_indices[method_pairs[[i]][2]]
    med_cor <- median(per_sample_mat[, i], na.rm = TRUE)
    cor_mat[m1_idx, m2_idx] <- med_cor
    cor_mat[m2_idx, m1_idx] <- med_cor
  }
  return(cor_mat)
}

median_spearman <- build_median_cor_matrix(spearman_per_sample)

.median_out <- data.frame(Method = rownames(median_spearman), median_spearman, check.names = FALSE)
if (.conc_use_dt) {
  data.table::fwrite(.median_out, file.path(TABLES_DIR, "median_correlation_matrix_spearman.csv"))
} else {
  write.csv(.median_out, file.path(TABLES_DIR, "median_correlation_matrix_spearman.csv"), row.names = FALSE)
}
rm(.median_out)

cat("  Spearman median correlations:\n")
print(round(median_spearman, 3))

# -----------------------------------------------
# 2.3 Median Spearman concordance heatmap
# -----------------------------------------------

cat("\n--- Generating Median Spearman concordance heatmap ---\n")

# Replace any NaN/Inf values in the correlation matrix with NA before visualization
# Compute off-diagonal mask once — avoids 2× O(M²) row()/col() allocations
.off_diag <- row(median_spearman) != col(median_spearman)
median_spearman[!is.finite(median_spearman) & .off_diag] <- NA
min_cor <- min(median_spearman[.off_diag], na.rm = TRUE)
if (!is.finite(min_cor) || min_cor >= 1) min_cor <- 0.9
if (all(is.na(median_spearman[.off_diag]))) {
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
  column_title = get_concordance_title(),
  column_title_gp = gpar(fontsize = 16, fontface = "bold"),
  heatmap_legend_param = list(
    title = "Median Spearman\nCorrelation\n(Cross-Method\nTPM)",
    legend_height = unit(5, "cm")
  ),
  width = unit(14, "cm"),
  height = unit(14, "cm")
)

.fig_layout <- calc_figure_layout(
  row_labels = rownames(median_spearman),
  col_labels = colnames(median_spearman),
  col_rot = 45, hm_body_cm = c(14, 14),
  has_dendro = TRUE, has_title = TRUE,
  legend_width_cm = 5, font_size = 13
)
.dev_open <- FALSE
tryCatch({
  png(file.path(FIGURES_DIR, "method_concordance_heatmap_spearman.png"),
      width = .fig_layout$width, height = .fig_layout$height, res = FIGURE_DPI)
  .dev_open <- TRUE
  draw(ht, padding = .fig_layout$padding)
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
