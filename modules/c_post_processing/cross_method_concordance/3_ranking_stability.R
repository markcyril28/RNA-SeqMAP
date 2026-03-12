#!/usr/bin/env Rscript

# ===============================================
# CROSS-METHOD CONCORDANCE - RANKING STABILITY
# ===============================================
# For gene groups of interest (SmelDMPs, SmelGIF, SmelGRF), compares
# expression rankings across methods. Flags genes whose ranking changes
# drastically between methods.
#
# Outputs:
#   - tables/ranking_stability_{gene_group}.csv
#   - tables/ranking_instability_flagged.csv
#   - figures/ranking_bump_chart_{gene_group}.png
#   - figures/ranking_heatmap_{gene_group}.png

source(file.path(Sys.getenv("CONCORDANCE_SCRIPT_DIR", "."), "0_concordance_config.R"))

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

cat("\n=== STEP 3: Ranking Stability Analysis ===\n\n")

# Load harmonized data
data <- readRDS(HARMONIZED_RDS)
tpm_matrices <- data$tpm_matrices
common_genes <- data$common_genes
common_samples <- data$common_samples
methods <- names(tpm_matrices)

# -----------------------------------------------
# Load gene group definitions
# -----------------------------------------------

cat("--- Loading gene groups ---\n")
cat("  Gene groups dir:", GENE_GROUPS_DIR, "\n")

all_ranking_results <- list()
all_flagged_list <- list()  # Collect flagged DFs in list; rbind once at end

for (gene_group in CONCORDANCE_GENE_GROUPS) {
  cat("\n=== Gene group:", gene_group, "===\n")

  # Find CSV file
  csv_file <- file.path(GENE_GROUPS_DIR, paste0(gene_group, ".csv"))
  if (!file.exists(csv_file)) {
    cat("  [WARN] Gene group CSV not found:", csv_file, "\n")
    next
  }

  gene_df <- read.csv(csv_file, stringsAsFactors = FALSE, header = TRUE)
  gene_ids <- trimws(gene_df$Gene_ID)
  gene_names <- if ("Shortened_Name" %in% colnames(gene_df)) {
    setNames(trimws(gene_df$Shortened_Name), gene_ids)
  } else {
    setNames(gene_ids, gene_ids)
  }

  cat("  Genes in group:", length(gene_ids), "\n")

  # Match gene IDs to harmonized matrix (handles suffix differences)
  matched_genes <- match_gene_ids(gene_ids, common_genes)
  cat("  Matched in harmonized data:", length(matched_genes), "\n")

  if (length(matched_genes) < 2) {
    cat("  [SKIP] Too few matched genes\n")
    next
  }

  # -----------------------------------------------
  # Compute mean expression rank per method
  # -----------------------------------------------
  # For each method, compute mean TPM per gene across all samples,
  # then rank genes (1 = highest expressed).

  rank_matrix <- matrix(NA, nrow = length(matched_genes), ncol = length(methods),
                         dimnames = list(matched_genes, sapply(methods, get_short_name)))
  mean_tpm_matrix <- matrix(NA, nrow = length(matched_genes), ncol = length(methods),
                             dimnames = list(matched_genes, sapply(methods, get_short_name)))

  for (j in seq_along(methods)) {
    m <- methods[j]
    mean_tpm <- rowMeans(tpm_matrices[[m]][matched_genes, , drop = FALSE])
    mean_tpm_matrix[, j] <- mean_tpm
    # Rank: 1 = highest expression (ties = average)
    rank_matrix[, j] <- rank(-mean_tpm, ties.method = "average")
  }

  # -----------------------------------------------
  # Compute ranking stability metrics
  # -----------------------------------------------

  n_genes <- length(matched_genes)
  # Vectorized: avoid 3 apply() calls over rank_matrix rows
  rank_row_means <- rowMeans(rank_matrix)
  n_meth <- ncol(rank_matrix)
  rank_range <- matrixStats::rowMaxs(rank_matrix) - matrixStats::rowMins(rank_matrix)
  rank_sd <- sqrt(rowSums((rank_matrix - rank_row_means)^2) / (n_meth - 1))
  rank_cv <- rank_sd / rank_row_means
  median_rank <- matrixStats::rowMedians(rank_matrix)

  # Fractional rank change: max rank shift / total genes
  frac_rank_change <- rank_range / n_genes

  # Flag genes with drastic ranking changes
  flagged <- frac_rank_change > RANKING_CHANGE_THRESHOLD

  # Map to shortened names for display
  display_names <- gene_names[matched_genes]
  # For genes not in the original mapping (matched via suffix stripping), use the gene ID
  na_mask <- is.na(display_names)
  if (any(na_mask)) {
    # Try matching base IDs
    base_matched <- sub("\\.[0-9]+$", "", matched_genes[na_mask])
    for (k in which(na_mask)) {
      base_id <- sub("\\.[0-9]+$", "", matched_genes[k])
      if (base_id %in% names(gene_names)) {
        display_names[k] <- gene_names[base_id]
      } else {
        display_names[k] <- matched_genes[k]
      }
    }
  }

  # Build results table
  ranking_df <- data.frame(
    Gene_ID = matched_genes,
    Shortened_Name = display_names,
    Median_Rank = round(median_rank, 1),
    Rank_Range = rank_range,
    Rank_SD = round(rank_sd, 2),
    Frac_Rank_Change = round(frac_rank_change, 3),
    Flagged = flagged,
    rank_matrix,
    Mean_TPM = round(rowMeans(mean_tpm_matrix), 2),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  ranking_df <- ranking_df[order(ranking_df$Median_Rank), ]

  write.csv(ranking_df,
            file.path(TABLES_DIR, paste0("ranking_stability_", gene_group, ".csv")),
            row.names = FALSE)
  cat("  Flagged genes (rank change >", RANKING_CHANGE_THRESHOLD * 100, "%):",
      sum(flagged), "/", n_genes, "\n")

  all_ranking_results[[gene_group]] <- ranking_df

  # Collect flagged genes across groups (append to list; single rbind at end)
  if (sum(flagged) > 0) {
    flagged_rows <- ranking_df[ranking_df$Flagged, ]
    flagged_rows$Gene_Group <- gene_group
    all_flagged_list[[length(all_flagged_list) + 1]] <- flagged_rows
  }

  # -----------------------------------------------
  # Ranking heatmap
  # -----------------------------------------------

  cat("  Generating ranking heatmap...\n")

  # Use shortened names as row labels
  display_rank_mat <- rank_matrix
  rownames(display_rank_mat) <- display_names

  if (n_genes == 2) {
    col_fun <- colorRamp2(c(1, 2), c("#4A148C", "#F3E5F5"))
  } else {
    col_fun <- colorRamp2(
      c(1, ceiling(n_genes / 2), n_genes),
      c("#4A148C", "#CE93D8", "#F3E5F5")
    )
  }

  # Annotation: flag column
  flag_colors <- ifelse(flagged, "#EF8A62", "#CCCCCC")
  row_ha <- rowAnnotation(
    Stability = anno_simple(ifelse(flagged, "Unstable", "Stable"),
                            col = c(Unstable = "#EF8A62", Stable = "#81C784"),
                            width = unit(1, "cm")),
    Rank_Range = anno_barplot(rank_range,
                              gp = gpar(fill = ifelse(flagged, "#EF8A62", "#81C784")),
                              width = unit(2.5, "cm")),
    annotation_name_gp = gpar(fontsize = 9)
  )

  ht <- Heatmap(display_rank_mat,
    name = "Rank",
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.text(display_rank_mat[i, j], x, y,
                gp = gpar(fontsize = if (n_genes > 15) 9 else 11, fontface = "bold"))
    },
    row_names_gp = gpar(fontsize = if (n_genes > 15) 9 else 11),
    row_names_max_width = unit(8, "cm"),
    column_names_gp = gpar(fontsize = 12),
    column_names_rot = 45,
    column_title = paste0("Expression Ranking Across Methods: ", gene_group),
    column_title_gp = gpar(fontsize = 14, fontface = "bold"),
    right_annotation = row_ha,
    heatmap_legend_param = list(
      title = "Rank\n(1=highest)",
      legend_height = unit(4, "cm")
    )
  )

  fig_height <- max(700, 200 + n_genes * 35)
  png(file.path(FIGURES_DIR, paste0("ranking_heatmap_", gene_group, ".png")),
      width = 1600, height = fig_height, res = 150)
  on.exit(try(dev.off(), silent = TRUE), add = TRUE)
  draw(ht, padding = unit(c(30, 30, 25, 40), "mm"))
  dev.off()
  on.exit(NULL)
  cat("  Saved: ranking_heatmap_", gene_group, ".png\n", sep = "")

  # -----------------------------------------------
  # Z-score scaled (0–10) heatmap for gene group
  # -----------------------------------------------
  # Normalizes each gene's expression across methods to a common 0–10 scale
  # so all genes are directly comparable regardless of absolute expression.

  cat("  Generating Z-score scaled heatmap...\n")

  log2_mean_tpm <- log2(mean_tpm_matrix + 1)
  zscore_grp <- t(scale(t(log2_mean_tpm)))
  zscore_grp[is.nan(zscore_grp)] <- 0
  grp_row_mins <- apply(zscore_grp, 1, min, na.rm = TRUE)
  grp_row_maxs <- apply(zscore_grp, 1, max, na.rm = TRUE)
  grp_row_range <- grp_row_maxs - grp_row_mins
  zscore_grp_scaled <- (zscore_grp - grp_row_mins) / ifelse(grp_row_range == 0, 1, grp_row_range) * 10
  zscore_grp_scaled[grp_row_range == 0, ] <- 5.0

  # Use display names and same row order as ranking heatmap (by median rank)
  display_zscore_mat <- zscore_grp_scaled
  rownames(display_zscore_mat) <- display_names
  display_zscore_mat <- display_zscore_mat[order(median_rank), , drop = FALSE]
  ordered_rank_range <- rank_range[order(median_rank)]
  ordered_flagged <- flagged[order(median_rank)]

  col_fun_z <- colorRamp2(c(0, 5, 10), c("#2166AC", "#F7F7F7", "#B2182B"))

  cell_fun_z <- function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.1f", display_zscore_mat[i, j]), x, y,
              gp = gpar(fontsize = if (n_genes > 15) 9 else 11, fontface = "bold"))
  }

  row_ha_z <- rowAnnotation(
    Rank_Range = anno_barplot(ordered_rank_range,
                              gp = gpar(fill = ifelse(ordered_flagged, "#EF8A62", "#81C784")),
                              width = unit(2.5, "cm")),
    annotation_name_gp = gpar(fontsize = 9)
  )

  ht_z <- Heatmap(display_zscore_mat,
    name = "Z-Score\n(0-10)",
    col = col_fun_z,
    cell_fun = cell_fun_z,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_names_gp = gpar(fontsize = if (n_genes > 15) 9 else 11),
    row_names_max_width = unit(8, "cm"),
    column_names_gp = gpar(fontsize = 12),
    column_names_rot = 45,
    column_title = paste0("Cross-Method Expression (Z-Score 0\u201310): ", gene_group),
    column_title_gp = gpar(fontsize = 14, fontface = "bold"),
    right_annotation = row_ha_z,
    heatmap_legend_param = list(
      title = "Z-Score\n(0-10)",
      at = c(0, 2.5, 5, 7.5, 10),
      labels = c("0 (low)", "2.5", "5 (avg)", "7.5", "10 (high)"),
      legend_height = unit(4, "cm")
    )
  )

  png(file.path(FIGURES_DIR, paste0("zscore_heatmap_", gene_group, ".png")),
      width = 1600, height = fig_height, res = 150)
  on.exit(try(dev.off(), silent = TRUE), add = TRUE)
  draw(ht_z, padding = unit(c(30, 30, 25, 40), "mm"))
  dev.off()
  on.exit(NULL)
  cat("  Saved: zscore_heatmap_", gene_group, ".png\n", sep = "")

  # Save Z-score table
  zscore_grp_df <- data.frame(
    Gene_ID = matched_genes[order(median_rank)],
    Shortened_Name = display_names[order(median_rank)],
    zscore_grp_scaled[order(median_rank), , drop = FALSE],
    Mean_TPM = round(rowMeans(mean_tpm_matrix)[order(median_rank)], 2),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  write.csv(zscore_grp_df,
            file.path(TABLES_DIR, paste0("zscore_expression_", gene_group, ".csv")),
            row.names = FALSE)
  cat("  Saved: zscore_expression_", gene_group, ".csv\n", sep = "")

  # -----------------------------------------------
  # Bump chart (ranking across methods)
  # -----------------------------------------------

  cat("  Generating bump chart...\n")

  n_methods_plot <- ncol(rank_matrix)
  method_labels <- colnames(rank_matrix)

  # Assign colors per gene
  gene_colors <- colorRampPalette(c("#4A148C", "#7B1FA2", "#AB47BC",
                                     "#CE93D8", "#2196F3", "#4CAF50",
                                     "#FF9800", "#F44336", "#795548"))(n_genes)

  png(file.path(FIGURES_DIR, paste0("ranking_bump_chart_", gene_group, ".png")),
      width = max(1200, 260 * n_methods_plot), height = max(850, 120 + n_genes * 40), res = 150)
  on.exit(try(dev.off(), silent = TRUE), add = TRUE)

  par(mar = c(10, 10, 6, 20), xpd = TRUE)
  plot(1, type = "n",
       xlim = c(0.5, n_methods_plot + 0.5),
       ylim = c(n_genes + 0.5, 0.5),
       xlab = "", ylab = "Rank (1 = highest)",
       xaxt = "n", yaxt = "n",
       main = paste0("Expression Ranking Stability: ", gene_group),
       cex.main = 1.3, cex.lab = 1.2)

  axis(1, at = seq_len(n_methods_plot), labels = method_labels, las = 2, cex.axis = 1.0)
  axis(2, at = seq_len(n_genes), las = 1, cex.axis = 0.9)

  for (i in seq_len(n_genes)) {
    ranks <- rank_matrix[i, ]
    lwd_val <- if (flagged[i]) 2.5 else 1.2
    lty_val <- if (flagged[i]) 1 else 2

    lines(seq_len(n_methods_plot), ranks, col = gene_colors[i],
          lwd = lwd_val, lty = lty_val)
    points(seq_len(n_methods_plot), ranks, col = gene_colors[i],
           pch = 16, cex = 1.2)
  }

  # Legend outside plot
  legend("right", inset = c(-0.35, 0),
         legend = display_names,
         col = gene_colors, lwd = 2, pch = 16,
         cex = if (n_genes > 12) 0.75 else 0.95,
         ncol = if (n_genes > 20) 2 else 1,
         bg = "white")

  dev.off()
  on.exit(NULL)
  cat("  Saved: ranking_bump_chart_", gene_group, ".png\n", sep = "")
}

# -----------------------------------------------
# Save combined flagged genes
# -----------------------------------------------

all_flagged_genes <- if (length(all_flagged_list) > 0) do.call(rbind, all_flagged_list) else data.frame()
if (nrow(all_flagged_genes) > 0) {
  write.csv(all_flagged_genes,
            file.path(TABLES_DIR, "ranking_instability_flagged.csv"),
            row.names = FALSE)
  cat("\n  Total flagged genes across all groups:", nrow(all_flagged_genes), "\n")
} else {
  cat("\n  No genes flagged for ranking instability\n")
}

# Save ranking results for report
saveRDS(all_ranking_results, file.path(OUTPUT_DIR, "ranking_results.rds"))
cat("\n[DONE] Ranking stability analysis complete\n")
