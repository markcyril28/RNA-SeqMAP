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

# Ranking stability compares gene rankings across alignment methods — only meaningful
# in cross_method mode where tpm_matrices keys are method names.
if (CONCORDANCE_MODE != "cross_method") {
  cat("[SKIP] Ranking stability analysis only applies to cross_method mode",
      "(current:", CONCORDANCE_MODE, ")\n")
  quit(save = "no", status = 0)
}

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

# matrixStats provides fast row-wise min/max/median; fall back to base R apply()
.HAS_MATRIXSTATS <- requireNamespace("matrixStats", quietly = TRUE)
.rowMaxs <- if (.HAS_MATRIXSTATS) matrixStats::rowMaxs else function(x, ...) apply(x, 1, max, na.rm = TRUE)
.rowMins <- if (.HAS_MATRIXSTATS) matrixStats::rowMins else function(x, ...) apply(x, 1, min, na.rm = TRUE)
.rowMedians <- if (.HAS_MATRIXSTATS) matrixStats::rowMedians else function(x, ...) apply(x, 1, median, na.rm = TRUE)

cat("\n=== STEP 3: Ranking Stability Analysis ===\n\n")

# Load harmonized data
data <- readRDS(HARMONIZED_RDS)
tpm_matrices <- data$tpm_matrices
common_genes <- data$common_genes
common_samples <- data$common_samples
methods <- names(tpm_matrices)

if (length(methods) < 2) {
  stop("Ranking stability requires >= 2 methods, but only found: ",
       paste(methods, collapse = ", "))
}

# -----------------------------------------------
# Load gene group definitions
# -----------------------------------------------

cat("--- Loading gene groups ---\n")
cat("  Gene groups dir:", GENE_GROUPS_DIR, "\n")

all_ranking_results <- list()
# Pre-allocate flagged list to avoid O(n) reallocation per append
all_flagged_list <- vector("list", length(CONCORDANCE_GENE_GROUPS))
.flagged_idx <- 0L

# Cache directory listing once — O(D) scan reused across all gene groups
# Avoids O(G × D) repeated list.files() calls for G gene groups
.gene_group_csv_cache <- list.files(GENE_GROUPS_DIR, pattern = "\\.csv$",
                                     recursive = TRUE, full.names = TRUE)
# Pre-compute basenames once (O(C)) instead of per-group (O(G×C))
.gene_group_csv_basenames <- basename(.gene_group_csv_cache)

for (gene_group in CONCORDANCE_GENE_GROUPS) {
  cat("\n=== Gene group:", gene_group, "===\n")

  # Find CSV file (search top-level first, then cached subdirectory listing)
  csv_file <- file.path(GENE_GROUPS_DIR, paste0(gene_group, ".csv"))
  if (!file.exists(csv_file)) {
    # O(C) vectorized match on pre-computed basenames instead of per-group basename() call
    candidates <- .gene_group_csv_cache[.gene_group_csv_basenames == paste0(gene_group, ".csv")]
    if (length(candidates) > 0) {
      csv_file <- candidates[1]
    } else {
      cat("  [WARN] Gene group CSV not found in", GENE_GROUPS_DIR, "or subdirectories:", gene_group, "\n")
      next
    }
  }

  gene_df <- .fast_read_csv(csv_file)
  if (!"Gene_ID" %in% colnames(gene_df)) {
    cat("  [WARN] Gene group CSV lacks 'Gene_ID' column:", csv_file, "\n")
    next
  }
  gene_ids <- trimws(gene_df$Gene_ID)
  # Filter out empty/NA gene IDs that would create invalid lookup keys
  valid_mask <- nzchar(gene_ids) & !is.na(gene_ids)
  gene_ids <- gene_ids[valid_mask]
  gene_names <- if ("Shortened_Name" %in% colnames(gene_df)) {
    setNames(trimws(gene_df$Shortened_Name[valid_mask]), gene_ids)
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

  # Pre-compute short names once — avoids O(M) sapply per matrix creation. O(M) total.
  .short_methods <- vapply(methods, get_short_name, character(1))
  rank_matrix <- matrix(NA, nrow = length(matched_genes), ncol = length(methods),
                         dimnames = list(matched_genes, .short_methods))
  mean_tpm_matrix <- matrix(NA, nrow = length(matched_genes), ncol = length(methods),
                             dimnames = list(matched_genes, .short_methods))

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
  rank_range <- .rowMaxs(rank_matrix) - .rowMins(rank_matrix)
  rank_sd <- sqrt(rowSums((rank_matrix - rank_row_means)^2) / (n_meth - 1))
  median_rank <- .rowMedians(rank_matrix)

  # Fractional rank change: max rank shift / total genes
  frac_rank_change <- rank_range / n_genes

  # Pre-compute overall mean TPM (mean-of-means across methods) — O(G×M)
  # Reused in both the summary table and zscore table, avoiding redundant rowMeans()
  .overall_mean_tpm <- rowMeans(mean_tpm_matrix)

  # Flag genes with drastic ranking changes
  flagged <- frac_rank_change > RANKING_CHANGE_THRESHOLD

  # Map to shortened names for display
  display_names <- gene_names[matched_genes]
  # For genes not in the original mapping (matched via suffix stripping), use the gene ID
  na_mask <- is.na(display_names)
  if (any(na_mask)) {
    # Vectorized base ID matching — O(N) vs O(N × lookup) loop
    base_ids_na <- sub("\\.[0-9]+$", "", matched_genes[na_mask])
    in_names <- base_ids_na %in% names(gene_names)
    display_names[na_mask][in_names] <- gene_names[base_ids_na[in_names]]
    display_names[na_mask][!in_names] <- matched_genes[na_mask][!in_names]
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
    Mean_TPM = round(.overall_mean_tpm, 2),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  ranking_df <- ranking_df[order(ranking_df$Median_Rank), ]

  if (.conc_use_dt) {
    data.table::fwrite(ranking_df, file.path(TABLES_DIR, paste0("ranking_stability_", gene_group, ".csv")))
  } else {
    write.csv(ranking_df, file.path(TABLES_DIR, paste0("ranking_stability_", gene_group, ".csv")), row.names = FALSE)
  }
  cat("  Flagged genes (rank change >", RANKING_CHANGE_THRESHOLD * 100, "%):",
      sum(flagged), "/", n_genes, "\n")

  all_ranking_results[[gene_group]] <- ranking_df

  # Collect flagged genes into pre-allocated list; O(1) indexed insert
  if (sum(flagged) > 0) {
    flagged_rows <- ranking_df[ranking_df$Flagged, ]
    flagged_rows$Gene_Group <- gene_group
    .flagged_idx <- .flagged_idx + 1L
    all_flagged_list[[.flagged_idx]] <- flagged_rows
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

  .rank_font <- if (n_genes > 15) 9 else 11
  .rank_body_h <- max(8, n_genes * 0.8)
  .rank_body_w <- max(10, ncol(display_rank_mat) * 2.5)
  .fig_layout <- calc_figure_layout(
    row_labels = rownames(display_rank_mat),
    col_labels = colnames(display_rank_mat),
    col_rot = 45, hm_body_cm = c(.rank_body_w, .rank_body_h),
    has_dendro = FALSE, has_title = TRUE,
    legend_width_cm = 6, font_size = .rank_font
  )
  .dev_open <- FALSE
  tryCatch({
    png(file.path(FIGURES_DIR, paste0("ranking_heatmap_", gene_group, ".png")),
        width = .fig_layout$width, height = .fig_layout$height, res = FIGURE_DPI)
    .dev_open <- TRUE
    draw(ht, padding = .fig_layout$padding)
    dev.off()
    .dev_open <- FALSE
    cat("  Saved: ranking_heatmap_", gene_group, ".png\n", sep = "")
  }, error = function(e) {
    if (.dev_open) try(dev.off(), silent = TRUE)
    cat("  Error generating ranking heatmap:", e$message, "\n")
  })

  # -----------------------------------------------
  # Z-score scaled (0–10) heatmap for gene group
  # -----------------------------------------------
  # Normalizes each gene's expression across methods to a common 0–10 scale
  # so all genes are directly comparable regardless of absolute expression.

  cat("  Generating Z-score scaled heatmap...\n")

  log2_mean_tpm <- log2(mean_tpm_matrix + 1)
  # Row-wise z-score: vectorized rowMeans/rowSums avoids 3 G×M matrix allocations
  # from t(scale(t(...))). O(G×M) arithmetic with zero transpositions.
  .rm <- rowMeans(log2_mean_tpm, na.rm = TRUE)
  .rsd <- sqrt(rowSums((log2_mean_tpm - .rm)^2, na.rm = TRUE) / max(ncol(log2_mean_tpm) - 1, 1))
  .rsd[.rsd == 0] <- 1  # prevent division by zero (constant rows → z-score = 0)
  zscore_grp <- (log2_mean_tpm - .rm) / .rsd
  zscore_grp[is.nan(zscore_grp)] <- 0
  grp_row_mins <- .rowMins(zscore_grp)
  grp_row_maxs <- .rowMaxs(zscore_grp)
  grp_row_range <- grp_row_maxs - grp_row_mins
  zscore_grp_scaled <- (zscore_grp - grp_row_mins) / ifelse(grp_row_range == 0, 1, grp_row_range) * 10
  zscore_grp_scaled[grp_row_range == 0, ] <- 5.0

  # Pre-compute sort order once — reused for zscore matrix, rank_range, flagged, and CSV export
  # O(G log G) sort done once instead of 4× redundant sorts
  .rank_order <- order(median_rank)

  # Use display names and same row order as ranking heatmap (by median rank)
  display_zscore_mat <- zscore_grp_scaled
  rownames(display_zscore_mat) <- display_names
  display_zscore_mat <- display_zscore_mat[.rank_order, , drop = FALSE]
  ordered_rank_range <- rank_range[.rank_order]
  ordered_flagged <- flagged[.rank_order]

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

  .zscore_layout <- calc_figure_layout(
    row_labels = rownames(display_zscore_mat),
    col_labels = colnames(display_zscore_mat),
    col_rot = 45, hm_body_cm = c(.rank_body_w, .rank_body_h),
    has_dendro = FALSE, has_title = TRUE,
    legend_width_cm = 6, font_size = .rank_font
  )
  .dev_open <- FALSE
  tryCatch({
    png(file.path(FIGURES_DIR, paste0("zscore_heatmap_", gene_group, ".png")),
        width = .zscore_layout$width, height = .zscore_layout$height, res = FIGURE_DPI)
    .dev_open <- TRUE
    draw(ht_z, padding = .zscore_layout$padding)
    dev.off()
    .dev_open <- FALSE
    cat("  Saved: zscore_heatmap_", gene_group, ".png\n", sep = "")
  }, error = function(e) {
    if (.dev_open) try(dev.off(), silent = TRUE)
    cat("  Error generating zscore heatmap:", e$message, "\n")
  })

  # Save Z-score table (reuse pre-computed .rank_order)
  zscore_grp_df <- data.frame(
    Gene_ID = matched_genes[.rank_order],
    Shortened_Name = display_names[.rank_order],
    zscore_grp_scaled[.rank_order, , drop = FALSE],
    Mean_TPM = round(.overall_mean_tpm[.rank_order], 2),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  if (.conc_use_dt) {
    data.table::fwrite(zscore_grp_df, file.path(TABLES_DIR, paste0("zscore_expression_", gene_group, ".csv")))
  } else {
    write.csv(zscore_grp_df, file.path(TABLES_DIR, paste0("zscore_expression_", gene_group, ".csv")), row.names = FALSE)
  }
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

  # Auto-calculate bump chart dimensions from content
  .bump_max_method_chars <- max(nchar(method_labels))
  .bump_max_gene_chars   <- max(nchar(rownames(rank_matrix)))
  .bump_mar_bottom <- max(8, ceiling(.bump_max_method_chars * 0.7))
  .bump_mar_right  <- max(16, ceiling(.bump_max_gene_chars * 1.2))
  .bump_w <- max(1400, 300 * n_methods_plot + .bump_mar_right * 12) / 300 * FIGURE_DPI
  .bump_h <- max(1050, 200 + n_genes * 40 + .bump_mar_bottom * 12) / 300 * FIGURE_DPI

  .dev_open <- FALSE
  tryCatch({
    png(file.path(FIGURES_DIR, paste0("ranking_bump_chart_", gene_group, ".png")),
        width = .bump_w, height = .bump_h, res = FIGURE_DPI)
    .dev_open <- TRUE

    par(mar = c(.bump_mar_bottom, 12, 8, .bump_mar_right), xpd = TRUE)
    plot(1, type = "n",
         xlim = c(0.5, n_methods_plot + 0.5),
         ylim = c(n_genes + 0.5, 0.5),
         xlab = "", ylab = "Rank (1 = highest)",
         xaxt = "n", yaxt = "n",
         main = paste0("Expression Ranking Stability: ", gene_group),
         cex.main = 1.3, cex.lab = 1.2)

    axis(1, at = seq_len(n_methods_plot), labels = method_labels, las = 2, cex.axis = 1.0)
    axis(2, at = seq_len(n_genes), las = 1, cex.axis = 0.9)

    # Pre-compute loop constants — avoids G redundant seq_len() calls
    # and G conditional branches for lwd/lty values. O(G) vectorized ifelse.
    x_seq <- seq_len(n_methods_plot)
    lwd_vals <- ifelse(flagged, 2.5, 1.2)
    lty_vals <- ifelse(flagged, 1L, 2L)
    # matlines/matpoints: single C-level call replaces O(G) R-level lines()/points()
    matlines(x_seq, t(rank_matrix), col = gene_colors, lwd = lwd_vals, lty = lty_vals)
    matpoints(x_seq, t(rank_matrix), col = gene_colors, pch = 16, cex = 1.2)

    # Legend outside plot
    legend("right", inset = c(-0.35, 0),
           legend = display_names,
           col = gene_colors, lwd = 2, pch = 16,
           cex = if (n_genes > 12) 0.75 else 0.95,
           ncol = if (n_genes > 20) 2 else 1,
           bg = "white")

    dev.off()
    .dev_open <- FALSE
    cat("  Saved: ranking_bump_chart_", gene_group, ".png\n", sep = "")
  }, error = function(e) {
    if (.dev_open) try(dev.off(), silent = TRUE)
    cat("  Error generating bump chart:", e$message, "\n")
  })
}

# -----------------------------------------------
# Save combined flagged genes
# -----------------------------------------------

# Trim unused pre-allocated slots before rbind
all_flagged_list <- all_flagged_list[seq_len(.flagged_idx)]
# Use rbindlist when data.table available — avoids O(R²) copy overhead of do.call(rbind).
# Reuse .conc_use_dt from 0_concordance_config.R (avoids redundant requireNamespace probe).
all_flagged_genes <- if (.flagged_idx > 0L) {
  if (.conc_use_dt) data.table::setDF(data.table::rbindlist(all_flagged_list, use.names = TRUE, fill = TRUE))
  else do.call(rbind, all_flagged_list)
} else data.frame()
if (nrow(all_flagged_genes) > 0) {
  if (.conc_use_dt) {
    data.table::fwrite(all_flagged_genes, file.path(TABLES_DIR, "ranking_instability_flagged.csv"))
  } else {
    write.csv(all_flagged_genes, file.path(TABLES_DIR, "ranking_instability_flagged.csv"), row.names = FALSE)
  }
  cat("\n  Total flagged genes across all groups:", nrow(all_flagged_genes), "\n")
} else {
  cat("\n  No genes flagged for ranking instability\n")
}

# Save ranking results for report
saveRDS(all_ranking_results, file.path(OUTPUT_DIR, "ranking_results.rds"))
cat("\n[DONE] Ranking stability analysis complete\n")
