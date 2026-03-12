#!/usr/bin/env Rscript

# ===============================================
# CROSS-METHOD CONCORDANCE - REPORT GENERATOR
# ===============================================
# Generates a unified Markdown report with embedded figure references,
# correlation matrices, gene lists, and analysis summaries.
#
# Output: 3_POST_PROC/cross_method_concordance_report.md

source(file.path(Sys.getenv("CONCORDANCE_SCRIPT_DIR", "."), "0_concordance_config.R"))

cat("\n=== STEP 4: Generating Concordance Report ===\n\n")

# Load all results
data <- readRDS(HARMONIZED_RDS)
concordance <- readRDS(file.path(OUTPUT_DIR, "concordance_results.rds"))
ranking_results <- if (file.exists(file.path(OUTPUT_DIR, "ranking_results.rds"))) {
  readRDS(file.path(OUTPUT_DIR, "ranking_results.rds"))
} else {
  list()
}

methods <- names(data$tpm_matrices)
short_names <- sapply(methods, get_short_name)

# Report output path
report_path <- file.path(POST_PROC_BASE, "cross_method_concordance_report.md")

# -----------------------------------------------
# Build report
# -----------------------------------------------

lines <- character()
add <- function(...) {
  lines <<- c(lines, paste0(...))
}

add("# Cross-Method Concordance Report")
add("")
add("**Reference genome:** ", MASTER_REFERENCE)
add("**Date generated:** ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
add("**Methods compared:** ", paste(short_names, collapse = ", "))
add("**Common genes:** ", length(data$common_genes))
add("**Common samples:** ", length(data$common_samples))
add("")

# -----------------------------------------------
# Section 1: Method Overview
# -----------------------------------------------

add("## 1. Method Overview")
add("")
add("| Method | Short Name | Genes (Raw) | Samples (Raw) | Genes (Harmonized) | Samples (Harmonized) |")
add("|--------|-----------|-------------|--------------|-------------------|---------------------|")
for (i in seq_len(nrow(data$method_stats))) {
  row <- data$method_stats[i, ]
  add("| ", row$method, " | ", row$short_name, " | ",
      format(row$n_genes_raw, big.mark = ","), " | ", row$n_samples_raw, " | ",
      format(row$n_genes_harmonized, big.mark = ","), " | ", row$n_samples_harmonized, " |")
}
add("")

# Gene set overlap summary
add("### Gene Set Overlaps")
add("")
add("| Method Pair | Shared Genes |")
add("|------------|-------------|")
gene_sets <- data$gene_sets_raw
method_names <- names(gene_sets)
for (i in seq_len(length(method_names) - 1)) {
  for (j in seq(i + 1, length(method_names))) {
    overlap <- length(intersect(gene_sets[[method_names[i]]], gene_sets[[method_names[j]]]))
    add("| ", get_short_name(method_names[i]), " & ", get_short_name(method_names[j]),
        " | ", format(overlap, big.mark = ","), " |")
  }
}
add("")

# -----------------------------------------------
# Section 2: Quantification Concordance
# -----------------------------------------------

add("## 2. Quantification Concordance")
add("")

add("### 2.1 Median Spearman Correlation Matrix")
add("")
add("Pairwise median Spearman correlations across all samples (computed on log2(TPM+1) of expressed genes).")
add("")

# Format correlation matrix as markdown table
sp_mat <- concordance$median_spearman
add("| Method |", paste(colnames(sp_mat), collapse = " | "), " |")
add("|", paste(rep("------", ncol(sp_mat) + 1), collapse = "|"), "|")
for (i in seq_len(nrow(sp_mat))) {
  vals <- sprintf("%.3f", sp_mat[i, ])
  add("| **", rownames(sp_mat)[i], "** | ", paste(vals, collapse = " | "), " |")
}
add("")

add("### 2.2 Median Pearson Correlation Matrix")
add("")

pe_mat <- concordance$median_pearson
add("| Method |", paste(colnames(pe_mat), collapse = " | "), " |")
add("|", paste(rep("------", ncol(pe_mat) + 1), collapse = "|"), "|")
for (i in seq_len(nrow(pe_mat))) {
  vals <- sprintf("%.3f", pe_mat[i, ])
  add("| **", rownames(pe_mat)[i], "** | ", paste(vals, collapse = " | "), " |")
}
add("")

add("### 2.3 Concordance Heatmaps")
add("")
add("![Spearman Concordance](cross_method_concordance/figures/method_concordance_heatmap_spearman.png)")
add("")
add("![Pearson Concordance](cross_method_concordance/figures/method_concordance_heatmap_pearson.png)")
add("")

# ---- Data-driven interpretation of concordance heatmaps ----
sp_upper_vals <- concordance$median_spearman[upper.tri(concordance$median_spearman)]
pe_upper_vals <- concordance$median_pearson[upper.tri(concordance$median_pearson)]
sp_overall <- median(sp_upper_vals)
pe_overall <- median(pe_upper_vals)
sp_min <- min(sp_upper_vals)
sp_max <- max(sp_upper_vals)
pe_min <- min(pe_upper_vals)
pe_max <- max(pe_upper_vals)

# Find best/worst pairs
sp_full <- concordance$median_spearman
diag(sp_full) <- NA
sp_full[lower.tri(sp_full)] <- NA
pe_full <- concordance$median_pearson
diag(pe_full) <- NA
pe_full[lower.tri(pe_full)] <- NA

sp_best_idx <- which(sp_full == max(sp_full, na.rm = TRUE), arr.ind = TRUE)
sp_worst_idx <- which(sp_full == min(sp_full, na.rm = TRUE), arr.ind = TRUE)
pe_best_idx <- which(pe_full == max(pe_full, na.rm = TRUE), arr.ind = TRUE)
pe_worst_idx <- which(pe_full == min(pe_full, na.rm = TRUE), arr.ind = TRUE)

sp_best_name <- paste0(rownames(sp_full)[sp_best_idx[1,1]], " & ", colnames(sp_full)[sp_best_idx[1,2]])
sp_worst_name <- paste0(rownames(sp_full)[sp_worst_idx[1,1]], " & ", colnames(sp_full)[sp_worst_idx[1,2]])
pe_best_name <- paste0(rownames(pe_full)[pe_best_idx[1,1]], " & ", colnames(pe_full)[pe_best_idx[1,2]])
pe_worst_name <- paste0(rownames(pe_full)[pe_worst_idx[1,1]], " & ", colnames(pe_full)[pe_worst_idx[1,2]])

add("#### Interpretation")
add("")

# Overall assessment
if (sp_overall >= 0.95) {
  add("**Overall concordance is excellent.** The median pairwise Spearman correlation is ",
      sprintf("%.3f", sp_overall), " (range: ", sprintf("%.3f", sp_min), "\u2013",
      sprintf("%.3f", sp_max), "), indicating that all methods produce highly consistent ",
      "expression rankings. The Pearson median is ", sprintf("%.3f", pe_overall),
      " (range: ", sprintf("%.3f", pe_min), "\u2013", sprintf("%.3f", pe_max),
      "), confirming strong agreement in absolute expression magnitudes as well.")
} else if (sp_overall >= 0.90) {
  add("**Overall concordance is good.** The median pairwise Spearman correlation is ",
      sprintf("%.3f", sp_overall), " (range: ", sprintf("%.3f", sp_min), "\u2013",
      sprintf("%.3f", sp_max), "). Most method pairs agree well on gene expression rankings, ",
      "though some pairs show moderate divergence. The Pearson median is ",
      sprintf("%.3f", pe_overall), " (range: ", sprintf("%.3f", pe_min), "\u2013",
      sprintf("%.3f", pe_max), ").")
} else if (sp_overall >= 0.85) {
  add("**Overall concordance is moderate.** The median pairwise Spearman correlation is ",
      sprintf("%.3f", sp_overall), " (range: ", sprintf("%.3f", sp_min), "\u2013",
      sprintf("%.3f", sp_max), "). Notable method-specific biases are present. The Pearson ",
      "median is ", sprintf("%.3f", pe_overall), " (range: ", sprintf("%.3f", pe_min), "\u2013",
      sprintf("%.3f", pe_max), "). Downstream results should be validated across multiple methods.")
} else {
  add("**Overall concordance is low.** The median pairwise Spearman correlation is only ",
      sprintf("%.3f", sp_overall), " (range: ", sprintf("%.3f", sp_min), "\u2013",
      sprintf("%.3f", sp_max), "). Substantial disagreement exists between methods. The Pearson ",
      "median is ", sprintf("%.3f", pe_overall), " (range: ", sprintf("%.3f", pe_min), "\u2013",
      sprintf("%.3f", pe_max), "). Results should be interpreted with caution and validated ",
      "with orthogonal data.")
}
add("")

# Best/worst pairs
add("- **Most concordant (Spearman):** ", sp_best_name, " (\u03C1 = ",
    sprintf("%.3f", max(sp_full, na.rm = TRUE)), ")")
add("- **Least concordant (Spearman):** ", sp_worst_name, " (\u03C1 = ",
    sprintf("%.3f", min(sp_full, na.rm = TRUE)), ")")
add("- **Most concordant (Pearson):** ", pe_best_name, " (r = ",
    sprintf("%.3f", max(pe_full, na.rm = TRUE)), ")")
add("- **Least concordant (Pearson):** ", pe_worst_name, " (r = ",
    sprintf("%.3f", min(pe_full, na.rm = TRUE)), ")")
add("")

# Spearman vs Pearson comparison
sp_pe_diff <- abs(sp_overall - pe_overall)
if (sp_pe_diff < 0.01) {
  add("Spearman and Pearson correlations are nearly identical, indicating that methods agree ",
      "on both the rank ordering and absolute magnitudes of gene expression. No systematic ",
      "non-linear distortions are present between methods.")
} else if (pe_overall > sp_overall + 0.02) {
  add("Pearson correlations are notably higher than Spearman correlations (",
      sprintf("%.3f", pe_overall), " vs ", sprintf("%.3f", sp_overall),
      "). This suggests methods agree well on highly-expressed genes (which dominate Pearson) ",
      "but diverge more in the rank ordering of moderate- and low-expression genes. ",
      "Downstream analyses focused on low-expression genes should use caution.")
} else if (sp_overall > pe_overall + 0.02) {
  add("Spearman correlations are notably higher than Pearson correlations (",
      sprintf("%.3f", sp_overall), " vs ", sprintf("%.3f", pe_overall),
      "). Methods preserve relative gene ordering well, but absolute expression magnitudes ",
      "differ. This may reflect differences in normalization or quantification scale between ",
      "methods. Rank-based downstream analyses (e.g., GSEA) are more robust here than those ",
      "relying on absolute expression values.")
} else {
  add("Spearman and Pearson correlations are similar (median \u03C1 = ", sprintf("%.3f", sp_overall),
      ", median r = ", sprintf("%.3f", pe_overall),
      "), indicating consistent agreement in both rank ordering and expression magnitudes.")
}
add("")

# Clustering pattern analysis
# Perform hierarchical clustering to identify method groupings
sp_dist <- as.dist(1 - concordance$median_spearman)
sp_hclust <- hclust(sp_dist, method = "complete")
sp_clusters <- cutree(sp_hclust, h = median(sp_dist))
n_clusters <- length(unique(sp_clusters))

if (n_clusters > 1) {
  add("**Clustering pattern:** Methods form ", n_clusters, " distinct clusters based on ",
      "Spearman correlation similarity:")
  add("")
  for (cl in sort(unique(sp_clusters))) {
    cl_methods <- names(sp_clusters)[sp_clusters == cl]
    add("- Cluster ", cl, ": ", paste(cl_methods, collapse = ", "))
  }
  add("")
  add("Methods within the same cluster produce more similar quantification results ",
      "and can be considered more interchangeable for this dataset.")
} else {
  add("**Clustering pattern:** All methods form a single tight cluster, indicating ",
      "uniformly high agreement across all method pairs. No method stands out as ",
      "systematically different from the others.")
}
add("")

# Count how many pairs exceed key thresholds
n_pairs_total <- length(sp_upper_vals)
n_above_95 <- sum(sp_upper_vals >= 0.95)
n_above_90 <- sum(sp_upper_vals >= 0.90)
n_below_90 <- sum(sp_upper_vals < 0.90)

add("**Threshold summary (Spearman):**")
add("")
add("- ", n_above_95, "/", n_pairs_total, " pairs with \u03C1 \u2265 0.95 (excellent agreement)")
add("- ", n_above_90 - n_above_95, "/", n_pairs_total, " pairs with 0.90 \u2264 \u03C1 < 0.95 (good agreement)")
if (n_below_90 > 0) {
  add("- ", n_below_90, "/", n_pairs_total, " pairs with \u03C1 < 0.90 (needs investigation)")
}
add("")

add("### 2.4 Per-Sample Correlation Distribution")
add("")
add("Box plots showing the distribution of Spearman correlations per sample for each method pair.")
add("")
add("![Per-sample Correlations](cross_method_concordance/figures/per_sample_correlation_boxplot.png)")
add("")

# Summary statistics for per-sample correlations
add("**Summary of per-sample Spearman correlations:**")
add("")
sp_summary <- apply(concordance$spearman_per_sample, 2, function(x) {
  c(median = median(x, na.rm = TRUE),
    min = min(x, na.rm = TRUE),
    max = max(x, na.rm = TRUE),
    IQR = IQR(x, na.rm = TRUE))
})
add("| Method Pair | Median rho | Min | Max | IQR |")
add("|------------|-----------|-----|-----|-----|")
for (j in seq_len(ncol(sp_summary))) {
  add("| ", colnames(sp_summary)[j], " | ",
      sprintf("%.3f", sp_summary["median", j]), " | ",
      sprintf("%.3f", sp_summary["min", j]), " | ",
      sprintf("%.3f", sp_summary["max", j]), " | ",
      sprintf("%.3f", sp_summary["IQR", j]), " |")
}
add("")

# ---- Data-driven interpretation of per-sample correlations ----
add("#### Interpretation")
add("")

# Identify problematic pairs (any sample below 0.90)
min_per_pair <- sp_summary["min", ]
pairs_below_90 <- names(min_per_pair)[min_per_pair < 0.90]
pairs_below_95 <- names(min_per_pair)[min_per_pair < 0.95]

# IQR analysis - tight vs variable
iqr_vals <- sp_summary["IQR", ]
median_iqr <- median(iqr_vals)
widest_pair <- names(which.max(iqr_vals))
tightest_pair <- names(which.min(iqr_vals))

if (all(min_per_pair >= 0.95)) {
  add("All method pairs maintain Spearman \u03C1 \u2265 0.95 across every sample, indicating ",
      "uniformly strong agreement. No individual samples show concerning divergence between any methods.")
} else if (all(min_per_pair >= 0.90)) {
  add("All method pairs maintain Spearman \u03C1 \u2265 0.90 across every sample. ",
      "While most pairs show excellent agreement (\u03C1 > 0.95 median), some individual ",
      "samples in certain pairs dip below 0.95, suggesting sample-specific variability ",
      "in method agreement.")
} else {
  add("**Some method pairs show concerning per-sample variability.** The following pairs ",
      "have at least one sample with Spearman \u03C1 < 0.90:")
  add("")
  for (p in pairs_below_90) {
    add("- **", p, "**: minimum \u03C1 = ", sprintf("%.3f", sp_summary["min", p]),
        " (median = ", sprintf("%.3f", sp_summary["median", p]), ")")
  }
  add("")
  add("Samples with low correlations should be investigated for sequencing quality issues, ",
      "low read depth, or method-specific mapping artifacts.")
}
add("")

add("- **Most consistent pair** (smallest IQR): ", tightest_pair,
    " (IQR = ", sprintf("%.4f", min(iqr_vals)), ")")
add("- **Most variable pair** (largest IQR): ", widest_pair,
    " (IQR = ", sprintf("%.4f", max(iqr_vals)), ")")
add("")

if (median_iqr < 0.005) {
  add("The interquartile ranges are very narrow across all pairs (median IQR = ",
      sprintf("%.4f", median_iqr), "), confirming that method agreement is highly consistent ",
      "from sample to sample with minimal sample-dependent variability.")
} else if (median_iqr < 0.02) {
  add("The interquartile ranges are modest (median IQR = ",
      sprintf("%.4f", median_iqr), "), indicating reasonably consistent agreement ",
      "across samples, with minor sample-to-sample variability.")
} else {
  add("The interquartile ranges are relatively wide (median IQR = ",
      sprintf("%.4f", median_iqr), "), indicating notable sample-to-sample variability ",
      "in method agreement. Some samples produce much more concordant results than others, ",
      "possibly reflecting differences in sequencing depth or library complexity.")
}
add("")

add("### 2.5 Pairwise Scatter Plots")
add("")
add("Mean log2(TPM+1) per gene, comparing all method pairs. Red dashed line = y=x (perfect agreement).")
add("")
add("![Pairwise Scatters](cross_method_concordance/figures/pairwise_scatter_plots.png)")
add("")

# ---- Data-driven interpretation of scatter plots ----
add("#### Interpretation")
add("")

# Compute per-pair Pearson vs Spearman on mean expression to detect non-linearity
pe_vs_sp_diffs <- sp_summary["median", ] - apply(concordance$pearson_per_sample, 2, median, na.rm = TRUE)
pairs_with_divergence <- names(pe_vs_sp_diffs[abs(pe_vs_sp_diffs) > 0.02])

add("The scatter plots visualize the gene-by-gene relationship between each pair of methods. ",
    "Points tightly along the diagonal (red dashed line) indicate strong quantitative agreement. ",
    "Systematic offsets or fan-shaped spread at low expression reveal method-specific biases.")
add("")

# Check for pairs where Pearson >> Spearman (suggesting outlier influence)
if (length(pairs_with_divergence) > 0) {
  add("**Non-linear patterns detected** in the following pairs (Pearson and Spearman ",
      "diverge by > 0.02):")
  add("")
  for (p in pairs_with_divergence) {
    sp_val <- sp_summary["median", p]
    pe_val <- median(concordance$pearson_per_sample[, p], na.rm = TRUE)
    if (pe_val > sp_val) {
      add("- **", p, "**: Pearson (", sprintf("%.3f", pe_val), ") > Spearman (",
          sprintf("%.3f", sp_val), ") \u2014 suggests agreement driven by highly-expressed ",
          "genes; low-expression genes may diverge.")
    } else {
      add("- **", p, "**: Spearman (", sprintf("%.3f", sp_val), ") > Pearson (",
          sprintf("%.3f", pe_val), ") \u2014 rank ordering is preserved but absolute ",
          "magnitudes differ, possibly due to normalization differences.")
    }
  }
  add("")
} else {
  add("Pearson and Spearman correlations are consistent across all method pairs (",
      "difference < 0.02), indicating no systematic non-linear distortions between methods.")
  add("")
}

add("")

# -----------------------------------------------
# Section 3: Discordant Genes
# -----------------------------------------------

add("## 3. Discordant Genes")
add("")
add("Genes with high coefficient of variation (CV > ", DISCORDANCE_CV_THRESHOLD,
    ") in mean expression across methods.")
add("These genes show substantially different quantification depending on the method used.")
add("")
add("**Total discordant genes:** ", concordance$n_discordant,
    " / ", length(data$common_genes), " (",
    sprintf("%.1f%%", 100 * concordance$n_discordant / length(data$common_genes)), ")")
add("")

if (concordance$n_discordant > 0) {
  add("### Top 20 Most Discordant Genes")
  add("")
  top_disc <- head(concordance$discordant_df, 20)
  method_cols <- sapply(methods, get_short_name)
  add("| Gene ID | CV | Mean log2TPM | ", paste(method_cols, collapse = " | "), " |")
  add("|---------|---:|------------:|", paste(rep("---:", length(method_cols)), collapse = "|"), "|")
  for (i in seq_len(nrow(top_disc))) {
    row <- top_disc[i, ]
    method_vals <- sprintf("%.2f", as.numeric(row[method_cols]))
    add("| ", row$Gene_ID, " | ", sprintf("%.2f", row$CV_across_methods),
        " | ", sprintf("%.2f", row$Mean_log2TPM),
        " | ", paste(method_vals, collapse = " | "), " |")
  }
  add("")
  add("![Discordant Genes Heatmap](cross_method_concordance/figures/discordant_genes_heatmap.png)")
  add("")
  add("#### Z-Score Normalized View (0\u201310 Scale)")
  add("")
  add("Each gene's expression is normalized across methods to a 0\u201310 scale (per-gene Z-score rescaled).")
  add("0 = lowest method estimate, 5 = average, 10 = highest. This removes absolute expression differences")
  add("and highlights *which methods* over- or under-estimate each gene relative to the consensus.")
  add("")
  add("![Discordant Genes Z-Score Heatmap](cross_method_concordance/figures/discordant_genes_zscore_heatmap.png)")
  add("")
  add("Full list: [`tables/discordant_genes.csv`](cross_method_concordance/tables/discordant_genes.csv)")
  add("")
  add("Z-score table: [`tables/discordant_genes_zscore.csv`](cross_method_concordance/tables/discordant_genes_zscore.csv)")
  add("")

  # ---- Data-driven interpretation of discordant genes ----
  add("#### Interpretation")
  add("")

  disc_df <- concordance$discordant_df
  disc_pct <- 100 * concordance$n_discordant / length(data$common_genes)

  if (disc_pct < 1) {
    add("Only ", sprintf("%.1f%%", disc_pct), " of genes are discordant (CV > ",
        DISCORDANCE_CV_THRESHOLD, "), indicating **very high method agreement**. ",
        "The vast majority of genes are quantified consistently regardless of method.")
  } else if (disc_pct < 5) {
    add("A small fraction (", sprintf("%.1f%%", disc_pct), ") of genes show discordant ",
        "quantification across methods. This is typical for RNA-seq and most of these ",
        "are likely low-expression genes where small absolute differences produce large CVs.")
  } else {
    add("A notable fraction (", sprintf("%.1f%%", disc_pct), ") of genes are discordant. ",
        "This suggests meaningful method-dependent biases for a subset of genes. ",
        "Downstream results for these genes should be validated across methods.")
  }
  add("")

  # Analyze which methods tend to be outliers among discordant genes
  disc_method_cols <- sapply(methods, get_short_name)
  disc_method_cols <- intersect(disc_method_cols, colnames(disc_df))
  if (length(disc_method_cols) >= 2 && nrow(disc_df) >= 5) {
    # For top discordant genes, identify which method(s) tend to be highest/lowest
    top_n <- min(50, nrow(disc_df))
    top_disc_mat <- as.matrix(disc_df[seq_len(top_n), disc_method_cols, drop = FALSE])
    # Count how often each method is the max or min per gene
    max_counts <- table(factor(disc_method_cols[apply(top_disc_mat, 1, which.max)],
                               levels = disc_method_cols))
    min_counts <- table(factor(disc_method_cols[apply(top_disc_mat, 1, which.min)],
                               levels = disc_method_cols))

    add("Among the top ", top_n, " discordant genes:")
    add("")
    add("| Method | Highest Estimate | Lowest Estimate |")
    add("|--------|:----------------:|:---------------:|")
    for (mc in disc_method_cols) {
      add("| ", mc, " | ", as.integer(max_counts[mc]), " genes | ",
          as.integer(min_counts[mc]), " genes |")
    }
    add("")

    most_high_method <- names(which.max(max_counts))
    most_low_method <- names(which.max(min_counts))
    add("**", most_high_method, "** most frequently gives the highest expression estimate ",
        "for discordant genes (", as.integer(max_counts[most_high_method]), "/", top_n, "), ",
        "while **", most_low_method, "** most frequently gives the lowest estimate (",
        as.integer(min_counts[most_low_method]), "/", top_n, "). ",
        "This pattern may reflect differences in how these methods handle multi-mapped reads, ",
        "transcript isoform resolution, or gene boundary definitions.")
    add("")

    # Expression level analysis of discordant genes
    mean_expr <- disc_df$Mean_log2TPM[seq_len(top_n)]
    n_low_expr <- sum(mean_expr < 2)
    n_mid_expr <- sum(mean_expr >= 2 & mean_expr < 6)
    n_high_expr <- sum(mean_expr >= 6)
    add("**Expression level distribution of discordant genes:**")
    add("")
    add("- Low expression (log2TPM < 2): ", n_low_expr, "/", top_n, " genes")
    add("- Moderate expression (2\u20136): ", n_mid_expr, "/", top_n, " genes")
    add("- High expression (> 6): ", n_high_expr, "/", top_n, " genes")
    add("")
    if (n_high_expr > top_n * 0.3) {
      add("A substantial number of discordant genes have high mean expression, ",
          "suggesting the discordance is not merely a low-expression noise artifact. ",
          "These genes warrant closer investigation.")
    } else if (n_low_expr > top_n * 0.6) {
      add("Most discordant genes have low expression, suggesting the disagreement ",
          "is largely driven by noise at low expression levels. These discordances ",
          "are less biologically concerning.")
    } else {
      add("Discordant genes span a range of expression levels, indicating both ",
          "noise-driven (low-expression) and genuine method-dependent differences.")
    }
    add("")
  }
}

# -----------------------------------------------
# Section 4: Ranking Stability
# -----------------------------------------------

add("## 4. Ranking Stability for Gene Groups of Interest")
add("")
add("For each gene group, genes are ranked by mean TPM per method (rank 1 = highest).")
add("Genes are flagged if their fractional rank change exceeds ",
    RANKING_CHANGE_THRESHOLD * 100, "% of the group size.")
add("")

for (gene_group in names(ranking_results)) {
  rdf <- ranking_results[[gene_group]]
  add("### 4.", which(names(ranking_results) == gene_group), " ", gene_group)
  add("")

  n_total <- nrow(rdf)
  n_flagged <- sum(rdf$Flagged)
  add("**Genes:** ", n_total, " | **Flagged (unstable):** ", n_flagged)
  add("")

  # Ranking table
  method_cols <- sapply(methods, get_short_name)
  display_cols <- intersect(method_cols, colnames(rdf))

  add("| Gene | Short Name | Median Rank | Range | SD | Flagged |",
      paste(display_cols, collapse = " | "), " | Mean TPM |")
  sep_cols <- paste(rep("---:", length(display_cols)), collapse = "|")
  add("|------|-----------|----------:|-----:|---:|---------|", sep_cols, "|--------:|")

  for (i in seq_len(nrow(rdf))) {
    row <- rdf[i, ]
    method_ranks <- sprintf("%.0f", as.numeric(row[display_cols]))
    flag_str <- if (row$Flagged) "**YES**" else ""
    add("| ", row$Gene_ID, " | ", row$Shortened_Name,
        " | ", row$Median_Rank,
        " | ", row$Rank_Range,
        " | ", row$Rank_SD,
        " | ", flag_str,
        " | ", paste(method_ranks, collapse = " | "),
        " | ", row$Mean_TPM, " |")
  }
  add("")

  add("![Ranking Heatmap](cross_method_concordance/figures/ranking_heatmap_", gene_group, ".png)")
  add("")
  add("#### Z-Score Normalized Expression (0\u201310 Scale)")
  add("")
  add("Per-gene Z-score scaled to 0\u201310 across methods. Blue (0) = method gives lowest estimate,")
  add("white (5) = average, red (10) = highest. Uniform rows indicate strong method agreement.")
  add("")
  add("![Z-Score Heatmap](cross_method_concordance/figures/zscore_heatmap_", gene_group, ".png)")
  add("")
  add("![Bump Chart](cross_method_concordance/figures/ranking_bump_chart_", gene_group, ".png)")
  add("")

  # ---- Data-driven interpretation of ranking stability per gene group ----
  add("#### Interpretation: ", gene_group)
  add("")

  pct_flagged <- 100 * n_flagged / n_total

  if (n_flagged == 0) {
    add("**All ", n_total, " genes have stable rankings across methods.** No gene's rank ",
        "changes by more than ", RANKING_CHANGE_THRESHOLD * 100, "% of the group size. ",
        "This indicates robust cross-method consistency for this gene group.")
  } else if (pct_flagged < 20) {
    add("**", n_flagged, "/", n_total, " genes (", sprintf("%.0f%%", pct_flagged),
        ") show unstable rankings.** Most genes in this group are consistently ranked ",
        "across methods, but a few show notable rank shifts.")
  } else if (pct_flagged < 50) {
    add("**", n_flagged, "/", n_total, " genes (", sprintf("%.0f%%", pct_flagged),
        ") show unstable rankings.** A substantial portion of this gene group is ",
        "ranked differently depending on which method is used. Biological conclusions ",
        "about relative expression levels within this group should be verified.")
  } else {
    add("**", n_flagged, "/", n_total, " genes (", sprintf("%.0f%%", pct_flagged),
        ") show unstable rankings.** The majority of genes in this group are ranked ",
        "differently across methods, suggesting high method sensitivity for this gene set. ",
        "Use caution when interpreting expression hierarchies within this group.")
  }
  add("")

  # Identify most/least stable genes
  if (n_total >= 3) {
    most_stable <- rdf[which.min(rdf$Rank_SD), ]
    least_stable <- rdf[which.max(rdf$Rank_SD), ]

    add("- **Most stable gene:** ", most_stable$Shortened_Name, " (",
        most_stable$Gene_ID, ") \u2014 median rank ", most_stable$Median_Rank,
        ", range ", most_stable$Rank_Range, ", SD ", most_stable$Rank_SD)
    add("- **Least stable gene:** ", least_stable$Shortened_Name, " (",
        least_stable$Gene_ID, ") \u2014 median rank ", least_stable$Median_Rank,
        ", range ", least_stable$Rank_Range, ", SD ", least_stable$Rank_SD)
    add("")

    # Analyze whether instability correlates with expression level
    if ("Mean_TPM" %in% colnames(rdf)) {
      mean_tpm_vals <- as.numeric(rdf$Mean_TPM)
      rank_sd_vals <- as.numeric(rdf$Rank_SD)
      if (length(mean_tpm_vals) >= 5 && sd(rank_sd_vals) > 0) {
        cor_test <- tryCatch(
          cor.test(log2(mean_tpm_vals + 1), rank_sd_vals, method = "spearman"),
          error = function(e) NULL
        )
        if (!is.null(cor_test)) {
          if (cor_test$p.value < 0.05 && cor_test$estimate < -0.3) {
            add("Rank instability is significantly correlated with lower expression levels ",
                "(Spearman \u03C1 = ", sprintf("%.2f", cor_test$estimate),
                ", p = ", sprintf("%.3g", cor_test$p.value),
                "), confirming that low-expression genes are harder to rank consistently.")
          } else if (cor_test$p.value < 0.05 && cor_test$estimate > 0.3) {
            add("Unexpectedly, rank instability is correlated with *higher* expression levels ",
                "(Spearman \u03C1 = ", sprintf("%.2f", cor_test$estimate),
                ", p = ", sprintf("%.3g", cor_test$p.value),
                "). This suggests method-specific biases affect highly-expressed genes in this group.")
          } else {
            add("Rank instability does not strongly correlate with expression level (",
                "Spearman \u03C1 = ", sprintf("%.2f", cor_test$estimate),
                ", p = ", sprintf("%.2g", cor_test$p.value),
                "), indicating that method disagreement on ranking is not simply a function of ",
                "expression magnitude.")
          }
          add("")
        }
      }
    }
  }
}

# -----------------------------------------------
# Section 5: Conclusions
# -----------------------------------------------

add("## 5. Summary and Recommendations")
add("")

# Compute overall concordance assessment
all_medians <- concordance$median_spearman[upper.tri(concordance$median_spearman)]
overall_median <- median(all_medians)

add("### Overall Concordance")
add("")
add("- **Median pairwise Spearman correlation:** ", sprintf("%.3f", overall_median))
if (overall_median > 0.9) {
  add("- **Assessment:** High overall concordance across methods.")
} else if (overall_median > 0.8) {
  add("- **Assessment:** Moderate concordance. Some method-specific biases observed.")
} else {
  add("- **Assessment:** Low concordance. Significant method-dependent differences in quantification.")
}
add("")

# Most/least concordant pairs
sp_upper <- concordance$median_spearman
diag(sp_upper) <- NA
sp_upper[lower.tri(sp_upper)] <- NA
best_pair <- which(sp_upper == max(sp_upper, na.rm = TRUE), arr.ind = TRUE)
worst_pair <- which(sp_upper == min(sp_upper, na.rm = TRUE), arr.ind = TRUE)

add("- **Most concordant pair:** ", rownames(sp_upper)[best_pair[1, 1]], " & ",
    colnames(sp_upper)[best_pair[1, 2]],
    " (rho = ", sprintf("%.3f", max(sp_upper, na.rm = TRUE)), ")")
add("- **Least concordant pair:** ", rownames(sp_upper)[worst_pair[1, 1]], " & ",
    colnames(sp_upper)[worst_pair[1, 2]],
    " (rho = ", sprintf("%.3f", min(sp_upper, na.rm = TRUE)), ")")
add("")

add("### Discordance Summary")
add("")
add("- ", concordance$n_discordant, " genes (",
    sprintf("%.1f%%", 100 * concordance$n_discordant / length(data$common_genes)),
    ") show high cross-method variability (CV > ", DISCORDANCE_CV_THRESHOLD, ").")
add("- These genes should be interpreted with caution in downstream analyses.")
add("- Full discordant gene list: [`tables/discordant_genes.csv`](cross_method_concordance/tables/discordant_genes.csv)")
add("- Full concordance scores: [`tables/all_genes_concordance.csv`](cross_method_concordance/tables/all_genes_concordance.csv)")
add("")

if (length(ranking_results) > 0) {
  add("### Ranking Stability Summary")
  add("")
  total_flagged <- sum(sapply(ranking_results, function(r) sum(r$Flagged)))
  total_genes <- sum(sapply(ranking_results, nrow))
  add("- Across all gene groups: ", total_flagged, "/", total_genes,
      " genes show unstable rankings across methods.")
  if (total_flagged > 0) {
    add("- Flagged genes: [`tables/ranking_instability_flagged.csv`](cross_method_concordance/tables/ranking_instability_flagged.csv)")
  }
  add("")
}

add("---")
add("*Generated by HeatSeq Cross-Method Concordance Analysis*")

# -----------------------------------------------
# Write report
# -----------------------------------------------

writeLines(lines, report_path)
cat("Report saved to:", report_path, "\n")

# -----------------------------------------------
# Write figure interpretation guide
# -----------------------------------------------

guide_path <- file.path(FIGURES_DIR, "figure_interpretation_guide.txt")
guide <- c(
"================================================================================",
" CROSS-METHOD CONCORDANCE FIGURES - INTERPRETATION GUIDE",
"================================================================================",
"",
"This guide explains how to read and interpret each figure produced by the",
"cross-method concordance analysis. All figures compare quantification results",
"across RNA-seq alignment/quantification methods:",
"",
paste0("  Methods compared: ", paste(short_names, collapse = ", ")),
paste0("  Reference genome: ", MASTER_REFERENCE),
paste0("  Common genes: ", length(data$common_genes)),
paste0("  Common samples: ", length(data$common_samples)),
"",
"",
"================================================================================",
"1. METHOD CONCORDANCE HEATMAPS",
"   Files: method_concordance_heatmap_spearman.png",
"          method_concordance_heatmap_pearson.png",
"================================================================================",
"",
"WHAT IT SHOWS:",
"  A symmetric matrix where each cell contains the median correlation between",
"  two methods, computed across all samples.",
"",
"HOW TO READ IT:",
"  - Each row/column is one method (short name labels on axes).",
"  - Cell values range from 0 to 1 (the diagonal is always 1.00).",
"  - Colors: blue = low correlation, red = high correlation.",
"  - Dendrograms on the top and left show hierarchical clustering of methods",
"    by similarity -- methods that cluster together produce the most similar",
"    quantification results.",
"",
"WHAT TO LOOK FOR:",
"  - Values > 0.95: Excellent agreement -- the two methods are nearly",
"    interchangeable for this dataset.",
"  - Values 0.85-0.95: Good agreement -- minor differences, usually at low-",
"    expression genes.",
"  - Values < 0.85: Notable divergence -- investigate which genes differ",
"    (see discordant genes heatmap).",
"  - Clustering patterns: genome-level methods (M1, M2, M3) often cluster",
"    together, separate from transcript-level methods (M4, M5).",
"",
"SPEARMAN vs PEARSON:",
"  - Spearman (rank-based): robust to outliers and non-linearity. Measures",
"    whether methods agree on the relative ordering of gene expression.",
"  - Pearson (linear): sensitive to absolute magnitude. High Pearson with",
"    low Spearman suggests agreement on highly-expressed genes but divergence",
"    in rank ordering of moderate/low-expression genes.",
"",
"",
"================================================================================",
"2. PER-SAMPLE CORRELATION BOXPLOT",
"   File: per_sample_correlation_boxplot.png",
"================================================================================",
"",
"WHAT IT SHOWS:",
"  Distribution of Spearman correlations for each method pair, where each",
"  data point is one sample's correlation between two methods.",
"",
"HOW TO READ IT:",
"  - X-axis: method pair labels (e.g., 'M1:HISAT2-Ref vs M3:STAR-Salmon').",
"  - Y-axis: Spearman rho for that sample.",
"  - Each box shows the interquartile range (IQR) of correlations across",
"    all samples for that method pair.",
"  - Red dashed line: rho = 0.90 (minimum acceptable agreement).",
"  - Green dotted line: rho = 0.95 (strong agreement threshold).",
"",
"WHAT TO LOOK FOR:",
"  - Boxes entirely above the green line (0.95): strong, consistent agreement.",
"  - Boxes between red and green lines: acceptable but variable agreement.",
"  - Any box or whisker extending below the red line (0.90): some samples",
"    show poor agreement -- investigate those samples for quality issues.",
"  - Outlier points below the boxes: individual samples where the two methods",
"    disagree significantly, possibly due to low sequencing depth, mapping",
"    artifacts, or method-specific biases.",
"  - Tight boxes (small IQR): consistent agreement across all samples.",
"  - Wide boxes: variable agreement -- some samples agree well, others poorly.",
"",
"",
"================================================================================",
"3. DISCORDANT GENES HEATMAP",
"   File: discordant_genes_heatmap.png",
"================================================================================",
"",
"WHAT IT SHOWS:",
"  The top 50 most discordant genes -- genes whose mean expression (log2 TPM+1)",
"  varies the most across methods, measured by the coefficient of variation (CV).",
"",
"HOW TO READ IT:",
"  - Rows: gene IDs (clustered by expression pattern similarity).",
"  - Columns: methods (fixed order, not clustered).",
"  - Cell color intensity: mean log2(TPM+1) expression for that gene in that",
"    method. Dark purple = high expression, light/white = low expression.",
"  - Right-side bar plot ('CV'): the CV of each gene's expression across",
"    methods. Taller orange bars = more disagreement between methods.",
"",
"WHAT TO LOOK FOR:",
"  - Genes with high CV (tall bars) but high expression: these are the most",
"    concerning -- methods disagree substantially on genes that should be",
"    well-quantified.",
"  - Genes with high CV but low expression: less concerning -- low-expression",
"    genes are inherently noisier and small absolute differences produce",
"    large CVs.",
"  - Patterns across methods: if a gene is consistently dark in genome-level",
"    methods (M1-M3) but light in transcript-level methods (M4-M5), this",
"    suggests a systematic quantification difference (e.g., multi-mapping",
"    handling, transcript vs gene-level aggregation).",
paste0("  - A gene is flagged as discordant when CV > ", DISCORDANCE_CV_THRESHOLD,
       " and mean expression above log2(", CONCORDANCE_MIN_EXPR, " + 1)."),
"",
"",
"================================================================================",
"3b. DISCORDANT GENES Z-SCORE HEATMAP",
"    File: discordant_genes_zscore_heatmap.png",
"================================================================================",
"",
"WHAT IT SHOWS:",
"  The same top 50 discordant genes as section 3, but with each gene's",
"  expression normalized to a 0-10 scale across methods using Z-scores.",
"  This removes absolute expression differences so you can directly compare",
"  how methods agree or disagree for every gene on a common scale.",
"",
"HOW THE SCALE WORKS:",
"  For each gene (row): compute Z = (value - row_mean) / row_sd across",
"  methods, then rescale to 0-10 where 0 = lowest method estimate and",
"  10 = highest method estimate. A score of 5.0 means that method's",
"  estimate equals the average across all methods.",
"",
"HOW TO READ IT:",
"  - Blue cells (0-3): that method produces a relatively LOW estimate for",
"    this gene compared to other methods.",
"  - White cells (~5): that method agrees with the cross-method average.",
"  - Red cells (7-10): that method produces a relatively HIGH estimate.",
"  - Each cell shows the 0-10 score value.",
"  - Right-side CV barplot: same as the raw heatmap (higher = more discord).",
"",
"WHAT TO LOOK FOR:",
"  - Columns that are systematically blue or red: that method consistently",
"    under- or over-estimates expression for discordant genes.",
"  - Rows with extreme blue-red contrast: genes where methods disagree most.",
"  - Rows that are mostly white (all ~5.0): genes where the discordance",
"    is relatively mild despite crossing the CV threshold.",
"  - Method family patterns: e.g., genome-level methods (M1-M3) grouping",
"    blue while transcript-level methods (M4-M5) group red, or vice versa.",
"",
"",
"================================================================================",
"4. PAIRWISE SCATTER PLOTS",
"   File: pairwise_scatter_plots.png",
"================================================================================",
"",
"WHAT IT SHOWS:",
"  For every pair of methods, a scatter plot of mean gene expression (log2 TPM+1)",
"  in method A (x-axis) vs method B (y-axis).",
"",
"HOW TO READ IT:",
"  - Each point is one gene. Purple points have nonzero expression; grey",
"    points are unexpressed in both methods.",
"  - Red dashed line: y = x (perfect agreement). Points on this line mean",
"    both methods produce identical expression estimates.",
"  - Title shows both Pearson (r) and Spearman (rho) correlations.",
"",
"WHAT TO LOOK FOR:",
"  - Points tightly along the diagonal: strong quantitative agreement.",
"  - Systematic offsets (cloud shifted above or below the diagonal): one",
"    method consistently estimates higher or lower expression. This can",
"    indicate differences in normalization or quantification approach.",
"  - Fan-shaped spread at low expression: expected -- low-expression genes",
"    are inherently noisy.",
"  - Outlier genes far from diagonal: candidates for method-specific",
"    artifacts. Cross-reference these with the discordant genes table.",
"  - Differences between method families:",
"      * M1 vs M2: tests the impact of reference-guided vs de novo assembly.",
"      * M1/M2 vs M3: tests HISAT2+StringTie vs STAR+Salmon approaches.",
"      * M3 vs M4: tests genome-aligned Salmon vs pseudo-alignment Salmon.",
"      * M4 vs M5: tests Salmon vs RSEM at the transcript level.",
"",
"",
"================================================================================",
"5. RANKING HEATMAPS (per gene group)",
"   Files: ranking_heatmap_<gene_group>.png",
"================================================================================",
"",
"WHAT IT SHOWS:",
"  How each method ranks the genes in a gene group by mean expression. Rank 1",
"  is the highest expressed gene.",
"",
"HOW TO READ IT:",
"  - Rows: genes (using shortened/display names), ordered by median rank.",
"  - Columns: methods.",
"  - Cell values: the rank assigned by that method (1 = highest expression).",
"  - Cell color: dark purple = high rank (highly expressed), light pink =",
"    low rank (lowly expressed).",
"  - Right-side annotations:",
"      * 'Stability' column: green = stable ranking across methods,",
"        orange = unstable (flagged for rank instability).",
"      * 'Rank_Range' bar plot: the difference between the highest and",
"        lowest rank across methods. Larger bars = more disagreement in",
"        where a gene falls in the expression hierarchy.",
"",
"WHAT TO LOOK FOR:",
"  - Uniform rows (same rank across all columns): all methods agree on",
"    this gene's relative expression level -- high confidence result.",
"  - Rows with large rank jumps (e.g., rank 2 in one method, rank 15 in",
"    another): methods disagree on how expressed this gene is relative to",
"    others in the group.",
paste0("  - A gene is flagged as unstable when its fractional rank change exceeds ",
       RANKING_CHANGE_THRESHOLD * 100, "% (rank range / total genes > ",
       RANKING_CHANGE_THRESHOLD, ")."),
"  - Consistently top-ranked genes (dark purple across all methods): the",
"    most reliably highly-expressed genes in the group.",
"  - Consistently bottom-ranked genes (light across all methods): reliably",
"    lowly-expressed, regardless of method.",
"",
"",
"================================================================================",
"5b. Z-SCORE EXPRESSION HEATMAPS (per gene group)",
"    Files: zscore_heatmap_<gene_group>.png",
"================================================================================",
"",
"WHAT IT SHOWS:",
"  Each gene's mean expression normalized to a 0-10 scale across methods.",
"  Unlike the ranking heatmap (which shows ordinal ranks), this shows the",
"  *magnitude* of cross-method differences on a common per-gene scale.",
"",
"HOW THE SCALE WORKS:",
"  For each gene: log2(mean_TPM + 1) is computed per method, then Z-scored",
"  across methods and rescaled so 0 = the method with the lowest estimate",
"  and 10 = the method with the highest estimate. A score of 5.0 means",
"  that method matches the cross-method average.",
"",
"HOW TO READ IT:",
"  - Blue cells (0-3): method gives a relatively low expression estimate.",
"  - White cells (~5): method matches the cross-method average.",
"  - Red cells (7-10): method gives a relatively high expression estimate.",
"  - Rows are in the same order as the ranking heatmap (by median rank).",
"  - Right-side Rank_Range barplot: same as the ranking heatmap.",
"  - Cell values show the 0-10 score to one decimal place.",
"",
"WHAT TO LOOK FOR:",
"  - Mostly white rows: strong agreement -- all methods give similar estimates.",
"  - Rows with strong blue-red contrast: methods disagree on this gene.",
"  - A column that is consistently blue or red across many genes: that method",
"    systematically under/overestimates this gene group.",
"  - Compare with the ranking heatmap: a rank change may reflect a tiny",
"    expression difference (narrow Z-range) or a large one (extreme blue-red).",
"",
"",
"================================================================================",
"6. RANKING BUMP CHARTS (per gene group)",
"   Files: ranking_bump_chart_<gene_group>.png",
"================================================================================",
"",
"WHAT IT SHOWS:",
"  A bump chart / slope graph tracking each gene's rank across methods.",
"  Each line is one gene, and horizontal position shows the method.",
"",
"HOW TO READ IT:",
"  - X-axis: methods (left to right).",
"  - Y-axis: rank (1 at top = highest expression, N at bottom = lowest).",
"  - Each colored line represents one gene across methods.",
"  - Solid thick lines: flagged genes (unstable ranking).",
"  - Dashed thin lines: stable genes (consistent ranking).",
"  - Legend (right side): maps line colors to gene names.",
"",
"WHAT TO LOOK FOR:",
"  - Horizontal lines: the gene maintains the same rank across all methods",
"    -- strong cross-method consistency.",
"  - Lines with steep slopes: the gene's ranking changes substantially",
"    between methods -- interpret with caution.",
"  - Crossing lines: genes that 'swap' positions between methods. If many",
"    lines cross between the same two methods, those methods disagree on",
"    the relative expression ordering within this gene group.",
"  - Clusters of stable lines at the top: genes that are reliably the",
"    highest expressed in the group, regardless of method.",
"  - One method with many crossings: that method may handle this gene",
"    group differently (e.g., multi-mapped reads, transcript isoforms).",
"",
"",
"================================================================================",
" GENERAL NOTES",
"================================================================================",
"",
"EXPRESSION UNITS:",
"  All comparisons use TPM (Transcripts Per Million), log2-transformed as",
"  log2(TPM + 1). This normalization puts all methods on a comparable scale",
"  and reduces the impact of highly-expressed outlier genes.",
"",
"THRESHOLDS USED:",
paste0("  - Minimum expression: ", CONCORDANCE_MIN_EXPR, " TPM (genes below this are treated as unexpressed)"),
paste0("  - Minimum samples: ", CONCORDANCE_MIN_SAMPLES, " (gene must be expressed in >= this many samples)"),
paste0("  - Discordance CV threshold: ", DISCORDANCE_CV_THRESHOLD, " (CV across methods above this flags a gene)"),
paste0("  - Ranking instability: ", RANKING_CHANGE_THRESHOLD * 100,
       "% (rank range / total genes > ", RANKING_CHANGE_THRESHOLD, " flags a gene)"),
"",
"WHAT DRIVES METHOD DISAGREEMENT:",
"  - Multi-mapping reads: methods handle ambiguous reads differently.",
"  - Gene vs transcript quantification: M1/M2 quantify at gene level,",
"    M4/M5 at transcript level (aggregated to gene level).",
"  - Assembly approach: M2 uses de novo assembly (may discover novel",
"    transcripts not in the reference GTF).",
"  - Alignment vs pseudo-alignment: M4 (Salmon) skips full alignment;",
"    M1/M2/M3/M5 perform read alignment to the genome or transcriptome.",
"  - Low-expression genes: inherently noisy in all methods -- small absolute",
"    differences produce large relative disagreements.",
"",
"RECOMMENDED ACTIONS:",
"  - If overall concordance is high (>0.95 Spearman): methods are reliable",
"    for this dataset. Use any method with confidence.",
"  - If specific genes are discordant: validate with orthogonal data (qPCR,",
"    proteomics) before drawing biological conclusions about those genes.",
"  - If one method consistently diverges: consider whether its assumptions",
"    match your experimental design (e.g., de novo assembly may not be",
"    appropriate if a high-quality reference genome is available).",
"  - For downstream differential expression: use methods with high pairwise",
"    concordance to ensure reproducibility of DEG lists.",
"================================================================================",
paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
)

writeLines(guide, guide_path)
cat("Figure interpretation guide saved to:", guide_path, "\n")

cat("\n[DONE] Report generation complete\n")
