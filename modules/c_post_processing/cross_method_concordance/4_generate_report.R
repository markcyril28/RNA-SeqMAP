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

add("### 2.5 Pairwise Scatter Plots")
add("")
add("Mean log2(TPM+1) per gene, comparing all method pairs. Red dashed line = y=x (perfect agreement).")
add("")
add("![Pairwise Scatters](cross_method_concordance/figures/pairwise_scatter_plots.png)")
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
  add("Full list: [`tables/discordant_genes.csv`](cross_method_concordance/tables/discordant_genes.csv)")
  add("")
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
  add("![Bump Chart](cross_method_concordance/figures/ranking_bump_chart_", gene_group, ".png)")
  add("")
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
cat("\n[DONE] Report generation complete\n")
