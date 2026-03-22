#!/usr/bin/env Rscript

# ===============================================
# CROSS-METHOD CONCORDANCE - REPORT GENERATOR
# ===============================================
# Generates a unified Markdown report with embedded figure references,
# correlation matrices, gene lists, and analysis summaries.
#
# Output: REPORT_BASE/concordance_report.md (or POST_PROC_BASE fallback)

source(file.path(Sys.getenv("CONCORDANCE_SCRIPT_DIR", "."), "0_concordance_config.R"))

cat("\n=== STEP 4: Generating Concordance Report ===\n\n")

# Load all results
if (!file.exists(HARMONIZED_RDS)) {
  stop("Harmonized data not found: ", HARMONIZED_RDS, "\n  Run 1_load_matrices.R first.")
}
concordance_rds <- file.path(OUTPUT_DIR, "concordance_results.rds")
if (!file.exists(concordance_rds)) {
  stop("Concordance results not found: ", concordance_rds, "\n  Run 2_quantification_concordance.R first.")
}
data <- readRDS(HARMONIZED_RDS)
concordance <- readRDS(concordance_rds)
# Pre-compute path to avoid redundant file.path() call
ranking_rds_path <- file.path(OUTPUT_DIR, "ranking_results.rds")
ranking_results <- if (file.exists(ranking_rds_path)) {
  readRDS(ranking_rds_path)
} else {
  list()
}

methods <- names(data$tpm_matrices)
short_names <- vapply(methods, get_short_name, character(1))

# Report output path (separate from POST_PROC_BASE inputs)
report_base <- Sys.getenv("REPORT_BASE", POST_PROC_BASE)
report_path <- file.path(report_base, "concordance_report.md")

# Mode-aware labels for report headers
.item_label <- switch(CONCORDANCE_MODE,
  cross_method          = "Method",
  cross_genome          = "Genome",
  cross_gene_group      = "Gene Group",
  cross_equivalent_gene = "Equivalent Gene",
  "Item"
)

# -----------------------------------------------
# Build report
# -----------------------------------------------

# Use list accumulation instead of c() concatenation to avoid O(n²) reallocation.
# Dynamic growth: doubles capacity when full (amortized O(1) per append).
.line_capacity <- 500L
lines <- vector("list", .line_capacity)
.line_idx <- 0L
add <- function(...) {
  .line_idx <<- .line_idx + 1L
  if (.line_idx > .line_capacity) {
    .line_capacity <<- .line_capacity * 2L
    length(lines) <<- .line_capacity
  }
  lines[[.line_idx]] <<- paste0(...)
}

add("# Cross-", .item_label, " Concordance Report")
add("")
add("**Reference genome:** ", MASTER_REFERENCE)
add("**Date generated:** ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
add("**", .item_label, "s compared:** ", paste(short_names, collapse = ", "))
if (length(data$common_genes) > 0) add("**Common genes:** ", length(data$common_genes))
add("**Common samples:** ", length(data$common_samples))
add("")

# -----------------------------------------------
# Section 1: Method Overview
# -----------------------------------------------

add("## 1. ", .item_label, " Overview")
add("")
add("| ", .item_label, " | Short Name | Genes (Raw) | Samples (Raw) | Genes (Harmonized) | Samples (Harmonized) |")
add("|--------|-----------|-------------|--------------|-------------------|---------------------|")
# Vectorized row construction: O(R) paste0 + single lapply instead of per-row add()
.ms <- data$method_stats
.ms_rows <- paste0("| ", .ms$method, " | ", .ms$short_name, " | ",
                   format(.ms$n_genes_raw, big.mark = ","), " | ", .ms$n_samples_raw, " | ",
                   format(.ms$n_genes_harmonized, big.mark = ","), " | ", .ms$n_samples_harmonized, " |")
invisible(lapply(.ms_rows, add))
add("")

# Gene set overlap summary (only when gene_sets_raw is available)
gene_sets <- data$gene_sets_raw
if (!is.null(gene_sets)) {
method_names <- names(gene_sets)
add("### Gene Set Overlaps")
add("")
add("| ", .item_label, " Pair | Shared Genes |")
add("|------------|-------------|")
# O(M²) pairwise method comparison where M = number of methods.
# Pre-compute short name map to avoid O(M²) get_short_name() calls.
.sn_map <- setNames(vapply(method_names, get_short_name, character(1)), method_names)
for (i in seq_len(length(method_names) - 1)) {
  for (j in seq(i + 1, length(method_names))) {
    overlap <- length(intersect(gene_sets[[method_names[i]]], gene_sets[[method_names[j]]]))
    add("| ", .sn_map[method_names[i]], " & ", .sn_map[method_names[j]],
        " | ", format(overlap, big.mark = ","), " |")
  }
}
add("")
}  # end gene_sets_raw guard

# -----------------------------------------------
# Section 2: Quantification Concordance
# -----------------------------------------------

# Sections 2-4 require median_spearman from 2_quantification_concordance.R.
# cross_gene_group mode uses 2_gene_group_concordance.R which produces
# gene_cor_matrices instead — skip the method-pair concordance sections.
.has_spearman <- !is.null(concordance$median_spearman)

if (.has_spearman) {

add("## 2. Quantification Concordance")
add("")

add("### 2.1 Median Spearman Correlation Matrix")
add("")
add("Pairwise median Spearman correlations across all samples (computed on log2(TPM+1) of expressed genes).")
add("")

# Format correlation matrix as markdown table
sp_mat <- concordance$median_spearman
add("| Method |", paste(colnames(sp_mat), collapse = " | "), " |")
add("|", strrep("------|", ncol(sp_mat) + 1), "|")
# Pre-format entire matrix at once: O(R×C) total, then vectorized row construction
.sp_mat_str <- matrix(sprintf("%.3f", sp_mat), nrow = nrow(sp_mat), dimnames = dimnames(sp_mat))
# Vectorized row-paste: single C-level call vs O(R) per-row paste dispatch
.sp_rows <- paste0("| **", rownames(sp_mat), "** | ",
                   do.call(paste, c(as.data.frame(.sp_mat_str, stringsAsFactors = FALSE), sep = " | ")), " |")
invisible(lapply(.sp_rows, add))
add("")

add("### 2.2 Concordance Heatmap")
add("")
add("![Spearman Concordance](figures/method_concordance_heatmap_spearman.png)")
add("")

# ---- Data-driven interpretation of concordance heatmap ----
sp_upper_vals <- concordance$median_spearman[upper.tri(concordance$median_spearman)]
sp_overall <- median(sp_upper_vals, na.rm = TRUE)
sp_min <- min(sp_upper_vals, na.rm = TRUE)
sp_max <- max(sp_upper_vals, na.rm = TRUE)

# Find best/worst pairs — upper-triangle mask computed once and reused in Section 4
# (avoids redundant O(M²) copy + mask at lines ~399-401)
.sp_upper_masked <- concordance$median_spearman
diag(.sp_upper_masked) <- NA
.sp_upper_masked[lower.tri(.sp_upper_masked)] <- NA

sp_full_max <- max(.sp_upper_masked, na.rm = TRUE)
sp_full_min <- min(.sp_upper_masked, na.rm = TRUE)

sp_best_name <- "N/A"
sp_worst_name <- "N/A"
if (is.finite(sp_full_max)) {
  sp_best_idx <- which(.sp_upper_masked == sp_full_max, arr.ind = TRUE)
  if (nrow(sp_best_idx) > 0) {
    sp_best_name <- paste0(rownames(.sp_upper_masked)[sp_best_idx[1,1]], " & ", colnames(.sp_upper_masked)[sp_best_idx[1,2]])
  }
}
if (is.finite(sp_full_min)) {
  sp_worst_idx <- which(.sp_upper_masked == sp_full_min, arr.ind = TRUE)
  if (nrow(sp_worst_idx) > 0) {
    sp_worst_name <- paste0(rownames(.sp_upper_masked)[sp_worst_idx[1,1]], " & ", colnames(.sp_upper_masked)[sp_worst_idx[1,2]])
  }
}

add("#### Interpretation")
add("")

# Overall assessment (guard against NA from all-NA pairwise correlations)
if (!is.finite(sp_overall)) {
  add("**Overall concordance could not be assessed** — all pairwise correlations are NA. ",
      "Check that at least two methods have overlapping expressed genes.")
} else if (sp_overall >= 0.95) {
  add("**Overall concordance is excellent.** The median pairwise Spearman correlation is ",
      sprintf("%.3f", sp_overall), " (range: ", sprintf("%.3f", sp_min), "\u2013",
      sprintf("%.3f", sp_max), "), indicating that all methods produce highly consistent ",
      "expression rankings.")
} else if (sp_overall >= 0.90) {
  add("**Overall concordance is good.** The median pairwise Spearman correlation is ",
      sprintf("%.3f", sp_overall), " (range: ", sprintf("%.3f", sp_min), "\u2013",
      sprintf("%.3f", sp_max), "). Most method pairs agree well on gene expression rankings, ",
      "though some pairs show moderate divergence.")
} else if (sp_overall >= 0.85) {
  add("**Overall concordance is moderate.** The median pairwise Spearman correlation is ",
      sprintf("%.3f", sp_overall), " (range: ", sprintf("%.3f", sp_min), "\u2013",
      sprintf("%.3f", sp_max), "). Notable method-specific biases are present. ",
      "Downstream results should be validated across multiple methods.")
} else {
  add("**Overall concordance is low.** The median pairwise Spearman correlation is only ",
      sprintf("%.3f", sp_overall), " (range: ", sprintf("%.3f", sp_min), "\u2013",
      sprintf("%.3f", sp_max), "). Substantial disagreement exists between methods. ",
      "Results should be interpreted with caution and validated with orthogonal data.")
}
add("")

# Best/worst pairs
if (is.finite(sp_overall)) {
  add("- **Most concordant:** ", sp_best_name, " (\u03C1 = ",
      sprintf("%.3f", sp_full_max), ")")
  add("- **Least concordant:** ", sp_worst_name, " (\u03C1 = ",
      sprintf("%.3f", sp_full_min), ")")
}
add("")

# Clustering pattern analysis
# Replace any remaining NAs in the distance matrix with 1.0 (max dissimilarity)
# to prevent hclust from crashing when pairwise correlations are unavailable.
sp_dist_mat <- 1 - concordance$median_spearman
sp_dist_mat[is.na(sp_dist_mat)] <- 1.0
sp_dist <- as.dist(sp_dist_mat)
sp_hclust <- hclust(sp_dist, method = "complete")
sp_clusters <- cutree(sp_hclust, h = median(sp_dist))
n_clusters <- length(unique(sp_clusters))

if (n_clusters > 1) {
  add("**Clustering pattern:** Methods form ", n_clusters, " distinct clusters based on ",
      "Spearman correlation similarity:")
  add("")
  # O(K) where K = number of clusters (K ≤ M methods)
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

# Count how many pairs exceed key thresholds (NA values excluded from counts)
n_pairs_total <- sum(!is.na(sp_upper_vals))
n_above_95 <- sum(sp_upper_vals >= 0.95, na.rm = TRUE)
n_above_90 <- sum(sp_upper_vals >= 0.90, na.rm = TRUE)
n_below_90 <- sum(sp_upper_vals < 0.90, na.rm = TRUE)

add("**Threshold summary (Spearman):**")
add("")
add("- ", n_above_95, "/", n_pairs_total, " pairs with \u03C1 \u2265 0.95 (excellent agreement)")
add("- ", n_above_90 - n_above_95, "/", n_pairs_total, " pairs with 0.90 \u2264 \u03C1 < 0.95 (good agreement)")
if (n_below_90 > 0) {
  add("- ", n_below_90, "/", n_pairs_total, " pairs with \u03C1 < 0.90 (needs investigation)")
}
add("")

# -----------------------------------------------
# Section 3: Ranking Stability
# -----------------------------------------------

add("## 3. Ranking Stability for Gene Groups of Interest")
add("")
add("For each gene group, genes are ranked by mean TPM per method (rank 1 = highest).")
add("Genes are flagged if their fractional rank change exceeds ",
    RANKING_CHANGE_THRESHOLD * 100, "% of the group size.")
add("")

for (.gi in seq_along(ranking_results)) {
  gene_group <- names(ranking_results)[.gi]
  rdf <- ranking_results[[.gi]]
  add("### 3.", .gi, " ", gene_group)
  add("")

  n_total <- nrow(rdf)
  n_flagged <- sum(rdf$Flagged)
  add("**Genes:** ", n_total, " | **Flagged (unstable):** ", n_flagged)
  add("")

  # Ranking table
  method_cols <- vapply(methods, get_short_name, character(1))
  display_cols <- intersect(method_cols, colnames(rdf))

  add("| Gene | Short Name | Median Rank | Range | SD | Flagged | ",
      paste(display_cols, collapse = " | "), " | Mean TPM |")
  sep_cols <- paste(rep.int("---:", length(display_cols)), collapse = " | ")
  add("|------|-----------|----------:|-----:|---:|---------| ", sep_cols, " |--------:|")

  # Pre-extract columns as vectors: O(C + R) vs O(R × C) for rdf[i,] subsetting
  .gene_ids <- rdf$Gene_ID
  .short_names <- rdf$Shortened_Name
  .med_ranks <- rdf$Median_Rank
  .rank_ranges <- rdf$Rank_Range
  .rank_sds <- rdf$Rank_SD
  .flagged <- ifelse(rdf$Flagged, "**YES**", "")
  .mean_tpms <- rdf$Mean_TPM
  .rank_mat <- as.matrix(rdf[, display_cols, drop = FALSE])
  # Pre-format entire matrix: vectorized row construction avoids O(G) per-row paste
  .rank_mat_str <- matrix(sprintf("%.0f", .rank_mat), nrow = nrow(.rank_mat),
                          dimnames = dimnames(.rank_mat))
  # Vectorized row-paste: single C-level call vs O(G) per-row paste dispatch
  .rank_col_str <- do.call(paste, c(as.data.frame(.rank_mat_str, stringsAsFactors = FALSE), sep = " | "))
  .table_rows <- paste0("| ", .gene_ids, " | ", .short_names,
                        " | ", .med_ranks,
                        " | ", .rank_ranges,
                        " | ", .rank_sds,
                        " | ", .flagged,
                        " | ", .rank_col_str,
                        " | ", .mean_tpms, " |")
  invisible(lapply(.table_rows, add))
  add("")

  add("![Ranking Heatmap](figures/ranking_heatmap_", gene_group, ".png)")
  add("")
  add("#### Z-Score Normalized Expression (0\u201310 Scale)")
  add("")
  add("Per-gene Z-score scaled to 0\u201310 across methods. Blue (0) = method gives lowest estimate,")
  add("white (5) = average, red (10) = highest. Uniform rows indicate strong method agreement.")
  add("")
  add("![Z-Score Heatmap](figures/zscore_heatmap_", gene_group, ".png)")
  add("")
  add("![Bump Chart](figures/ranking_bump_chart_", gene_group, ".png)")
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
      # Columns are already numeric from round() in 3_ranking_stability.R — skip coercion
      mean_tpm_vals <- rdf$Mean_TPM
      rank_sd_vals <- rdf$Rank_SD
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
# Section 4: Conclusions
# -----------------------------------------------

add("## 4. Summary and Recommendations")
add("")

# Compute overall concordance assessment
all_medians <- concordance$median_spearman[upper.tri(concordance$median_spearman)]
overall_median <- median(all_medians, na.rm = TRUE)

add("### Overall Concordance")
add("")
if (!is.finite(overall_median)) {
  add("- **Median pairwise Spearman correlation:** N/A (all correlations unavailable)")
  add("- **Assessment:** Could not assess concordance.")
} else {
  add("- **Median pairwise Spearman correlation:** ", sprintf("%.3f", overall_median))
  if (overall_median > 0.9) {
    add("- **Assessment:** High overall concordance across methods.")
  } else if (overall_median > 0.8) {
    add("- **Assessment:** Moderate concordance. Some method-specific biases observed.")
  } else {
    add("- **Assessment:** Low concordance. Significant method-dependent differences in quantification.")
  }
}
add("")

# Most/least concordant pairs — reuse .sp_upper_masked computed in Section 2
if (is.finite(sp_full_max) && is.finite(sp_full_min)) {
  best_pair <- which(.sp_upper_masked == sp_full_max, arr.ind = TRUE)
  worst_pair <- which(.sp_upper_masked == sp_full_min, arr.ind = TRUE)
  add("- **Most concordant pair:** ", rownames(.sp_upper_masked)[best_pair[1, 1]], " & ",
      colnames(.sp_upper_masked)[best_pair[1, 2]],
      " (rho = ", sprintf("%.3f", sp_full_max), ")")
  add("- **Least concordant pair:** ", rownames(.sp_upper_masked)[worst_pair[1, 1]], " & ",
      colnames(.sp_upper_masked)[worst_pair[1, 2]],
      " (rho = ", sprintf("%.3f", sp_full_min), ")")
} else {
  add("- Best/worst concordant pairs could not be determined (insufficient valid correlations)")
}
add("")

add("### Recommendations")
add("")
if (is.finite(overall_median) && overall_median > 0.9) {
  add("- All methods produce largely consistent results; any single method can be used with confidence.")
} else {
  add("- Consider using consensus results from multiple methods for higher confidence.")
}
add("")

if (length(ranking_results) > 0) {
  add("### Ranking Stability Summary")
  add("")
  total_flagged <- sum(vapply(ranking_results, function(r) sum(r$Flagged), integer(1)))
  total_genes <- sum(vapply(ranking_results, nrow, integer(1)))
  add("- Across all gene groups: ", total_flagged, "/", total_genes,
      " genes show unstable rankings across methods.")
  if (total_flagged > 0) {
    add("- Flagged genes: [`tables/ranking_instability_flagged.csv`](tables/ranking_instability_flagged.csv)")
  }
  add("")
}

} else if (CONCORDANCE_MODE == "cross_equivalent_gene") {
  # cross_equivalent_gene mode: per-gene correlations across genome pairs
  add("## 2. Equivalent-Gene Concordance")
  add("")
  gene_cors <- concordance$gene_correlations
  pair_labels <- concordance$pair_labels
  if (!is.null(gene_cors) && length(gene_cors) > 0) {
    add("Per-gene Spearman correlations between equivalent genes across genome pairs.")
    add("")
    for (pi in seq_along(gene_cors)) {
      .label <- if (!is.null(pair_labels) && pi <= length(pair_labels)) pair_labels[pi] else paste("Pair", pi)
      .cors <- gene_cors[[pi]]
      .finite <- .cors[is.finite(.cors)]
      if (length(.finite) == 0) next
      add("### ", .label)
      add("")
      add("- **Genes compared:** ", length(.cors))
      add("- **Median Spearman:** ", sprintf("%.3f", median(.finite)))
      add("- **Range:** ", sprintf("%.3f", min(.finite)), " to ", sprintf("%.3f", max(.finite)))
      add("")
    }
    # Link to heatmap figure if it exists
    add("![Equivalent Gene Concordance](cross_method_concordance/figures/equivalent_gene_heatmap.png)")
    add("")
  } else {
    add("No equivalent-gene concordance results available.")
    add("")
  }

} else {
  # cross_gene_group mode: report on gene-vs-gene correlation matrices
  add("## 2. Gene-vs-Gene Concordance")
  add("")
  gene_cor_matrices <- concordance$gene_cor_matrices
  if (!is.null(gene_cor_matrices) && length(gene_cor_matrices) > 0) {
    for (gg_name in names(gene_cor_matrices)) {
      safe_name <- gsub("[^[:alnum:]_.-]", "_", gg_name)
      cor_mat <- gene_cor_matrices[[gg_name]]
      n_genes <- nrow(cor_mat)
      off_diag <- cor_mat[row(cor_mat) != col(cor_mat)]
      add("### ", gg_name)
      add("")
      add("- **Genes:** ", n_genes)
      add("- **Median pairwise Spearman:** ", sprintf("%.3f", median(off_diag, na.rm = TRUE)))
      add("- **Range:** ", sprintf("%.3f", min(off_diag, na.rm = TRUE)),
          " to ", sprintf("%.3f", max(off_diag, na.rm = TRUE)))
      add("")
      add("![Gene Concordance](cross_method_concordance/figures/gene_concordance_heatmap_", safe_name, ".png)")
      add("")
    }
  } else {
    add("No gene-vs-gene concordance results available.")
    add("")
  }
}  # end .has_spearman guard

add("---")
add("*Generated by HeatSeq Cross-", .item_label, " Concordance Analysis*")

# -----------------------------------------------
# Write report
# -----------------------------------------------

# Flatten list to character vector (discard unused pre-allocated slots)
lines <- unlist(lines[seq_len(.line_idx)])
writeLines(lines, report_path)
cat("Report saved to:", report_path, "\n")

# -----------------------------------------------
# Write figure interpretation guide
# -----------------------------------------------

# Figure interpretation guide is cross_method-specific (references method-specific figures)
if (.has_spearman) {
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
"1. MEDIAN SPEARMAN CONCORDANCE HEATMAP",
"   File: method_concordance_heatmap_spearman.png",
"================================================================================",
"",
"WHAT IT SHOWS:",
"  A symmetric matrix where each cell contains the median Spearman correlation",
"  between two methods, computed across all samples on log2(TPM+1) of expressed",
"  genes. Higher values indicate better agreement in gene expression rankings.",
"",
"HOW TO READ IT:",
"  - Each row/column is one method (short name labels on axes).",
"  - Cell values range from 0 to 1 (the diagonal is always 1.00).",
"  - Colors: warm (orange/yellow) = lower correlation, cool (blue) = higher.",
"  - Dendrograms on the top and left show hierarchical clustering of methods",
"    by similarity -- methods that cluster together produce the most similar",
"    quantification results.",
"",
"WHAT TO LOOK FOR:",
"  - Values > 0.95: Excellent agreement -- the two methods are nearly",
"    interchangeable for this dataset.",
"  - Values 0.85-0.95: Good agreement -- minor differences, usually at low-",
"    expression genes.",
"  - Values < 0.85: Notable divergence -- investigate which genes differ.",
"  - Clustering patterns: genome-level methods (M1, M2, M3) often cluster",
"    together, separate from transcript-level methods (M4, M5).",
"",
"",
"================================================================================",
"2. RANKING HEATMAPS (per gene group)",
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
"2b. Z-SCORE EXPRESSION HEATMAPS (per gene group)",
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
"3. RANKING BUMP CHARTS (per gene group)",
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
}  # end .has_spearman guard for figure guide

cat("\n[DONE] Report generation complete\n")
