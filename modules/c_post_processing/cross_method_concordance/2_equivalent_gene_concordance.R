#!/usr/bin/env Rscript

# ===============================================
# EQUIVALENT-GENE CONCORDANCE ACROSS GENOMES
# ===============================================
# For a fixed method, compares expression of positionally-equivalent genes
# across two (or more) reference genomes. Each gene's expression vector across
# samples is correlated between genomes to assess per-gene quantification
# agreement.
#
# Requires: Step 1 (1_load_matrices_cross_genome.R) to have built the
# positional orthology mapping and produced HARMONIZED_RDS.
#
# Big O: O(G x S) where G=genes, S=samples.
#
# Outputs:
#   - tables/per_gene_cross_genome_correlation.csv
#   - tables/genome_expression_comparison.csv
#   - figures/equivalent_gene_concordance_heatmap.png
#   - figures/equivalent_gene_scatter.png

# Skip re-sourcing when running under concordance_batch_dispatcher.R (already loaded)
if (!exists(".CONC_BATCH_CONFIG_LOADED") || !isTRUE(.CONC_BATCH_CONFIG_LOADED)) {
  source(file.path(Sys.getenv("CONCORDANCE_SCRIPT_DIR", "."), "0_concordance_config.R"))

  suppressPackageStartupMessages({
    library(ComplexHeatmap)
    library(circlize)
    library(grid)
  })
}

# Reuse cached probe from 1_utility_functions.R when available; fall back to requireNamespace
if (!exists(".HAS_MATRIXSTATS")) .HAS_MATRIXSTATS <- requireNamespace("matrixStats", quietly = TRUE)

cat("\n=== STEP 2: Equivalent-Gene Concordance Across Genomes ===\n\n")

# Load harmonized data (produced by 1_load_matrices_cross_genome.R with positional mapping)
if (!file.exists(HARMONIZED_RDS)) {
  stop("Harmonized data not found. Run 1_load_matrices_cross_genome.R first.")
}
data <- readRDS(HARMONIZED_RDS)
tpm_matrices <- data$tpm_matrices
common_genes <- data$common_genes
common_samples <- data$common_samples
ortho_gene_ids <- data$ortho_gene_ids          # per-genome original Gene_IDs (may be NULL)
ortho_short_names <- data$ortho_short_names    # per-genome Shortened_Names (may be NULL)
ortho_gene_groups <- data$ortho_gene_groups    # common_label -> gene group name (may be NULL)
genomes <- names(tpm_matrices)
n_genomes <- length(genomes)

short_names <- vapply(genomes, get_short_name, character(1))
cat("Genomes:", paste(short_names, collapse = ", "), "\n")
cat("Equivalent genes:", length(common_genes), "| Samples:", length(common_samples), "\n\n")

if (n_genomes < 2) {
  stop("Need at least 2 genomes for equivalent-gene concordance. Found: ", n_genomes)
}
if (length(common_genes) < 2) {
  stop("Need at least 2 common genes for equivalent-gene concordance. Found: ", length(common_genes))
}

# Build per-genome display labels using Shortened_Names from each genome's CSV.
# Format: "ShortGenome:ShortenedName" (e.g., "GPE001970:SmelDMP5_01.730")
genome_gene_labels <- list()
for (g in genomes) {
  if (!is.null(ortho_short_names) && g %in% names(ortho_short_names)) {
    sn <- ortho_short_names[[g]]
    genome_gene_labels[[g]] <- setNames(sn[common_genes], common_genes)
  } else {
    genome_gene_labels[[g]] <- setNames(common_genes, common_genes)
  }
}

# -----------------------------------------------
# 2.1 Per-gene Spearman correlation across genomes
# -----------------------------------------------
# For each gene, correlate its expression vector across samples between
# each pair of genomes. This measures how consistently the gene is
# quantified across different reference assemblies.

cat("--- Computing per-gene cross-genome correlations ---\n")

genome_pairs <- combn(genomes, 2, simplify = FALSE)
pair_labels <- vapply(genome_pairs, function(p) {
  paste(short_names[p[1]], "vs", short_names[p[2]])
}, character(1))

# Pre-compute log2(TPM+1) for correlation stability
log2_matrices <- lapply(tpm_matrices, function(mat) log2(mat + 1))

# Per-gene correlation: for each gene, correlate its sample vector between genomes
gene_cor_list <- vector("list", length(genome_pairs))

for (pi in seq_along(genome_pairs)) {
  g1 <- genome_pairs[[pi]][1]
  g2 <- genome_pairs[[pi]][2]
  mat1 <- log2_matrices[[g1]]
  mat2 <- log2_matrices[[g2]]

  # Vectorized per-gene Spearman: row-wise rank transform + vectorized Pearson.
  # Avoids G sequential cor() calls (each with R-level overhead) by computing
  # rank transforms once then using matrix algebra for Pearson-on-ranks.
  # O(G × S × log(S)) for rank, O(G × S) for correlation — same asymptotic
  # complexity but ~10-50x fewer R function call overheads.
  m1_sub <- mat1[common_genes, , drop = FALSE]
  m2_sub <- mat2[common_genes, , drop = FALSE]
  nonzero_mask <- (m1_sub > 0) | (m2_sub > 0)  # logical G×S
  too_few <- rowSums(nonzero_mask) < 3
  # Mask non-expressed entries as NA for rank computation
  m1_sub[!nonzero_mask] <- NA
  m2_sub[!nonzero_mask] <- NA
  # Row-wise rank transform (each gene ranked across samples)
  # matrixStats::rowRanks is a C-level loop — avoids G R-level apply() calls.
  # .HAS_MATRIXSTATS is set at module load (line 31) — no exists() guard needed
  if (.HAS_MATRIXSTATS) {
    r1 <- matrixStats::rowRanks(m1_sub, ties.method = "average")
    r2 <- matrixStats::rowRanks(m2_sub, ties.method = "average")
    # rowRanks drops NAs differently — re-apply NA mask from source matrices
    r1[is.na(m1_sub)] <- NA
    r2[is.na(m2_sub)] <- NA
    dimnames(r1) <- dimnames(m1_sub)
    dimnames(r2) <- dimnames(m2_sub)
  } else {
    r1 <- t(apply(m1_sub, 1, rank, na.last = "keep"))
    r2 <- t(apply(m2_sub, 1, rank, na.last = "keep"))
  }
  # Vectorized Pearson on ranks = Spearman (row-wise centering + dot product)
  r1_mean <- rowMeans(r1, na.rm = TRUE)
  r2_mean <- rowMeans(r2, na.rm = TRUE)
  # Guard: all-NA rows produce NaN means → replace with 0 to prevent NaN propagation
  r1_mean[!is.finite(r1_mean)] <- 0
  r2_mean[!is.finite(r2_mean)] <- 0
  r1c <- r1 - r1_mean; r1c[is.na(r1c)] <- 0
  r2c <- r2 - r2_mean; r2c[is.na(r2c)] <- 0
  num <- rowSums(r1c * r2c)
  # r*r avoids ^ S3 method dispatch on matrix
  den <- sqrt(rowSums(r1c*r1c) * rowSums(r2c*r2c))
  den[den == 0] <- 1
  cors <- setNames(num / den, common_genes)
  cors[too_few] <- NA_real_

  gene_cor_list[[pi]] <- cors
}

# Build per-gene results table in one cbind — avoids O(P+G) copy-on-modify from
# repeated gene_results[[ ]] <- assignments (each copies the entire data.frame)
.spearman_cols <- setNames(
  lapply(seq_along(genome_pairs), function(pi) gene_cor_list[[pi]]),
  paste0("Spearman_", pair_labels)
)
.mean_tpm_cols <- setNames(
  lapply(genomes, function(g) rowMeans(tpm_matrices[[g]][common_genes, , drop = FALSE])),
  paste0("Mean_TPM_", short_names[genomes])
)
gene_results <- cbind(
  data.frame(Gene = common_genes, stringsAsFactors = FALSE, check.names = FALSE),
  as.data.frame(.spearman_cols, stringsAsFactors = FALSE, check.names = FALSE),
  as.data.frame(.mean_tpm_cols, stringsAsFactors = FALSE, check.names = FALSE)
)

# Add fold change (genome2 / genome1 mean TPM, log2)
if (n_genomes == 2) {
  mean1 <- gene_results[[paste0("Mean_TPM_", short_names[genomes[1]])]]
  mean2 <- gene_results[[paste0("Mean_TPM_", short_names[genomes[2]])]]
  gene_results$Log2_FoldChange <- log2((mean2 + 0.01) / (mean1 + 0.01))
}

# Sort by correlation (ascending — worst agreement first)
cor_col <- paste0("Spearman_", pair_labels[1])
gene_results <- gene_results[order(gene_results[[cor_col]], na.last = TRUE), ]

if (.conc_use_dt) {
  data.table::fwrite(gene_results, file.path(TABLES_DIR, "per_gene_cross_genome_correlation.csv"))
} else {
  write.csv(gene_results, file.path(TABLES_DIR, "per_gene_cross_genome_correlation.csv"), row.names = FALSE)
}

# Print summary (cap at 50 genes to avoid flooding console for large gene sets)
# O(min(G, 50)) output lines instead of unbounded O(G)
.n_print <- min(nrow(gene_results), 50L)
cat("  Per-gene correlations")
if (nrow(gene_results) > .n_print) cat(" (showing first ", .n_print, " of ", nrow(gene_results), ")")
cat(":\n")
.genes_print <- gene_results$Gene[seq_len(.n_print)]
.cors_print <- gene_results[[cor_col]][seq_len(.n_print)]
.cors_str <- ifelse(is.finite(.cors_print), sprintf("%.3f", .cors_print), "N/A")
# Vectorized paste + single cat call: O(G) string ops but 1 I/O call instead of G
cat(paste0("    ", .genes_print, " : ", .cors_str, "\n"), sep = "")
cat("  Saved: per_gene_cross_genome_correlation.csv\n\n")

# -----------------------------------------------
# 2.2 Expression comparison table (gene x sample for each genome)
# -----------------------------------------------

cat("--- Building expression comparison table ---\n")

# Vectorized construction: O(n_genomes) vectorized ops instead of O(G × S) R-level loop.
# Replaces per-cell data.frame assignment + O(n²) list accumulation pattern.
comparison_dfs <- lapply(genomes, function(g) {
  mat <- tpm_matrices[[g]][common_genes, common_samples, drop = FALSE]
  # Single-step construction avoids intermediate as.data.frame copy
  data.frame(Gene = common_genes, Genome = short_names[g], mat,
             check.names = FALSE, stringsAsFactors = FALSE)
})
comparison_df <- if (.conc_use_dt) {
  data.table::setDF(data.table::rbindlist(comparison_dfs, use.names = TRUE))
} else {
  do.call(rbind, comparison_dfs)
}

if (.conc_use_dt) {
  data.table::fwrite(comparison_df, file.path(TABLES_DIR, "genome_expression_comparison.csv"))
} else {
  write.csv(comparison_df, file.path(TABLES_DIR, "genome_expression_comparison.csv"), row.names = FALSE)
}
cat("  Saved: genome_expression_comparison.csv\n\n")

# -----------------------------------------------
# 2.3 & 2.4  Per-gene-group heatmaps and scatter plots
# -----------------------------------------------
# Generate separate figures for each gene group.
# Heatmap: rows = genome A genes (Shortened_Names), cols = genome B genes
# Scatter: one panel per equivalent gene pair

if (n_genomes == 2) {
  g1 <- genomes[1]; g2 <- genomes[2]
  precomputed_cors <- gene_cor_list[[1]]

  # Determine gene groups to iterate over
  if (!is.null(ortho_gene_groups) && length(ortho_gene_groups) > 0) {
    group_names <- unique(ortho_gene_groups[common_genes])
    group_names <- group_names[!is.na(group_names)]
  } else {
    # No group info — treat all genes as one group
    group_names <- "all_genes"
    ortho_gene_groups <- setNames(rep("all_genes", length(common_genes)), common_genes)
  }

  cat("--- Generating per-gene-group figures (", length(group_names), "groups ) ---\n")

  # Hoist color palette outside loop — colorRampPalette() does interpolation
  # setup once; palette(100) generates 100 colors. Both are loop-invariant.
  .pal_100 <- colorRampPalette(c("#D73027", "#F46D43", "#FDAE61", "#FEE090",
                                   "#FFFFBF",
                                   "#E0F3F8", "#ABD9E9", "#74ADD1", "#4575B4"))(100)

  # Cache FIXED_METHOD short name — loop-invariant, avoids G get_short_name() calls
  .fixed_method_short <- get_short_name(FIXED_METHOD)

  # Shared rank helper — .HAS_MATRIXSTATS set at module load (line 31)
  .rank_rows <- function(m) {
    if (.HAS_MATRIXSTATS) {
      r <- matrixStats::rowRanks(m, ties.method = "average")
      dimnames(r) <- dimnames(m)
    } else {
      r <- t(apply(m, 1, rank))
    }
    r - rowMeans(r)
  }

  for (.grp in group_names) {
    .grp_genes <- common_genes[ortho_gene_groups[common_genes] == .grp]
    if (length(.grp_genes) < 2) {
      cat("  [INFO] Skipping group '", .grp, "' — fewer than 2 genes\n")
      next
    }

    # Sanitize group name for filenames
    .grp_tag <- gsub("[^[:alnum:]_.-]", "_", .grp)
    .n_grp <- length(.grp_genes)
    cat("\n  Gene group:", .grp, "(", .n_grp, "genes )\n")

    # --- 2.3  Gene-vs-gene correlation heatmap ---

    mat1 <- log2_matrices[[g1]][.grp_genes, , drop = FALSE]
    mat2 <- log2_matrices[[g2]][.grp_genes, , drop = FALSE]
    row_labels <- genome_gene_labels[[g1]][.grp_genes]
    col_labels <- genome_gene_labels[[g2]][.grp_genes]

    r1 <- .rank_rows(mat1)
    r2 <- .rank_rows(mat2)
    num <- gpu_matmult(r1, t(r2))
    # r*r avoids ^ S3 method dispatch on matrix
    den <- sqrt(rowSums(r1*r1)) %o% sqrt(rowSums(r2*r2))
    den[den == 0] <- 1
    gene_gene_cor <- num / den
    rownames(gene_gene_cor) <- row_labels
    colnames(gene_gene_cor) <- col_labels

    # Save table
    .gg_out <- data.frame(Gene = rownames(gene_gene_cor), gene_gene_cor, check.names = FALSE)
    .tbl_file <- paste0("gene_vs_gene_correlation_", .grp_tag, ".csv")
    if (.conc_use_dt) {
      data.table::fwrite(.gg_out, file.path(TABLES_DIR, .tbl_file))
    } else {
      write.csv(.gg_out, file.path(TABLES_DIR, .tbl_file), row.names = FALSE)
    }
    rm(.gg_out)

    # Color scale — uses pre-computed .pal_100 (hoisted above loop)
    .abs_lim <- max(abs(range(gene_gene_cor, na.rm = TRUE)), 0.5)
    col_fun <- colorRamp2(seq(-.abs_lim, .abs_lim, length.out = 100), .pal_100)

    # Capture gene_gene_cor in cell_fun closure
    .local_cor <- gene_gene_cor
    .local_n <- .n_grp
    cell_fun <- function(j, i, x, y, width, height, fill) {
      val <- .local_cor[i, j]
      label <- if (is.finite(val)) sprintf("%.2f", val) else ""
      grid.text(label, x, y, gp = gpar(fontsize = max(6, min(11, 120 / .local_n)),
                                         fontface = if (i == j) "bold" else "plain"))
    }

    .hm_sz <- min(20, max(8, .n_grp * 1.8))
    ht <- Heatmap(gene_gene_cor,
      name = "Spearman",
      col = col_fun,
      cell_fun = cell_fun,
      cluster_rows = FALSE, cluster_columns = FALSE,
      show_row_dend = FALSE, show_column_dend = FALSE,
      row_names_gp = gpar(fontsize = 11),
      column_names_gp = gpar(fontsize = 11),
      column_names_rot = 45,
      row_title = short_names[g1],
      row_title_gp = gpar(fontsize = 14, fontface = "bold"),
      column_title = short_names[g2],
      column_title_gp = gpar(fontsize = 14, fontface = "bold"),
      column_title_side = "bottom",
      heatmap_legend_param = list(
        title = paste0("Spearman\nCorrelation\n(Cross-Genome\nExpression)"),
        legend_height = unit(5, "cm")
      ),
      width = unit(.hm_sz, "cm"), height = unit(.hm_sz, "cm")
    )

    .fig_layout <- calc_figure_layout(
      row_labels = row_labels, col_labels = col_labels,
      col_rot = 45, hm_body_cm = c(.hm_sz, .hm_sz),
      has_dendro = FALSE, has_title = TRUE,
      legend_width_cm = 5, font_size = 11
    )
    .pad <- as.numeric(.fig_layout$padding)
    .pad[2] <- .pad[2] + 15  # left +15mm for row title
    .pad[3] <- .pad[3] + 20  # top +10mm base + 10mm for figure title

    # Figure title: gene group + method
    .grp_display <- gsub("_", " ", .grp)
    .fig_title <- paste0(.grp_display, "\nEquivalent Gene Concordance (",
                         .fixed_method_short, ")")

    # Compute final pixel dimensions from total padding
    .extra_w_px <- ceiling(15 / 10 * FIGURE_DPI / 2.54)
    .extra_h_px <- ceiling(20 / 10 * FIGURE_DPI / 2.54)

    .hm_file <- paste0("equivalent_gene_heatmap_", .grp_tag, ".png")
    .dev_open <- FALSE
    tryCatch({
      png(file.path(FIGURES_DIR, .hm_file),
          width = .fig_layout$width + .extra_w_px,
          height = .fig_layout$height + .extra_h_px,
          res = FIGURE_DPI)
      .dev_open <- TRUE
      draw(ht, padding = unit(.pad, "mm"),
           column_title = .fig_title,
           column_title_gp = gpar(fontsize = 13, fontface = "bold"))
      dev.off(); .dev_open <- FALSE
      cat("    Saved:", .hm_file, "\n")
    }, error = function(e) {
      if (.dev_open) try(dev.off(), silent = TRUE)
      cat("    Error generating heatmap:", e$message, "\n")
    })

    # --- 2.4  Per-gene scatter plots ---

    n_cols_sc <- min(4, .n_grp)
    n_rows_sc <- ceiling(.n_grp / n_cols_sc)
    panel_w <- 3; panel_h <- 3

    .sc_file <- paste0("equivalent_gene_scatter_", .grp_tag, ".png")
    .dev_open <- FALSE
    tryCatch({
      png(file.path(FIGURES_DIR, .sc_file),
          width = n_cols_sc * panel_w, height = n_rows_sc * panel_h + 1,
          units = "in", res = FIGURE_DPI)
      .dev_open <- TRUE

      par(mfrow = c(n_rows_sc, n_cols_sc),
          mar = c(4, 4, 3, 1), oma = c(0, 0, 3, 0))

      .xlab <- paste0(short_names[g1], " log2(TPM+1)")
      .ylab <- paste0(short_names[g2], " log2(TPM+1)")
      .pt_col <- adjustcolor("#4575B4", 0.7)

      for (gene in .grp_genes) {
        # Reuse pre-computed log2_matrices (line 87) instead of redundant log2(tpm+1)
        lx <- log2_matrices[[g1]][gene, common_samples]
        ly <- log2_matrices[[g2]][gene, common_samples]
        rho <- precomputed_cors[gene]
        rho_str <- if (is.finite(rho)) sprintf("rho=%.2f", rho) else "rho=N/A"

        # Use Shortened_Names from both genomes for the panel title
        .lbl_g1 <- genome_gene_labels[[g1]][gene]
        .lbl_g2 <- genome_gene_labels[[g2]][gene]
        .panel_title <- if (!is.na(.lbl_g1) && !is.na(.lbl_g2) && .lbl_g1 != .lbl_g2) {
          paste0(.lbl_g1, " / ", .lbl_g2, "\n", rho_str)
        } else {
          paste0(ifelse(!is.na(.lbl_g1), .lbl_g1, gene), "\n", rho_str)
        }

        plot(lx, ly, xlab = .xlab, ylab = .ylab,
             main = .panel_title,
             pch = 19, col = .pt_col, cex = 1.2, cex.main = 0.9)
        abline(0, 1, col = "grey50", lty = 2)
        .ok <- is.finite(lx) & is.finite(ly)
        if (sum(.ok) >= 3) {
          .lx <- lx[.ok]; .ly <- ly[.ok]
          .slope <- cov(.lx, .ly) / var(.lx)
          .intercept <- mean(.ly) - .slope * mean(.lx)
          abline(.intercept, .slope, col = "#D73027", lwd = 1.5)
        }
      }

      remaining <- (n_cols_sc * n_rows_sc) - .n_grp
      if (remaining > 0) for (i in seq_len(remaining)) plot.new()

      mtext(paste0(.grp_display, ": ", short_names[g1], " vs ", short_names[g2],
                    " (", .fixed_method_short, ")"),
            outer = TRUE, cex = 1.1, font = 2)

      dev.off(); .dev_open <- FALSE
      cat("    Saved:", .sc_file, "\n")
    }, error = function(e) {
      if (.dev_open) try(dev.off(), silent = TRUE)
      cat("    Error generating scatter:", e$message, "\n")
    })
  }
} else {
  cat("  [INFO] Per-gene-group figures require exactly 2 genomes; skipping\n")
}

# -----------------------------------------------
# Save concordance results for report generation
# -----------------------------------------------

concordance_results <- list(
  gene_correlations = gene_results,
  per_pair_cors = gene_cor_list,
  pair_labels = pair_labels,
  genomes = genomes,
  short_names = short_names,
  common_genes = common_genes,
  common_samples = common_samples
)

saveRDS(concordance_results, file.path(OUTPUT_DIR, "concordance_results.rds"))
cat("\n[DONE] Equivalent-gene concordance results saved\n")
