#!/usr/bin/env Rscript

# ===============================================
# GENE-VS-GENE CONCORDANCE (WITHIN GENE GROUPS)
# ===============================================
# For each gene group, computes a gene x gene Spearman correlation matrix
# using expression profiles across samples. This shows which genes within
# a group have concordant expression patterns.
#
# Each gene's "profile" is its expression vector across all samples.
# The gene x gene correlation heatmap reveals co-expression structure
# within each gene group.
#
# Outputs (per gene group):
#   - tables/gene_correlation_matrix_<group>.csv
#   - figures/gene_concordance_heatmap_<group>.png

# Skip re-sourcing when running under concordance_batch_dispatcher.R (already loaded)
if (!exists(".CONC_BATCH_CONFIG_LOADED") || !isTRUE(.CONC_BATCH_CONFIG_LOADED)) {
  source(file.path(Sys.getenv("CONCORDANCE_SCRIPT_DIR", "."), "0_concordance_config.R"))

  suppressPackageStartupMessages({
    library(ComplexHeatmap)
    library(circlize)
    library(grid)
  })
}

cat("\n=== STEP 2: Gene-vs-Gene Concordance (Within Gene Groups) ===\n\n")

# Load data
if (!file.exists(HARMONIZED_RDS)) {
  stop("Gene group data not found. Run 1_load_matrices_cross_gene_group.R first.")
}
data <- readRDS(HARMONIZED_RDS)

tpm_matrices <- data$tpm_matrices
groups <- names(tpm_matrices)
common_samples <- data$common_samples

cat("Gene groups:", paste(groups, collapse = ", "), "\n")
cat("Samples:", length(common_samples), "\n\n")

if (length(groups) == 0) {
  stop("No gene groups with data found.")
}

# Store per-group results
all_cor_matrices <- list()

# Cache directory listing once — O(D) scan reused across all gene groups.
# Avoids O(G × D) repeated list.files() calls for G gene groups.
.gene_group_csv_cache <- if (dir.exists(GENE_GROUPS_DIR)) {
  list.files(GENE_GROUPS_DIR, pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)
} else character(0)
# Pre-compute basenames once (O(C)) and build a name→paths hash map via split()
# for O(1) lookup per gene group instead of O(C) vectorized == scan.
.gene_group_csv_basenames <- basename(.gene_group_csv_cache)
.gene_group_csv_by_name <- split(.gene_group_csv_cache, .gene_group_csv_basenames)

# Hoist palette outside loop — the 100-color interpolation is constant across gene groups
.gg_palette <- colorRampPalette(c("#D73027", "#FEE090", "#E0F3F8", "#91BFDB", "#4575B4"))(100)

for (gg_name in groups) {
  cat("\n--- Processing gene group:", gg_name, "---\n")

  sub_mat <- tpm_matrices[[gg_name]]  # genes x samples
  n_genes <- nrow(sub_mat)
  cat("  Genes:", n_genes, "| Samples:", ncol(sub_mat), "\n")

  if (n_genes < 3) {
    cat("  [WARN] Too few genes (", n_genes, ") for meaningful correlation, skipping\n")
    next
  }

  # Guard: cor(t()) produces a G×G matrix — O(G²×S) time and O(G²) memory.
  # For gene groups > 5000, this becomes prohibitively expensive (25M+ cells).
  if (n_genes > 5000) {
    cat("  [WARN] Gene group '", gg_name, "' has ", n_genes,
        " genes — skipping full G×G correlation (would allocate ",
        round(n_genes^2 * 8 / 1e6, 1), " MB)\n")
    next
  }
  # Soft warning for 3000-5000 genes: O(G²) matrix is 72-200 MB.
  # Proceeding, but user should be aware of memory cost.
  if (n_genes > 3000) {
    cat("  [NOTE] Large gene group (", n_genes, " genes): G×G correlation will allocate ~",
        round(n_genes^2 * 8 / 1e6, 0), " MB\n")
  }

  # Use log2(TPM+1) to reduce skewness from highly expressed genes
  log2_mat <- log2(sub_mat + 1)

  # Compute gene x gene Spearman correlation matrix — O(G²×S)
  # Each row is a gene, each column is a sample — cor() correlates columns,
  # so transpose to get gene-vs-gene correlations across samples.
  # gpu_cor() rank-transforms on CPU then offloads the O(G²×S) matmul to GPU.
  cor_mat <- gpu_cor(t(log2_mat), method = "spearman")

  # Use Shortened_Name labels if available from the gene group CSV
  gene_labels <- rownames(cor_mat)
  # Check hash map first (O(1)) before falling back to file.exists probes.
  # On WSL2, each file.exists() costs 5-20ms due to cross-filesystem stat().
  gg_file <- NULL
  found <- .gene_group_csv_by_name[[paste0(gg_name, ".csv")]]
  if (!is.null(found) && length(found) > 0L && file.exists(found[1L])) {
    gg_file <- found[1L]
  } else {
    for (cand in c(file.path(GENE_GROUPS_DIR, paste0(gg_name, ".csv")),
                   file.path(GENE_GROUPS_DIR, gg_name))) {
      if (file.exists(cand)) { gg_file <- cand; break }
    }
  }
  if (!is.null(gg_file)) {
    gg_df <- .fast_read_csv(gg_file)
    if (all(c("Gene_ID", "Shortened_Name") %in% colnames(gg_df))) {
      # Strip suffixes from Gene_ID to match matrix rownames
      gg_df$Gene_ID_clean <- sub("(\\.[0-9]+){1,2}$", "", trimws(gg_df$Gene_ID))
      name_map <- setNames(trimws(gg_df$Shortened_Name), gg_df$Gene_ID_clean)
      mapped <- name_map[gene_labels]
      # Use short name where available, fall back to gene ID.
      # Direct index assignment avoids ifelse() full-vector allocation.
      valid <- !is.na(mapped) & nzchar(mapped)
      gene_labels[valid] <- mapped[valid]
    }
  }

  rownames(cor_mat) <- gene_labels
  colnames(cor_mat) <- gene_labels

  all_cor_matrices[[gg_name]] <- cor_mat

  # Save correlation matrix
  safe_name <- gsub("[^[:alnum:]_.-]", "_", gg_name)
  # Use as.data.table(keep.rownames=) when available — avoids O(G²) data.frame() copy
  .cor_csv_path <- file.path(TABLES_DIR, paste0("gene_correlation_matrix_", safe_name, ".csv"))
  if (.conc_use_dt) {
    .cor_out <- data.table::as.data.table(cor_mat, keep.rownames = "Gene")
    data.table::fwrite(.cor_out, .cor_csv_path)
  } else {
    .cor_out <- data.frame(Gene = rownames(cor_mat), cor_mat, check.names = FALSE)
    write.csv(.cor_out, .cor_csv_path, row.names = FALSE)
  }
  rm(.cor_out, .cor_csv_path)

  # Compute off-diagonal mask once — reused for range stats, NA cleanup, and min_cor.
  # Avoids 4 redundant O(G²) row()/col() calls (was computed 4× independently).
  .off_diag_mask <- row(cor_mat) != col(cor_mat)

  cat("  Correlation matrix: ", n_genes, "x", n_genes, "\n")
  cat("  Range:", round(min(cor_mat[.off_diag_mask], na.rm = TRUE), 3),
      "to", round(max(cor_mat[.off_diag_mask], na.rm = TRUE), 3), "\n")

  # -----------------------------------------------
  # Generate concordance heatmap
  # -----------------------------------------------

  cat("  Generating heatmap...\n")

  cor_mat[!is.finite(cor_mat) & .off_diag_mask] <- NA
  off_diag <- cor_mat[.off_diag_mask]
  min_cor <- min(off_diag, na.rm = TRUE)
  if (!is.finite(min_cor) || min_cor >= 1) min_cor <- -1

  col_fun <- colorRamp2(
    seq(min_cor, 1, length.out = 100),
    .gg_palette
  )

  # Only show cell values if the matrix is small enough to be readable
  cell_fun <- if (n_genes <= 30) {
    function(j, i, x, y, width, height, fill) {
      val <- cor_mat[i, j]
      label <- if (is.finite(val)) sprintf("%.2f", val) else ""
      grid.text(label, x, y, gp = gpar(fontsize = max(6, 12 - n_genes / 5), fontface = "bold"))
    }
  } else {
    NULL
  }

  # Scale heatmap size with gene count — compute once, reuse for unit() and layout
  .hm_body_val <- max(10, min(30, n_genes * 0.6))
  hm_size <- unit(.hm_body_val, "cm")
  font_size <- max(6, 13 - n_genes / 5)

  # Pre-compute distance matrix for clustering — gpu_dist() offloads the O(G²×G)
  # Euclidean distance computation to GPU when available, falls back to CPU otherwise.
  # cor_mat is symmetric so row and column distances are identical.
  .cl_dist <- gpu_dist(cor_mat)

  ht <- Heatmap(cor_mat,
    name = "Spearman",
    col = col_fun,
    cell_fun = cell_fun,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    clustering_distance_rows = .cl_dist,
    clustering_distance_columns = .cl_dist,
    show_row_dend = (n_genes <= 50),
    show_column_dend = (n_genes <= 50),
    row_names_gp = gpar(fontsize = font_size),
    column_names_gp = gpar(fontsize = font_size),
    column_names_rot = 45,
    column_title = paste0("Gene-vs-Gene Concordance: ", gg_name,
                          "\n(", FIXED_METHOD, " | ", MASTER_REFERENCE, ")"),
    column_title_gp = gpar(fontsize = 14, fontface = "bold"),
    heatmap_legend_param = list(
      title = "Spearman\nCorrelation\n(Within-Group\nGene Co-expression)",
      legend_height = unit(5, "cm")
    ),
    width = hm_size,
    height = hm_size
  )

  # Auto-calculate figure dimensions from content (reuses .hm_body_val computed above)
  .fig_layout <- calc_figure_layout(
    row_labels = rownames(cor_mat),
    col_labels = colnames(cor_mat),
    col_rot = 45, hm_body_cm = c(.hm_body_val, .hm_body_val),
    has_dendro = (n_genes <= 50), has_title = TRUE,
    legend_width_cm = 5, font_size = font_size
  )

  .dev_open <- FALSE
  tryCatch({
    png(file.path(FIGURES_DIR, paste0("gene_concordance_heatmap_", safe_name, ".png")),
        width = .fig_layout$width, height = .fig_layout$height, res = FIGURE_DPI)
    .dev_open <- TRUE
    draw(ht, padding = .fig_layout$padding)
    dev.off()
    .dev_open <- FALSE
    cat("  Saved: gene_concordance_heatmap_", safe_name, ".png\n")
  }, error = function(e) {
    if (.dev_open) try(dev.off(), silent = TRUE)
    cat("  Error generating heatmap:", e$message, "\n")
  })
}

# Save results for report
concordance_results <- list(
  gene_cor_matrices = all_cor_matrices,
  gene_group_sizes = data$gene_group_sizes,
  tpm_matrices = tpm_matrices
)

saveRDS(concordance_results, file.path(OUTPUT_DIR, "concordance_results.rds"))
cat("\n[DONE] Gene-vs-gene concordance results saved\n")
