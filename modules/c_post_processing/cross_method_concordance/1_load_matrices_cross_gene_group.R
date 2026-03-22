#!/usr/bin/env Rscript

# ===============================================
# CROSS-GENE-GROUP CONCORDANCE - LOAD & SUBSET MATRICES
# ===============================================
# For a fixed method+genome, loads the full TPM matrix, then for each gene
# group produces a submatrix and per-sample median expression summary.
# This enables gene-group-vs-gene-group concordance analysis.
#
# Output: HARMONIZED_RDS containing:
#   $tpm_matrices    - named list: one genes x samples submatrix per gene group
#   $summary_vectors - matrix: gene_groups x samples (per-sample median TPM)
#   $common_samples  - sample IDs
#   $method_stats    - data.frame with per-group stats
#   $gene_group_sizes - named integer vector of gene counts per group

source(file.path(Sys.getenv("CONCORDANCE_SCRIPT_DIR", "."), "0_concordance_config.R"))
source(file.path(Sys.getenv("CONCORDANCE_SCRIPT_DIR", "."), "0_method_loaders.R"))

cat("\n=== STEP 1: Loading Expression Matrices (Cross-Gene-Group Mode) ===\n\n")

if (!nzchar(FIXED_METHOD)) {
  stop("FIXED_METHOD must be set for cross-gene-group mode")
}
if (length(CONCORDANCE_GENE_GROUPS) < 2) {
  stop("Need at least 2 gene groups for cross-gene-group concordance. Found: ",
       length(CONCORDANCE_GENE_GROUPS))
}

cat("Fixed method:", FIXED_METHOD, "\n")
cat("Reference:", MASTER_REFERENCE, "\n")
cat("Gene groups to compare:", paste(CONCORDANCE_GENE_GROUPS, collapse = ", "), "\n\n")

# Cache directory listing once — O(D) scan reused across all gene groups.
# Avoids O(G × D) repeated list.files() calls for G gene groups.
.gene_group_csv_cache <- if (dir.exists(GENE_GROUPS_DIR)) {
  list.files(GENE_GROUPS_DIR, pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)
} else character(0)
# Pre-compute basenames once (O(C)) instead of per-group basename() call (O(G×C))
.gene_group_csv_basenames <- basename(.gene_group_csv_cache)

# matrixStats::colMedians is a C-level column-wise median — ~3x faster than apply(x, 2, median)
.HAS_MATRIXSTATS <- requireNamespace("matrixStats", quietly = TRUE)

# -----------------------------------------------
# Load the full TPM matrix for the fixed method
# -----------------------------------------------

ref_dir <- get_method_ref_dir(FIXED_METHOD)
cat("--- Loading full TPM matrix for", FIXED_METHOD, "(ref:", ref_dir, ") ---\n")

full_matrix <- tryCatch(load_method_for_ref(FIXED_METHOD, ref_dir), error = function(e) {
  cat("  [ERROR]", e$message, "\n")
  NULL
})

if (is.null(full_matrix) || nrow(full_matrix) == 0 || ncol(full_matrix) == 0) {
  stop("Failed to load TPM matrix for ", FIXED_METHOD, " with reference ", ref_dir)
}

# Strip transcript suffixes
rn <- rownames(full_matrix)
new_rn <- sub("(\\.[0-9]+){1,2}$", "", rn)
if (!identical(new_rn, rn)) {
  if (any(duplicated(new_rn))) {
    full_matrix <- rowsum(full_matrix, group = new_rn, reorder = FALSE)
  } else {
    rownames(full_matrix) <- new_rn
  }
}

cat("  Full matrix:", nrow(full_matrix), "genes x", ncol(full_matrix), "samples\n")

# -----------------------------------------------
# Load gene group definitions and subset
# -----------------------------------------------

cat("\n--- Loading gene group definitions ---\n")
cat("  Gene groups dir:", GENE_GROUPS_DIR, "\n")

tpm_matrices <- list()
summary_mat <- matrix(NA_real_, nrow = length(CONCORDANCE_GENE_GROUPS),
                      ncol = ncol(full_matrix),
                      dimnames = list(CONCORDANCE_GENE_GROUPS, colnames(full_matrix)))
gene_group_sizes <- setNames(integer(length(CONCORDANCE_GENE_GROUPS)), CONCORDANCE_GENE_GROUPS)
stats_list <- list()

for (gg_name in CONCORDANCE_GENE_GROUPS) {
  cat("\n  Processing gene group:", gg_name, "\n")

  # Find the gene group CSV file using cached directory listing (O(G) basename match
  # on cached listing instead of O(D) filesystem rescan per group)
  gg_file <- NULL
  candidates <- c(
    file.path(GENE_GROUPS_DIR, paste0(gg_name, ".csv")),
    file.path(GENE_GROUPS_DIR, gg_name)
  )
  # Also search subdirectories via cached listing
  found <- .gene_group_csv_cache[.gene_group_csv_basenames == paste0(gg_name, ".csv")]
  candidates <- c(candidates, found)
  for (cand in candidates) {
    if (file.exists(cand)) { gg_file <- cand; break }
  }

  if (is.null(gg_file)) {
    cat("    [WARN] Gene group file not found for:", gg_name, "\n")
    next
  }

  cat("    File:", gg_file, "\n")
  gg_df <- .fast_read_csv(gg_file)

  # Get gene IDs (support Gene_ID and Shortened_Name columns)
  if ("Gene_ID" %in% colnames(gg_df)) {
    gene_ids <- unique(trimws(gg_df$Gene_ID))
  } else {
    gene_ids <- unique(trimws(gg_df[[1]]))
  }
  # Strip suffixes to match matrix
  gene_ids <- sub("(\\.[0-9]+){1,2}$", "", gene_ids)
  gene_ids <- gene_ids[nzchar(gene_ids)]

  # Match against full matrix
  matched <- intersect(gene_ids, rownames(full_matrix))
  cat("    Genes defined:", length(gene_ids), "| Matched in matrix:", length(matched), "\n")

  if (length(matched) < 3) {
    cat("    [WARN] Too few matched genes (", length(matched), "), skipping\n")
    next
  }

  # Subset matrix
  sub_mat <- full_matrix[matched, , drop = FALSE]
  tpm_matrices[[gg_name]] <- sub_mat
  gene_group_sizes[gg_name] <- length(matched)

  # Per-sample median expression (summary vector for correlation)
  # matrixStats::colMedians uses C-level implementation (~3x faster than apply+median)
  summary_mat[gg_name, ] <- if (.HAS_MATRIXSTATS) {
    matrixStats::colMedians(sub_mat, na.rm = TRUE)
  } else {
    apply(sub_mat, 2, median, na.rm = TRUE)
  }

  stats_list[[gg_name]] <- data.frame(
    method = gg_name,
    short_name = gg_name,
    n_genes_raw = length(gene_ids),
    n_samples_raw = ncol(sub_mat),
    n_genes_harmonized = length(matched),
    n_samples_harmonized = ncol(sub_mat),
    stringsAsFactors = FALSE
  )
}

# Remove gene groups with no data
valid_groups <- names(tpm_matrices)
summary_mat <- summary_mat[valid_groups, , drop = FALSE]
gene_group_sizes <- gene_group_sizes[valid_groups]

if (length(tpm_matrices) < 2) {
  stop("Need at least 2 gene groups with data. Found: ", length(tpm_matrices))
}

# Use rbindlist when data.table available — avoids O(R²) copy overhead of do.call(rbind)
method_stats <- if (exists(".use_dt") && .use_dt) {
  data.table::setDF(data.table::rbindlist(stats_list, use.names = TRUE, fill = TRUE))
} else {
  do.call(rbind, stats_list)
}

cat("\n--- Loaded", length(tpm_matrices), "gene groups ---\n")
for (gg in names(tpm_matrices)) {
  cat("  ", gg, ":", nrow(tpm_matrices[[gg]]), "genes x", ncol(tpm_matrices[[gg]]), "samples\n")
}

# -----------------------------------------------
# Save data
# -----------------------------------------------

result <- list(
  tpm_matrices = tpm_matrices,
  summary_vectors = summary_mat,
  common_samples = colnames(full_matrix),
  common_genes = character(0),  # gene groups don't share genes
  method_stats = method_stats,
  gene_group_sizes = gene_group_sizes
)

saveRDS(result, HARMONIZED_RDS)
cat("\n[DONE] Gene group data saved to:", HARMONIZED_RDS, "\n")
if (.conc_use_dt) {
  data.table::fwrite(method_stats, file.path(TABLES_DIR, "gene_group_stats.csv"))
} else {
  write.csv(method_stats, file.path(TABLES_DIR, "gene_group_stats.csv"), row.names = FALSE)
}
