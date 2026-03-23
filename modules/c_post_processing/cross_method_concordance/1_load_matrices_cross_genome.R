#!/usr/bin/env Rscript

# ===============================================
# CROSS-GENOME CONCORDANCE - LOAD & HARMONIZE MATRICES
# ===============================================
# For a fixed alignment method, loads TPM matrices from multiple reference
# genomes and harmonizes gene IDs and samples. This enables genome-vs-genome
# concordance analysis (how do expression estimates agree across references?).
#
# Output: HARMONIZED_RDS containing:
#   $tpm_matrices  - named list of gene x sample TPM matrices (one per genome)
#   $common_genes  - character vector of genes present in all genomes
#   $common_samples - character vector of samples present in all genomes
#   $method_stats  - data.frame with per-genome stats

source(file.path(Sys.getenv("CONCORDANCE_SCRIPT_DIR", "."), "0_concordance_config.R"))
source(file.path(Sys.getenv("CONCORDANCE_SCRIPT_DIR", "."), "0_method_loaders.R"))

cat("\n=== STEP 1: Loading Expression Matrices (Cross-Genome Mode) ===\n\n")

if (!nzchar(FIXED_METHOD)) {
  stop("FIXED_METHOD must be set for cross-genome mode")
}
if (length(CONCORDANCE_GENOMES) < 2) {
  stop("Need at least 2 genomes for cross-genome concordance. Found: ",
       length(CONCORDANCE_GENOMES))
}

cat("Fixed method:", FIXED_METHOD, "\n")
cat("Genomes to compare:", paste(CONCORDANCE_GENOMES, collapse = ", "), "\n\n")

# -----------------------------------------------
# Load fixed method across all genomes
# -----------------------------------------------

tpm_matrices <- list()
stats_list <- vector("list", length(CONCORDANCE_GENOMES))
stats_idx <- 0L

for (genome in CONCORDANCE_GENOMES) {
  cat("\nLoading", FIXED_METHOD, "for genome:", genome, "...\n")

  mat <- tryCatch(load_method_for_ref(FIXED_METHOD, genome), error = function(e) {
    cat("  [ERROR]", e$message, "\n")
    NULL
  })

  if (!is.null(mat) && nrow(mat) > 0 && ncol(mat) > 0) {
    tpm_matrices[[genome]] <- mat
    stats_idx <- stats_idx + 1L
    stats_list[[stats_idx]] <- data.frame(
      method = genome,
      short_name = get_short_name(genome),
      n_genes_raw = nrow(mat),
      n_samples_raw = ncol(mat),
      stringsAsFactors = FALSE
    )
    cat(" [", genome, "] Loaded:", nrow(mat), "genes x", ncol(mat), "samples\n")
  } else {
    cat("  [WARN] No data loaded for genome:", genome, "\n")
  }
}
stats_list <- stats_list[seq_len(stats_idx)]
method_stats <- if (.use_dt) {
  data.table::setDF(data.table::rbindlist(stats_list, use.names = TRUE, fill = TRUE))
} else {
  do.call(rbind, stats_list)
}

if (length(tpm_matrices) < 2) {
  stop("Need at least 2 genomes with data for cross-genome concordance. Found: ",
       length(tpm_matrices))
}

cat("\n--- Loaded", length(tpm_matrices), "genomes ---\n")

# -----------------------------------------------
# Harmonize gene IDs across genomes
# -----------------------------------------------

cat("\n--- Harmonizing gene IDs ---\n")

# Strip transcript suffixes from gene IDs in each genome's matrix.
# lapply avoids copy-on-modify: for-loop `tpm_matrices[[g]] <- mat` triggers a
# full list copy on each assignment; lapply builds the new list in one pass.
tpm_matrices <- lapply(tpm_matrices, function(mat) {
  rn <- rownames(mat)
  new_rn <- sub("(\\.[0-9]+){1,2}$", "", rn)
  if (identical(new_rn, rn)) return(mat)
  if (any(duplicated(new_rn))) return(rowsum(mat, group = new_rn, reorder = FALSE))
  rownames(mat) <- new_rn
  mat
})

gene_sets <- lapply(tpm_matrices, rownames)
n_genomes <- length(gene_sets)

# -----------------------------------------------
# Build positional orthology mapping from gene group CSVs
# -----------------------------------------------
# When genomes use different gene ID schemas (e.g., SMEL4.1_* vs SMEL5_*),
# direct ID intersection yields zero common genes. Instead, use per-genome
# gene group CSVs where row N in genome A corresponds to row N in genome B
# (positional equivalence). This builds a mapping to common labels.

.use_positional_mapping <- FALSE

if (length(GENOME_GENE_GROUPS_MAP) >= 2) {
  # Verify all comparison genomes have a mapping entry
  .mapped_genomes <- intersect(names(tpm_matrices), names(GENOME_GENE_GROUPS_MAP))
  if (length(.mapped_genomes) >= 2) {
    cat("  Building positional orthology mapping from gene group CSVs...\n")

    # For each gene group slot (by position), load CSVs from each genome
    # and build a cross-genome gene ID mapping using row position.
    .ref_genome <- .mapped_genomes[1]  # first genome provides common labels
    .n_groups <- length(GENOME_GENE_GROUPS_MAP[[.ref_genome]])

    # Resolve gene groups directory per genome
    .gg_base <- file.path(BASE_DIR, "inputs", "3_post_proc_inputs", "gene_groups_csv")

    .resolve_gg_dir <- function(genome) {
      tag <- sub("_genome$", "", sub("_transcripts.*$", "", genome))
      d <- file.path(.gg_base, "experimental", tag)
      if (dir.exists(d)) return(d)
      return(.gg_base)
    }

    # Build mapping: for each genome, collect gene IDs from its CSVs (in order)
    # ortho_map[[genome]] = character vector of gene IDs (suffix-stripped, ordered)
    ortho_map <- list()
    # Pre-allocate list collectors to avoid O(G²) growing-vector copies
    .label_chunks <- vector("list", .n_groups)
    .gene_chunks <- setNames(
      lapply(.mapped_genomes, function(g) vector("list", .n_groups)),
      .mapped_genomes)
    .shortname_chunks <- setNames(
      lapply(.mapped_genomes, function(g) vector("list", .n_groups)),
      .mapped_genomes)
    .group_name_chunks <- vector("list", .n_groups)  # gene group name per gene
    .mapping_ok <- TRUE

    for (.gi in seq_len(.n_groups)) {
      .per_genome_genes <- list()
      .per_genome_dfs <- list()  # Cache full data frames to avoid re-reading for Shortened_Name

      for (.genome in .mapped_genomes) {
        .csv_name <- GENOME_GENE_GROUPS_MAP[[.genome]][.gi]
        .gg_dir <- .resolve_gg_dir(.genome)
        .csv_path <- file.path(.gg_dir, paste0(.csv_name, ".csv"))

        if (!file.exists(.csv_path)) {
          cat("  [WARN] Gene group CSV not found:", .csv_path, "\n")
          .mapping_ok <- FALSE
          break
        }

        .df <- .fast_read_csv(.csv_path)
        if (!"Gene_ID" %in% colnames(.df)) {
          cat("  [WARN] No Gene_ID column in:", .csv_path, "\n")
          .mapping_ok <- FALSE
          break
        }

        .per_genome_dfs[[.genome]] <- .df  # Cache for Shortened_Name lookup below
        # Strip suffixes to match the matrix rownames
        .ids <- sub("(\\.[0-9]+){1,2}$", "", .df$Gene_ID)
        .per_genome_genes[[.genome]] <- .ids
      }
      if (!.mapping_ok) break

      # Verify all genomes have the same number of rows for this gene group
      .row_counts <- vapply(.per_genome_genes, length, integer(1))
      if (length(unique(.row_counts)) != 1) {
        cat("  [WARN] Row count mismatch for gene group slot", .gi, ":",
            paste(paste(names(.row_counts), .row_counts, sep = "="), collapse = ", "), "\n")
        .mapping_ok <- FALSE
        break
      }

      # Use Shortened_Name from reference genome's CSV as common label (fallback: positional index)
      # Reuse cached data frame from inner loop — avoids redundant O(rows) fread call.
      .ref_df <- .per_genome_dfs[[.ref_genome]]
      if ("Shortened_Name" %in% colnames(.ref_df)) {
        .labels <- .ref_df$Shortened_Name
      } else {
        .prev_count <- if (.gi == 1L) 0L else sum(lengths(.label_chunks[seq_len(.gi - 1L)]))
        .labels <- paste0("gene_", seq_len(.row_counts[1]) + .prev_count)
      }

      .label_chunks[[.gi]] <- .labels
      # Tag each gene with its gene group name (from reference genome's CSV basename)
      .group_name_chunks[[.gi]] <- rep(
        GENOME_GENE_GROUPS_MAP[[.ref_genome]][.gi],
        length(.labels)
      )
      for (.genome in .mapped_genomes) {
        .gene_chunks[[.genome]][[.gi]] <- .per_genome_genes[[.genome]]
        # Collect per-genome Shortened_Names for display in cross-equivalent-gene analysis
        .gdf <- .per_genome_dfs[[.genome]]
        if ("Shortened_Name" %in% colnames(.gdf)) {
          .shortname_chunks[[.genome]][[.gi]] <- .gdf$Shortened_Name
        } else {
          .shortname_chunks[[.genome]][[.gi]] <- .labels
        }
      }
    }

    # Single concatenation: O(N_total) vs O(G² / 2) from incremental c()
    .common_labels <- unlist(.label_chunks, use.names = FALSE)
    for (.genome in .mapped_genomes) {
      ortho_map[[.genome]] <- unlist(.gene_chunks[[.genome]], use.names = FALSE)
    }

    if (.mapping_ok && length(.common_labels) > 0) {
      .use_positional_mapping <- TRUE
      cat("  Positional orthology mapping built:", length(.common_labels), "genes across",
          length(.mapped_genomes), "genomes\n")

      # Remap each genome's matrix to common labels (subset to mapped genes only)
      for (.genome in .mapped_genomes) {
        mat <- tpm_matrices[[.genome]]
        .genome_ids <- ortho_map[[.genome]]
        # Keep only genes that exist in the matrix
        .present <- .genome_ids %in% rownames(mat)
        if (sum(.present) == 0) {
          cat("  [WARN] No mapped genes found in matrix for:", .genome, "\n")
          .use_positional_mapping <- FALSE
          break
        }
        if (sum(!.present) > 0) {
          cat("  [INFO]", sum(!.present), "mapped gene(s) missing from", get_short_name(.genome),
              "matrix — excluded from analysis\n")
        }
        # Subset matrix and rename rows to common labels
        mat <- mat[.genome_ids[.present], , drop = FALSE]
        rownames(mat) <- .common_labels[.present]
        tpm_matrices[[.genome]] <- mat
      }

      if (.use_positional_mapping) {
        # Recompute gene_sets and common_genes from remapped matrices
        gene_sets <- lapply(tpm_matrices, rownames)
        n_genomes <- length(gene_sets)
        .gene_counts <- table(unlist(gene_sets, use.names = FALSE))
        common_genes <- names(.gene_counts[.gene_counts == n_genomes])
        rm(.gene_counts)
        cat("  Common genes after positional mapping:", length(common_genes), "\n")
      }
    } else {
      cat("  [WARN] Positional mapping failed — falling back to direct ID matching\n")
    }

    # Guard: .common_labels is only assigned if the mapping loop completes successfully (line 210)
    # Keep .common_labels alive — it is reused at line 372 to build ortho_gene_ids/ortho_short_names
    suppressWarnings(rm("ortho_map"))
  }
}

# Fallback: direct gene ID intersection (works when genomes share the same ID schema)
if (!.use_positional_mapping) {
  if (n_genomes == 1L) {
    common_genes <- gene_sets[[1L]]
  } else {
    .gene_counts <- table(unlist(gene_sets, use.names = FALSE))
    common_genes <- names(.gene_counts[.gene_counts == n_genomes])
    rm(.gene_counts)
  }
  cat("  Common genes across all genomes:", length(common_genes), "\n")
}

if (length(common_genes) == 0) {
  cat("\n  [WARN] Zero common genes across genomes (different annotations?).\n")
  cat("  Per-genome gene counts:\n")
  for (g in names(gene_sets)) {
    cat("    ", get_short_name(g), ":", length(gene_sets[[g]]), "genes",
        "(example:", paste(head(gene_sets[[g]], 3), collapse = ", "), "...)\n")
  }
  cat("\n  [SKIP] Cross-genome concordance requires shared gene IDs across references.\n")
  cat("  Different reference assemblies use different gene ID schemas.\n")
  cat("  Provide a genome_gene_groups_map in the TOML config to enable positional mapping.\n")
  # Write a sentinel so downstream steps know to skip
  writeLines("skipped:no_common_genes", file.path(OUTPUT_DIR, ".skip_sentinel"))
  quit(save = "no", status = 0)
}

# Pairwise gene overlaps
for (i in seq_len(n_genomes - 1)) {
  for (j in seq(i + 1, n_genomes)) {
    overlap <- length(intersect(gene_sets[[names(gene_sets)[i]]], gene_sets[[names(gene_sets)[j]]]))
    cat("  ", get_short_name(names(gene_sets)[i]), " & ",
        get_short_name(names(gene_sets)[j]), ":", overlap, "shared genes\n")
  }
}

# -----------------------------------------------
# Harmonize sample IDs across genomes
# -----------------------------------------------

cat("\n--- Harmonizing sample IDs ---\n")

sample_sets <- lapply(tpm_matrices, colnames)
if (length(sample_sets) == 1L) {
  common_samples <- sample_sets[[1L]]
} else {
  .sample_counts <- table(unlist(sample_sets, use.names = FALSE))
  common_samples <- names(.sample_counts[.sample_counts == length(sample_sets)])
  rm(.sample_counts)
}
cat("  Common samples across all genomes:", length(common_samples), "\n")

if (length(common_samples) < CONCORDANCE_MIN_SAMPLES) {
  stop("Too few common samples (", length(common_samples), "). Need at least ",
       CONCORDANCE_MIN_SAMPLES)
}

# -----------------------------------------------
# Subset to common genes and samples
# -----------------------------------------------

cat("\n--- Subsetting to common features ---\n")

for (genome in names(tpm_matrices)) {
  tpm_matrices[[genome]] <- tpm_matrices[[genome]][common_genes, common_samples, drop = FALSE]
  cat("  ", get_short_name(genome), ":", nrow(tpm_matrices[[genome]]), "x",
      ncol(tpm_matrices[[genome]]), "\n")
}

.nrow_map <- vapply(tpm_matrices, nrow, integer(1))
.ncol_map <- vapply(tpm_matrices, ncol, integer(1))
method_stats$n_genes_harmonized <- .nrow_map[method_stats$method]
method_stats$n_samples_harmonized <- .ncol_map[method_stats$method]

# -----------------------------------------------
# Filter lowly-expressed genes
# -----------------------------------------------

cat("\n--- Filtering lowly-expressed genes ---\n")
cat("  Criteria: TPM >=", CONCORDANCE_MIN_EXPR, "in at least", CONCORDANCE_MIN_SAMPLES, "samples\n")

expressed_genes_per_genome <- lapply(tpm_matrices, function(mat) {
  rownames(mat)[rowSums(mat >= CONCORDANCE_MIN_EXPR) >= CONCORDANCE_MIN_SAMPLES]
})
expressed_union <- unique(unlist(expressed_genes_per_genome))
.expr_counts <- table(unlist(expressed_genes_per_genome, use.names = FALSE))
expressed_intersect <- names(.expr_counts[.expr_counts == length(expressed_genes_per_genome)])
rm(.expr_counts)

cat("  Genes expressed in ANY genome:", length(expressed_union), "\n")
cat("  Genes expressed in ALL genomes:", length(expressed_intersect), "\n")

# lapply avoids copy-on-modify: for-loop list assignment copies the list spine on each iteration
filtered_genes <- intersect(common_genes, expressed_union)
tpm_matrices <- lapply(tpm_matrices, function(mat) mat[filtered_genes, , drop = FALSE])
cat("  Final gene count:", length(filtered_genes), "\n")

rm(expressed_genes_per_genome, expressed_union)
gc(verbose = FALSE)

# -----------------------------------------------
# Save harmonized data
# -----------------------------------------------

# Build per-genome original gene ID lookup for cross-equivalent-gene analysis.
# Maps common labels back to original per-genome Gene_IDs from the CSV files.
ortho_gene_ids <- NULL
ortho_short_names <- NULL
ortho_gene_groups <- NULL
if (.use_positional_mapping && exists(".gene_chunks") && exists(".label_chunks")) {
  ortho_gene_ids <- list()
  ortho_short_names <- list()
  # Reuse .common_labels (computed at line 210) — avoids redundant O(N) unlist()
  .all_labels <- .common_labels
  for (.genome in names(.gene_chunks)) {
    .all_ids <- unlist(.gene_chunks[[.genome]], use.names = FALSE)
    ortho_gene_ids[[.genome]] <- setNames(.all_ids, .all_labels)
    .all_short <- unlist(.shortname_chunks[[.genome]], use.names = FALSE)
    ortho_short_names[[.genome]] <- setNames(.all_short, .all_labels)
  }
  # Subset to filtered genes
  ortho_gene_ids <- lapply(ortho_gene_ids, function(ids) {
    ids[names(ids) %in% filtered_genes]
  })
  ortho_short_names <- lapply(ortho_short_names, function(sn) {
    sn[names(sn) %in% filtered_genes]
  })
  # Gene group membership: common_label -> gene group CSV basename
  .all_groups <- unlist(.group_name_chunks, use.names = FALSE)
  ortho_gene_groups <- setNames(.all_groups, .all_labels)
  ortho_gene_groups <- ortho_gene_groups[names(ortho_gene_groups) %in% filtered_genes]
}

result <- list(
  tpm_matrices = tpm_matrices,
  common_genes = filtered_genes,
  common_samples = common_samples,
  expressed_in_all = expressed_intersect,
  method_stats = method_stats,
  gene_sets_raw = gene_sets,
  sample_sets_raw = sample_sets,
  ortho_gene_ids = ortho_gene_ids,
  ortho_short_names = ortho_short_names,
  ortho_gene_groups = ortho_gene_groups
)

saveRDS(result, HARMONIZED_RDS)
cat("\n[DONE] Harmonized data saved to:", HARMONIZED_RDS, "\n")
if (.conc_use_dt) {
  data.table::fwrite(method_stats, file.path(TABLES_DIR, "genome_stats.csv"))
} else {
  write.csv(method_stats, file.path(TABLES_DIR, "genome_stats.csv"), row.names = FALSE)
}
