#!/usr/bin/env Rscript

# ===============================================
# MATRIX CREATION SHARED UTILITIES
# ===============================================
# Shared functions used by method-specific Matrix Creation scripts.
# Source this AFTER 0_shared_config.R and 1_utility_functions.R.
# Functions here rely on: ensure_output_dir, convert_to_organ_labels,
# get_output_folder_name, load_runtime_config (all from 0_shared_config.R /
# 1_utility_functions.R).

# ===============================================
# MATRIX SAVING
# ===============================================

# Save count (and optionally TPM) matrices with the naming convention expected
# by build_input_path():
#   {prefix}_{count_type}_{gene_type}_from_{master_ref}_{level}.csv
save_count_matrices <- function(counts, output_dir, prefix, master_ref, level,
                                tpm = NULL, count_type_label = NULL) {
  # Auto-detect count label from method type when not explicitly provided.
  # Salmon/STAR use "NumReads"; RSEM uses "expected_count"; StringTie/prepDE uses "counts".
  if (is.null(count_type_label)) {
    method_type <- get_method_type(CURRENT_METHOD)
    count_type_label <- switch(method_type,
      "salmon" = , "star" = "NumReads",
      "stringtie" = "counts",
      "expected_count"  # default (RSEM and others)
    )
  }
  ensure_output_dir(output_dir)

  save_matrix <- function(matrix_data, count_type, gene_type) {
    output_file <- file.path(output_dir,
      paste0(prefix, "_", count_type, "_", gene_type, "_from_", master_ref, "_", level, ".csv"))
    # Single data.frame() call avoids intermediate O(G×S) copy from as.data.frame + cbind
    matrix_df <- data.frame(GeneID = rownames(matrix_data), matrix_data,
                            check.names = FALSE, stringsAsFactors = FALSE)
    rownames(matrix_df) <- NULL
    if (.HAS_DATATABLE) {
      data.table::fwrite(matrix_df, output_file, sep = ",", quote = FALSE)
    } else {
      write.table(matrix_df, output_file, sep = ",", quote = FALSE, row.names = FALSE)
    }
    cat("Saved:", basename(output_file), "\n")
  }

  # Always save Gene_ID (needed as source for Shortened_Name and by DESeq2)
  save_matrix(counts, count_type_label, "Gene_ID")

  # Hoist gene_group_for_map computation — used by both counts and TPM Shortened_Name paths.
  # Do NOT pre-apply convert_to_organ_labels() here — column label conversion
  # is handled at runtime by apply_labels() in the processing engine.
  # Pre-applying Organ labels bakes them into the file, making it impossible
  # for the processing engine to produce SRR_ID-labeled output from this file.
  # Row name conversion (Gene_ID -> Shortened_Name) is safe to bake in since
  # apply_labels() detects already-shortened names and skips re-conversion.
  gene_group_for_map <- if (nzchar(CURRENT_DATASET)) {
    sub(paste0("_in_", CURRENT_DATASET, "$"), "", prefix)
  } else prefix

  if ("Shortened_Name" %in% GENE_TYPES) {
    counts_short <- tryCatch(
      convert_to_shortened_names(counts, gene_group_for_map),
      error = function(e) counts
    )
    save_matrix(counts_short, count_type_label, "Shortened_Name")
  }

  if (!is.null(tpm)) {
    save_matrix(tpm, "tpm", "Gene_ID")
    if ("Shortened_Name" %in% GENE_TYPES) {
      tpm_short <- tryCatch(
        convert_to_shortened_names(tpm, gene_group_for_map),
        error = function(e) tpm
      )
      save_matrix(tpm_short, "tpm", "Shortened_Name")
    }
  }
}

# ===============================================
# GENE GROUP FILTERING
# ===============================================

# In-memory cache for parsed gene lists — avoids re-reading the same CSV/TXT
# file when filter_by_gene_group() is called across multiple levels in
# run_matrix_saving(). O(1) hash lookup per subsequent call vs O(parse_time).
.gene_list_cache <- new.env(hash = TRUE, parent = emptyenv())

# Filter a count matrix to rows matching genes listed in a CSV or plain-text file.
# Supports prefix matching for isoform IDs (e.g., "SMEL4.1_01g005840" matches
# "SMEL4.1_01g005840.1.01").
filter_by_gene_group <- function(counts_matrix, gene_list_file) {
  if (is.null(counts_matrix) || nrow(counts_matrix) == 0) {
    cat("Empty counts matrix — skipping gene group filtering\n")
    return(NULL)
  }
  if (!file.exists(gene_list_file)) {
    cat("Gene list file not found:", gene_list_file, "\n")
    return(NULL)
  }

  # Check in-memory cache first — O(1) hash lookup
  .cache_key <- gene_list_file
  if (exists(.cache_key, envir = .gene_list_cache, inherits = FALSE)) {
    gene_list <- get(.cache_key, envir = .gene_list_cache, inherits = FALSE)
  } else {
    file_ext <- tools::file_ext(gene_list_file)

    if (tolower(file_ext) == "csv") {
      # Use data.table::fread when available (10-50x faster for large CSVs)
      gene_df <- tryCatch(
        if (.HAS_DATATABLE) {
          data.table::fread(gene_list_file, header = TRUE, data.table = FALSE)
        } else {
          read.csv(gene_list_file, stringsAsFactors = FALSE, header = TRUE)
        },
        error = function(e) NULL)
      if (is.null(gene_df) || nrow(gene_df) == 0) {
        cat("Failed to read CSV gene list\n")
        return(NULL)
      }
      if (!"Gene_ID" %in% colnames(gene_df)) {
        cat("  Warning: 'Gene_ID' column not found in", basename(gene_list_file),
            "— using first column ('", colnames(gene_df)[1], "') as gene IDs\n")
      }
      gene_list <- trimws(
        if ("Gene_ID" %in% colnames(gene_df)) gene_df$Gene_ID else gene_df[[1]])
    } else {
      gene_list <- suppressWarnings(readLines(gene_list_file))
      gene_list <- gene_list[!grepl("^#|^Gene_ID(\\t|$)|^Gene(\\t|$)", gene_list, ignore.case = TRUE) & nzchar(gene_list)]
      gene_list <- trimws(sub("\t.*", "", gene_list))  # strip tab-delimited extra fields (e.g. gene names)
    }

    gene_list <- gene_list[nzchar(gene_list)]
    # Store in cache for subsequent calls with the same file
    assign(.cache_key, gene_list, envir = .gene_list_cache)
  }

  # Use shared gene ID matching from 1_utility_functions.R
  matched_genes <- match_gene_ids(gene_list, rownames(counts_matrix))

  if (length(matched_genes) == 0) {
    cat("No genes matched from list\n")
    return(NULL)
  }

  cat("Matched", length(matched_genes), "of", length(gene_list), "genes\n")
  return(counts_matrix[matched_genes, , drop = FALSE])
}

# ===============================================
# SHARED MATRIX SAVING LOOP
# ===============================================
# Iterates over levels in `results`, saves full-genome and gene-group-filtered
# matrices for each level, using the dataset-namespaced folder structure expected
# by build_input_path().
#
# Arguments:
#   results        - named list; main keys are level names ("gene_level",
#                    "isoform_level"); companion TPM keys end in "_tpm"
#   output_dir     - base output directory (relative to method CWD)
#   master_ref     - MASTER_REFERENCE string
#   count_label    - raw-count type label ("expected_count" or "NumReads")
#   gene_groups_dir - path to gene_groups directory (NULL to skip filtering)

run_matrix_saving <- function(results, output_dir, master_ref,
                               count_label, gene_groups_dir = NULL) {
  if (length(results) == 0) {
    cat("No results to save.\n")
    return(invisible(NULL))
  }

  method_base_dir <- Sys.getenv("METHOD_BASE_DIR", unset = ".")
  config <- load_runtime_config(method_base_dir)

  # Skip companion "_tpm" keys — they are processed alongside their parent level
  level_names <- names(results)[!grepl("_tpm$", names(results))]

  # Pre-build gene group file lookup (single list.files call instead of one per group per level)
  .gg_file_map <- NULL
  if (!is.null(gene_groups_dir) && dir.exists(gene_groups_dir)) {
    .all_gg_files <- list.files(gene_groups_dir, pattern = "\\.(csv|txt|tsv)$",
                                recursive = TRUE, full.names = TRUE)
    .gg_names <- tools::file_path_sans_ext(basename(.all_gg_files))
    .gg_first <- !duplicated(.gg_names)  # keep first occurrence (consistent with tximport script)
    .gg_file_map <- setNames(.all_gg_files[.gg_first], .gg_names[.gg_first])
  }

  for (level in level_names) {
    level_output <- file.path(output_dir, master_ref, level)
    ensure_output_dir(level_output)

    # Full-genome matrices in dataset-namespaced subfolder (avoids overwrite across datasets)
    full_ref_folder <- get_output_folder_name(master_ref)
    full_ref_output <- file.path(level_output, full_ref_folder)
    ensure_output_dir(full_ref_output)

    tpm_data <- results[[paste0(level, "_tpm")]]   # NULL when TPM not captured
    save_count_matrices(results[[level]], full_ref_output, full_ref_folder,
                        master_ref, level,
                        tpm = tpm_data, count_type_label = count_label)

    # Gene-group-filtered matrices
    if (!is.null(.gg_file_map)) {
      # O(G × F) where G = gene groups, F = filter + save cost per group
      for (gene_group in config$gene_groups) {
        # Check pre-built map first (O(1) hash lookup) before filesystem probes
        if (gene_group %in% names(.gg_file_map)) {
          gf <- .gg_file_map[[gene_group]]
        } else {
          .gf_candidates <- file.path(gene_groups_dir, paste0(gene_group, c(".csv", ".txt")))
          .gf_exist <- file.exists(.gf_candidates)
          gf <- if (any(.gf_exist)) .gf_candidates[which.max(.gf_exist)] else ""
        }
        if (!nzchar(gf) || !file.exists(gf)) {
          cat("Gene group file not found:", gene_group, "\n")
          next
        }

        filtered <- filter_by_gene_group(results[[level]], gf)
        if (!is.null(filtered)) {
          folder_name  <- get_output_folder_name(gene_group)
          group_output <- file.path(level_output, folder_name)
          ensure_output_dir(group_output)
          filtered_tpm <- if (!is.null(tpm_data)) {
            common_rows <- intersect(rownames(filtered), rownames(tpm_data))
            if (length(common_rows) > 0) tpm_data[common_rows, , drop = FALSE] else NULL
          } else NULL
          save_count_matrices(filtered, group_output, folder_name,
                              master_ref, level,
                              tpm = filtered_tpm, count_type_label = count_label)
        }
      }
    }
  }

  cat("\nMatrix saving complete.\n")
  invisible(results)
}
