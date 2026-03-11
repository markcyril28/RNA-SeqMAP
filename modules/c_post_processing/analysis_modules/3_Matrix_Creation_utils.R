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
#   {prefix}_{count_type}_{gene_type}_from_{master_ref}_{level}.tsv
save_count_matrices <- function(counts, output_dir, prefix, master_ref, level,
                                tpm = NULL, count_type_label = "expected_count") {
  ensure_output_dir(output_dir)

  save_matrix <- function(matrix_data, count_type, gene_type) {
    output_file <- file.path(output_dir,
      paste0(prefix, "_", count_type, "_", gene_type, "_from_", master_ref, "_", level, ".tsv"))
    matrix_df <- as.data.frame(matrix_data, check.names = FALSE)
    matrix_df <- cbind(GeneID = rownames(matrix_data), matrix_df)
    rownames(matrix_df) <- NULL
    data.table::fwrite(matrix_df, output_file, sep = "\t", quote = FALSE)
    cat("Saved:", basename(output_file), "\n")
  }

  # Always save Gene_ID (needed as source for Shortened_Name and by DESeq2)
  save_matrix(counts, count_type_label, "Gene_ID")
  if ("Shortened_Name" %in% GENE_TYPES) {
    # Only convert column headers from SRR IDs to organ labels when "Organ" is
    # an active label type; otherwise keep original SRR column names.
    counts_short <- if ("Organ" %in% LABEL_TYPES) convert_to_organ_labels(counts) else counts
    # Also convert row names from Gene_ID to Shortened_Name using the gene group mapping.
    # prefix is typically the gene_group or folder_name; extract gene_group portion
    # (strip "_in_<dataset>" suffix if present) for the mapping lookup.
    gene_group_for_map <- sub("_in_.*$", "", prefix)
    counts_short <- tryCatch(
      convert_to_shortened_names(counts_short, gene_group_for_map),
      error = function(e) counts_short
    )
    save_matrix(counts_short, count_type_label, "Shortened_Name")
  }

  if (!is.null(tpm)) {
    save_matrix(tpm, "tpm", "Gene_ID")
    if ("Shortened_Name" %in% GENE_TYPES) {
      tpm_short <- if ("Organ" %in% LABEL_TYPES) convert_to_organ_labels(tpm) else tpm
      gene_group_for_map <- sub("_in_.*$", "", prefix)
      tpm_short <- tryCatch(
        convert_to_shortened_names(tpm_short, gene_group_for_map),
        error = function(e) tpm_short
      )
      save_matrix(tpm_short, "tpm", "Shortened_Name")
    }
  }
}

# ===============================================
# GENE GROUP FILTERING
# ===============================================

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

  file_ext <- tools::file_ext(gene_list_file)

  if (tolower(file_ext) == "csv") {
    gene_df <- tryCatch(
      read.csv(gene_list_file, stringsAsFactors = FALSE, header = TRUE),
      error = function(e) NULL)
    if (is.null(gene_df) || nrow(gene_df) == 0) {
      cat("Failed to read CSV gene list\n")
      return(NULL)
    }
    gene_list <- trimws(
      if ("Gene_ID" %in% colnames(gene_df)) gene_df$Gene_ID else gene_df[[1]])
  } else {
    gene_list <- suppressWarnings(readLines(gene_list_file))
    gene_list <- gene_list[!grepl("^#|^Gene", gene_list, ignore.case = TRUE) & nzchar(gene_list)]
    gene_list <- trimws(sub("\t.*", "", gene_list))  # strip tab-delimited extra fields (e.g. gene names)
  }

  gene_list <- gene_list[nzchar(gene_list)]

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
    if (!is.null(gene_groups_dir) && dir.exists(gene_groups_dir)) {
      for (gene_group in config$gene_groups) {
        gf <- file.path(gene_groups_dir, paste0(gene_group, ".csv"))
        if (!file.exists(gf)) gf <- file.path(gene_groups_dir, paste0(gene_group, ".txt"))
        # Search subdirectories if not found at top level
        if (!file.exists(gf)) {
          hits <- list.files(gene_groups_dir, pattern = paste0("^", gene_group, "\\.(csv|txt)$"),
                             recursive = TRUE, full.names = TRUE)
          if (length(hits) > 0) gf <- hits[1]
        }
        if (!file.exists(gf)) {
          cat("Gene group file not found:", gene_group, "\n")
          next
        }

        filtered <- filter_by_gene_group(results[[level]], gf)
        if (!is.null(filtered)) {
          folder_name  <- get_output_folder_name(gene_group)
          group_output <- file.path(level_output, folder_name)
          ensure_output_dir(group_output)
          filtered_tpm <- if (!is.null(tpm_data)) tpm_data[rownames(filtered), , drop = FALSE] else NULL
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
