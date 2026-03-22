#!/usr/bin/env Rscript

# ===============================================
# TXIMPORT: RSEM TO DESEQ2 MATRICES - GENERIC
# ===============================================
# Processes RSEM quantification using tximport
# for statistically sound gene expression analysis

suppressPackageStartupMessages({
  library(tximport)
})

# Source shared utilities (provides ensure_output_dir, convert_to_organ_labels, etc.)
SCRIPT_DIR <- Sys.getenv("ANALYSIS_MODULES_DIR", ".")
source(file.path(SCRIPT_DIR, "0_shared_config.R"))
source(file.path(SCRIPT_DIR, "1_utility_functions.R"))
source(file.path(SCRIPT_DIR, "3_Matrix_Creation_utils.R"))

# Ensure match_gene_ids is available (fallback to standalone utility if not in 1_utility_functions.R)
if (!exists("match_gene_ids", mode = "function")) {
  .match_ids_path <- file.path(dirname(SCRIPT_DIR), "utilities", "match_gene_ids.R")
  if (file.exists(.match_ids_path)) source(.match_ids_path)
}

# Thin wrapper around shared save_count_matrices() for backward compatibility.
# Converts the legacy level_suffix (e.g., "_gene_level") to a plain level name.
save_count_matrix <- function(counts_matrix, output_dir, base_name, master_ref, level_suffix,
                              tpm_matrix = NULL) {
  level <- gsub("^_", "", level_suffix)
  save_count_matrices(counts_matrix, output_dir, base_name, master_ref, level,
                      tpm = tpm_matrix, count_type_label = "expected_count")
}

# ===============================================
# CONFIGURATION
# ===============================================

base_dir <- Sys.getenv("BASE_DIR", "")
# Use RSEM_QUANT_ROOT env var (set by setup_method_env) when available;
# it already includes the fasta_tag as a subdir — no MASTER_REFERENCE suffix needed.
RSEM_QUANT_ROOT_ENV <- Sys.getenv("RSEM_QUANT_ROOT", unset = "")
QUANT_DIR_INCLUDES_REF <- nzchar(RSEM_QUANT_ROOT_ENV)
QUANT_DIR <- if (QUANT_DIR_INCLUDES_REF) {
  RSEM_QUANT_ROOT_ENV
} else if (nzchar(base_dir)) {
  file.path(base_dir, "2_ALIGNMENT_RESULTs", "M5_RSEM_Bowtie2", "RSEM_Quant_WD")
} else {
  "RSEM_Quant_WD"  # fallback for standalone execution
}
MATRICES_OUTPUT_DIR <- if (nzchar(base_dir)) {
  file.path(base_dir, "3_POST_PROC", "M5_RSEM_Bowtie2", "count_matrices_from_RSEM_Quant")
} else {
  "count_matrices_from_RSEM_Quant"  # relative fallback when called from pushd context
}
# Use shared GENE_GROUPS_DIR from 0_shared_config.R (already sourced)

# Toggle to generate both gene-level and isoform-level matrices
GENERATE_GENE_LEVEL    <- isTRUE(as.logical(Sys.getenv("RSEM_GENERATE_GENE_LEVEL",    "TRUE")))  # Use .genes.results (aggregated)
GENERATE_ISOFORM_LEVEL <- isTRUE(as.logical(Sys.getenv("RSEM_GENERATE_ISOFORM_LEVEL", "TRUE")))  # Use .isoforms.results (transcript-specific)

# Use SAMPLE_IDS and SAMPLE_LABELS from shared config (0_shared_config.R)
# Already sourced above - no duplication needed

# Create output directories
output_dir <- file.path(MATRICES_OUTPUT_DIR, MASTER_REFERENCE)
ensure_output_dir(output_dir)

cat("\n", strrep("=", 70), "\n")
cat("TXIMPORT: RSEM QUANTIFICATION TO MATRICES\n")
cat(strrep("=", 70), "\n\n")

cat("Configuration:\n")
cat("  • Master Reference:", MASTER_REFERENCE, "\n")
cat("  • Generate gene-level matrices:", GENERATE_GENE_LEVEL, "\n")
cat("  • Generate isoform-level matrices:", GENERATE_ISOFORM_LEVEL, "\n\n")

# Define processing levels
processing_levels <- list()
if (GENERATE_GENE_LEVEL) {
  processing_levels[["gene_level"]] <- list(
    file_type = ".genes.results",
    tx_in = FALSE,
    tx_out = FALSE,
    label = "Gene-Level",
    output_suffix = "_gene_level"
  )
}
if (GENERATE_ISOFORM_LEVEL) {
  processing_levels[["isoform_level"]] <- list(
    file_type = ".isoforms.results",
    tx_in = TRUE,
    tx_out = TRUE,
    label = "Isoform-Level",
    output_suffix = "_isoform_level"
  )
}

if (length(processing_levels) == 0) {
  stop("ERROR: At least one of GENERATE_GENE_LEVEL or GENERATE_ISOFORM_LEVEL must be TRUE")
}

# ===============================================
# PROCESS EACH LEVEL
# ===============================================

for (level_name in names(processing_levels)) {
  level_config <- processing_levels[[level_name]]
  
  cat(strrep("=", 70), "\n")
  cat("PROCESSING:", level_config$label, "\n")
  cat(strrep("=", 70), "\n\n")
  
  # ===============================================
  # STEP 1: LOCATE RSEM OUTPUT FILES
  # ===============================================
  
  cat("Step 1: Locating RSEM", level_config$label, "output files...\n")
  
  rsem_quant_dir <- if (QUANT_DIR_INCLUDES_REF) QUANT_DIR else file.path(QUANT_DIR, MASTER_REFERENCE)

  if (!dir.exists(rsem_quant_dir)) {
    cat("ERROR: RSEM quantification directory not found:", rsem_quant_dir, "\n")
    cat("  Check RSEM_QUANT_ROOT env var or alignment output at:\n")
    cat("  ", file.path(base_dir, "2_ALIGNMENT_RESULTs", "M5_RSEM_Bowtie2", "RSEM_Quant_WD"), "\n")
    next
  }

  # Build paths to RSEM files
  files <- file.path(rsem_quant_dir, SAMPLE_IDS, paste0(SAMPLE_IDS, level_config$file_type))
  names(files) <- SAMPLE_IDS
  
  # Check which files exist
  files_exist <- file.exists(files)
  if (sum(files_exist) == 0) {
    cat("ERROR: No", level_config$label, "quantification files found in", rsem_quant_dir, "\n")
    cat("Skipping this level...\n\n")
    next
  }
  
  current_sample_ids <- SAMPLE_IDS
  if (sum(files_exist) < length(files)) {
    missing <- SAMPLE_IDS[!files_exist]
    cat("WARNING: Missing quantification for samples:", paste(missing, collapse = ", "), "\n")
    files <- files[files_exist]
    current_sample_ids <- SAMPLE_IDS[files_exist]
  }
  
  cat("Found", length(files), "RSEM quantification files\n")
  cat("Samples:", paste(current_sample_ids, collapse = ", "), "\n\n")
  
  # ===============================================
  # STEP 2: IMPORT WITH TXIMPORT
  # ===============================================
  
  cat("Step 2: Importing RSEM data with tximport...\n")
  
  # tximport for RSEM
  txi <- tryCatch(
    tximport(files, type = "rsem", txIn = level_config$tx_in, txOut = level_config$tx_out),
    error = function(e) {
      message("ERROR: tximport failed for ", level_config$label, ": ", e$message)
      message("  Check RSEM output files for corruption or format issues")
      return(NULL)
    }
  )
  if (is.null(txi)) {
    cat("Skipping", level_config$label, "due to tximport failure\n\n")
    next
  }

  entity_type <- if (level_config$tx_out) "transcripts" else "genes"
  cat("Successfully imported data for", ncol(txi$counts), "samples\n")
  cat("Total", entity_type, ":", nrow(txi$counts), "\n\n")

  # ===============================================
  # STEP 3: CREATE SAMPLE METADATA
  # ===============================================
  
  cat("Step 3: Creating sample metadata...\n")
  
  # Get conditions from SAMPLE_LABELS, fallback to sample ID if not found
  conditions <- SAMPLE_LABELS[current_sample_ids]
  missing_labels <- is.na(conditions)
  if (any(missing_labels)) {
    cat("WARNING: No labels found for samples:", paste(current_sample_ids[missing_labels], collapse = ", "), "\n")
    cat("Using sample IDs as condition labels for these samples\n")
    conditions[missing_labels] <- current_sample_ids[missing_labels]
  }
  
  sample_data <- data.frame(
    SampleID = current_sample_ids,
    Condition = conditions,
    row.names = current_sample_ids,
    stringsAsFactors = FALSE
  )
  
  cat("Sample metadata:\n")
  print(sample_data)
  cat("\n")
  
  # ===============================================
  # STEP 4: CHECK FOR REPLICATES
  # ===============================================
  
  cat("Step 4: Checking for biological replicates...\n")
  
  condition_counts <- table(sample_data$Condition)
  has_replicates <- any(condition_counts > 1)
  
  cat("Samples per condition:\n")
  print(condition_counts)
  cat("\n")
  
  if (!has_replicates) {
    cat("WARNING: No biological replicates detected!\n")
    cat("Using TPM normalization for visualization purposes\n\n")
  }
  
  # ===============================================
  # STEP 5: FILTER AND EXTRACT MATRICES
  # ===============================================

  cat("Step 5: Extracting raw expected counts and TPM values...\n")

  # Filter entries with zero effective length (RSEM produces these for unaligned transcripts)
  if (!is.null(txi$length)) {
    zero_length_mask <- rowSums(txi$length == 0) > 0
    n_zero_length <- sum(zero_length_mask)
    if (n_zero_length > 0) {
      cat("  Filtering out", n_zero_length, entity_type, "with zero effective length...\n")
      txi$counts    <- txi$counts[!zero_length_mask, , drop = FALSE]
      txi$abundance <- txi$abundance[!zero_length_mask, , drop = FALSE]
      txi$length    <- txi$length[!zero_length_mask, , drop = FALSE]
      cat("  Remaining", entity_type, ":", nrow(txi$counts), "\n")
      if (nrow(txi$counts) == 0) {
        cat("  ERROR: All", entity_type, "removed after zero-length filtering — skipping level\n")
        next
      }
    }
  }

  # Save full tximport object for DESeq2 (preserves transcript-length offsets)
  # Saved AFTER zero-length filtering so DESeq2 won't encounter NaN/Inf from division by zero
  # Save for both gene-level and isoform-level so DESeq2 can use either
  {
    txi_rds_dir <- file.path(output_dir, level_name)
    dir.create(txi_rds_dir, recursive = TRUE, showWarnings = FALSE)
    rds_name <- if (level_config$tx_out) "tximport_isoform_level.rds" else "tximport_gene_level.rds"
    saveRDS(txi, file.path(txi_rds_dir, rds_name))
    cat("Saved tximport RDS for DESeq2:", rds_name, "\n\n")
  }

  # raw_counts : RSEM posterior expected counts (fractional; correct input for tximport-aware DESeq2)
  # tpm_matrix : transcripts per million (for visualization / heatmaps)
  raw_counts <- txi$counts
  tpm_matrix <- txi$abundance

  cat("Matrix dimensions:", nrow(raw_counts), entity_type, "x", ncol(raw_counts), "samples\n\n")
  
  # ===============================================
  # STEP 6: SAVE RAW EXPECTED COUNT AND TPM MATRICES
  # ===============================================
  
  cat("Step 6: Saving raw expected count and TPM matrices...\n")

  # Dataset-namespaced folder (matches gene group structure; avoids overwrite across datasets)
  CURRENT_DATASET <- Sys.getenv("CURRENT_DATASET", unset = "")
  full_ref_folder <- if (nzchar(CURRENT_DATASET)) {
    paste0(MASTER_REFERENCE, "_in_", CURRENT_DATASET)
  } else {
    MASTER_REFERENCE
  }

  level_output_dir <- file.path(output_dir, level_name)
  full_ref_dir     <- file.path(level_output_dir, full_ref_folder)
  ensure_output_dir(full_ref_dir)

  save_count_matrix(raw_counts, full_ref_dir, full_ref_folder,
                   MASTER_REFERENCE, level_config$output_suffix,
                   tpm_matrix = tpm_matrix)
  cat("\n")
  
  # ===============================================
  # STEP 7: PROCESS GENE GROUPS
  # ===============================================
  
  cat("Step 7: Processing gene groups...\n")

  # Reuse pre-computed gene group file list (hoisted above level loop to avoid
  # repeated list.files() filesystem traversals — same result across levels)
  if (!exists(".gg_files_cached")) {
    .gg_files_cached <- list.files(GENE_GROUPS_DIR, pattern = "\\.(csv|txt|tsv)$", recursive = TRUE, full.names = TRUE)
    .gg_base_cached <- tools::file_path_sans_ext(basename(.gg_files_cached))
    if (length(.gg_files_cached) > 1) {
      dup_idx <- duplicated(.gg_base_cached)
      if (any(dup_idx)) {
        cat("  Note: removing", sum(dup_idx), "duplicate gene group file(s) by basename\n")
        .gg_files_cached <- .gg_files_cached[!dup_idx]
        .gg_base_cached <- .gg_base_cached[!dup_idx]
      }
    }
    gene_groups_str <- Sys.getenv("GENE_GROUPS_STR", unset = "")
    if (nzchar(gene_groups_str)) {
      enabled_groups <- trimws(strsplit(gene_groups_str, " ")[[1]])
      .gg_files_cached <- .gg_files_cached[.gg_base_cached %in% enabled_groups]
      cat("Filtering to configured gene groups:", paste(enabled_groups, collapse = ", "), "\n")
    }
  }
  gene_group_files <- .gg_files_cached
  
  if (length(gene_group_files) == 0) {
    cat("No gene group files found in", GENE_GROUPS_DIR, "\n")
  } else {
    cat("Found", length(gene_group_files), "gene group files\n\n")
    
    # CURRENT_DATASET already loaded in Step 6 above
    # Cache requireNamespace probe once before loop (avoids per-iteration PATH scan)
    .use_dt <- requireNamespace("data.table", quietly = TRUE)

    for (gene_group_file in gene_group_files) {
      gene_group_name <- tools::file_path_sans_ext(basename(gene_group_file))
      # Generate combined output name with dataset suffix
      output_folder_name <- if (nzchar(CURRENT_DATASET)) {
        paste0(gene_group_name, "_in_", CURRENT_DATASET)
      } else {
        gene_group_name
      }
      cat("  Processing:", gene_group_name, "-> Output:", output_folder_name, "\n")
  
      # Read gene list from CSV (first column is Gene_ID)
      # tryCatch returns its value to the outer assignment — fixes scoping issue
      gene_list <- tryCatch({
        if (grepl("\\.csv$", gene_group_file, ignore.case = TRUE)) {
          gene_df <- if (.use_dt) {
            data.table::fread(gene_group_file, header = TRUE, data.table = FALSE)
          } else {
            read.csv(gene_group_file, stringsAsFactors = FALSE, header = TRUE)
          }
          gl <- if ("Gene_ID" %in% colnames(gene_df)) gene_df$Gene_ID else gene_df[[1]]
        } else {
          gl <- suppressWarnings(readLines(gene_group_file))
          gl <- gl[!grepl("^#|^Gene_ID", gl, ignore.case = TRUE) & nzchar(gl)]
          gl <- sub("\t.*", "", gl)  # strip tab-delimited extra fields (e.g. gene names, descriptions)
        }
        trimws(gl)
      }, error = function(e) {
        cat("    Error reading file:", e$message, "\n")
        character(0)
      })
  
      if (length(gene_list) == 0) {
        cat("    Skipping: No genes in list\n")
        next
      }
  
      # Match genes using shared utility (handles exact, prefix, and reverse suffix matching)
      genes_in_data <- match_gene_ids(gene_list, rownames(raw_counts))
  
      if (length(genes_in_data) == 0) {
        cat("    Skipping: No matching", entity_type, "found\n")
        next
      }
  
      cat("    Found", length(genes_in_data), "/", length(gene_list), entity_type, "\n")
  
      gene_group_dir <- file.path(level_output_dir, output_folder_name)
      ensure_output_dir(gene_group_dir)
  
      subset_counts <- raw_counts[genes_in_data, , drop = FALSE]
      subset_tpm <- tpm_matrix[genes_in_data, , drop = FALSE]
      save_count_matrix(subset_counts, gene_group_dir, output_folder_name,
                       MASTER_REFERENCE, level_config$output_suffix,
                       tpm_matrix = subset_tpm)
      cat("    Saved subset matrices\n")
    }
  }
  
  # ===============================================
  # LEVEL SUMMARY
  # ===============================================
  
  cat("\n", strrep("=", 70), "\n")
  cat(level_config$label, "PROCESSING COMPLETE\n")
  cat(strrep("=", 70), "\n\n")
}

# ===============================================
# FINAL SUMMARY
# ===============================================

cat(strrep("=", 70), "\n")
cat("ALL PROCESSING COMPLETE\n")
cat(strrep("=", 70), "\n\n")

cat("Generated matrices for", length(processing_levels), "level(s)\n")
cat("Output directory:", output_dir, "\n")
cat("Matrices ready for heatmap generation\n\n")
