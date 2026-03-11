#!/usr/bin/env Rscript

# ===============================================
# TXIMPORT: SALMON TO DESEQ2 MATRICES - GENERIC
# ===============================================
# Processes Salmon quantification using tximport
# for statistically sound gene expression analysis
# GENERIC VERSION - adaptable for any method

suppressPackageStartupMessages({
  library(tximport)
  # dplyr and tibble are NOT loaded here — all operations use base R.
  # They are available via 1_utility_functions.R if sourced scripts need them.
  # DESeq2 is NOT used here for normalization — raw tximport counts are saved
  # directly so downstream analysis scripts can run DESeq2 without double-normalization.
})

# ===============================================
# CONFIGURATION
# ===============================================

# Source shared config and utilities (DRY principle)
SCRIPT_DIR <- Sys.getenv("ANALYSIS_MODULES_DIR", ".")
source(file.path(SCRIPT_DIR, "0_shared_config.R"))
source(file.path(SCRIPT_DIR, "1_utility_functions.R"))
source(file.path(SCRIPT_DIR, "3_Matrix_Creation_utils.R"))

# Wrapper to adapt the shared save_count_matrices() signature to the call sites
# below, which pass a level_suffix like "_gene_level".
save_count_matrix <- function(counts_matrix, output_dir, base_name, master_ref, level_suffix,
                              tpm_matrix = NULL) {
  level <- gsub("^_", "", level_suffix)
  save_count_matrices(counts_matrix, output_dir, base_name, master_ref, level,
                      tpm = tpm_matrix, count_type_label = "NumReads")
}

# Use absolute paths from env vars to avoid working-directory dependency.
# BASE_DIR is exported by run_all_post_processing.sh.
# SALMON_QUANT_ROOT is optionally exported by the alignment pipeline after
# set_fasta_output_dirs() — it already includes the fasta_tag as a subdir.
BASE_DIR <- Sys.getenv("BASE_DIR", unset = "")
SALMON_QUANT_ROOT_ENV <- Sys.getenv("SALMON_QUANT_ROOT", unset = "")

# When SALMON_QUANT_ROOT is set it already ends with the fasta_tag, so use it
# directly as salmon_quant_dir (no MASTER_REFERENCE suffix needed).
# When it is not set, build the parent path from BASE_DIR and let the loop
# below append MASTER_REFERENCE.
QUANT_DIR_INCLUDES_REF <- nzchar(SALMON_QUANT_ROOT_ENV)
QUANT_DIR <- if (QUANT_DIR_INCLUDES_REF) {
  SALMON_QUANT_ROOT_ENV
} else if (nzchar(BASE_DIR)) {
  file.path(BASE_DIR, "2_ALIGNMENT_RESULTs", "M4_Salmon_Saf", "Salmon_Quant")
} else {
  "Salmon_Quant"  # last-resort relative fallback
}
MASTER_REFERENCE <- Sys.getenv("MASTER_REFERENCE", unset = "")
if (!nzchar(MASTER_REFERENCE)) {
  warning("MASTER_REFERENCE env var is not set — defaulting to 'All_Smel_Genes'. ",
          "Export MASTER_REFERENCE before running this script to suppress this warning.")
  MASTER_REFERENCE <- "All_Smel_Genes"
}
# Use an absolute path when BASE_DIR is available (run_method_analysis does pushd, so the
# working directory is correct, but an absolute path allows the script to be run standalone).
MATRICES_OUTPUT_DIR <- if (nzchar(BASE_DIR)) {
  file.path(BASE_DIR, "3_POST_PROC", "M4_Salmon_Saf", "count_matrices_from_Salmon_Quant")
} else {
  "count_matrices_from_Salmon_Quant"  # relative fallback when called from pushd context
}
# Use shared GENE_GROUPS_DIR from 0_shared_config.R (already sourced)

# Toggle to generate both gene-level and isoform-level matrices
GENERATE_GENE_LEVEL <- TRUE      # Summarize transcripts to genes
GENERATE_ISOFORM_LEVEL <- TRUE   # Keep transcript-level data

# Use SAMPLE_IDS from shared config (0_shared_config.R)
# Override here if needed for method-specific samples

# Create output directories
output_dir <- file.path(MATRICES_OUTPUT_DIR, MASTER_REFERENCE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("TXIMPORT: SALMON QUANTIFICATION TO MATRICES\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

cat("Configuration:\n")
cat("  • Master Reference:", MASTER_REFERENCE, "\n")
cat("  • Generate gene-level matrices:", GENERATE_GENE_LEVEL, "\n")
cat("  • Generate isoform-level matrices:", GENERATE_ISOFORM_LEVEL, "\n\n")

# Define processing levels
processing_levels <- list()
if (GENERATE_GENE_LEVEL) {
  processing_levels[["gene_level"]] <- list(
    tx_out = FALSE,
    label = "Gene-Level",
    output_suffix = "_gene_level"
  )
}
if (GENERATE_ISOFORM_LEVEL) {
  processing_levels[["isoform_level"]] <- list(
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

  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat("PROCESSING:", level_config$label, "\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")

  # ===============================================
  # STEP 1: LOCATE SALMON OUTPUT FILES
  # ===============================================

  cat("Step 1: Locating Salmon", level_config$label, "output files...\n")

  # If QUANT_DIR already includes the fasta_tag (from SALMON_QUANT_ROOT), use as-is.
  # Otherwise append MASTER_REFERENCE (which equals the fasta_tag in the normal run).
  salmon_quant_dir <- if (QUANT_DIR_INCLUDES_REF) QUANT_DIR else file.path(QUANT_DIR, MASTER_REFERENCE)

  # Build paths to Salmon quant.sf files
  files <- file.path(salmon_quant_dir, SAMPLE_IDS, "quant.sf")
  names(files) <- SAMPLE_IDS

  # Check which files exist
  files_exist <- file.exists(files)
  if (sum(files_exist) == 0) {
    cat("ERROR: No", level_config$label, "quantification files found in", salmon_quant_dir, "\n")
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

  cat("Found", length(files), "Salmon quantification files\n")
  cat("Samples:", paste(current_sample_ids, collapse = ", "), "\n\n")

  # ===============================================
  # STEP 2: IMPORT WITH TXIMPORT
  # ===============================================

  cat("Step 2: Importing Salmon data with tximport...\n")

  # For Salmon, we need a tx2gene mapping if summarizing to genes
  # Use the gene_trans_map file corresponding to the MASTER_REFERENCE
  INPUT_FASTAS_DIR <- Sys.getenv("INPUT_FASTAS_DIR", unset = "")
  if (!nzchar(INPUT_FASTAS_DIR)) {
    INPUT_FASTAS_DIR <- if (nzchar(BASE_DIR)) file.path(BASE_DIR, "inputs") else "../../inputs"
  }

  # Search for the tx2gene mapping file produced by create_gene_trans_map()
  # (saved as <fasta>.gene_trans_map next to the FASTA, or exportable via GENE_TRANS_MAP_FILE)
  tx2gene_candidates <- c(
    Sys.getenv("GENE_TRANS_MAP_FILE", unset = ""),
    file.path(INPUT_FASTAS_DIR, "mapping", paste0(MASTER_REFERENCE, ".fa.gene_trans_map")),
    file.path(INPUT_FASTAS_DIR, "mapping", paste0(MASTER_REFERENCE, ".fasta.gene_trans_map")),
    file.path(INPUT_FASTAS_DIR, "fasta", paste0(MASTER_REFERENCE, ".fa.gene_trans_map"))
  )
  # Glob-based fallback: find any .gene_trans_map containing MASTER_REFERENCE
  all_maps <- list.files(INPUT_FASTAS_DIR, pattern = "\\.gene_trans_map$",
                         recursive = TRUE, full.names = TRUE)
  ref_maps <- all_maps[grepl(MASTER_REFERENCE, all_maps, fixed = TRUE)]
  tx2gene_candidates <- c(tx2gene_candidates, ref_maps)

  tx2gene_file <- NULL
  for (.cand in tx2gene_candidates) {
    if (nzchar(.cand) && file.exists(.cand)) { tx2gene_file <- .cand; break }
  }
  rm(.cand)

  if (level_config$tx_out == FALSE) {
    # Gene-level: need tx2gene mapping
    if (is.null(tx2gene_file)) {
      cat("ERROR: tx2gene mapping file not found for", MASTER_REFERENCE, "\n")
      cat("  Searched in:", INPUT_FASTAS_DIR, "\n")
      cat("  Tip: export GENE_TRANS_MAP_FILE=/path/to/<ref>.fa.gene_trans_map\n")
      cat("Skipping gene-level processing...\n\n")
      next
    }

    # Read tx2gene mapping (columns: GENEID, TXNAME → reorder to TXNAME, GENEID for tximport)
    tx2gene <- read.table(tx2gene_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                         colClasses = c("character", "character"))
    colnames(tx2gene) <- c("GENEID", "TXNAME")
    tx2gene$GENEID <- trimws(tx2gene$GENEID)
    tx2gene$TXNAME <- trimws(tx2gene$TXNAME)
    tx2gene <- tx2gene[, c("TXNAME", "GENEID")]
    cat("Loaded tx2gene mapping:", nrow(tx2gene), "entries\n")
    cat("  Using:", tx2gene_file, "\n")

    # Aggregate to gene level using tximport's built-in summarization.
    # This correctly recomputes gene-level TPM (not a simple rowsum of per-transcript TPMs)
    # and computes the isoform-usage-weighted effective length required by DESeq2 offsets.
    cat("  Importing at transcript level and aggregating to genes via tximport...\n")
    txi <- tximport(files, type = "salmon",
                    txIn = TRUE, txOut = FALSE,
                    tx2gene = tx2gene,
                    ignoreTxVersion = TRUE, ignoreAfterBar = FALSE)
  } else {
    # Isoform-level: no tx2gene needed
    cat("  Importing at transcript level...\n")
    txi <- tximport(files, type = "salmon", txIn = TRUE, txOut = TRUE,
                   ignoreTxVersion = TRUE, ignoreAfterBar = FALSE)
  }

  entity_type <- if (level_config$tx_out) "transcripts" else "genes"
  cat("Successfully imported data for", ncol(txi$counts), "samples\n")
  cat("Total", entity_type, ":", nrow(txi$counts), "\n\n")

  # Save full tximport object for DESeq2 (preserves transcript-length offsets)
  if (!level_config$tx_out) {
    txi_rds_dir <- file.path(output_dir, level_name)
    dir.create(txi_rds_dir, recursive = TRUE, showWarnings = FALSE)
    saveRDS(txi, file.path(txi_rds_dir, "tximport_gene_level.rds"))
    cat("Saved tximport RDS for DESeq2: tximport_gene_level.rds\n\n")
  }

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
  # STEP 5: NORMALIZATION
  # ===============================================

  # Extract raw counts (Salmon NumReads) and TPM for output.
  # Raw counts are saved as-is for DESeq2 input (downstream scripts normalise
  # internally); TPM is saved for visualisation (heatmaps, PCA, etc.).
  raw_counts <- txi$counts
  tpm_matrix <- txi$abundance

  cat("Step 5: Matrices ready (raw NumReads + TPM from tximport).\n")
  if (!has_replicates) {
    cat("WARNING: No biological replicates — DESeq2 analyses will not be meaningful.\n")
    cat("         Use TPM matrices for visualisation only.\n")
  }

  cat("Matrix dimensions:", nrow(raw_counts), entity_type, "x", ncol(raw_counts), "samples\n\n")

  # ===============================================
  # STEP 6: SAVE FULL COUNT MATRIX
  # ===============================================

  cat("Step 6: Saving count matrices...\n")

  level_output_dir <- file.path(output_dir, level_name)
  dir.create(level_output_dir, recursive = TRUE, showWarnings = FALSE)

  # Dataset-namespaced folder (matches gene group structure; avoids overwrite across datasets)
  full_ref_folder <- if (nzchar(CURRENT_DATASET)) {
    paste0(MASTER_REFERENCE, "_in_", CURRENT_DATASET)
  } else {
    MASTER_REFERENCE
  }
  # Save full-genome matrix in a named subdirectory so the path matches
  # build_input_path() which always expects: {level}/{folder_name}/{folder_name}_...tsv
  full_ref_dir <- file.path(level_output_dir, full_ref_folder)
  dir.create(full_ref_dir, recursive = TRUE, showWarnings = FALSE)

  # Save raw_counts as "NumReads" — downstream DESeq2 scripts require raw integer-like
  # counts as input and perform their own normalization internally. Saving DESeq2-normalized
  # values here would cause double-normalization. TPM is saved alongside for visualization.
  save_count_matrix(raw_counts, full_ref_dir, full_ref_folder,
                   MASTER_REFERENCE, level_config$output_suffix,
                   tpm_matrix = tpm_matrix)
  cat("\n")

  # ===============================================
  # STEP 7: PROCESS GENE GROUPS
  # ===============================================

  cat("Step 7: Processing gene groups...\n")

  gene_group_files <- list.files(GENE_GROUPS_DIR, pattern = "\\.(csv|txt|tsv)$", recursive = TRUE, full.names = TRUE)

  # Filter to only process gene groups specified in GENE_GROUPS_STR (from bash config)
  gene_groups_str <- Sys.getenv("GENE_GROUPS_STR", unset = "")
  if (nzchar(gene_groups_str)) {
    enabled_groups <- trimws(strsplit(gene_groups_str, " ")[[1]])
    gene_group_files <- gene_group_files[tools::file_path_sans_ext(basename(gene_group_files)) %in% enabled_groups]
    cat("Filtering to configured gene groups:", paste(enabled_groups, collapse = ", "), "\n")
  }

  if (length(gene_group_files) == 0) {
    cat("No gene group files found in", GENE_GROUPS_DIR, "\n")
  } else {
    cat("Found", length(gene_group_files), "gene group files\n\n")

    successful_groups <- 0

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
      # Use tryCatch return value so the assignment is visible in this scope
      gene_list <- tryCatch({
        if (grepl("\\.csv$", gene_group_file, ignore.case = TRUE)) {
          gene_df <- read.csv(gene_group_file, stringsAsFactors = FALSE, header = TRUE)
          gl <- if ("Gene_ID" %in% colnames(gene_df)) gene_df$Gene_ID else gene_df[[1]]
        } else {
          gl <- suppressWarnings(readLines(gene_group_file))
          gl <- gl[!grepl("^#|^Gene_ID", gl, ignore.case = TRUE) & nzchar(gl)]
          # For multi-column TSV/TXT files, take only the first tab-delimited field
          # so that gene names, descriptions, etc. are not concatenated into the ID.
          gl <- trimws(sub("\t.*", "", gl))
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

      # Match gene list to row names using shared helper (handles version suffixes)
      genes_in_data <- match_gene_ids(gene_list, rownames(raw_counts))

      if (length(genes_in_data) == 0) {
        cat("    Skipping: No matching", entity_type, "found\n")
        next
      }

      cat("    Found", length(genes_in_data), "/", length(gene_list), entity_type, "\n")

      # Create gene group subdirectory with dataset suffix (consistent naming)
      gene_group_dir <- file.path(level_output_dir, output_folder_name)
      dir.create(gene_group_dir, recursive = TRUE, showWarnings = FALSE)

      # Save gene group subset using helper (saves NumReads, TPM, Gene_ID, Shortened_Name variants)
      subset_counts <- raw_counts[genes_in_data, , drop = FALSE]
      subset_tpm <- tpm_matrix[genes_in_data, , drop = FALSE]
      save_count_matrix(subset_counts, gene_group_dir, output_folder_name,
                       MASTER_REFERENCE, level_config$output_suffix,
                       tpm_matrix = subset_tpm)
      cat("    Saved subset matrices\n")
      successful_groups <- successful_groups + 1
    }

    cat("\n  Summary: Processed", successful_groups, "/", length(gene_group_files), "gene groups\n")
  }

  # ===============================================
  # LEVEL SUMMARY
  # ===============================================

  cat("\n", paste(rep("=", 70), collapse = ""), "\n")
  cat(level_config$label, "PROCESSING COMPLETE\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")
}

# ===============================================
# FINAL SUMMARY
# ===============================================

cat(paste(rep("=", 70), collapse = ""), "\n")
cat("ALL PROCESSING COMPLETE\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

cat("Generated matrices for", length(processing_levels), "level(s)\n")
cat("Output directory:", output_dir, "\n")
cat("Matrices ready for heatmap generation\n\n")
