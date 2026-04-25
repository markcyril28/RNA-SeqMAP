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

# Default CURRENT_METHOD for this Salmon-specific script before shared config
# sets a generic fallback (M5). Prevents wrong-directory lookups in standalone mode.
if (!nzchar(Sys.getenv("CURRENT_METHOD", unset = ""))) {
  Sys.setenv(CURRENT_METHOD = "M4_Salmon_Saf")
}

# Source shared config and utilities (DRY principle)
SCRIPT_DIR <- Sys.getenv("ANALYSIS_MODULES_DIR", {
  if (nzchar(Sys.getenv("WF_MANAGED_ENV", "")))
    stop("[SALMON_TXIMPORT] ANALYSIS_MODULES_DIR is required under workflow manager (WF_MANAGED_ENV is set).")
  "."
})
source(file.path(SCRIPT_DIR, "0_shared_config.R"))
source(file.path(SCRIPT_DIR, "1_utility_functions.R"))
source(file.path(SCRIPT_DIR, "3_Matrix_Creation_utils.R"))

# Ensure match_gene_ids is available (fallback to standalone utility if not in 1_utility_functions.R)
if (!exists("match_gene_ids", mode = "function")) {
  .match_ids_path <- file.path(dirname(SCRIPT_DIR), "utilities", "match_gene_ids.R")
  if (file.exists(.match_ids_path)) {
    source(.match_ids_path)
  } else if (nzchar(Sys.getenv("WF_MANAGED_ENV", ""))) {
    stop("[SALMON TXIMPORT] match_gene_ids.R not found at ", .match_ids_path,
         ". Ensure ANALYSIS_MODULES_DIR points to a directory whose sibling 'utilities/' contains match_gene_ids.R.")
  }
}

# Wrapper to adapt the shared save_count_matrices() signature to the call sites
# below, which pass a level_suffix like "_gene_level".
save_count_matrix <- function(counts_matrix, output_dir, base_name, master_ref, level_suffix,
                              tpm_matrix = NULL) {
  level <- gsub("^_", "", level_suffix)
  save_count_matrices(counts_matrix, output_dir, base_name, master_ref, level,
                      tpm = tpm_matrix, count_type_label = "NumReads")
}

# Use absolute paths from env vars to avoid working-directory dependency.
# BASE_DIR is exported by run_post_processing.sh.
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
  file.path(BASE_DIR, "II_RESULTS/2_ALIGNMENT_RESULTs", CURRENT_METHOD, "Salmon_Quant")
} else if (nzchar(Sys.getenv("WF_MANAGED_ENV", ""))) {
  stop("[SALMON TXIMPORT] BASE_DIR or SALMON_QUANT_ROOT is required under workflow manager (WF_MANAGED_ENV is set).")
} else {
  message("[SALMON TXIMPORT] WARN: BASE_DIR and SALMON_QUANT_ROOT not set; using relative 'Salmon_Quant/'. ",
          "Set BASE_DIR for orchestrated execution (Nextflow/Snakemake).")
  "Salmon_Quant"  # last-resort relative fallback
}
# MASTER_REFERENCE is already set by 0_shared_config.R (sourced above).
# Do not re-read from env here — the defaults would diverge.
# Use an absolute path when BASE_DIR is available (run_method_analysis does pushd, so the
# working directory is correct, but an absolute path allows the script to be run standalone).
MATRICES_OUTPUT_DIR <- if (nzchar(BASE_DIR)) {
  file.path(BASE_DIR, "II_RESULTS", "3_POST_PROC", Sys.getenv("CURRENT_GENE_GROUP", "_active"), CURRENT_METHOD, "count_matrices_from_Salmon_Quant")
} else if (nzchar(Sys.getenv("WF_MANAGED_ENV", ""))) {
  stop("[SALMON TXIMPORT] BASE_DIR is required for output directory under workflow manager.")
} else {
  message("[SALMON TXIMPORT] WARN: BASE_DIR not set; using relative 'count_matrices_from_Salmon_Quant/'. ",
          "Set BASE_DIR for orchestrated execution (Nextflow/Snakemake).")
  "count_matrices_from_Salmon_Quant"  # relative fallback when called from pushd context
}
# Use shared GENE_GROUPS_DIR from 0_shared_config.R (already sourced)

# Toggle to generate both gene-level and isoform-level matrices (overridable via env vars)
GENERATE_GENE_LEVEL <- isTRUE(as.logical(Sys.getenv("SALMON_GENERATE_GENE_LEVEL", unset = "TRUE")))
GENERATE_ISOFORM_LEVEL <- isTRUE(as.logical(Sys.getenv("SALMON_GENERATE_ISOFORM_LEVEL", unset = "TRUE")))

# Reuse .HAS_DATATABLE from 0_shared_config.R (sourced above) — avoids redundant PATH scan
.use_dt <- .HAS_DATATABLE

# Use SAMPLE_IDS from shared config (0_shared_config.R)
# Override here if needed for method-specific samples
if (length(SAMPLE_IDS) == 0) {
  stop("No samples loaded. Check SRR_COMBINED_LIST_STR and SRR_csv files.")
}

# Create output directories
output_dir <- file.path(MATRICES_OUTPUT_DIR, MASTER_REFERENCE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("\n", strrep("=", 70), "\n")
cat("TXIMPORT: SALMON QUANTIFICATION TO MATRICES\n")
cat(strrep("=", 70), "\n\n")

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

  cat(strrep("=", 70), "\n")
  cat("PROCESSING:", level_config$label, "\n")
  cat(strrep("=", 70), "\n\n")

  # ===============================================
  # STEP 1: LOCATE SALMON OUTPUT FILES
  # ===============================================

  cat("Step 1: Locating Salmon", level_config$label, "output files...\n")

  # If QUANT_DIR already includes the fasta_tag (from SALMON_QUANT_ROOT), use as-is.
  # Otherwise append MASTER_REFERENCE (which equals the fasta_tag in the normal run).
  salmon_quant_dir <- if (QUANT_DIR_INCLUDES_REF) QUANT_DIR else file.path(QUANT_DIR, MASTER_REFERENCE)

  if (!dir.exists(salmon_quant_dir)) {
    cat("ERROR: Salmon quantification directory not found:", salmon_quant_dir, "\n")
    cat("Run the M4 Salmon SAF alignment stage first, or check SALMON_QUANT_ROOT / MASTER_REFERENCE.\n")
    cat("Skipping this level...\n\n")
    next
  }

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
    INPUT_FASTAS_DIR <- if (nzchar(BASE_DIR)) {
      file.path(BASE_DIR, "inputs")
    } else if (nzchar(Sys.getenv("WF_MANAGED_ENV", ""))) {
      stop("[SALMON TXIMPORT] INPUT_FASTAS_DIR or BASE_DIR is required under workflow manager.")
    } else {
      message("[SALMON TXIMPORT] WARN: INPUT_FASTAS_DIR and BASE_DIR not set; using relative '../../inputs'. ",
              "Set INPUT_FASTAS_DIR or BASE_DIR for orchestrated execution.")
      "../../inputs"
    }
  }

  # Search for the tx2gene mapping file produced by create_gene_trans_map()
  # (saved as <fasta>.gene_trans_map next to the FASTA, or exportable via GENE_TRANS_MAP_FILE)
  tx2gene_candidates <- c(
    Sys.getenv("GENE_TRANS_MAP_FILE", unset = ""),
    file.path(INPUT_FASTAS_DIR, "mapping", paste0(MASTER_REFERENCE, ".fa.gene_trans_map")),
    file.path(INPUT_FASTAS_DIR, "mapping", paste0(MASTER_REFERENCE, ".fasta.gene_trans_map")),
    file.path(INPUT_FASTAS_DIR, "fasta", paste0(MASTER_REFERENCE, ".fa.gene_trans_map")),
    file.path(INPUT_FASTAS_DIR, "fasta", "reference_genome", paste0(MASTER_REFERENCE, ".fa.gene_trans_map")),
    file.path(INPUT_FASTAS_DIR, "fasta", "reference_genome", paste0(MASTER_REFERENCE, ".fasta.gene_trans_map"))
  )

  # Vectorized candidate check: single batch stat instead of sequential loop  O(1) syscall batch
  tx2gene_file <- NULL
  .nonempty <- nzchar(tx2gene_candidates)
  if (any(.nonempty)) {
    .found <- .nonempty & file.exists(tx2gene_candidates)
    if (any(.found)) tx2gene_file <- tx2gene_candidates[which(.found)[1L]]
  }
  # Lazy fallback: only recurse directory tree if direct paths failed
  if (is.null(tx2gene_file) && nzchar(INPUT_FASTAS_DIR) && dir.exists(INPUT_FASTAS_DIR)) {
    all_maps <- list.files(INPUT_FASTAS_DIR, pattern = "\\.gene_trans_map$",
                           recursive = TRUE, full.names = TRUE)
    # Match on basename to avoid substring false positives (e.g., "V4" matching "V4.1")
    ref_maps <- all_maps[grepl(paste0("(^|[/\\\\])", MASTER_REFERENCE, "\\."), all_maps)]
    if (length(ref_maps) > 0) tx2gene_file <- ref_maps[1]
  }
  suppressWarnings(rm(".nonempty", ".found"))

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
    # fread fast path: 5-10x faster for large tx2gene files (50K+ transcripts)
    tx2gene <- if (.use_dt) {
      data.table::fread(tx2gene_file, header = FALSE, sep = "\t",
                        strip.white = TRUE, showProgress = FALSE, data.table = FALSE)
    } else {
      read.table(tx2gene_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                 strip.white = TRUE)
    }
    if (nrow(tx2gene) == 0) {
      cat("ERROR: tx2gene file is empty (0 rows):", tx2gene_file, "\n")
      cat("Skipping gene-level processing...\n\n")
      next
    }
    if (ncol(tx2gene) < 2) {
      cat("ERROR: tx2gene file must have at least 2 tab-separated columns, found", ncol(tx2gene), "\n")
      cat("  File:", tx2gene_file, "\n")
      cat("Skipping gene-level processing...\n\n")
      next
    }
    # Keep only first 2 columns (gene_id, transcript_id)
    tx2gene <- tx2gene[, 1:2, drop = FALSE]
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
    txi <- tryCatch(
      tximport(files, type = "salmon",
               txIn = TRUE, txOut = FALSE,
               tx2gene = tx2gene,
               ignoreTxVersion = FALSE, ignoreAfterBar = FALSE),
      error = function(e) {
        message("ERROR: tximport failed for ", level_config$label, ": ", conditionMessage(e))
        NULL
      }
    )
  } else {
    # Isoform-level: no tx2gene needed
    cat("  Importing at transcript level...\n")
    txi <- tryCatch(
      tximport(files, type = "salmon", txIn = TRUE, txOut = TRUE,
               ignoreTxVersion = FALSE, ignoreAfterBar = FALSE),
      error = function(e) {
        message("ERROR: tximport failed for ", level_config$label, ": ", conditionMessage(e))
        NULL
      }
    )
  }

  if (is.null(txi)) {
    cat("Skipping level:", level_name, "\n\n")
    next
  }
  if (nrow(txi$counts) == 0 || ncol(txi$counts) == 0) {
    message("ERROR: tximport returned empty counts matrix (",
            nrow(txi$counts), " rows x ", ncol(txi$counts), " cols)")
    cat("Skipping level:", level_name, "\n\n")
    next
  }

  # Filter entries with zero effective length (prevents NaN/Inf in DESeq2 normalization)
  if (!is.null(txi$length)) {
    zero_mask <- rowSums(txi$length == 0) > 0
    if (any(zero_mask)) {
      entity_type_tmp <- if (level_config$tx_out) "transcripts" else "genes"
      cat("  Filtering", sum(zero_mask), entity_type_tmp, "with zero effective length\n")
      txi$counts    <- txi$counts[!zero_mask, , drop = FALSE]
      txi$abundance <- txi$abundance[!zero_mask, , drop = FALSE]
      txi$length    <- txi$length[!zero_mask, , drop = FALSE]
      if (nrow(txi$counts) == 0) {
        cat("  ERROR: All entries removed after zero-length filtering — skipping level\n")
        next
      }
    }
  }

  entity_type <- if (level_config$tx_out) "transcripts" else "genes"
  cat("Successfully imported data for", ncol(txi$counts), "samples\n")
  cat("Total", entity_type, ":", nrow(txi$counts), "\n\n")

  # Save full tximport object for DESeq2 (preserves transcript-length offsets)
  {
    txi_rds_dir <- file.path(output_dir, level_name)
    dir.create(txi_rds_dir, recursive = TRUE, showWarnings = FALSE)
    rds_name <- if (level_config$tx_out) "tximport_isoform_level.rds" else "tximport_gene_level.rds"
    tryCatch({
      saveRDS(txi, file.path(txi_rds_dir, rds_name))
      cat("Saved tximport RDS for DESeq2:", rds_name, "\n\n")
    }, error = function(e) cat("  Warning: Failed to save tximport RDS:", e$message, "\n"))
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
  # build_input_path() which always expects: {level}/{folder_name}/{folder_name}_...csv
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

  # Reuse pre-computed gene group file list (hoisted above level loop to avoid
  # repeated list.files() filesystem traversals — same result across levels)
  if (!exists(".gg_files_cached")) {
    .gg_files_cached <- if (dir.exists(GENE_GROUPS_DIR)) {
      list.files(GENE_GROUPS_DIR, pattern = "\\.(csv|txt|tsv)$", recursive = TRUE, full.names = TRUE)
    } else character(0)
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

    successful_groups <- 0
    # .use_dt already cached at top of script (avoids per-iteration PATH scan)

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
          gene_df <- if (.use_dt) {
            data.table::fread(gene_group_file, header = TRUE, data.table = FALSE)
          } else {
            read.csv(gene_group_file, stringsAsFactors = FALSE, header = TRUE)
          }
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
        cat("    Error reading file:", conditionMessage(e), "\n")
        character(0)
      })
      # Filter empty strings/NAs after trimming (matches STAR pattern at tximport_star_to_matrices.R)
      gene_list <- gene_list[nzchar(gene_list) & !is.na(gene_list)]

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
