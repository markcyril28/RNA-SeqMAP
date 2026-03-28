#!/usr/bin/env Rscript

# ===============================================
# TXIMPORT: STAR+SALMON TO DESEQ2 MATRICES
# ===============================================
# Processes Salmon quantification output from STAR alignment (M3) using tximport.
# Produces standardized count matrices consumed by all downstream analysis modules.
#
# Runs from: 3_POST_PROC/M3_STAR_Align/   (via pushd in run_method_analysis)
# Quant files: ../../2_ALIGNMENT_RESULTs/M3_STAR_Align/{MASTER_REFERENCE}/6_salmon/quant/{SRR_ID}/quant.sf
# tx2gene:     count_matrices_from_STAR/{MASTER_REFERENCE}/tx2gene_{MASTER_REFERENCE}.tsv
# Output:      count_matrices_from_STAR/{MASTER_REFERENCE}/{level}/{gene_group}/
#
# Output naming convention (matches build_input_path() in 0_shared_config.R):
#   {gene_group}_{count_type}_{gene_type}_from_{master_ref}_{processing_level}.csv
# ===============================================

suppressPackageStartupMessages({
  library(tximport)
})

# ===============================================
# CONFIGURATION
# ===============================================

SCRIPT_DIR <- Sys.getenv("ANALYSIS_MODULES_DIR", {
  if (nzchar(Sys.getenv("WF_MANAGED_ENV", "")))
    stop("[STAR_TXIMPORT] ANALYSIS_MODULES_DIR is required under workflow manager (WF_MANAGED_ENV is set).")
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
    stop("[STAR TXIMPORT] match_gene_ids.R not found at ", .match_ids_path,
         ". Ensure ANALYSIS_MODULES_DIR points to a directory whose sibling 'utilities/' contains match_gene_ids.R.")
  }
}

# Salmon quant output is in the alignment results directory, NOT post-proc.
# Path includes MASTER_REFERENCE (= fasta_tag) to isolate per-reference outputs.
# Use BASE_DIR (absolute path set by run_post_processing.sh) when available;
# fall back to relative path for standalone usage (script runs from 3_POST_PROC/M3_STAR_Align/).
.base_dir     <- Sys.getenv("BASE_DIR", unset = "")
QUANT_DIR     <- if (nzchar(.base_dir)) {
  file.path(.base_dir, "2_ALIGNMENT_RESULTs", "M3_STAR_Align",
            MASTER_REFERENCE, "6_salmon", "quant")
} else if (nzchar(Sys.getenv("WF_MANAGED_ENV", ""))) {
  stop("[STAR TXIMPORT] BASE_DIR is required when running under a workflow manager (WF_MANAGED_ENV is set).")
} else {
  message("[STAR TXIMPORT] WARN: BASE_DIR not set; using relative path '../..'. ",
          "Set BASE_DIR for orchestrated execution (Nextflow/Snakemake).")
  file.path("..", "..", "2_ALIGNMENT_RESULTs", "M3_STAR_Align",
            MASTER_REFERENCE, "6_salmon", "quant")
}
MATRICES_DIR  <- if (nzchar(.base_dir)) {
  file.path(.base_dir, "3_POST_PROC", "M3_STAR_Align", "count_matrices_from_STAR")
} else if (nzchar(Sys.getenv("WF_MANAGED_ENV", ""))) {
  stop("[STAR TXIMPORT] BASE_DIR is required for output directory under workflow manager (WF_MANAGED_ENV is set).")
} else {
  "count_matrices_from_STAR"  # relative fallback when called from pushd context
}
# Note: MASTER_REFERENCE is already set by 0_shared_config.R above; no re-assignment needed.

GENERATE_GENE_LEVEL     <- isTRUE(as.logical(Sys.getenv("STAR_GENERATE_GENE_LEVEL",    "TRUE")))
GENERATE_ISOFORM_LEVEL  <- isTRUE(as.logical(Sys.getenv("STAR_GENERATE_ISOFORM_LEVEL", "TRUE")))

# ===============================================
# HELPER: SAVE MATRICES WITH STANDARD NAMING
# ===============================================
# Delegates to save_count_matrices() from 3_Matrix_Creation_utils.R to avoid
# duplicate matrix-saving logic that can drift out of sync.

save_count_matrix <- function(counts_matrix, output_dir, base_name, master_ref,
                               level_suffix, tpm_matrix = NULL) {
  processing_level <- gsub("^_", "", level_suffix)
  save_count_matrices(counts_matrix, output_dir, prefix = base_name,
                      master_ref = master_ref, level = processing_level,
                      tpm = tpm_matrix, count_type_label = "NumReads")
}

# ===============================================
# LOCATE TX2GENE MAPPING
# ===============================================

find_tx2gene <- function(matrices_dir, master_ref) {
  # Created by star_alignment_pipeline() in m3_star_alignment.sh
  candidate <- file.path(matrices_dir, master_ref, paste0("tx2gene_", master_ref, ".tsv"))
  if (file.exists(candidate)) return(candidate)

  # Fallback: search for any tx2gene file under the master_ref directory
  search_dir <- file.path(matrices_dir, master_ref)
  if (dir.exists(search_dir)) {
    hits <- list.files(search_dir, pattern = "^tx2gene.*\\.tsv$", full.names = TRUE)
    if (length(hits) > 0) {
      cat("  Using tx2gene:", hits[1], "\n")
      return(hits[1])
    }
  }

  # Legacy fallback: gene_trans_map files (created by M4 Salmon, not M3 STAR).
  # Kept for cross-method compatibility when tx2gene was not generated by STAR pipeline.
  .input_env <- Sys.getenv("INPUT_FASTAS_DIR", "")
  input_dir <- if (nzchar(.input_env)) .input_env else if (nzchar(Sys.getenv("WF_MANAGED_ENV", ""))) {
    stop("[STAR TXIMPORT] INPUT_FASTAS_DIR is required under workflow manager.")
  } else {
    message("[STAR TXIMPORT] WARN: INPUT_FASTAS_DIR not set; using relative '../inputs'. ",
            "Set INPUT_FASTAS_DIR or BASE_DIR for orchestrated execution.")
    file.path("..", "..", "inputs")
  }
  # Vectorized batch stat: O(1) syscall batch instead of sequential checks
  .alts <- file.path(input_dir, "mapping",
                     paste0(master_ref, c(".fa.gene_trans_map", ".fasta.gene_trans_map")))
  .hit <- which(file.exists(.alts))[1L]
  if (!is.na(.hit)) return(.alts[.hit])

  return(NULL)
}

# ===============================================
# BANNER
# ===============================================

cat("\n", strrep("=", 70), "\n")
cat("TXIMPORT: STAR+SALMON QUANTIFICATION TO MATRICES\n")
cat(strrep("=", 70), "\n\n")
cat("Master Reference:", MASTER_REFERENCE, "\n")
cat("Quant directory: ", QUANT_DIR, "\n")
cat("Output directory:", MATRICES_DIR, "\n\n")

# ===============================================
# DEFINE PROCESSING LEVELS
# ===============================================

processing_levels <- list()
if (GENERATE_GENE_LEVEL) {
  processing_levels[["gene_level"]] <- list(
    tx_out        = FALSE,
    label         = "Gene-Level",
    output_suffix = "_gene_level"
  )
}
if (GENERATE_ISOFORM_LEVEL) {
  processing_levels[["isoform_level"]] <- list(
    tx_out        = TRUE,
    label         = "Isoform-Level",
    output_suffix = "_isoform_level"
  )
}

if (length(processing_levels) == 0) {
  stop("At least one of GENERATE_GENE_LEVEL or GENERATE_ISOFORM_LEVEL must be TRUE")
}

# ===============================================
# PROCESS EACH LEVEL
# ===============================================

for (level_name in names(processing_levels)) {
  level_config <- processing_levels[[level_name]]

  cat(strrep("=", 70), "\n")
  cat("PROCESSING:", level_config$label, "\n")
  cat(strrep("=", 70), "\n\n")

  # -------------------------------------------------
  # STEP 1: LOCATE QUANT.SF FILES
  # -------------------------------------------------
  cat("Step 1: Locating Salmon quant.sf files...\n")

  # Try flat layout first: QUANT_DIR/{SRR}/quant.sf
  # If not found, check tissue-specific layout: QUANT_DIR/{tissue}/{SRR}/quant.sf
  files <- file.path(QUANT_DIR, SAMPLE_IDS, "quant.sf")
  names(files) <- SAMPLE_IDS

  if (sum(file.exists(files)) == 0) {
    cat("  No quant.sf in flat layout; scanning tissue subdirectories...\n")
    tissue_dirs <- list.dirs(QUANT_DIR, recursive = FALSE, full.names = TRUE)
    if (length(tissue_dirs) > 0) {
      # Vectorized file.exists: build all candidate paths at once via outer(),
      # then batch-check existence. O(T × S) file.exists calls in one vector
      # vs O(T × S) individual calls in nested loop. Matches .resolve_quant_files().
      all_cands <- outer(tissue_dirs, SAMPLE_IDS, function(td, sid) file.path(td, sid, "quant.sf"))
      all_exist <- matrix(file.exists(all_cands), nrow = length(tissue_dirs))
      # Vectorized: max.col() finds first TRUE per column in one C-level pass
      any_found <- colSums(all_exist) > 0
      if (any(any_found)) {
        hits <- max.col(t(all_exist), ties.method = "first")
        idx <- which(any_found)
        files[SAMPLE_IDS[idx]] <- all_cands[cbind(hits[idx], idx)]
      }
    }
  }

  files_exist    <- file.exists(files)
  current_ids    <- SAMPLE_IDS

  if (sum(files_exist) == 0) {
    cat("ERROR: No quant.sf files found under", QUANT_DIR, "\n")
    cat("Ensure STAR+Salmon alignment completed successfully.\n\n")
    next
  }
  if (sum(files_exist) < length(files)) {
    missing <- SAMPLE_IDS[!files_exist]
    cat("WARNING: Missing quantification for:", paste(missing, collapse = ", "), "\n")
    files       <- files[files_exist]
    current_ids <- SAMPLE_IDS[files_exist]
  }
  if (length(files) < 2) {
    cat("ERROR: Need >= 2 quant.sf files for tximport (found", length(files), ")\n\n")
    next
  }
  cat("Found", length(files), "quant.sf files\n\n")

  # -------------------------------------------------
  # STEP 2: LOAD TX2GENE MAPPING (gene-level only)
  # -------------------------------------------------
  tx2gene <- NULL
  if (!level_config$tx_out) {
    cat("Step 2: Loading tx2gene mapping...\n")
    tx2gene_file <- find_tx2gene(MATRICES_DIR, MASTER_REFERENCE)

    if (is.null(tx2gene_file)) {
      cat("ERROR: tx2gene mapping not found. Cannot summarize to gene level.\n")
      cat("Run STAR+Salmon alignment first to generate tx2gene mapping.\n\n")
      next
    }

    # Detect column order: tximport needs c(TXNAME, GENEID)
    # O(N) where N = tx2gene rows (50K-200K for full transcriptomes); fread is 10-50x faster
    # Reuse .HAS_DATATABLE from 0_shared_config.R (sourced at line 27) — avoids redundant PATH scan
    .use_dt_tx2gene <- .HAS_DATATABLE
    raw <- if (.use_dt_tx2gene) {
      data.table::fread(tx2gene_file, header = FALSE, colClasses = "character",
                        data.table = FALSE)
    } else {
      read.delim(tx2gene_file, header = FALSE, stringsAsFactors = FALSE,
                 colClasses = "character")
    }
    if (nrow(raw) == 0) {
      cat("ERROR: tx2gene file is empty (0 rows):", tx2gene_file, "\n")
      next
    }
    if (ncol(raw) < 2) {
      cat("ERROR: tx2gene file has", ncol(raw), "column(s), expected >= 2:", tx2gene_file, "\n")
      next
    }
    # star_alignment_pipeline writes: transcript_id TAB gene_id  (col1=TX, col2=GENE)
    # gene_trans_map fallback writes: gene_id TAB transcript_id  (col1=GENE, col2=TX)
    # Detect: if file is a gene_trans_map, swap columns
    if (grepl("gene_trans_map$", tx2gene_file)) {
      tx2gene <- raw[, c(2, 1), drop = FALSE]
    } else {
      tx2gene <- raw[, 1:2, drop = FALSE]
    }
    colnames(tx2gene) <- c("TXNAME", "GENEID")
    tx2gene$TXNAME <- trimws(tx2gene$TXNAME)
    tx2gene$GENEID <- trimws(tx2gene$GENEID)
    cat("Loaded tx2gene:", nrow(tx2gene), "entries\n\n")

    # Validate tx2gene transcript IDs match quant.sf transcript IDs
    # Use fread when available (5-10x faster for header-only sampling)
    sample_qf <- if (.use_dt_tx2gene) {
      data.table::fread(files[1], header = TRUE, nrows = 100, data.table = FALSE)
    } else {
      read.delim(files[1], header = TRUE, nrows = 100, stringsAsFactors = FALSE)
    }
    qf_ids <- sample_qf$Name
    tx_ids <- tx2gene$TXNAME
    # Use %in% instead of intersect(): avoids allocating the intersection vector,
    # calling unique(), and length(). Both use match() internally with same hash cost.
    overlap <- sum(qf_ids %in% tx_ids)
    match_rate <- if (length(qf_ids) > 0) overlap / length(qf_ids) else 0
    if (match_rate < 0.5) {
      cat("ERROR: tx2gene transcript IDs poorly match quant.sf IDs!\n")
      cat("  Match rate:", round(match_rate * 100), "% (", overlap, "/", length(qf_ids), "sampled)\n")
      cat("  tx2gene IDs (first 3):", paste(head(tx2gene$TXNAME, 3), collapse = ", "), "\n")
      cat("  quant.sf IDs (first 3):", paste(head(sample_qf$Name, 3), collapse = ", "), "\n")
      cat("  This usually means tx2gene was generated from the wrong GTF.\n")
      cat("  Re-run STAR+Salmon alignment to regenerate tx2gene.\n")
      cat("  Skipping gene-level import.\n\n")
      next
    }
    cat("  tx2gene/quant.sf ID check: PASS (", round(match_rate * 100), "%, ", overlap, "/", length(qf_ids), " sampled IDs match)\n\n")
  }

  # -------------------------------------------------
  # STEP 3: TXIMPORT
  # -------------------------------------------------
  cat("Step 3: Importing with tximport...\n")

  txi <- tryCatch({
    if (!level_config$tx_out) {
      tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = FALSE)
    } else {
      tximport(files, type = "salmon", txIn = TRUE, txOut = TRUE,
               ignoreTxVersion = FALSE, ignoreAfterBar = FALSE)
    }
  }, error = function(e) {
    message("ERROR in tximport: ", conditionMessage(e))
    NULL
  })

  if (is.null(txi)) {
    cat("Skipping level:", level_name, "\n\n")
    next
  }

  # Filter entries with zero effective length (prevents NaN/Inf in DESeq2 normalization)
  if (!is.null(txi$length)) {
    zero_mask <- rowSums(txi$length == 0) > 0
    if (any(zero_mask)) {
      cat("  Filtering", sum(zero_mask), "entries with zero effective length\n")
      txi$counts    <- txi$counts[!zero_mask, , drop = FALSE]
      txi$abundance <- txi$abundance[!zero_mask, , drop = FALSE]
      txi$length    <- txi$length[!zero_mask, , drop = FALSE]
    }
  }

  if (nrow(txi$counts) == 0 || ncol(txi$counts) == 0) {
    message("ERROR: tximport returned empty counts matrix (",
            nrow(txi$counts), " rows x ", ncol(txi$counts), " cols)")
    cat("Skipping level:", level_name, "\n\n")
    next
  }

  entity_type <- if (level_config$tx_out) "transcripts" else "genes"
  cat("Imported:", ncol(txi$counts), "samples,", nrow(txi$counts), entity_type, "\n\n")

  # Save full tximport object for DESeq2 (preserves transcript-length offsets)
  txi_rds_dir <- file.path(MATRICES_DIR, MASTER_REFERENCE, level_name)
  dir.create(txi_rds_dir, recursive = TRUE, showWarnings = FALSE)
  if (!level_config$tx_out) {
    saveRDS(txi, file.path(txi_rds_dir, "tximport_gene_level.rds"))
    cat("Saved tximport RDS for DESeq2: tximport_gene_level.rds\n\n")
  } else {
    saveRDS(txi, file.path(txi_rds_dir, "tximport_isoform_level.rds"))
    cat("Saved tximport RDS for isoform-level DESeq2: tximport_isoform_level.rds\n\n")
  }

  # -------------------------------------------------
  # STEP 4: SAMPLE METADATA
  # -------------------------------------------------
  cat("Step 4: Building sample metadata...\n")
  conditions <- SAMPLE_LABELS[current_ids]
  missing_labels <- is.na(conditions)
  if (any(missing_labels)) {
    cat("WARNING: No labels for:", paste(current_ids[missing_labels], collapse = ", "), "\n")
    conditions[missing_labels] <- current_ids[missing_labels]
  }
  sample_data <- data.frame(
    SampleID  = current_ids,
    Condition = unname(conditions),
    row.names = current_ids,
    stringsAsFactors = FALSE
  )
  cat("Conditions per tissue:\n")
  print(table(sample_data$Condition))
  cat("\n")

  # -------------------------------------------------
  # STEP 5: EXTRACT MATRICES
  # -------------------------------------------------
  raw_counts <- txi$counts
  tpm_matrix <- txi$abundance

  # -------------------------------------------------
  # STEP 6: SAVE FULL REFERENCE MATRIX
  # -------------------------------------------------
  cat("Step 6: Saving full reference matrices...\n")
  level_output_dir <- file.path(MATRICES_DIR, MASTER_REFERENCE, level_name)
  dir.create(level_output_dir, recursive = TRUE, showWarnings = FALSE)

  # Dataset-namespaced folder (matches gene group structure; avoids overwrite across datasets)
  # CURRENT_DATASET is already set by 0_shared_config.R from the same env var
  full_ref_folder <- if (nzchar(CURRENT_DATASET)) {
    paste0(MASTER_REFERENCE, "_in_", CURRENT_DATASET)
  } else {
    MASTER_REFERENCE
  }
  full_ref_dir <- file.path(level_output_dir, full_ref_folder)
  dir.create(full_ref_dir, recursive = TRUE, showWarnings = FALSE)

  save_count_matrix(raw_counts, full_ref_dir,
                    base_name    = full_ref_folder,
                    master_ref   = MASTER_REFERENCE,
                    level_suffix = level_config$output_suffix,
                    tpm_matrix   = tpm_matrix)
  cat("\n")

  # -------------------------------------------------
  # STEP 7: PROCESS GENE GROUPS
  # -------------------------------------------------
  cat("Step 7: Processing gene groups...\n")

  # Reuse pre-computed gene group file list (hoisted above level loop to avoid
  # repeated list.files() filesystem traversals — same result across levels)
  if (!exists(".gg_files_cached")) {
    .gg_files_cached <- if (dir.exists(GENE_GROUPS_DIR)) {
      list.files(GENE_GROUPS_DIR, pattern = "\\.(csv|txt|tsv)$",
                 recursive = TRUE, full.names = TRUE)
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
      cat("Gene groups to process:", paste(enabled_groups, collapse = ", "), "\n")
    }
  }
  gene_group_files <- .gg_files_cached

  if (length(gene_group_files) == 0) {
    cat("No matching gene group files found in", GENE_GROUPS_DIR, "\n")
  } else {
    cat("Found", length(gene_group_files), "gene group file(s)\n\n")

    successful_groups <- 0
    # Reuse .HAS_DATATABLE from 0_shared_config.R — avoids redundant PATH scan
    .use_dt <- .HAS_DATATABLE

    for (gene_group_file in gene_group_files) {
      gene_group_name  <- tools::file_path_sans_ext(basename(gene_group_file))
      output_folder_name <- if (nzchar(CURRENT_DATASET)) {
        paste0(gene_group_name, "_in_", CURRENT_DATASET)
      } else {
        gene_group_name
      }
      cat("  Processing:", gene_group_name, "->", output_folder_name, "\n")

      # Read gene list
      gene_list <- tryCatch({
        if (grepl("\\.csv$", gene_group_file, ignore.case = TRUE)) {
          gdf <- if (.use_dt) {
            data.table::fread(gene_group_file, header = TRUE, data.table = FALSE)
          } else {
            read.csv(gene_group_file, stringsAsFactors = FALSE, header = TRUE)
          }
          if ("Gene_ID" %in% colnames(gdf)) gdf$Gene_ID else gdf[[1]]
        } else {
          raw_lines <- suppressWarnings(readLines(gene_group_file))
          raw_lines <- raw_lines[
            !grepl("^#|^Gene_ID", raw_lines, ignore.case = TRUE) & nzchar(raw_lines)
          ]
          sub("\t.*", "", raw_lines)  # strip tab-delimited extra fields (e.g. gene names, descriptions)
        }
      }, error = function(e) {
        cat("    Error reading gene list:", conditionMessage(e), "\n")
        character(0)
      })
      gene_list <- trimws(gene_list[nzchar(gene_list)])

      if (length(gene_list) == 0) {
        cat("    Skipping: empty gene list\n")
        next
      }

      # Match genes using shared helper (handles version suffixes like GENE.1, GENE.1.01)
      genes_in_data <- match_gene_ids(gene_list, rownames(raw_counts))

      if (length(genes_in_data) == 0) {
        cat("    No matching", entity_type, "found\n")
        next
      }
      cat("    Matched", length(genes_in_data), "/", length(gene_list), entity_type, "\n")

      # Save gene-group subset
      gene_group_dir <- file.path(level_output_dir, output_folder_name)
      dir.create(gene_group_dir, recursive = TRUE, showWarnings = FALSE)

      subset_counts <- raw_counts[genes_in_data, , drop = FALSE]
      subset_tpm    <- tpm_matrix[genes_in_data, , drop = FALSE]

      save_count_matrix(subset_counts, gene_group_dir,
                        base_name    = output_folder_name,
                        master_ref   = MASTER_REFERENCE,
                        level_suffix = level_config$output_suffix,
                        tpm_matrix   = subset_tpm)
      successful_groups <- successful_groups + 1
    }

    cat("\n  Summary:", successful_groups, "/", length(gene_group_files),
        "gene groups processed\n")
  }

  cat("\n", strrep("=", 70), "\n")
  cat(level_config$label, "COMPLETE\n")
  cat(strrep("=", 70), "\n\n")
}

cat(strrep("=", 70), "\n")
cat("ALL LEVELS COMPLETE\n")
cat("Output directory:", file.path(MATRICES_DIR, MASTER_REFERENCE), "\n")
cat(strrep("=", 70), "\n\n")
