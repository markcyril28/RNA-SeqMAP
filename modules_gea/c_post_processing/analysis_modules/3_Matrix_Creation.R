#!/usr/bin/env Rscript

# ===============================================
# MATRIX CREATION MODULE (fallback dispatcher)
# ===============================================
# Creates count matrices from quantification outputs.
# For M3/M4/M5, pipeline_utils.sh dispatches to method-specific scripts
# (3_Matrix_Creation_STAR.R / _Salmon.R / _RSEM.R).
# This script handles unknown/legacy methods and the M1/M2 StringTie fallback.

# ===============================================
# CONFIGURATION
# ===============================================

# Detect method-appropriate env var for generate level control.
# Method-specific scripts (3_Matrix_Creation_STAR.R, _Salmon.R, _RSEM.R) use
# their own prefixed env vars; this fallback dispatcher checks all of them.
.method_prefix <- switch(
  sub("^(M[0-9]+).*", "\\1", Sys.getenv("CURRENT_METHOD", "M5")),
  "M3" = "STAR", "M4" = "SALMON", "M5" = "RSEM", "RSEM"
)
GENERATE_GENE_LEVEL    <- isTRUE(as.logical(Sys.getenv(paste0(.method_prefix, "_GENERATE_GENE_LEVEL"),    "TRUE")))
GENERATE_ISOFORM_LEVEL <- isTRUE(as.logical(Sys.getenv(paste0(.method_prefix, "_GENERATE_ISOFORM_LEVEL"), "TRUE")))

# tximport is only needed for M3/M4/M5 (salmon/RSEM/STAR), not for M1/M2 (StringTie)
.HAS_TXIMPORT <- requireNamespace("tximport", quietly = TRUE)
if (.HAS_TXIMPORT) suppressPackageStartupMessages(library(tximport))

SCRIPT_DIR <- Sys.getenv("ANALYSIS_MODULES_DIR", {
  if (nzchar(Sys.getenv("WF_MANAGED_ENV", "")))
    stop("[MATRIX_CREATION] ANALYSIS_MODULES_DIR is required under workflow manager (WF_MANAGED_ENV is set).")
  "."
})
source(file.path(SCRIPT_DIR, "0_shared_config.R"))
source(file.path(SCRIPT_DIR, "1_utility_functions.R"))
source(file.path(SCRIPT_DIR, "3_Matrix_Creation_utils.R"))

# ===============================================
# STRINGTIE IMPORT (M1/M2 fallback)
# ===============================================

import_stringtie <- function(ballgown_dir, sample_ids) {
  gene_count_file       <- file.path(ballgown_dir, "gene_count_matrix.csv")
  transcript_count_file <- file.path(ballgown_dir, "transcript_count_matrix.csv")

  result <- list()
  .read_count_csv <- function(path) {
    df <- if (.HAS_DATATABLE) {
      data.table::fread(path, header = TRUE, data.table = FALSE)
    } else {
      read.csv(path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
    }
    rownames(df) <- df[[1]]; df[[1]] <- NULL
    as.matrix(df)
  }
  if (file.exists(gene_count_file)) {
    result$gene <- .read_count_csv(gene_count_file)
  }
  if (file.exists(transcript_count_file)) {
    result$transcript <- .read_count_csv(transcript_count_file)
  }
  return(result)
}

# ===============================================
# QUANT FILE RESOLUTION HELPER
# ===============================================

# Resolve quant.sf paths with tissue-specific directory fallback.
# Big O: O(T × S) where T = tissue dirs, S = samples; file.exists is batched.
# Caches result to avoid redundant filesystem scans when called for both
# gene-level and isoform-level imports within the same run.
.resolve_quant_files <- function(quant_dir, sample_ids) {
  files <- setNames(file.path(quant_dir, sample_ids, "quant.sf"), sample_ids)
  if (sum(file.exists(files)) == 0 && dir.exists(quant_dir)) {
    tissue_dirs <- list.dirs(quant_dir, recursive = FALSE, full.names = TRUE)
    if (length(tissue_dirs) > 0) {
      all_cands <- outer(tissue_dirs, sample_ids, function(td, sid) file.path(td, sid, "quant.sf"))
      all_exist <- matrix(file.exists(all_cands), nrow = length(tissue_dirs))
      # Vectorized: max.col() finds first TRUE per column in one C-level pass
      # O(T×S) single pass vs O(S) which() calls
      any_found <- colSums(all_exist) > 0
      if (any(any_found)) {
        hits <- max.col(t(all_exist), ties.method = "first")
        idx <- which(any_found)
        files[sample_ids[idx]] <- all_cands[cbind(hits[idx], idx)]
      }
    }
  }
  files[file.exists(files)]
}

# ===============================================
# MAIN PROCESSING FUNCTION
# ===============================================

run_matrix_creation <- function(method, quant_dir, output_dir, master_ref,
                                 sample_ids, gene_groups_dir = NULL) {
  cat("\n", strrep("=", 60), "\n")
  cat("MATRIX CREATION -", method, "\n")
  cat(strrep("=", 60), "\n\n")

  results     <- list()
  # Salmon/STAR use "NumReads"; RSEM uses "expected_count"; StringTie/prepDE uses "counts"
  count_label <- if (grepl("Salmon|M4|STAR|M3", method, ignore.case = TRUE)) {
    "NumReads"
  } else if (grepl("HISAT|StringTie|M1|M2", method, ignore.case = TRUE)) {
    "counts"
  } else {
    "expected_count"
  }

  if (grepl("STAR|M3", method, ignore.case = TRUE)) {
    if (!.HAS_TXIMPORT) { cat("WARNING: tximport not installed, skipping M3 matrix creation.\n"); return(results) }
    if (!dir.exists(quant_dir)) {
      cat("ERROR: M3 STAR+Salmon quantification directory not found:", quant_dir,
          "\n  Run the M3 STAR+Salmon alignment stage first, or check BASE_DIR / MASTER_REFERENCE.\n")
      return(results)
    }
    # Prefer 3_Matrix_Creation_STAR.R; this branch is a fallback only.
    # Resolve quant.sf paths once — reused by both gene-level and isoform-level imports.
    # Avoids duplicate filesystem scan when both GENERATE_GENE_LEVEL and GENERATE_ISOFORM_LEVEL are TRUE.
    .m3_quant_files <- if (GENERATE_GENE_LEVEL || GENERATE_ISOFORM_LEVEL) {
      .resolve_quant_files(quant_dir, sample_ids)
    } else {
      character(0)
    }
    if (GENERATE_GENE_LEVEL) {
      .tx2gene_search_dir <- file.path(output_dir, master_ref)
      tx2gene_file <- if (dir.exists(.tx2gene_search_dir)) {
        list.files(.tx2gene_search_dir, pattern = "^tx2gene.*\\.tsv$", full.names = TRUE)
      } else character(0)
      if (length(tx2gene_file) > 0) {
        # Read without col.names to detect actual column count and order.
        # Matches the approach in 3_Matrix_Creation_STAR.R / tximport_star_to_matrices.R.
        # O(N) where N = tx2gene rows (50K-200K); fread is 10-50x faster
        raw_tx2gene <- if (.HAS_DATATABLE) {
          data.table::fread(tx2gene_file[1], header = FALSE, colClasses = "character",
                            data.table = FALSE)
        } else {
          read.delim(tx2gene_file[1], header = FALSE,
                     stringsAsFactors = FALSE, colClasses = "character")
        }
        if (ncol(raw_tx2gene) >= 2) {
          # Detect column order: tximport needs c(TXNAME, GENEID)
          # star_alignment_pipeline writes: transcript_id TAB gene_id (col1=TX, col2=GENE)
          # gene_trans_map fallback writes: gene_id TAB transcript_id (col1=GENE, col2=TX)
          if (grepl("gene_trans_map$", tx2gene_file[1])) {
            tx2gene <- raw_tx2gene[, c(2, 1), drop = FALSE]
          } else {
            tx2gene <- raw_tx2gene[, 1:2, drop = FALSE]
          }
          colnames(tx2gene) <- c("TXNAME", "GENEID")
          tx2gene$TXNAME <- trimws(tx2gene$TXNAME)
          tx2gene$GENEID <- trimws(tx2gene$GENEID)
        } else {
          cat("  Error: tx2gene file has", ncol(raw_tx2gene), "column(s), expected >= 2:",
              tx2gene_file[1], "\n")
          tx2gene <- NULL
        }
        if (is.null(tx2gene)) {
          cat("  Skipping M3 gene-level import (malformed tx2gene)\n")
        } else {
        # quant.sf paths already resolved above (hoisted before gene/isoform blocks)
        if (length(.m3_quant_files) < 2) {
          cat("  Error: Need >= 2 quant.sf files for M3 gene-level import, found",
              length(.m3_quant_files), "\n")
        } else {
        # Validate tx2gene transcript IDs match quant.sf IDs
        .sample_qf <- read.delim(.m3_quant_files[1], header = TRUE, nrows = 100,
                                  stringsAsFactors = FALSE)
        .qf_ids <- .sample_qf$Name
        .overlap <- length(intersect(.qf_ids, tx2gene$TXNAME))
        .match_rate <- if (length(.qf_ids) > 0) .overlap / length(.qf_ids) else 0
        if (.match_rate < 0.5) {
          cat("  ERROR: tx2gene transcript IDs poorly match quant.sf IDs (match rate:",
              round(.match_rate * 100), "%) — skipping M3 gene-level import\n")
        } else {
        txi_gene <- tryCatch(
          tximport(.m3_quant_files,
                   type = "salmon", tx2gene = tx2gene,
                   ignoreTxVersion = FALSE, ignoreAfterBar = FALSE),
          error = function(e) { cat("  M3 gene-level import error:", e$message, "\n"); NULL })
        if (!is.null(txi_gene)) {
          if (nrow(txi_gene$counts) == 0 || ncol(txi_gene$counts) == 0) {
            cat("  WARNING: M3 gene-level import produced empty matrix — skipping\n")
          } else {
            results$gene_level     <- txi_gene$counts
            results$gene_level_tpm <- txi_gene$abundance
            # Save tximport RDS for DESeq2 (preserves transcript-length offsets)
            txi_rds_dir <- file.path(output_dir, master_ref, "gene_level")
            ensure_output_dir(txi_rds_dir)
            saveRDS(txi_gene, file.path(txi_rds_dir, "tximport_gene_level.rds"))
            cat("  Saved tximport RDS for DESeq2: tximport_gene_level.rds\n")
          }
        }
        }  # end if (match_rate >= 0.5)
        }  # end if (length >= 2)
        }  # end if (!is.null(tx2gene))
      } else {
        cat("  Warning: tx2gene file not found for M3 gene-level import\n")
      }
    }
    if (GENERATE_ISOFORM_LEVEL) {
      # quant.sf paths already resolved above (hoisted before gene/isoform blocks)
      if (length(.m3_quant_files) < 2) {
        cat("  Error: Need >= 2 quant.sf files for M3 isoform-level import, found",
            length(.m3_quant_files), "\n")
      } else {
      txi_iso <- tryCatch(
        tximport(.m3_quant_files,
                 type = "salmon", txIn = TRUE, txOut = TRUE,
                 ignoreTxVersion = FALSE, ignoreAfterBar = FALSE),
        error = function(e) { cat("  M3 isoform-level import error:", e$message, "\n"); NULL })
      if (!is.null(txi_iso)) {
        if (nrow(txi_iso$counts) == 0 || ncol(txi_iso$counts) == 0) {
          cat("  WARNING: M3 isoform-level import produced empty matrix — skipping\n")
        } else {
          results$isoform_level     <- txi_iso$counts
          results$isoform_level_tpm <- txi_iso$abundance
          # Save tximport RDS for isoform-level DESeq2
          txi_iso_rds_dir <- file.path(output_dir, master_ref, "isoform_level")
          ensure_output_dir(txi_iso_rds_dir)
          saveRDS(txi_iso, file.path(txi_iso_rds_dir, "tximport_isoform_level.rds"))
          cat("  Saved tximport RDS for isoform-level DESeq2: tximport_isoform_level.rds\n")
        }
      }
      }  # end if (length >= 2)
    }

  } else if (grepl("RSEM|M5", method, ignore.case = TRUE)) {
    if (!.HAS_TXIMPORT) { cat("WARNING: tximport not installed, skipping M5 matrix creation.\n"); return(results) }
    # Prefer 3_Matrix_Creation_RSEM.R; this branch is a fallback only.
    if (GENERATE_GENE_LEVEL) {
      txi_gene <- tryCatch(
        tximport(setNames(file.path(quant_dir, sample_ids,
                                    paste0(sample_ids, ".genes.results")), sample_ids),
                 type = "rsem", txIn = FALSE, txOut = FALSE),
        error = function(e) { cat("  M5 gene-level import error:", e$message, "\n"); NULL })
      if (!is.null(txi_gene)) {
        # Filter entries with zero effective length (RSEM produces these for unaligned genes)
        if (!is.null(txi_gene$length)) {
          zero_mask <- rowSums(txi_gene$length == 0) > 0
          if (any(zero_mask)) {
            cat("  Filtering", sum(zero_mask), "genes with zero effective length\n")
            txi_gene$counts    <- txi_gene$counts[!zero_mask, , drop = FALSE]
            txi_gene$abundance <- txi_gene$abundance[!zero_mask, , drop = FALSE]
            txi_gene$length    <- txi_gene$length[!zero_mask, , drop = FALSE]
          }
        }
        if (nrow(txi_gene$counts) == 0 || ncol(txi_gene$counts) == 0) {
          cat("  WARNING: M5 gene-level import produced empty matrix (possibly all zero-length) — skipping\n")
        } else {
          results$gene_level     <- txi_gene$counts
          results$gene_level_tpm <- txi_gene$abundance
          # Save tximport RDS for DESeq2 (preserves transcript-length offsets)
          txi_rds_dir <- file.path(output_dir, master_ref, "gene_level")
          ensure_output_dir(txi_rds_dir)
          saveRDS(txi_gene, file.path(txi_rds_dir, "tximport_gene_level.rds"))
          cat("  Saved tximport RDS for DESeq2: tximport_gene_level.rds\n")
        }
      }
    }
    if (GENERATE_ISOFORM_LEVEL) {
      txi_iso <- tryCatch(
        tximport(setNames(file.path(quant_dir, sample_ids,
                                    paste0(sample_ids, ".isoforms.results")), sample_ids),
                 type = "rsem", txIn = TRUE, txOut = TRUE),
        error = function(e) { cat("  M5 isoform-level import error:", e$message, "\n"); NULL })
      if (!is.null(txi_iso)) {
        # Filter entries with zero effective length (RSEM produces these for unaligned transcripts)
        if (!is.null(txi_iso$length)) {
          zero_mask <- rowSums(txi_iso$length == 0) > 0
          if (any(zero_mask)) {
            cat("  Filtering", sum(zero_mask), "isoforms with zero effective length\n")
            txi_iso$counts    <- txi_iso$counts[!zero_mask, , drop = FALSE]
            txi_iso$abundance <- txi_iso$abundance[!zero_mask, , drop = FALSE]
            txi_iso$length    <- txi_iso$length[!zero_mask, , drop = FALSE]
          }
        }
        if (nrow(txi_iso$counts) == 0 || ncol(txi_iso$counts) == 0) {
          cat("  WARNING: M5 isoform-level import produced empty matrix (possibly all zero-length) — skipping\n")
        } else {
          results$isoform_level     <- txi_iso$counts
          results$isoform_level_tpm <- txi_iso$abundance
          # Save tximport RDS for isoform-level DESeq2
          txi_iso_rds_dir <- file.path(output_dir, master_ref, "isoform_level")
          ensure_output_dir(txi_iso_rds_dir)
          saveRDS(txi_iso, file.path(txi_iso_rds_dir, "tximport_isoform_level.rds"))
          cat("  Saved tximport RDS for isoform-level DESeq2: tximport_isoform_level.rds\n")
        }
      }
    }

  } else if (grepl("Salmon|M4", method, ignore.case = TRUE)) {
    if (!.HAS_TXIMPORT) { cat("WARNING: tximport not installed, skipping M4 matrix creation.\n"); return(results) }
    if (!dir.exists(quant_dir)) {
      cat("ERROR: M4 Salmon quantification directory not found:", quant_dir,
          "\n  Run the M4 Salmon alignment stage first, or check BASE_DIR / MASTER_REFERENCE.\n")
      return(results)
    }
    # Prefer 3_Matrix_Creation_Salmon.R; this branch is a fallback only.
    if (GENERATE_GENE_LEVEL) {
      # Look up tx2gene mapping (required for gene-level Salmon import)
      .input_fastas_dir <- Sys.getenv("INPUT_FASTAS_DIR", unset = "")
      if (!nzchar(.input_fastas_dir)) {
        .base_dir_fb <- Sys.getenv("BASE_DIR", "")
        if (nzchar(.base_dir_fb)) .input_fastas_dir <- file.path(.base_dir_fb, "I_INPUTS", "inputs", "eggplant")
      }
      .tx2gene_m4 <- NULL
      if (nzchar(.input_fastas_dir)) {
        .cands <- c(
          Sys.getenv("GENE_TRANS_MAP_FILE", unset = ""),
          file.path(.input_fastas_dir, "mapping", paste0(master_ref, ".fa.gene_trans_map")),
          file.path(.input_fastas_dir, "mapping", paste0(master_ref, ".fasta.gene_trans_map"))
        )
        # O(C) where C = candidate paths (≤5); breaks on first valid file
        for (.c in .cands) {
          if (nzchar(.c) && file.exists(.c)) {
            # .rds sidecar cache: 5-10x faster reload vs TSV re-parsing
            .m4_rds <- paste0(.c, ".tx2gene.rds")
            # Batch file.info() for both files: 1 syscall instead of 2 (stat .rds + stat .c)
            .both_info <- file.info(c(.m4_rds, .c))
            if (!is.na(.both_info$mtime[1]) && .both_info$mtime[1] >= .both_info$mtime[2]) {
              .raw <- readRDS(.m4_rds)
              cat("  Loaded tx2gene from .rds cache:", basename(.m4_rds), "\n")
            } else {
              # O(N) where N = gene_trans_map rows; fread is 10-50x faster
              .raw <- if (.HAS_DATATABLE) {
                data.table::fread(.c, header = FALSE, sep = "\t", strip.white = TRUE,
                                  data.table = FALSE)
              } else {
                read.table(.c, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                           strip.white = TRUE)
              }
              tryCatch(saveRDS(.raw, .m4_rds), error = function(e) NULL)
            }
            .raw <- .raw[, 1:2, drop = FALSE]
            colnames(.raw) <- c("GENEID", "TXNAME")
            .tx2gene_m4 <- .raw[, c("TXNAME", "GENEID")]
            cat("  Loaded tx2gene mapping:", nrow(.tx2gene_m4), "entries from", basename(.c), "\n")
            break
          }
        }
        if (is.null(.tx2gene_m4) && dir.exists(.input_fastas_dir)) {
          # Lazy fallback: recurse directory tree
          .all_maps <- list.files(.input_fastas_dir, pattern = "\\.gene_trans_map$",
                                  recursive = TRUE, full.names = TRUE)
          # Use regex with boundary anchors to prevent substring false positives
          # (e.g., "V4" matching "V4.1")
          .hits <- .all_maps[grepl(paste0("(^|[/\\\\])", master_ref, "\\."), .all_maps)]
          if (length(.hits) > 0) {
            # O(N) where N = gene_trans_map rows; fread is 10-50x faster
            .raw <- if (.HAS_DATATABLE) {
              data.table::fread(.hits[1], header = FALSE, sep = "\t", strip.white = TRUE,
                                data.table = FALSE)
            } else {
              read.table(.hits[1], header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                         strip.white = TRUE)
            }
            .raw <- .raw[, 1:2, drop = FALSE]
            colnames(.raw) <- c("GENEID", "TXNAME")
            .tx2gene_m4 <- .raw[, c("TXNAME", "GENEID")]
            cat("  Loaded tx2gene mapping:", nrow(.tx2gene_m4), "entries from", basename(.hits[1]), "\n")
          }
        }
      }
      if (is.null(.tx2gene_m4)) {
        cat("  Warning: tx2gene mapping not found for M4 — skipping gene-level import\n")
        cat("  Use 3_Matrix_Creation_Salmon.R for proper gene-level aggregation.\n")
      } else {
        txi_gene <- tryCatch(
          tximport(setNames(file.path(quant_dir, sample_ids, "quant.sf"), sample_ids),
                   type = "salmon", txIn = TRUE, txOut = FALSE,
                   tx2gene = .tx2gene_m4, ignoreTxVersion = FALSE, ignoreAfterBar = FALSE),
          error = function(e) {
            cat("  Warning: M4 gene-level import failed:", e$message, "\n")
            cat("  Use 3_Matrix_Creation_Salmon.R for proper gene-level aggregation.\n")
            NULL
          })
        if (!is.null(txi_gene)) {
          if (nrow(txi_gene$counts) == 0 || ncol(txi_gene$counts) == 0) {
            cat("  WARNING: M4 gene-level import produced empty matrix — skipping\n")
          } else {
            results$gene_level     <- txi_gene$counts
            results$gene_level_tpm <- txi_gene$abundance
            # Save tximport RDS for DESeq2 (preserves transcript-length offsets)
            txi_rds_dir <- file.path(output_dir, master_ref, "gene_level")
            ensure_output_dir(txi_rds_dir)
            saveRDS(txi_gene, file.path(txi_rds_dir, "tximport_gene_level.rds"))
            cat("  Saved tximport RDS for DESeq2: tximport_gene_level.rds\n")
          }
        }
      }
      # suppressWarnings avoids "object not found" if var was never assigned; skips ls() env scan
      suppressWarnings(rm(".input_fastas_dir", ".base_dir_fb", ".cands", ".c",
                          ".raw", ".all_maps", ".hits", ".tx2gene_m4"))
    }
    if (GENERATE_ISOFORM_LEVEL) {
      txi_iso <- tryCatch(
        tximport(setNames(file.path(quant_dir, sample_ids, "quant.sf"), sample_ids),
                 type = "salmon", txOut = TRUE,
                 ignoreTxVersion = FALSE, ignoreAfterBar = FALSE),
        error = function(e) { cat("  M4 isoform-level import error:", e$message, "\n"); NULL })
      if (!is.null(txi_iso)) {
        if (nrow(txi_iso$counts) == 0 || ncol(txi_iso$counts) == 0) {
          cat("  WARNING: M4 isoform-level import produced empty matrix — skipping\n")
        } else {
          results$isoform_level     <- txi_iso$counts
          results$isoform_level_tpm <- txi_iso$abundance
          # Save tximport RDS for isoform-level DESeq2
          txi_iso_rds_dir <- file.path(output_dir, master_ref, "isoform_level")
          ensure_output_dir(txi_iso_rds_dir)
          saveRDS(txi_iso, file.path(txi_iso_rds_dir, "tximport_isoform_level.rds"))
          cat("  Saved tximport RDS for isoform-level DESeq2: tximport_isoform_level.rds\n")
        }
      }
    }

  } else if (grepl("HISAT|StringTie|M1|M2", method, ignore.case = TRUE)) {
    st_results <- import_stringtie(quant_dir, sample_ids)
    if (!is.null(st_results$gene))       results$gene_level    <- st_results$gene
    if (!is.null(st_results$transcript)) results$isoform_level <- st_results$transcript

  } else {
    cat("WARNING: Unknown method '", method, "' — no import performed.\n", sep = "")
  }

  run_matrix_saving(results, output_dir, master_ref, count_label, gene_groups_dir)
  return(results)
}

# ===============================================
# MAIN EXECUTION
# ===============================================

# Zero-argument entry point for batch_dispatcher.R (M1/M2 fallback path)
run_matrix_creation_main <- function() {
  cat("\n", strrep("=", 60), "\n")
  cat("MATRIX CREATION MODULE\n")
  cat(strrep("=", 60), "\n\n")
  cat("Method:           ", CURRENT_METHOD, "\n")
  cat("Master Reference: ", MASTER_REFERENCE, "\n")
  cat("Samples:          ", length(SAMPLE_IDS), "\n\n")

  if (length(SAMPLE_IDS) == 0) {
    stop("No samples loaded. Check SRR_COMBINED_LIST_STR and SRR_csv files.")
  }

  base_dir    <- Sys.getenv("BASE_DIR", "")
  method_type <- get_method_type(CURRENT_METHOD)
  quant_dir <- if (nzchar(base_dir)) {
    switch(method_type,
      "rsem"   = file.path(base_dir, "II_RESULTS/2_ALIGNMENT_RESULTs", CURRENT_METHOD,
                           "RSEM_Quant_WD", MASTER_REFERENCE),
      "salmon" = file.path(base_dir, "II_RESULTS/2_ALIGNMENT_RESULTs", CURRENT_METHOD,
                           "Salmon_Quant", MASTER_REFERENCE),
      "star"   = file.path(base_dir, "II_RESULTS/2_ALIGNMENT_RESULTs", CURRENT_METHOD,
                           MASTER_REFERENCE, "6_salmon", "quant"),
      get_quant_dir(CURRENT_METHOD)  # fallback: relative path for other methods
    )
  } else {
    get_quant_dir(CURRENT_METHOD)  # fallback for standalone execution
  }

  output_dir <- if (nzchar(base_dir)) {
    file.path(base_dir, "II_RESULTS", "3_POST_PROC", CURRENT_GENE_GROUP, CURRENT_METHOD, get_matrices_dir(CURRENT_METHOD))
  } else if (nzchar(Sys.getenv("WF_MANAGED_ENV", ""))) {
    stop("[MATRIX CREATION] BASE_DIR is required when running under a workflow manager ",
         "(WF_MANAGED_ENV is set). Export BASE_DIR pointing to the project root.")
  } else {
    get_matrices_dir(CURRENT_METHOD)  # relative fallback for standalone pushd context
  }

  cat("Quantification directory:", quant_dir, "\n")
  cat("Output directory:        ", output_dir, "\n\n")

  run_matrix_creation(
    method          = CURRENT_METHOD,
    quant_dir       = quant_dir,
    output_dir      = output_dir,
    master_ref      = MASTER_REFERENCE,
    sample_ids      = SAMPLE_IDS,
    gene_groups_dir = GENE_GROUPS_DIR
  )
}

# Run if executed directly (not sourced by batch_dispatcher.R)
if (!interactive() && identical(environment(), globalenv()) &&
    !isTRUE(get0(".BATCH_DISPATCHER_ACTIVE"))) {
  run_matrix_creation_main()
}
