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
GENERATE_GENE_LEVEL    <- as.logical(Sys.getenv(paste0(.method_prefix, "_GENERATE_GENE_LEVEL"),    "TRUE"))
GENERATE_ISOFORM_LEVEL <- as.logical(Sys.getenv(paste0(.method_prefix, "_GENERATE_ISOFORM_LEVEL"), "TRUE"))

# tximport is only needed for M3/M4/M5 (salmon/RSEM/STAR), not for M1/M2 (StringTie)
.HAS_TXIMPORT <- requireNamespace("tximport", quietly = TRUE)
if (.HAS_TXIMPORT) suppressPackageStartupMessages(library(tximport))

SCRIPT_DIR <- Sys.getenv("ANALYSIS_MODULES_DIR", ".")
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
# MAIN PROCESSING FUNCTION
# ===============================================

run_matrix_creation <- function(method, quant_dir, output_dir, master_ref,
                                 sample_ids, gene_groups_dir = NULL) {
  cat("\n", paste(rep("=", 60), collapse = ""), "\n")
  cat("MATRIX CREATION -", method, "\n")
  cat(paste(rep("=", 60), collapse = ""), "\n\n")

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
    # Prefer 3_Matrix_Creation_STAR.R; this branch is a fallback only.
    if (GENERATE_GENE_LEVEL) {
      tx2gene_file <- list.files(file.path("count_matrices_from_STAR", master_ref),
                                  pattern = "^tx2gene.*\\.tsv$", full.names = TRUE)
      if (length(tx2gene_file) > 0) {
        # Read without col.names to detect actual column count and order.
        # Matches the approach in 3_Matrix_Creation_STAR.R / tximport_star_to_matrices.R.
        raw_tx2gene <- read.delim(tx2gene_file[1], header = FALSE,
                                  stringsAsFactors = FALSE, colClasses = "character")
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
          tx2gene <- raw_tx2gene
          colnames(tx2gene) <- c("TXNAME", "GENEID")[seq_len(ncol(tx2gene))]
        }
        txi_gene <- tryCatch(
          tximport(setNames(file.path(quant_dir, sample_ids, "quant.sf"), sample_ids),
                   type = "salmon", tx2gene = tx2gene, ignoreTxVersion = FALSE),
          error = function(e) NULL)
        if (!is.null(txi_gene)) {
          results$gene_level     <- txi_gene$counts
          results$gene_level_tpm <- txi_gene$abundance
          # Save tximport RDS for DESeq2 (preserves transcript-length offsets)
          txi_rds_dir <- file.path(output_dir, master_ref, "gene_level")
          ensure_output_dir(txi_rds_dir)
          saveRDS(txi_gene, file.path(txi_rds_dir, "tximport_gene_level.rds"))
          cat("  Saved tximport RDS for DESeq2: tximport_gene_level.rds\n")
        }
      } else {
        cat("  Warning: tx2gene file not found for M3 gene-level import\n")
      }
    }
    if (GENERATE_ISOFORM_LEVEL) {
      txi_iso <- tryCatch(
        tximport(setNames(file.path(quant_dir, sample_ids, "quant.sf"), sample_ids),
                 type = "salmon", txIn = TRUE, txOut = TRUE,
                 ignoreTxVersion = FALSE, ignoreAfterBar = FALSE),
        error = function(e) NULL)
      if (!is.null(txi_iso)) {
        results$isoform_level     <- txi_iso$counts
        results$isoform_level_tpm <- txi_iso$abundance
        # Save tximport RDS for isoform-level DESeq2
        txi_iso_rds_dir <- file.path(output_dir, master_ref, "isoform_level")
        ensure_output_dir(txi_iso_rds_dir)
        saveRDS(txi_iso, file.path(txi_iso_rds_dir, "tximport_isoform_level.rds"))
        cat("  Saved tximport RDS for isoform-level DESeq2: tximport_isoform_level.rds\n")
      }
    }

  } else if (grepl("RSEM|M5", method, ignore.case = TRUE)) {
    if (!.HAS_TXIMPORT) { cat("WARNING: tximport not installed, skipping M5 matrix creation.\n"); return(results) }
    # Prefer 3_Matrix_Creation_RSEM.R; this branch is a fallback only.
    if (GENERATE_GENE_LEVEL) {
      txi_gene <- tryCatch(
        tximport(setNames(file.path(quant_dir, sample_ids,
                                    paste0(sample_ids, ".genes.results")), sample_ids),
                 type = "rsem", txIn = FALSE, txOut = FALSE),
        error = function(e) NULL)
      if (!is.null(txi_gene)) {
        results$gene_level     <- txi_gene$counts
        results$gene_level_tpm <- txi_gene$abundance
      }
    }
    if (GENERATE_ISOFORM_LEVEL) {
      txi_iso <- tryCatch(
        tximport(setNames(file.path(quant_dir, sample_ids,
                                    paste0(sample_ids, ".isoforms.results")), sample_ids),
                 type = "rsem", txIn = TRUE, txOut = TRUE),
        error = function(e) NULL)
      if (!is.null(txi_iso)) {
        results$isoform_level     <- txi_iso$counts
        results$isoform_level_tpm <- txi_iso$abundance
      }
    }

  } else if (grepl("Salmon|M4", method, ignore.case = TRUE)) {
    if (!.HAS_TXIMPORT) { cat("WARNING: tximport not installed, skipping M4 matrix creation.\n"); return(results) }
    # Prefer 3_Matrix_Creation_Salmon.R; this branch is a fallback only.
    if (GENERATE_GENE_LEVEL) {
      txi_gene <- tryCatch(
        tximport(setNames(file.path(quant_dir, sample_ids, "quant.sf"), sample_ids),
                 type = "salmon"),
        error = function(e) {
          cat("  Warning: M4 gene-level import failed (tx2gene mapping required).", e$message, "\n")
          cat("  Use 3_Matrix_Creation_Salmon.R for proper gene-level aggregation.\n")
          NULL
        })
      if (!is.null(txi_gene)) {
        results$gene_level     <- txi_gene$counts
        results$gene_level_tpm <- txi_gene$abundance
      }
    }
    if (GENERATE_ISOFORM_LEVEL) {
      txi_iso <- tryCatch(
        tximport(setNames(file.path(quant_dir, sample_ids, "quant.sf"), sample_ids),
                 type = "salmon", txOut = TRUE),
        error = function(e) NULL)
      if (!is.null(txi_iso)) {
        results$isoform_level     <- txi_iso$counts
        results$isoform_level_tpm <- txi_iso$abundance
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

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("MATRIX CREATION MODULE\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")
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
    "rsem"   = file.path(base_dir, "2_ALIGNMENT_RESULTs", "M5_RSEM_Bowtie2",
                         "RSEM_Quant_WD", MASTER_REFERENCE),
    "salmon" = file.path(base_dir, "2_ALIGNMENT_RESULTs", "M4_Salmon_Saf",
                         "Salmon_Quant", MASTER_REFERENCE),
    "star"   = file.path(base_dir, "2_ALIGNMENT_RESULTs", "M3_STAR_Align",
                         MASTER_REFERENCE, "6_salmon", "quant"),
    get_quant_dir(CURRENT_METHOD)  # fallback: relative path for other methods
  )
} else {
  get_quant_dir(CURRENT_METHOD)  # fallback for standalone execution
}

output_dir <- get_matrices_dir(CURRENT_METHOD)

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
