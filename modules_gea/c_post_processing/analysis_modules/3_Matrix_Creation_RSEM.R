#!/usr/bin/env Rscript

# ===============================================
# MATRIX CREATION - M5 RSEM/Bowtie2
# ===============================================
# Creates count matrices from RSEM quantification outputs using tximport.
# Expects environment variables set by pipeline_utils.sh:
#   CURRENT_METHOD, MASTER_REFERENCE, BASE_DIR,
#   SRR_COMBINED_LIST_STR, GENE_GROUPS_STR, GENE_GROUPS_DIR,
#   RSEM_GENERATE_GENE_LEVEL, RSEM_GENERATE_ISOFORM_LEVEL

suppressPackageStartupMessages(library(tximport))

SCRIPT_DIR <- Sys.getenv("ANALYSIS_MODULES_DIR", {
  if (nzchar(Sys.getenv("WF_MANAGED_ENV", "")))
    stop("[MATRIX_CREATION_RSEM] ANALYSIS_MODULES_DIR is required under workflow manager (WF_MANAGED_ENV is set).")
  "."
})
source(file.path(SCRIPT_DIR, "0_shared_config.R"))
source(file.path(SCRIPT_DIR, "1_utility_functions.R"))
source(file.path(SCRIPT_DIR, "3_Matrix_Creation_utils.R"))

# ===============================================
# RSEM IMPORT
# ===============================================

import_rsem <- function(quant_dir, sample_ids, level = "gene") {
  file_type <- if (level == "gene") ".genes.results" else ".isoforms.results"
  files <- file.path(quant_dir, sample_ids, paste0(sample_ids, file_type))
  names(files) <- sample_ids

  files_exist <- file.exists(files)
  if (sum(files_exist) == 0) {
    cat("ERROR: No", level, "quantification files found in", quant_dir, "\n")
    return(NULL)
  }
  if (sum(files_exist) < length(files)) {
    cat("WARNING: Missing", level, "results for:",
        paste(names(files)[!files_exist], collapse = ", "), "\n")
    files <- files[files_exist]
  }

  txi <- tryCatch(
    tximport(files, type = "rsem", txIn = (level != "gene"), txOut = (level != "gene")),
    error = function(e) {
      cat("ERROR: tximport failed for", level, "level:", e$message, "\n")
      cat("  Check RSEM output files for corruption or format issues\n")
      return(NULL)
    }
  )
  if (is.null(txi)) return(NULL)

  # Filter entries with zero effective length (unaligned transcripts)
  if (!is.null(txi$length)) {
    zero_mask <- rowSums(txi$length == 0) > 0
    n_zero    <- sum(zero_mask)
    if (n_zero > 0) {
      cat("  Filtering", n_zero, level, "entries with zero effective length\n")
      txi$counts    <- txi$counts[!zero_mask, , drop = FALSE]
      txi$abundance <- txi$abundance[!zero_mask, , drop = FALSE]
      txi$length    <- txi$length[!zero_mask, , drop = FALSE]
      if (nrow(txi$counts) == 0) {
        cat("  ERROR: All", level, "entries removed after zero-length filtering\n")
        return(NULL)
      }
    }
  }

  return(txi)
}

# ===============================================
# MAIN
# ===============================================

run_rsem_matrix_creation <- function() {
GENERATE_GENE_LEVEL    <- isTRUE(as.logical(Sys.getenv("RSEM_GENERATE_GENE_LEVEL",    "TRUE")))
GENERATE_ISOFORM_LEVEL <- isTRUE(as.logical(Sys.getenv("RSEM_GENERATE_ISOFORM_LEVEL", "TRUE")))

cat("\n", strrep("=", 60), "\n")
cat("MATRIX CREATION - M5 RSEM/Bowtie2\n")
cat(strrep("=", 60), "\n\n")
cat("Master Reference:", MASTER_REFERENCE, "\n")
cat("Samples:         ", length(SAMPLE_IDS), "\n\n")

if (length(SAMPLE_IDS) == 0) stop("No samples loaded. Check SRR_COMBINED_LIST_STR and SRR_csv files.")

base_dir <- Sys.getenv("BASE_DIR", "")
rsem_quant_root_env <- Sys.getenv("RSEM_QUANT_ROOT", "")
quant_dir <- if (nzchar(rsem_quant_root_env)) {
  rsem_quant_root_env  # already includes fasta_tag
} else if (nzchar(base_dir)) {
  file.path(base_dir, "II_RESULTS/2_ALIGNMENT_RESULTs", "M5_RSEM_Bowtie2", "RSEM_Quant_WD", MASTER_REFERENCE)
} else if (nzchar(Sys.getenv("WF_MANAGED_ENV", ""))) {
  stop("[RSEM MATRIX] BASE_DIR is required when running under a workflow manager (WF_MANAGED_ENV is set). ",
       "Export BASE_DIR pointing to the project root.")
} else {
  message("[RSEM MATRIX] WARN: BASE_DIR and RSEM_QUANT_ROOT not set; using relative 'RSEM_Quant_WD'. ",
          "Set BASE_DIR for orchestrated execution (Nextflow/Snakemake).")
  "RSEM_Quant_WD"  # fallback for standalone execution
}
output_dir  <- if (nzchar(base_dir)) {
  file.path(base_dir, "II_RESULTS", "3_POST_PROC", CURRENT_GENE_GROUP, CURRENT_METHOD, "count_matrices_from_RSEM_Quant")
} else if (nzchar(Sys.getenv("WF_MANAGED_ENV", ""))) {
  stop("[RSEM MATRIX] BASE_DIR is required when running under a workflow manager (WF_MANAGED_ENV is set). ",
       "Export BASE_DIR pointing to the project root.")
} else {
  "count_matrices_from_RSEM_Quant"
}
count_label <- "expected_count"

cat("Quantification directory:", quant_dir, "\n")
cat("Output directory:        ", output_dir, "\n\n")

if (!dir.exists(quant_dir)) {
  stop("RSEM quantification directory not found: ", quant_dir,
       "\n  Run the M5 Bowtie2+RSEM alignment stage first, or check RSEM_QUANT_ROOT / MASTER_REFERENCE.")
}

results <- list()

if (GENERATE_GENE_LEVEL) {
  txi <- import_rsem(quant_dir, SAMPLE_IDS, "gene")
  if (!is.null(txi)) {
    results$gene_level     <- txi$counts
    results$gene_level_tpm <- txi$abundance
    # Save full tximport object for DESeq2 (preserves transcript-length offsets)
    txi_rds_dir <- file.path(output_dir, MASTER_REFERENCE, "gene_level")
    ensure_output_dir(txi_rds_dir)
    tryCatch({
      saveRDS(txi, file.path(txi_rds_dir, "tximport_gene_level.rds"))
      cat("Saved tximport RDS for DESeq2: tximport_gene_level.rds\n")
    }, error = function(e) cat("  Warning: Failed to save tximport RDS:", e$message, "\n"))
  }
}

if (GENERATE_ISOFORM_LEVEL) {
  txi <- import_rsem(quant_dir, SAMPLE_IDS, "isoform")
  if (!is.null(txi)) {
    results$isoform_level     <- txi$counts
    results$isoform_level_tpm <- txi$abundance
    # Save full tximport object for DESeq2 isoform-level analysis
    txi_iso_rds_dir <- file.path(output_dir, MASTER_REFERENCE, "isoform_level")
    ensure_output_dir(txi_iso_rds_dir)
    tryCatch({
      saveRDS(txi, file.path(txi_iso_rds_dir, "tximport_isoform_level.rds"))
      cat("Saved tximport RDS for isoform-level DESeq2: tximport_isoform_level.rds\n")
    }, error = function(e) cat("  Warning: Failed to save isoform tximport RDS:", e$message, "\n"))
  }
}

run_matrix_saving(results, output_dir, MASTER_REFERENCE, count_label, GENE_GROUPS_DIR)
}  # end run_rsem_matrix_creation

# Run if executed directly (not sourced by batch_dispatcher.R)
if (!interactive() && identical(environment(), globalenv()) &&
    !isTRUE(get0(".BATCH_DISPATCHER_ACTIVE"))) {
  run_rsem_matrix_creation()
}
