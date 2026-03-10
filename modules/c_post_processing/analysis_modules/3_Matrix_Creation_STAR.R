#!/usr/bin/env Rscript

# ===============================================
# MATRIX CREATION - M3 STAR + SALMON
# ===============================================
# Creates count matrices from STAR-aligned Salmon quantification outputs.
# Expects environment variables set by pipeline_utils.sh:
#   CURRENT_METHOD, MASTER_REFERENCE, BASE_DIR,
#   SRR_COMBINED_LIST_STR, GENE_GROUPS_STR, GENE_GROUPS_DIR,
#   RSEM_GENERATE_GENE_LEVEL, RSEM_GENERATE_ISOFORM_LEVEL

GENERATE_GENE_LEVEL    <- as.logical(Sys.getenv("RSEM_GENERATE_GENE_LEVEL",    "TRUE"))
GENERATE_ISOFORM_LEVEL <- as.logical(Sys.getenv("RSEM_GENERATE_ISOFORM_LEVEL", "TRUE"))

suppressPackageStartupMessages(library(tximport))

SCRIPT_DIR <- Sys.getenv("ANALYSIS_MODULES_DIR", ".")
source(file.path(SCRIPT_DIR, "0_shared_config.R"))
source(file.path(SCRIPT_DIR, "1_utility_functions.R"))
source(file.path(SCRIPT_DIR, "3_Matrix_Creation_utils.R"))

# ===============================================
# MAIN
# ===============================================

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("MATRIX CREATION - M3 STAR + Salmon\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")
cat("Master Reference:", MASTER_REFERENCE, "\n")
cat("Samples:         ", length(SAMPLE_IDS), "\n\n")

if (length(SAMPLE_IDS) == 0) stop("No samples loaded. Check SRR_COMBINED_LIST_STR and SRR_csv files.")

# Salmon quant outputs live under the MASTER_REFERENCE subdirectory created by STAR alignment
base_dir  <- Sys.getenv("BASE_DIR", "")
quant_dir <- if (nzchar(base_dir)) {
  file.path(base_dir, "2_ALIGNMENT_RESULTs", "M3_STAR_Align",
            MASTER_REFERENCE, "6_salmon", "quant")
} else {
  "STAR_alignment_WD"  # fallback for standalone execution
}
output_dir  <- "count_matrices_from_STAR"
count_label <- "NumReads"

cat("Quantification directory:", quant_dir, "\n")
cat("Output directory:        ", output_dir, "\n\n")

results <- list()

# ----- Gene-level (requires tx2gene mapping) --------------------------------
if (GENERATE_GENE_LEVEL) {
  tx2gene_files <- list.files(file.path(output_dir, MASTER_REFERENCE),
                               pattern = "^tx2gene.*\\.tsv$", full.names = TRUE)
  if (length(tx2gene_files) > 0) {
    tx2gene <- read.delim(tx2gene_files[1], header = FALSE,
                          col.names = c("TXNAME", "GENEID"),
                          stringsAsFactors = FALSE)
    quant_files <- file.path(quant_dir, SAMPLE_IDS, "quant.sf")
    names(quant_files) <- SAMPLE_IDS

    txi_gene <- tryCatch(
      tximport(quant_files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE),
      error = function(e) { cat("  Gene-level import error:", e$message, "\n"); NULL })

    if (!is.null(txi_gene)) {
      results$gene_level     <- txi_gene$counts
      results$gene_level_tpm <- txi_gene$abundance
    }
  } else {
    cat("  Warning: tx2gene file not found in", file.path(output_dir, MASTER_REFERENCE),
        "— skipping gene-level import\n")
  }
}

# ----- Isoform-level (transcript-level, no tx2gene needed) -------------------
if (GENERATE_ISOFORM_LEVEL) {
  quant_files <- file.path(quant_dir, SAMPLE_IDS, "quant.sf")
  names(quant_files) <- SAMPLE_IDS

  txi_iso <- tryCatch(
    tximport(quant_files, type = "salmon", txIn = TRUE, txOut = TRUE),
    error = function(e) { cat("  Isoform-level import error:", e$message, "\n"); NULL })

  if (!is.null(txi_iso)) {
    results$isoform_level     <- txi_iso$counts
    results$isoform_level_tpm <- txi_iso$abundance
  }
}

run_matrix_saving(results, output_dir, MASTER_REFERENCE, count_label, GENE_GROUPS_DIR)
