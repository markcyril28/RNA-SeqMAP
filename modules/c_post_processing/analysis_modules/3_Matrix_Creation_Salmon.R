#!/usr/bin/env Rscript

# ===============================================
# MATRIX CREATION - M4 SALMON
# ===============================================
# Creates count matrices from Salmon quantification outputs using tximport.
# Expects environment variables set by pipeline_utils.sh:
#   CURRENT_METHOD, MASTER_REFERENCE, BASE_DIR,
#   SRR_COMBINED_LIST_STR, GENE_GROUPS_STR, GENE_GROUPS_DIR,
#   SALMON_GENERATE_GENE_LEVEL, SALMON_GENERATE_ISOFORM_LEVEL

GENERATE_GENE_LEVEL    <- as.logical(Sys.getenv("SALMON_GENERATE_GENE_LEVEL",    "TRUE"))
GENERATE_ISOFORM_LEVEL <- as.logical(Sys.getenv("SALMON_GENERATE_ISOFORM_LEVEL", "TRUE"))

suppressPackageStartupMessages(library(tximport))

SCRIPT_DIR <- Sys.getenv("ANALYSIS_MODULES_DIR", ".")
source(file.path(SCRIPT_DIR, "0_shared_config.R"))
source(file.path(SCRIPT_DIR, "1_utility_functions.R"))
source(file.path(SCRIPT_DIR, "3_Matrix_Creation_utils.R"))

# ===============================================
# SALMON IMPORT
# ===============================================

import_salmon <- function(quant_dir, sample_ids, tx2gene = NULL) {
  files <- file.path(quant_dir, sample_ids, "quant.sf")
  names(files) <- sample_ids

  missing <- files[!file.exists(files)]
  if (length(missing) == length(files)) {
    cat("ERROR: No Salmon quant.sf files found in", quant_dir, "\n")
    return(NULL)
  }
  if (length(missing) > 0) {
    cat("WARNING: Missing quant.sf for:", paste(names(missing), collapse = ", "), "\n")
    files <- files[file.exists(files)]
  }

  if (!is.null(tx2gene)) {
    tximport(files, type = "salmon", txIn = TRUE, txOut = FALSE,
             tx2gene = tx2gene, ignoreTxVersion = TRUE, ignoreAfterBar = FALSE)
  } else {
    tximport(files, type = "salmon", txIn = TRUE, txOut = TRUE,
             ignoreTxVersion = TRUE, ignoreAfterBar = FALSE)
  }
}

# ===============================================
# TX2GENE LOOKUP (mirrors tximport_salmon_to_matrices.R search strategy)
# ===============================================

.base_dir         <- Sys.getenv("BASE_DIR", "")
.input_fastas_dir <- Sys.getenv("INPUT_FASTAS_DIR", unset = "")
if (!nzchar(.input_fastas_dir) && nzchar(.base_dir))
  .input_fastas_dir <- file.path(.base_dir, "inputs")

.candidates <- c(
  Sys.getenv("GENE_TRANS_MAP_FILE", unset = ""),
  file.path(.input_fastas_dir, "mapping", paste0(MASTER_REFERENCE, ".fa.gene_trans_map")),
  file.path(.input_fastas_dir, "mapping", paste0(MASTER_REFERENCE, ".fasta.gene_trans_map")),
  file.path(.input_fastas_dir, "fasta",   paste0(MASTER_REFERENCE, ".fa.gene_trans_map"))
)
if (nzchar(.input_fastas_dir)) {
  .all_maps   <- list.files(.input_fastas_dir, pattern = "\\.gene_trans_map$",
                            recursive = TRUE, full.names = TRUE)
  .candidates <- c(.candidates, .all_maps[grepl(MASTER_REFERENCE, .all_maps, fixed = TRUE)])
}

.tx2gene_file <- NULL
for (.cand in .candidates) {
  if (nzchar(.cand) && file.exists(.cand)) { .tx2gene_file <- .cand; break }
}

.tx2gene <- if (!is.null(.tx2gene_file)) {
  .t2g <- read.table(.tx2gene_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                     colClasses = c("character", "character"))
  colnames(.t2g) <- c("GENEID", "TXNAME")
  .t2g$GENEID <- trimws(.t2g$GENEID); .t2g$TXNAME <- trimws(.t2g$TXNAME)
  .t2g[, c("TXNAME", "GENEID")]
} else {
  cat("  Warning: tx2gene mapping not found for M4 — gene-level import will be transcript-level\n")
  NULL
}

rm(list = intersect(c(".base_dir", ".input_fastas_dir", ".candidates",
                       ".tx2gene_file", ".cand", ".all_maps"), ls()))

# ===============================================
# MAIN
# ===============================================

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("MATRIX CREATION - M4 Salmon\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")
cat("Master Reference:", MASTER_REFERENCE, "\n")
cat("Samples:         ", length(SAMPLE_IDS), "\n\n")

if (length(SAMPLE_IDS) == 0) stop("No samples loaded. Check SRR_COMBINED_LIST_STR and SRR_csv files.")

base_dir  <- Sys.getenv("BASE_DIR", "")
quant_dir <- if (nzchar(base_dir)) {
  file.path(base_dir, "2_ALIGNMENT_RESULTs", "M4_Salmon_Saf", "Salmon_Quant", MASTER_REFERENCE)
} else {
  "Salmon_Quant"  # fallback for standalone execution
}
output_dir  <- "count_matrices_from_Salmon_Quant"
count_label <- "NumReads"

cat("Quantification directory:", quant_dir, "\n")
cat("Output directory:        ", output_dir, "\n\n")

results <- list()

if (GENERATE_GENE_LEVEL) {
  txi <- import_salmon(quant_dir, SAMPLE_IDS, tx2gene = .tx2gene)
  if (!is.null(txi)) {
    results$gene_level     <- txi$counts
    results$gene_level_tpm <- txi$abundance
  }
}

if (GENERATE_ISOFORM_LEVEL) {
  txi <- import_salmon(quant_dir, SAMPLE_IDS)   # txOut = TRUE (no tx2gene)
  if (!is.null(txi)) {
    results$isoform_level     <- txi$counts
    results$isoform_level_tpm <- txi$abundance
  }
}

rm(.tx2gene)

run_matrix_saving(results, output_dir, MASTER_REFERENCE, count_label, GENE_GROUPS_DIR)
