#!/usr/bin/env Rscript

# ===============================================
# MATRIX CREATION - M4 SALMON
# ===============================================
# Creates count matrices from Salmon quantification outputs using tximport.
# Expects environment variables set by pipeline_utils.sh:
#   CURRENT_METHOD, MASTER_REFERENCE, BASE_DIR,
#   SRR_COMBINED_LIST_STR, GENE_GROUPS_STR, GENE_GROUPS_DIR,
#   SALMON_GENERATE_GENE_LEVEL, SALMON_GENERATE_ISOFORM_LEVEL

GENERATE_GENE_LEVEL    <- isTRUE(as.logical(Sys.getenv("SALMON_GENERATE_GENE_LEVEL",    "TRUE")))
GENERATE_ISOFORM_LEVEL <- isTRUE(as.logical(Sys.getenv("SALMON_GENERATE_ISOFORM_LEVEL", "TRUE")))

suppressPackageStartupMessages(library(tximport))

SCRIPT_DIR <- Sys.getenv("ANALYSIS_MODULES_DIR", ".")
source(file.path(SCRIPT_DIR, "0_shared_config.R"))
source(file.path(SCRIPT_DIR, "1_utility_functions.R"))
source(file.path(SCRIPT_DIR, "3_Matrix_Creation_utils.R"))

# Cache data.table availability for fread fast paths below
.use_dt_salmon <- requireNamespace("data.table", quietly = TRUE)

# ===============================================
# SALMON IMPORT
# ===============================================

import_salmon <- function(quant_dir, sample_ids, tx2gene = NULL) {
  files <- file.path(quant_dir, sample_ids, "quant.sf")
  names(files) <- sample_ids

  files_exist <- file.exists(files)
  if (!any(files_exist)) {
    cat("ERROR: No Salmon quant.sf files found in", quant_dir, "\n")
    return(NULL)
  }
  if (!all(files_exist)) {
    cat("WARNING: Missing quant.sf for:", paste(names(files)[!files_exist], collapse = ", "), "\n")
    files <- files[files_exist]
  }
  if (length(files) < 2) {
    cat("ERROR: Need >= 2 quant.sf files for import (found", length(files), ")\n")
    return(NULL)
  }

  if (!is.null(tx2gene)) {
    tximport(files, type = "salmon", txIn = TRUE, txOut = FALSE,
             tx2gene = tx2gene, ignoreTxVersion = FALSE, ignoreAfterBar = FALSE)
  } else {
    tximport(files, type = "salmon", txIn = TRUE, txOut = TRUE,
             ignoreTxVersion = FALSE, ignoreAfterBar = FALSE)
  }
}

# ===============================================
# TX2GENE LOOKUP (mirrors tximport_salmon_to_matrices.R search strategy)
# ===============================================

.base_dir         <- Sys.getenv("BASE_DIR", "")
.input_fastas_dir <- Sys.getenv("INPUT_FASTAS_DIR", unset = "")
if (!nzchar(.input_fastas_dir) && nzchar(.base_dir))
  .input_fastas_dir <- file.path(.base_dir, "inputs")

# Lazy tx2gene search: check known paths first, only recurse as fallback
.candidates <- c(
  Sys.getenv("GENE_TRANS_MAP_FILE", unset = ""),
  file.path(.input_fastas_dir, "mapping", paste0(MASTER_REFERENCE, ".fa.gene_trans_map")),
  file.path(.input_fastas_dir, "mapping", paste0(MASTER_REFERENCE, ".fasta.gene_trans_map")),
  file.path(.input_fastas_dir, "fasta",   paste0(MASTER_REFERENCE, ".fa.gene_trans_map")),
  file.path(.input_fastas_dir, "fasta",   "reference_genome", paste0(MASTER_REFERENCE, ".fa.gene_trans_map")),
  file.path(.input_fastas_dir, "fasta",   "reference_genome", paste0(MASTER_REFERENCE, ".fasta.gene_trans_map"))
)

.tx2gene_file <- NULL
for (.cand in .candidates) {
  if (nzchar(.cand) && file.exists(.cand)) { .tx2gene_file <- .cand; break }
}
# Only recurse directory tree if direct paths failed
if (is.null(.tx2gene_file) && nzchar(.input_fastas_dir)) {
  .all_maps <- list.files(.input_fastas_dir, pattern = "\\.gene_trans_map$",
                          recursive = TRUE, full.names = TRUE)
  # Match on basename to avoid substring false positives (e.g., "V4" matching "V4.1")
  .hits <- .all_maps[grepl(paste0("(^|[/\\\\])", MASTER_REFERENCE, "\\."), .all_maps)]
  if (length(.hits) > 0) .tx2gene_file <- .hits[1]
}

.tx2gene <- if (!is.null(.tx2gene_file)) {
  tryCatch({
    .t2g <- if (.use_dt_salmon) {
      data.table::fread(.tx2gene_file, header = FALSE, sep = "\t",
                        strip.white = TRUE, showProgress = FALSE, data.table = FALSE)
    } else {
      read.table(.tx2gene_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                 strip.white = TRUE)
    }
    .t2g <- .t2g[, 1:2, drop = FALSE]
    colnames(.t2g) <- c("GENEID", "TXNAME")
    .t2g[, c("TXNAME", "GENEID")]
  }, error = function(e) {
    cat("  Warning: Failed to read tx2gene file:", e$message, "\n")
    cat("  Gene-level import will be skipped\n")
    NULL
  })
} else {
  cat("  Warning: tx2gene mapping not found for M4 — gene-level import will be skipped\n")
  NULL
}

rm(list = intersect(c(".base_dir", ".input_fastas_dir", ".candidates",
                       ".tx2gene_file", ".cand", ".all_maps", ".hits"), ls(all.names = TRUE)))

# ===============================================
# MAIN
# ===============================================

cat("\n", strrep("=", 60), "\n")
cat("MATRIX CREATION - M4 Salmon\n")
cat(strrep("=", 60), "\n\n")
cat("Master Reference:", MASTER_REFERENCE, "\n")
cat("Samples:         ", length(SAMPLE_IDS), "\n\n")

if (length(SAMPLE_IDS) == 0) stop("No samples loaded. Check SRR_COMBINED_LIST_STR and SRR_csv files.")

base_dir  <- Sys.getenv("BASE_DIR", "")
salmon_quant_root_env <- Sys.getenv("SALMON_QUANT_ROOT", "")
quant_dir <- if (nzchar(salmon_quant_root_env)) {
  salmon_quant_root_env  # already includes fasta_tag
} else if (nzchar(base_dir)) {
  file.path(base_dir, "2_ALIGNMENT_RESULTs", "M4_Salmon_Saf", "Salmon_Quant", MASTER_REFERENCE)
} else {
  "Salmon_Quant"  # fallback for standalone execution
}
output_dir  <- if (nzchar(base_dir)) {
  file.path(base_dir, "3_POST_PROC", "M4_Salmon_Saf", "count_matrices_from_Salmon_Quant")
} else {
  "count_matrices_from_Salmon_Quant"  # relative fallback for pushd context
}
count_label <- "NumReads"

cat("Quantification directory:", quant_dir, "\n")
cat("Output directory:        ", output_dir, "\n\n")

if (!dir.exists(quant_dir)) {
  stop("Salmon quantification directory not found: ", quant_dir,
       "\n  Run the M4 Salmon SAF alignment stage first, or check SALMON_QUANT_ROOT / MASTER_REFERENCE.")
}

results <- list()

if (GENERATE_GENE_LEVEL) {
  if (is.null(.tx2gene)) {
    cat("WARNING: Skipping gene-level import — tx2gene mapping not found.\n")
    cat("  Gene-level matrices require a .gene_trans_map file to aggregate transcripts to genes.\n")
  } else {
    # Validate tx2gene IDs against a sample quant.sf before full import
    .sample_qsf <- file.path(quant_dir, SAMPLE_IDS[1], "quant.sf")
    if (file.exists(.sample_qsf)) {
      .qsf_ids <- tryCatch(read.table(.sample_qsf, header = TRUE, sep = "\t",
                            nrows = 200, stringsAsFactors = FALSE)$Name, error = function(e) NULL)
      if (!is.null(.qsf_ids)) {
        .match_rate <- mean(.qsf_ids %in% .tx2gene$TXNAME)
        if (.match_rate < 0.5) {
          cat("WARNING: tx2gene match rate against quant.sf is", round(.match_rate * 100, 1),
              "% — transcript IDs may not match. Check gene_trans_map file.\n")
        }
      }
    }
    rm(list = intersect(c(".sample_qsf", ".qsf_ids", ".match_rate"), ls(all.names = TRUE)))
    txi <- tryCatch(
      import_salmon(quant_dir, SAMPLE_IDS, tx2gene = .tx2gene),
      error = function(e) { cat("  Gene-level import error:", e$message, "\n"); NULL })
    if (!is.null(txi) && is.list(txi) && "counts" %in% names(txi)) {
      if (nrow(txi$counts) == 0 || ncol(txi$counts) == 0) {
        cat("WARNING: Gene-level import produced empty matrix — skipping\n")
      } else {
      results$gene_level     <- txi$counts
      results$gene_level_tpm <- txi$abundance
      # Save full tximport object for DESeq2 (preserves transcript-length offsets)
      txi_rds_dir <- file.path(output_dir, MASTER_REFERENCE, "gene_level")
      ensure_output_dir(txi_rds_dir)
      tryCatch(
        saveRDS(txi, file.path(txi_rds_dir, "tximport_gene_level.rds")),
        error = function(e) cat("  Warning: Failed to save tximport RDS:", e$message, "\n")
      )
      cat("Saved tximport RDS for DESeq2: tximport_gene_level.rds\n")
      }
    }
  }
}

if (GENERATE_ISOFORM_LEVEL) {
  txi <- tryCatch(
    import_salmon(quant_dir, SAMPLE_IDS),   # txOut = TRUE (no tx2gene)
    error = function(e) { cat("  Isoform-level import error:", e$message, "\n"); NULL })
  if (!is.null(txi) && is.list(txi) && "counts" %in% names(txi)) {
    if (nrow(txi$counts) == 0 || ncol(txi$counts) == 0) {
      cat("WARNING: Isoform-level import produced empty matrix — skipping\n")
    } else {
      results$isoform_level     <- txi$counts
      results$isoform_level_tpm <- txi$abundance
      # Save tximport RDS for isoform-level DESeq2 (preserves transcript-length offsets)
      txi_iso_rds_dir <- file.path(output_dir, MASTER_REFERENCE, "isoform_level")
      ensure_output_dir(txi_iso_rds_dir)
      tryCatch(
        saveRDS(txi, file.path(txi_iso_rds_dir, "tximport_isoform_level.rds")),
        error = function(e) cat("  Warning: Failed to save isoform tximport RDS:", e$message, "\n")
      )
      cat("Saved tximport RDS for isoform-level DESeq2: tximport_isoform_level.rds\n")
    }
  }
}

if (exists(".tx2gene")) rm(.tx2gene)

run_matrix_saving(results, output_dir, MASTER_REFERENCE, count_label, GENE_GROUPS_DIR)
