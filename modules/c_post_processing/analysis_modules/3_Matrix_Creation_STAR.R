#!/usr/bin/env Rscript

# ===============================================
# MATRIX CREATION - M3 STAR + SALMON
# ===============================================
# Creates count matrices from STAR-aligned Salmon quantification outputs.
# Expects environment variables set by pipeline_utils.sh:
#   CURRENT_METHOD, MASTER_REFERENCE, BASE_DIR,
#   SRR_COMBINED_LIST_STR, GENE_GROUPS_STR, GENE_GROUPS_DIR,
#   STAR_GENERATE_GENE_LEVEL, STAR_GENERATE_ISOFORM_LEVEL

GENERATE_GENE_LEVEL    <- as.logical(Sys.getenv("STAR_GENERATE_GENE_LEVEL",    "TRUE"))
GENERATE_ISOFORM_LEVEL <- as.logical(Sys.getenv("STAR_GENERATE_ISOFORM_LEVEL", "TRUE"))

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
  file.path("..", "..", "2_ALIGNMENT_RESULTs", "M3_STAR_Align",
            MASTER_REFERENCE, "6_salmon", "quant")
}
output_dir  <- "count_matrices_from_STAR"
count_label <- "NumReads"

cat("Quantification directory:", quant_dir, "\n")
cat("Output directory:        ", output_dir, "\n\n")

results <- list()

# ----- Gene-level (requires tx2gene mapping) --------------------------------
if (GENERATE_GENE_LEVEL) {
  tx2gene_dir <- file.path(output_dir, MASTER_REFERENCE)
  tx2gene_files <- list.files(tx2gene_dir, pattern = "^tx2gene.*\\.tsv$", full.names = TRUE)

  # Generate tx2gene from GTF if missing (alignment step may not have been run)
  if (length(tx2gene_files) == 0) {
    gtf_ref_dir <- if (nzchar(base_dir)) {
      file.path(base_dir, "inputs", "gtf", "reference")
    } else {
      file.path("..", "..", "inputs", "gtf", "reference")
    }
    gtf_candidates <- c(
      Sys.getenv("STAR_GTF_FILE", unset = ""),
      # Prefer stringtie GTF variant (used by STAR alignment) to ensure
      # transcript_id attributes match the Salmon index built during alignment
      file.path(gtf_ref_dir, paste0(MASTER_REFERENCE, "_function_IPR_final_stringtie.gtf")),
      file.path(gtf_ref_dir, paste0(MASTER_REFERENCE, "_function_IPR_final.gtf")),
      file.path(gtf_ref_dir, paste0(MASTER_REFERENCE, ".gtf"))
    )
    gtf_file <- Filter(file.exists, gtf_candidates)[1]
    if (!is.na(gtf_file) && nzchar(gtf_file)) {
      cat("  Generating tx2gene from GTF:", gtf_file, "\n")
      dir.create(tx2gene_dir, recursive = TRUE, showWarnings = FALSE)
      tx2gene_out <- file.path(tx2gene_dir, paste0("tx2gene_", MASTER_REFERENCE, ".tsv"))
      awk_cmd <- sprintf(
        "awk '$3==\"transcript\" { tid=\"\"; gid=\"\"; for(i=9;i<=NF;i++) { if($i==\"transcript_id\") { gsub(/[\";]/,\"\",$(i+1)); tid=$(i+1) } if($i==\"gene_id\") { gsub(/[\";]/,\"\",$(i+1)); gid=$(i+1) } } if(tid!=\"\" && gid!=\"\") print tid \"\\t\" gid }' '%s' | sort -u > '%s'",
        gtf_file, tx2gene_out)
      system(awk_cmd)
      if (file.exists(tx2gene_out) && file.size(tx2gene_out) > 0) {
        tx2gene_files <- tx2gene_out
        cat("  Created tx2gene mapping:", tx2gene_out, "\n")
      } else {
        cat("  Error: Failed to generate tx2gene from GTF\n")
      }
    }
  }

  if (length(tx2gene_files) > 0) {
    tx2gene <- read.delim(tx2gene_files[1], header = FALSE,
                          col.names = c("TXNAME", "GENEID"),
                          stringsAsFactors = FALSE)
    if (nrow(tx2gene) == 0 || ncol(tx2gene) < 2) {
      cat("  Error: tx2gene file is empty or malformed:", tx2gene_files[1], "\n")
    } else {
      quant_files <- file.path(quant_dir, SAMPLE_IDS, "quant.sf")
      names(quant_files) <- SAMPLE_IDS
      # Tissue-specific fallback: scan subdirectories for quant.sf files
      if (sum(file.exists(quant_files)) == 0 && dir.exists(quant_dir)) {
        cat("  No quant.sf in flat layout; scanning tissue subdirectories...\n")
        tissue_dirs <- list.dirs(quant_dir, recursive = FALSE, full.names = TRUE)
        for (sid in SAMPLE_IDS) {
          for (td in tissue_dirs) {
            candidate <- file.path(td, sid, "quant.sf")
            if (file.exists(candidate)) {
              quant_files[sid] <- candidate
              break
            }
          }
        }
      }
      missing_qf <- quant_files[!file.exists(quant_files)]
      if (length(missing_qf) > 0) {
        cat("  Warning: Missing quant.sf for", length(missing_qf), "samples:",
            paste(names(missing_qf), collapse = ", "), "\n")
        quant_files <- quant_files[file.exists(quant_files)]
      }
      if (length(quant_files) < 2) {
        cat("  Error: Need >= 2 quant.sf files for gene-level import\n")
      } else {
        # Validate tx2gene transcript IDs match quant.sf transcript IDs
        sample_qf <- read.delim(quant_files[1], header = TRUE, nrows = 100,
                                 stringsAsFactors = FALSE)
        qf_ids <- sample_qf$Name
        tx_ids <- tx2gene$TXNAME
        overlap <- length(intersect(qf_ids, tx_ids))
        match_rate <- overlap / length(qf_ids)
        if (match_rate < 0.5) {
          cat("  ERROR: tx2gene transcript IDs poorly match quant.sf IDs!\n")
          cat("    Match rate:", round(match_rate * 100), "% (", overlap, "/", length(qf_ids), "sampled)\n")
          cat("    tx2gene IDs (first 3):", paste(head(tx2gene$TXNAME, 3), collapse = ", "), "\n")
          cat("    quant.sf IDs (first 3):", paste(head(sample_qf$Name, 3), collapse = ", "), "\n")
          cat("    This usually means tx2gene was generated from the wrong GTF.\n")
          cat("    Re-run STAR+Salmon alignment to regenerate tx2gene.\n")
          cat("    Skipping gene-level import.\n\n")
        } else {
          txi_gene <- tryCatch(
            tximport(quant_files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = FALSE),
            error = function(e) { cat("  Gene-level import error:", e$message, "\n"); NULL })

          if (!is.null(txi_gene)) {
            results$gene_level     <- txi_gene$counts
            results$gene_level_tpm <- txi_gene$abundance
            # Save full tximport object for DESeq2 (preserves transcript-length offsets)
            txi_rds_dir <- file.path(output_dir, MASTER_REFERENCE, "gene_level")
            ensure_output_dir(txi_rds_dir)
            saveRDS(txi_gene, file.path(txi_rds_dir, "tximport_gene_level.rds"))
            cat("Saved tximport RDS for DESeq2: tximport_gene_level.rds\n")
          }
        }
      }
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
  # Tissue-specific fallback: scan subdirectories for quant.sf files
  if (sum(file.exists(quant_files)) == 0 && dir.exists(quant_dir)) {
    cat("  No quant.sf in flat layout; scanning tissue subdirectories...\n")
    tissue_dirs <- list.dirs(quant_dir, recursive = FALSE, full.names = TRUE)
    for (sid in SAMPLE_IDS) {
      for (td in tissue_dirs) {
        candidate <- file.path(td, sid, "quant.sf")
        if (file.exists(candidate)) {
          quant_files[sid] <- candidate
          break
        }
      }
    }
  }
  quant_files <- quant_files[file.exists(quant_files)]
  if (length(quant_files) < 2) {
    cat("  Error: Need >= 2 quant.sf files for isoform-level import\n")
  } else {
    txi_iso <- tryCatch(
      tximport(quant_files, type = "salmon", txIn = TRUE, txOut = TRUE,
               ignoreTxVersion = FALSE, ignoreAfterBar = FALSE),
      error = function(e) { cat("  Isoform-level import error:", e$message, "\n"); NULL })

    if (!is.null(txi_iso)) {
      results$isoform_level     <- txi_iso$counts
      results$isoform_level_tpm <- txi_iso$abundance
      # Save full tximport object for DESeq2 isoform-level analysis
      txi_iso_rds_dir <- file.path(output_dir, MASTER_REFERENCE, "isoform_level")
      ensure_output_dir(txi_iso_rds_dir)
      saveRDS(txi_iso, file.path(txi_iso_rds_dir, "tximport_isoform_level.rds"))
      cat("Saved tximport RDS for isoform-level DESeq2: tximport_isoform_level.rds\n")
    }
  }
}

run_matrix_saving(results, output_dir, MASTER_REFERENCE, count_label, GENE_GROUPS_DIR)
