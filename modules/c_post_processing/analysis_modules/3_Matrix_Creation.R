#!/usr/bin/env Rscript

# ===============================================
# MATRIX CREATION MODULE
# ===============================================
# Creates count matrices from quantification outputs (tximport)
# Supports: Salmon, RSEM, StringTie outputs

# ===============================================
# CONFIGURATION (can be overridden by environment)
# ===============================================

GENERATE_GENE_LEVEL    <- as.logical(Sys.getenv("RSEM_GENERATE_GENE_LEVEL",    "TRUE"))
GENERATE_ISOFORM_LEVEL <- as.logical(Sys.getenv("RSEM_GENERATE_ISOFORM_LEVEL", "TRUE"))

suppressPackageStartupMessages({
  library(tximport)
})

# Get script directory for sourcing
SCRIPT_DIR <- Sys.getenv("ANALYSIS_MODULES_DIR", ".")
source(file.path(SCRIPT_DIR, "0_shared_config.R"))
source(file.path(SCRIPT_DIR, "1_utility_functions.R"))

# ===============================================
# TXIMPORT FUNCTIONS BY METHOD
# ===============================================

# RSEM import (Method 5)
# Returns a list with $counts (expected_count) and $abundance (TPM),
# both filtered to remove genes/transcripts with zero effective length.
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
    missing <- names(files)[!files_exist]
    cat("WARNING: Missing", level, "results for:", paste(missing, collapse = ", "), "\n")
    cat("Continuing with", sum(files_exist), "available samples\n")
    files <- files[files_exist]
  }

  txi <- tximport(files, type = "rsem", txIn = (level != "gene"), txOut = (level != "gene"))

  # Filter entries with zero effective length (same as tximport_rsem_to_matrices.R)
  if (!is.null(txi$length)) {
    zero_length_mask <- rowSums(txi$length == 0) > 0
    n_zero <- sum(zero_length_mask)
    if (n_zero > 0) {
      cat("  Filtering", n_zero, level, "entries with zero effective length\n")
      txi$counts    <- txi$counts[!zero_length_mask, , drop = FALSE]
      txi$abundance <- txi$abundance[!zero_length_mask, , drop = FALSE]
      txi$length    <- txi$length[!zero_length_mask, , drop = FALSE]
    }
  }

  return(txi)
}

# Salmon import (Method 4)
import_salmon <- function(quant_dir, sample_ids, tx2gene = NULL) {
  files <- file.path(quant_dir, sample_ids, "quant.sf")
  names(files) <- sample_ids
  
  if (!all(file.exists(files))) {
    missing <- files[!file.exists(files)]
    cat("Missing files:", paste(missing, collapse = ", "), "\n")
    return(NULL)
  }
  
  if (!is.null(tx2gene)) {
    txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
  } else {
    txi <- tximport(files, type = "salmon", txOut = TRUE)
  }
  return(txi)
}

# StringTie import (Methods 1, 2)
import_stringtie <- function(ballgown_dir, sample_ids) {
  # StringTie with prepDE.py output
  gene_count_file <- file.path(ballgown_dir, "gene_count_matrix.csv")
  transcript_count_file <- file.path(ballgown_dir, "transcript_count_matrix.csv")
  
  result <- list()
  if (file.exists(gene_count_file)) {
    gene_counts <- read.csv(gene_count_file, row.names = 1)
    result$gene <- as.matrix(gene_counts)
  }
  if (file.exists(transcript_count_file)) {
    tx_counts <- read.csv(transcript_count_file, row.names = 1)
    result$transcript <- as.matrix(tx_counts)
  }
  return(result)
}

# ===============================================
# MATRIX SAVING FUNCTIONS
# ===============================================

save_count_matrices <- function(counts, output_dir, prefix, master_ref, level,
                                tpm = NULL) {
  ensure_output_dir(output_dir)

  save_matrix <- function(matrix_data, count_type, gene_type) {
    output_file <- file.path(output_dir,
      paste0(prefix, "_", count_type, "_", gene_type, "_from_", master_ref, "_", level, ".tsv"))
    matrix_df <- as.data.frame(matrix_data, check.names = FALSE)
    matrix_df <- cbind(GeneID = rownames(matrix_data), matrix_df)
    rownames(matrix_df) <- NULL
    write.table(matrix_df, output_file, sep = "\t", quote = FALSE, row.names = FALSE)
    cat("Saved:", basename(output_file), "\n")
  }

  # Save expected_count matrices
  save_matrix(counts, "expected_count", "Gene_ID")
  counts_organ <- convert_to_organ_labels(counts)
  save_matrix(counts_organ, "expected_count", "Shortened_Name")

  # Save TPM matrices (if provided)
  if (!is.null(tpm)) {
    save_matrix(tpm, "tpm", "Gene_ID")
    tpm_organ <- convert_to_organ_labels(tpm)
    save_matrix(tpm_organ, "tpm", "Shortened_Name")
  }
}

# ===============================================
# GENE GROUP FILTERING
# ===============================================

filter_by_gene_group <- function(counts_matrix, gene_list_file) {
  if (!file.exists(gene_list_file)) {
    cat("Gene list file not found:", gene_list_file, "\n")
    return(NULL)
  }
  
  # Read gene list - handle both CSV and plain text formats
  file_ext <- tools::file_ext(gene_list_file)
  
  if (tolower(file_ext) == "csv") {
    # CSV format: read Gene_ID column or first column as gene IDs
    gene_df <- tryCatch({
      read.csv(gene_list_file, stringsAsFactors = FALSE, header = TRUE)
    }, error = function(e) NULL)
    
    if (is.null(gene_df) || nrow(gene_df) == 0) {
      cat("Failed to read CSV gene list\n")
      return(NULL)
    }
    # Use Gene_ID column if present, otherwise first column
    if ("Gene_ID" %in% colnames(gene_df)) {
      gene_list <- trimws(gene_df$Gene_ID)
    } else {
      gene_list <- trimws(gene_df[[1]])  # Fallback to first column
    }
  } else {
    # Plain text: one gene per line
    gene_list <- suppressWarnings(readLines(gene_list_file))
    gene_list <- gene_list[!grepl("^#|^Gene", gene_list, ignore.case = TRUE) & nzchar(gene_list)]
    gene_list <- trimws(gene_list)
  }
  
  gene_list <- gene_list[nzchar(gene_list)]  # Remove empty entries

  # Match genes: try exact match first, then prefix match for isoform IDs
  # e.g., gene list entry "SMEL4.1_01g005840" matches row "SMEL4.1_01g005840.1.01"
  all_row_ids <- rownames(counts_matrix)
  matched_genes <- character(0)
  for (gene in gene_list) {
    if (gene %in% all_row_ids) {
      matched_genes <- c(matched_genes, gene)
    } else {
      pattern <- paste0("^", gsub("\\.", "\\\\.", gene), "(\\..*)?$")
      hits <- all_row_ids[grepl(pattern, all_row_ids)]
      if (length(hits) > 0) matched_genes <- c(matched_genes, hits)
    }
  }
  matched_genes <- unique(matched_genes)

  if (length(matched_genes) == 0) {
    cat("No genes matched from list\n")
    return(NULL)
  }

  cat("Matched", length(matched_genes), "of", length(gene_list), "genes\n")
  return(counts_matrix[matched_genes, , drop = FALSE])
}

# ===============================================
# MAIN PROCESSING FUNCTION
# ===============================================

run_matrix_creation <- function(method, quant_dir, output_dir, master_ref, 
                                 sample_ids, gene_groups_dir = NULL) {
  cat("\n", paste(rep("=", 60), collapse = ""), "\n")
  cat("MATRIX CREATION -", method, "\n")
  cat(paste(rep("=", 60), collapse = ""), "\n\n")
  
  results <- list()
  
  # Import based on method
  if (grepl("STAR|M3", method, ignore.case = TRUE)) {
    # M3: STAR+Salmon - use Salmon quant.sf output with tx2gene mapping
    salmon_quant_dir <- file.path("..", "..", "2_ALIGNMENT_RESULTs", "M3_STAR_Align",
                                   "6_salmon", "quant")
    if (GENERATE_GENE_LEVEL) {
      tx2gene_file <- list.files(file.path("count_matrices_from_STAR", master_ref),
                                  pattern = "^tx2gene.*\\.tsv$", full.names = TRUE)
      if (length(tx2gene_file) > 0) {
        tx2gene <- read.delim(tx2gene_file[1], header = FALSE,
                              col.names = c("TXNAME", "GENEID"),
                              stringsAsFactors = FALSE)
        txi_gene <- tryCatch(
          tximport(file.path(salmon_quant_dir, sample_ids, "quant.sf") |>
                     setNames(sample_ids),
                   type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE),
          error = function(e) NULL)
        if (!is.null(txi_gene)) results$gene_level <- txi_gene$counts
      } else {
        cat("  Warning: tx2gene file not found for M3 gene-level import\n")
      }
    }
    if (GENERATE_ISOFORM_LEVEL) {
      txi_iso <- tryCatch(
        tximport(file.path(salmon_quant_dir, sample_ids, "quant.sf") |>
                   setNames(sample_ids),
                 type = "salmon", txIn = TRUE, txOut = TRUE),
        error = function(e) NULL)
      if (!is.null(txi_iso)) results$isoform_level <- txi_iso$counts
    }
  } else if (grepl("RSEM|M5", method, ignore.case = TRUE)) {
    if (GENERATE_GENE_LEVEL) {
      txi_gene <- import_rsem(quant_dir, sample_ids, "gene")
      if (!is.null(txi_gene)) {
        results$gene_level         <- txi_gene$counts
        results$gene_level_tpm     <- txi_gene$abundance
      }
    }
    if (GENERATE_ISOFORM_LEVEL) {
      txi_iso <- import_rsem(quant_dir, sample_ids, "isoform")
      if (!is.null(txi_iso)) {
        results$isoform_level      <- txi_iso$counts
        results$isoform_level_tpm  <- txi_iso$abundance
      }
    }
  } else if (grepl("Salmon|M4", method, ignore.case = TRUE)) {
    txi <- import_salmon(quant_dir, sample_ids)
    if (!is.null(txi)) results$gene_level <- txi$counts
  } else if (grepl("HISAT|StringTie|M1|M2", method, ignore.case = TRUE)) {
    st_results <- import_stringtie(quant_dir, sample_ids)
    if (!is.null(st_results$gene)) results$gene_level <- st_results$gene
    if (!is.null(st_results$transcript)) results$isoform_level <- st_results$transcript
  }
  
  # Load config once before loop (gene_groups list doesn't change per level)
  config <- load_runtime_config()

  # Save matrices for each level (skip companion TPM keys — handled below)
  level_names <- names(results)[!grepl("_tpm$", names(results))]
  for (level in level_names) {
    level_output <- file.path(output_dir, master_ref, level)
    ensure_output_dir(level_output)

    tpm_data <- results[[paste0(level, "_tpm")]]  # NULL for non-RSEM methods
    save_count_matrices(results[[level]], level_output, master_ref, master_ref, level,
                        tpm = tpm_data)

    # Filter by configured gene groups (from config, not all files in directory)
    if (!is.null(gene_groups_dir) && dir.exists(gene_groups_dir)) {
      for (gene_group in config$gene_groups) {
        # Look for matching CSV file
        gf <- file.path(gene_groups_dir, paste0(gene_group, ".csv"))
        if (!file.exists(gf)) {
          gf <- file.path(gene_groups_dir, paste0(gene_group, ".txt"))
        }
        if (!file.exists(gf)) {
          cat("Gene group file not found:", gene_group, "\n")
          next
        }

        filtered <- filter_by_gene_group(results[[level]], gf)
        if (!is.null(filtered)) {
          # Use GeneGroup_in_Dataset folder naming to match build_input_path()
          folder_name <- get_output_folder_name(gene_group)
          group_output <- file.path(level_output, folder_name)
          ensure_output_dir(group_output)
          # Also filter TPM to matching rows
          filtered_tpm <- if (!is.null(tpm_data)) tpm_data[rownames(filtered), , drop = FALSE] else NULL
          save_count_matrices(filtered, group_output, folder_name, master_ref, level,
                              tpm = filtered_tpm)
        }
      }
    }
  }
  
  cat("\nMatrix creation complete.\n")
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

# Build absolute quantification directory using BASE_DIR (same pattern as tximport_rsem_to_matrices.R)
base_dir <- Sys.getenv("BASE_DIR", "")
method_type <- get_method_type(CURRENT_METHOD)
quant_dir <- if (nzchar(base_dir)) {
  switch(method_type,
    "rsem"   = file.path(base_dir, "2_ALIGNMENT_RESULTs", "M5_RSEM_Bowtie2",
                         "RSEM_Quant_WD", MASTER_REFERENCE),
    "salmon" = file.path(base_dir, "2_ALIGNMENT_RESULTs", "M4_Salmon_Saf",
                         "Salmon_Quant", MASTER_REFERENCE),
    "star"   = file.path(base_dir, "2_ALIGNMENT_RESULTs", "M3_STAR_Align",
                         "6_salmon", "quant"),
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
