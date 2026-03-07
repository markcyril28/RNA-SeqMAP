#!/usr/bin/env Rscript

# ===============================================
# MATRIX CREATION MODULE
# ===============================================
# Creates count matrices from quantification outputs (tximport)
# Supports: Salmon, RSEM, StringTie outputs

# ===============================================
# CONFIGURATION (can be overridden by environment)
# ===============================================

GENERATE_GENE_LEVEL <- TRUE
GENERATE_ISOFORM_LEVEL <- TRUE

suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
  library(dplyr)
  library(tibble)
})

# Get script directory for sourcing
SCRIPT_DIR <- Sys.getenv("ANALYSIS_MODULES_DIR", ".")
source(file.path(SCRIPT_DIR, "0_shared_config.R"))
source(file.path(SCRIPT_DIR, "1_utility_functions.R"))

# ===============================================
# TXIMPORT FUNCTIONS BY METHOD
# ===============================================

# RSEM import (Method 5)
import_rsem <- function(quant_dir, sample_ids, level = "gene") {
  file_type <- if (level == "gene") ".genes.results" else ".isoforms.results"
  files <- file.path(quant_dir, sample_ids, paste0(sample_ids, file_type))
  names(files) <- sample_ids
  
  if (!all(file.exists(files))) {
    missing <- files[!file.exists(files)]
    cat("Missing files:", paste(missing, collapse = ", "), "\n")
    return(NULL)
  }
  
  txi <- tximport(files, type = "rsem", txIn = (level != "gene"), txOut = (level != "gene"))
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

save_count_matrices <- function(counts, output_dir, prefix, master_ref, level) {
  ensure_output_dir(output_dir)
  
  # Save with gene IDs
  output_file_id <- file.path(output_dir, 
    paste0(prefix, "_expected_count_geneID_from_", master_ref, "_", level, ".tsv"))
  counts_df <- as.data.frame(counts)
  counts_df <- tibble::rownames_to_column(counts_df, "GeneID")
  write.table(counts_df, output_file_id, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("Saved:", output_file_id, "\n")
  
  # Save with organ labels
  output_file_name <- file.path(output_dir,
    paste0(prefix, "_expected_count_geneName_from_", master_ref, "_", level, ".tsv"))
  counts_organ <- convert_to_organ_labels(counts)
  counts_organ_df <- as.data.frame(counts_organ)
  counts_organ_df <- tibble::rownames_to_column(counts_organ_df, "GeneID")
  write.table(counts_organ_df, output_file_name, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("Saved:", output_file_name, "\n")
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
  
  # Match genes
  matched_genes <- rownames(counts_matrix)[rownames(counts_matrix) %in% gene_list]
  
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
  if (grepl("RSEM|M5", method, ignore.case = TRUE)) {
    if (GENERATE_GENE_LEVEL) {
      txi_gene <- import_rsem(quant_dir, sample_ids, "gene")
      if (!is.null(txi_gene)) results$gene_level <- txi_gene$counts
    }
    if (GENERATE_ISOFORM_LEVEL) {
      txi_iso <- import_rsem(quant_dir, sample_ids, "isoform")
      if (!is.null(txi_iso)) results$isoform_level <- txi_iso$counts
    }
  } else if (grepl("Salmon|M4", method, ignore.case = TRUE)) {
    txi <- import_salmon(quant_dir, sample_ids)
    if (!is.null(txi)) results$gene_level <- txi$counts
  } else if (grepl("HISAT|StringTie|M1|M2", method, ignore.case = TRUE)) {
    st_results <- import_stringtie(quant_dir, sample_ids)
    if (!is.null(st_results$gene)) results$gene_level <- st_results$gene
    if (!is.null(st_results$transcript)) results$isoform_level <- st_results$transcript
  }
  
  # Save matrices for each level
  for (level in names(results)) {
    level_output <- file.path(output_dir, master_ref, level)
    ensure_output_dir(level_output)
    
    save_count_matrices(results[[level]], level_output, master_ref, master_ref, level)
    
    # Filter by configured gene groups (from config, not all files in directory)
    config <- load_runtime_config()
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
          group_output <- file.path(level_output, gene_group)
          ensure_output_dir(group_output)
          save_count_matrices(filtered, group_output, gene_group, master_ref, level)
        }
      }
    }
  }
  
  cat("\nMatrix creation complete.\n")
  return(results)
}
