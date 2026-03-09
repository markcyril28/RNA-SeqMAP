#!/usr/bin/env Rscript

# ===============================================
# TXIMPORT: SALMON TO DESEQ2 MATRICES - GENERIC
# ===============================================
# Processes Salmon quantification using tximport
# for statistically sound gene expression analysis
# GENERIC VERSION - adaptable for any method

suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
  library(dplyr)
  library(tibble)
})

# ===============================================
# CONFIGURATION
# ===============================================

# Source shared config for SAMPLE_LABELS (DRY principle)
SCRIPT_DIR <- Sys.getenv("ANALYSIS_MODULES_DIR", ".")
source(file.path(SCRIPT_DIR, "0_shared_config.R"))
source(file.path(SCRIPT_DIR, "1_utility_functions.R"))

# ===============================================
# SALMON-SPECIFIC HELPER
# ===============================================

save_count_matrix <- function(counts_matrix, output_dir, base_name, master_ref, level_suffix, sample_labels,
                              tpm_matrix = NULL) {
  # Saves matrices with naming convention expected by build_input_path:
  # {gene_group}_{count_type}_{gene_type}_from_{master_ref}_{processing_level}.tsv
  # count_type: "tpm" or "NumReads"  
  # gene_type: "Gene_ID" or "Shortened_Name"
  # processing_level: "gene_level" or "isoform_level"
  
  # Extract processing level name (remove leading underscore)
  processing_level <- gsub("^_", "", level_suffix)
  
  # Helper to save a matrix with proper naming
  save_matrix <- function(matrix_data, count_type, gene_type) {
    output_file <- file.path(output_dir, 
      paste0(base_name, "_", count_type, "_", gene_type, "_from_", master_ref, "_", processing_level, ".tsv"))
    
    # Convert to data frame, handling potential duplicate column names
    matrix_df <- as.data.frame(matrix_data, check.names = FALSE)
    # Add GeneID column manually to avoid tibble issues with duplicate names
    matrix_df <- cbind(GeneID = rownames(matrix_data), matrix_df)
    rownames(matrix_df) <- NULL
    write.table(matrix_df, output_file, sep = "\t", quote = FALSE, row.names = FALSE)
    cat("Saved:", basename(output_file), "\n")
  }
  
  # Save NumReads matrices (Salmon's equivalent of raw counts)
  # With Gene_ID (SRR sample columns)
  save_matrix(counts_matrix, "NumReads", "Gene_ID")
  
  # With Shortened_Name (Organ labels for columns)
  counts_organ <- convert_to_organ_labels(counts_matrix)
  save_matrix(counts_organ, "NumReads", "Shortened_Name")
  
  # Save TPM matrices (if provided)
  if (!is.null(tpm_matrix)) {
    save_matrix(tpm_matrix, "tpm", "Gene_ID")
    tpm_organ <- convert_to_organ_labels(tpm_matrix)
    save_matrix(tpm_organ, "tpm", "Shortened_Name")
  }
}

# Debug: Show what SAMPLE_IDS we have after loading config
cat("[DEBUG] SRR_COMBINED_LIST_STR env var:", Sys.getenv("SRR_COMBINED_LIST_STR", unset = "<NOT SET>"), "\n")
cat("[DEBUG] Number of SAMPLE_IDS loaded:", length(SAMPLE_IDS), "\n")

QUANT_DIR <- "Salmon_Quant"
MASTER_REFERENCE <- Sys.getenv("MASTER_REFERENCE", "All_Smel_Genes")
MATRICES_OUTPUT_DIR <- "count_matrices_from_Salmon_Quant"
# Use shared GENE_GROUPS_DIR from 0_shared_config.R (already sourced)

# Toggle to generate both gene-level and isoform-level matrices
GENERATE_GENE_LEVEL <- TRUE      # Summarize transcripts to genes
GENERATE_ISOFORM_LEVEL <- TRUE   # Keep transcript-level data

# Use SAMPLE_IDS from shared config (0_shared_config.R)
# Override here if needed for method-specific samples

# Create output directories
output_dir <- file.path(MATRICES_OUTPUT_DIR, MASTER_REFERENCE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("TXIMPORT: SALMON QUANTIFICATION TO MATRICES\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

cat("Configuration:\n")
cat("  • Master Reference:", MASTER_REFERENCE, "\n")
cat("  • Generate gene-level matrices:", GENERATE_GENE_LEVEL, "\n")
cat("  • Generate isoform-level matrices:", GENERATE_ISOFORM_LEVEL, "\n\n")

# Define processing levels
processing_levels <- list()
if (GENERATE_GENE_LEVEL) {
  processing_levels[["gene_level"]] <- list(
    tx_in = TRUE,
    tx_out = FALSE,
    label = "Gene-Level",
    output_suffix = "_gene_level"
  )
}
if (GENERATE_ISOFORM_LEVEL) {
  processing_levels[["isoform_level"]] <- list(
    tx_in = TRUE,
    tx_out = TRUE,
    label = "Isoform-Level",
    output_suffix = "_isoform_level"
  )
}

if (length(processing_levels) == 0) {
  stop("ERROR: At least one of GENERATE_GENE_LEVEL or GENERATE_ISOFORM_LEVEL must be TRUE")
}

# ===============================================
# PROCESS EACH LEVEL
# ===============================================

for (level_name in names(processing_levels)) {
  level_config <- processing_levels[[level_name]]
  
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat("PROCESSING:", level_config$label, "\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")
  
  # ===============================================
  # STEP 1: LOCATE SALMON OUTPUT FILES
  # ===============================================
  
  cat("Step 1: Locating Salmon", level_config$label, "output files...\n")
  
  salmon_quant_dir <- file.path(QUANT_DIR, MASTER_REFERENCE)
  
  # Build paths to Salmon quant.sf files
  files <- file.path(salmon_quant_dir, SAMPLE_IDS, "quant.sf")
  names(files) <- SAMPLE_IDS
  
  # Check which files exist
  files_exist <- file.exists(files)
  if (sum(files_exist) == 0) {
    cat("ERROR: No", level_config$label, "quantification files found in", salmon_quant_dir, "\n")
    cat("Skipping this level...\n\n")
    next
  }
  
  current_sample_ids <- SAMPLE_IDS
  if (sum(files_exist) < length(files)) {
    missing <- SAMPLE_IDS[!files_exist]
    cat("WARNING: Missing quantification for samples:", paste(missing, collapse = ", "), "\n")
    files <- files[files_exist]
    current_sample_ids <- SAMPLE_IDS[files_exist]
  }
  
  cat("Found", length(files), "Salmon quantification files\n")
  cat("Samples:", paste(current_sample_ids, collapse = ", "), "\n\n")
  
  # ===============================================
  # STEP 2: IMPORT WITH TXIMPORT
  # ===============================================
  
  cat("Step 2: Importing Salmon data with tximport...\n")
  
  # For Salmon, we need a tx2gene mapping if summarizing to genes
  # Use the gene_trans_map file corresponding to the MASTER_REFERENCE
  INPUT_FASTAS_DIR <- Sys.getenv("INPUT_FASTAS_DIR", "../../0_INPUTs")
  tx2gene_file <- file.path(INPUT_FASTAS_DIR, "mapping", paste0(MASTER_REFERENCE, ".fa.gene_trans_map"))
  # Try alternative naming conventions if not found
  if (!file.exists(tx2gene_file)) {
    tx2gene_file <- file.path(INPUT_FASTAS_DIR, "mapping", paste0(MASTER_REFERENCE, ".fasta.gene_trans_map"))
  }
  
  if (level_config$tx_out == FALSE) {
    # Gene-level: need tx2gene mapping
    if (!file.exists(tx2gene_file)) {
      cat("ERROR: tx2gene mapping file not found:", tx2gene_file, "\n")
      cat("Skipping gene-level processing...\n\n")
      next
    }
    
    # Read tx2gene mapping
    tx2gene <- read.table(tx2gene_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE, 
                         colClasses = c("character", "character"))
    colnames(tx2gene) <- c("GENEID", "TXNAME")
    tx2gene$GENEID <- trimws(tx2gene$GENEID)
    tx2gene$TXNAME <- trimws(tx2gene$TXNAME)
    tx2gene <- tx2gene[, c("TXNAME", "GENEID")]
    cat("Loaded tx2gene mapping:", nrow(tx2gene), "entries\n")
    
    # Import at transcript level, then aggregate to gene level
    cat("  Importing at transcript level and aggregating to genes...\n")
    txi_tx <- tximport(files, type = "salmon", txIn = TRUE, txOut = TRUE,
                      ignoreTxVersion = FALSE, ignoreAfterBar = FALSE)
    
    # Match transcripts to genes and aggregate
    tx_names <- rownames(txi_tx$counts)
    gene_ids <- tx2gene$GENEID[match(tx_names, tx2gene$TXNAME)]
    
    # Handle unmapped transcripts - filter out NAs before rowsum
    unmapped_count <- sum(is.na(gene_ids))
    if (unmapped_count > 0) {
      cat("  WARNING:", unmapped_count, "transcripts not found in tx2gene mapping - excluding from gene-level counts\n")
      valid_idx <- !is.na(gene_ids)
      gene_ids <- gene_ids[valid_idx]
      txi_tx$counts <- txi_tx$counts[valid_idx, , drop = FALSE]
      txi_tx$abundance <- txi_tx$abundance[valid_idx, , drop = FALSE]
      txi_tx$length <- txi_tx$length[valid_idx, , drop = FALSE]
    }
    
    gene_counts <- rowsum(txi_tx$counts, gene_ids, na.rm = TRUE)
    gene_abundance <- rowsum(txi_tx$abundance, gene_ids, na.rm = TRUE)
    gene_length <- rowsum(txi_tx$length, gene_ids, na.rm = TRUE)
    
    txi <- list(
      abundance = gene_abundance,
      counts = gene_counts,
      length = gene_length,
      countsFromAbundance = "no"
    )
  } else {
    # Isoform-level: no tx2gene needed
    cat("  Importing at transcript level...\n")
    txi <- tximport(files, type = "salmon", txIn = TRUE, txOut = TRUE,
                   ignoreTxVersion = FALSE, ignoreAfterBar = FALSE)
  }
  
  entity_type <- if (level_config$tx_out) "transcripts" else "genes"
  cat("Successfully imported data for", ncol(txi$counts), "samples\n")
  cat("Total", entity_type, ":", nrow(txi$counts), "\n\n")

  # ===============================================
  # STEP 3: CREATE SAMPLE METADATA
  # ===============================================
  
  cat("Step 3: Creating sample metadata...\n")
  
  # Get conditions from SAMPLE_LABELS, fallback to sample ID if not found
  conditions <- SAMPLE_LABELS[current_sample_ids]
  missing_labels <- is.na(conditions)
  if (any(missing_labels)) {
    cat("WARNING: No labels found for samples:", paste(current_sample_ids[missing_labels], collapse = ", "), "\n")
    cat("Using sample IDs as condition labels for these samples\n")
    conditions[missing_labels] <- current_sample_ids[missing_labels]
  }
  
  sample_data <- data.frame(
    SampleID = current_sample_ids,
    Condition = conditions,
    row.names = current_sample_ids,
    stringsAsFactors = FALSE
  )
  
  cat("Sample metadata:\n")
  print(sample_data)
  cat("\n")
  
  # ===============================================
  # STEP 4: CHECK FOR REPLICATES
  # ===============================================
  
  cat("Step 4: Checking for biological replicates...\n")
  
  condition_counts <- table(sample_data$Condition)
  has_replicates <- any(condition_counts > 1)
  
  cat("Samples per condition:\n")
  print(condition_counts)
  cat("\n")
  
  if (!has_replicates) {
    cat("WARNING: No biological replicates detected!\n")
    cat("Using TPM normalization for visualization purposes\n\n")
  }
  
  # ===============================================
  # STEP 5: NORMALIZATION
  # ===============================================
  
  # Keep raw counts (NumReads) and TPM (abundance) separately
  raw_counts <- txi$counts
  tpm_matrix <- txi$abundance
  
  if (has_replicates) {
    cat("Step 5: Running DESeq2 normalization...\n")
    
    dds <- DESeqDataSetFromTximport(
      txi = txi,
      colData = sample_data,
      design = ~ Condition
    )
    
    dds <- DESeq(dds)
    normalized_counts <- counts(dds, normalized = TRUE)
    cat("DESeq2 normalization complete\n")
    
  } else {
    cat("Step 5: Computing TPM values...\n")
    # Use raw counts as the "normalized" counts (will also save TPM separately)
    normalized_counts <- raw_counts
    cat("TPM normalization complete (visualization only)\n")
  }
  
  rownames(normalized_counts) <- rownames(normalized_counts)
  cat("Matrix dimensions:", nrow(normalized_counts), entity_type, "x", ncol(normalized_counts), "samples\n\n")
  
  # ===============================================
  # STEP 6: SAVE FULL NORMALIZED MATRIX
  # ===============================================
  
  cat("Step 6: Saving normalized count matrix...\n")
  
  level_output_dir <- file.path(output_dir, level_name)
  dir.create(level_output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Use helper function to save all matrix variants (NumReads, TPM, Gene_ID, Shortened_Name)
  save_count_matrix(normalized_counts, level_output_dir, MASTER_REFERENCE,
                   MASTER_REFERENCE, level_config$output_suffix, SAMPLE_LABELS,
                   tpm_matrix = tpm_matrix)
  cat("\n")
  
  # ===============================================
  # STEP 7: PROCESS GENE GROUPS
  # ===============================================
  
  cat("Step 7: Processing gene groups...\n")
  
  gene_group_files <- list.files(GENE_GROUPS_DIR, pattern = "\\.(csv|txt|tsv)$", full.names = TRUE)
  
  # Filter to only process gene groups specified in GENE_GROUPS_STR (from bash config)
  gene_groups_str <- Sys.getenv("GENE_GROUPS_STR", unset = "")
  if (nzchar(gene_groups_str)) {
    enabled_groups <- trimws(strsplit(gene_groups_str, " ")[[1]])
    gene_group_files <- gene_group_files[tools::file_path_sans_ext(basename(gene_group_files)) %in% enabled_groups]
    cat("Filtering to configured gene groups:", paste(enabled_groups, collapse = ", "), "\n")
  }
  
  if (length(gene_group_files) == 0) {
    cat("No gene group files found in", GENE_GROUPS_DIR, "\n")
  } else {
    cat("Found", length(gene_group_files), "gene group files\n\n")
    
    successful_groups <- 0
    
    # Get combined output folder name (GeneGroup_in_Dataset)
    CURRENT_DATASET <- Sys.getenv("CURRENT_DATASET", unset = "")
  
    for (gene_group_file in gene_group_files) {
      gene_group_name <- tools::file_path_sans_ext(basename(gene_group_file))
      # Generate combined output name with dataset suffix
      output_folder_name <- if (nzchar(CURRENT_DATASET)) {
        paste0(gene_group_name, "_in_", CURRENT_DATASET)
      } else {
        gene_group_name
      }
      cat("  Processing:", gene_group_name, "-> Output:", output_folder_name, "\n")
  
      # Read gene list from CSV (first column is Gene_ID)
      tryCatch({
        if (grepl("\\.csv$", gene_group_file, ignore.case = TRUE)) {
          gene_df <- read.csv(gene_group_file, stringsAsFactors = FALSE, header = TRUE)
          if ("Gene_ID" %in% colnames(gene_df)) {
            gene_list <- gene_df$Gene_ID
          } else {
            gene_list <- gene_df[[1]]  # Fallback to first column
          }
        } else {
          gene_list <- suppressWarnings(readLines(gene_group_file))
          gene_list <- gene_list[!grepl("^#|^Gene_ID", gene_list, ignore.case = TRUE) & nzchar(gene_list)]
        }
        gene_list <- trimws(gene_list)
      }, error = function(e) {
        cat("    Error reading file:", e$message, "\n")
        gene_list <- character(0)
      })
  
      if (length(gene_list) == 0) {
        cat("    Skipping: No genes in list\n")
        next
      }
  
      # Match gene list to row names in data
      # Handle suffix mismatch: gene_list has "SMEL4.1_06g023900" but 
      # Salmon output has "SMEL4.1_06g023900.1" (gene-level) or "SMEL4.1_06g023900.1.01" (isoform-level)
      data_rownames <- rownames(normalized_counts)
      
      # Create mapping from base gene ID to full data ID
      # First strip .X.XX suffix (isoform), then .X suffix (gene-level)
      base_ids <- sub("\\.[0-9]+\\.[0-9]+$", "", data_rownames)  # Strip .X.XX
      base_ids <- sub("\\.[0-9]+$", "", base_ids)                  # Strip remaining .X
      base_to_full <- setNames(data_rownames, base_ids)
      
      # Match genes while preserving CSV order
      genes_in_data <- character(0)
      for (gene in gene_list) {
        # Try exact match first
        if (gene %in% data_rownames) {
          genes_in_data <- c(genes_in_data, gene)
        } else if (gene %in% names(base_to_full)) {
          # Match by base ID (gene in CSV -> full ID in data)
          genes_in_data <- c(genes_in_data, base_to_full[[gene]])
        }
      }
      genes_in_data <- unique(genes_in_data)
  
      if (length(genes_in_data) == 0) {
        cat("    Skipping: No matching", entity_type, "found\n")
        next
      }
  
      cat("    Found", length(genes_in_data), "/", length(gene_list), entity_type, "\n")
  
      # Create gene group subdirectory with dataset suffix (consistent naming)
      gene_group_dir <- file.path(level_output_dir, output_folder_name)
      dir.create(gene_group_dir, recursive = TRUE, showWarnings = FALSE)
      
      # Save gene group subset using helper (saves NumReads, TPM, Gene_ID, Shortened_Name variants)
      subset_counts <- normalized_counts[genes_in_data, , drop = FALSE]
      subset_tpm <- tpm_matrix[genes_in_data, , drop = FALSE]
      save_count_matrix(subset_counts, gene_group_dir, output_folder_name,
                       MASTER_REFERENCE, level_config$output_suffix, SAMPLE_LABELS,
                       tpm_matrix = subset_tpm)
      cat("    Saved subset matrices\n")
      successful_groups <- successful_groups + 1
    }
    
    cat("\n  Summary: Processed", successful_groups, "/", length(gene_group_files), "gene groups\n")
  }
  
  # ===============================================
  # LEVEL SUMMARY
  # ===============================================
  
  cat("\n", paste(rep("=", 70), collapse = ""), "\n")
  cat(level_config$label, "PROCESSING COMPLETE\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")
}

# ===============================================
# FINAL SUMMARY
# ===============================================

cat(paste(rep("=", 70), collapse = ""), "\n")
cat("ALL PROCESSING COMPLETE\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

cat("Generated matrices for", length(processing_levels), "level(s)\n")
cat("Output directory:", output_dir, "\n")
cat("Matrices ready for heatmap generation\n\n")
