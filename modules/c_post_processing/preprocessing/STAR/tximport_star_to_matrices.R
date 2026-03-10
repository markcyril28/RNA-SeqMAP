#!/usr/bin/env Rscript

# ===============================================
# TXIMPORT: STAR+SALMON TO DESEQ2 MATRICES
# ===============================================
# Processes Salmon quantification output from STAR alignment (M3) using tximport.
# Produces standardized count matrices consumed by all downstream analysis modules.
#
# Runs from: 3_POST_PROC/M3_STAR_Align/   (via pushd in run_method_analysis)
# Quant files: ../../2_ALIGNMENT_RESULTs/M3_STAR_Align/{MASTER_REFERENCE}/6_salmon/quant/{SRR_ID}/quant.sf
# tx2gene:     count_matrices_from_STAR/{MASTER_REFERENCE}/tx2gene_{MASTER_REFERENCE}.tsv
# Output:      count_matrices_from_STAR/{MASTER_REFERENCE}/{level}/{gene_group}/
#
# Output naming convention (matches build_input_path() in 0_shared_config.R):
#   {gene_group}_{count_type}_{gene_type}_from_{master_ref}_{processing_level}.tsv
# ===============================================

suppressPackageStartupMessages({
  library(tximport)
})

# ===============================================
# CONFIGURATION
# ===============================================

SCRIPT_DIR <- Sys.getenv("ANALYSIS_MODULES_DIR", ".")
source(file.path(SCRIPT_DIR, "0_shared_config.R"))
source(file.path(SCRIPT_DIR, "1_utility_functions.R"))

# Salmon quant output is in the alignment results directory, NOT post-proc.
# Path includes MASTER_REFERENCE (= fasta_tag) to isolate per-reference outputs.
QUANT_DIR     <- file.path("..", "..", "2_ALIGNMENT_RESULTs", "M3_STAR_Align",
                           MASTER_REFERENCE, "6_salmon", "quant")
MATRICES_DIR  <- "count_matrices_from_STAR"
# Note: MASTER_REFERENCE is already set by 0_shared_config.R above; no re-assignment needed.

GENERATE_GENE_LEVEL     <- TRUE
GENERATE_ISOFORM_LEVEL  <- TRUE

# ===============================================
# HELPER: SAVE MATRICES WITH STANDARD NAMING
# ===============================================

save_count_matrix <- function(counts_matrix, output_dir, base_name, master_ref,
                               level_suffix, sample_labels, tpm_matrix = NULL) {
  processing_level <- gsub("^_", "", level_suffix)

  save_matrix <- function(matrix_data, count_type, gene_type) {
    output_file <- file.path(output_dir,
      paste0(base_name, "_", count_type, "_", gene_type,
             "_from_", master_ref, "_", processing_level, ".tsv"))
    matrix_df <- as.data.frame(matrix_data, check.names = FALSE)
    matrix_df <- cbind(GeneID = rownames(matrix_data), matrix_df)
    rownames(matrix_df) <- NULL
    write.table(matrix_df, output_file, sep = "\t", quote = FALSE, row.names = FALSE)
    cat("  Saved:", basename(output_file), "\n")
  }

  # NumReads (raw counts) - with SRR IDs and with Organ labels
  save_matrix(counts_matrix, "NumReads", "Gene_ID")
  counts_organ <- convert_to_organ_labels(counts_matrix)
  save_matrix(counts_organ, "NumReads", "Shortened_Name")

  # TPM (normalized abundance)
  if (!is.null(tpm_matrix)) {
    save_matrix(tpm_matrix, "tpm", "Gene_ID")
    tpm_organ <- convert_to_organ_labels(tpm_matrix)
    save_matrix(tpm_organ, "tpm", "Shortened_Name")
  }
}

# ===============================================
# LOCATE TX2GENE MAPPING
# ===============================================

find_tx2gene <- function(matrices_dir, master_ref) {
  # Created by star_alignment_pipeline() in m3_star_alignment.sh
  candidate <- file.path(matrices_dir, master_ref, paste0("tx2gene_", master_ref, ".tsv"))
  if (file.exists(candidate)) return(candidate)

  # Fallback: search for any tx2gene file under the master_ref directory
  search_dir <- file.path(matrices_dir, master_ref)
  if (dir.exists(search_dir)) {
    hits <- list.files(search_dir, pattern = "^tx2gene.*\\.tsv$", full.names = TRUE)
    if (length(hits) > 0) {
      cat("  Using tx2gene:", hits[1], "\n")
      return(hits[1])
    }
  }

  # Fallback: use inputs/mapping gene_trans_map (same as M4 Salmon)
  input_dir <- Sys.getenv("INPUT_FASTAS_DIR", file.path("..", "..", "inputs"))
  alt <- file.path(input_dir, "mapping", paste0(master_ref, ".fa.gene_trans_map"))
  if (file.exists(alt)) return(alt)
  alt2 <- file.path(input_dir, "mapping", paste0(master_ref, ".fasta.gene_trans_map"))
  if (file.exists(alt2)) return(alt2)

  return(NULL)
}

# ===============================================
# BANNER
# ===============================================

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("TXIMPORT: STAR+SALMON QUANTIFICATION TO MATRICES\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")
cat("Master Reference:", MASTER_REFERENCE, "\n")
cat("Quant directory: ", QUANT_DIR, "\n")
cat("Output directory:", MATRICES_DIR, "\n\n")

# ===============================================
# DEFINE PROCESSING LEVELS
# ===============================================

processing_levels <- list()
if (GENERATE_GENE_LEVEL) {
  processing_levels[["gene_level"]] <- list(
    tx_out        = FALSE,
    label         = "Gene-Level",
    output_suffix = "_gene_level"
  )
}
if (GENERATE_ISOFORM_LEVEL) {
  processing_levels[["isoform_level"]] <- list(
    tx_out        = TRUE,
    label         = "Isoform-Level",
    output_suffix = "_isoform_level"
  )
}

if (length(processing_levels) == 0) {
  stop("At least one of GENERATE_GENE_LEVEL or GENERATE_ISOFORM_LEVEL must be TRUE")
}

# ===============================================
# PROCESS EACH LEVEL
# ===============================================

for (level_name in names(processing_levels)) {
  level_config <- processing_levels[[level_name]]

  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat("PROCESSING:", level_config$label, "\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")

  # -------------------------------------------------
  # STEP 1: LOCATE QUANT.SF FILES
  # -------------------------------------------------
  cat("Step 1: Locating Salmon quant.sf files...\n")

  files <- file.path(QUANT_DIR, SAMPLE_IDS, "quant.sf")
  names(files) <- SAMPLE_IDS

  files_exist    <- file.exists(files)
  current_ids    <- SAMPLE_IDS

  if (sum(files_exist) == 0) {
    cat("ERROR: No quant.sf files found under", QUANT_DIR, "\n")
    cat("Ensure STAR+Salmon alignment completed successfully.\n\n")
    next
  }
  if (sum(files_exist) < length(files)) {
    missing <- SAMPLE_IDS[!files_exist]
    cat("WARNING: Missing quantification for:", paste(missing, collapse = ", "), "\n")
    files       <- files[files_exist]
    current_ids <- SAMPLE_IDS[files_exist]
  }
  cat("Found", length(files), "quant.sf files\n\n")

  # -------------------------------------------------
  # STEP 2: LOAD TX2GENE MAPPING (gene-level only)
  # -------------------------------------------------
  tx2gene <- NULL
  if (!level_config$tx_out) {
    cat("Step 2: Loading tx2gene mapping...\n")
    tx2gene_file <- find_tx2gene(MATRICES_DIR, MASTER_REFERENCE)

    if (is.null(tx2gene_file)) {
      cat("ERROR: tx2gene mapping not found. Cannot summarize to gene level.\n")
      cat("Run STAR+Salmon alignment first to generate tx2gene mapping.\n\n")
      next
    }

    # Detect column order: tximport needs c(TXNAME, GENEID)
    raw <- read.delim(tx2gene_file, header = FALSE, stringsAsFactors = FALSE,
                      colClasses = "character")
    # star_alignment_pipeline writes: transcript_id TAB gene_id
    tx2gene <- raw[, 1:2, drop = FALSE]
    colnames(tx2gene) <- c("TXNAME", "GENEID")
    tx2gene$TXNAME <- trimws(tx2gene$TXNAME)
    tx2gene$GENEID <- trimws(tx2gene$GENEID)
    cat("Loaded tx2gene:", nrow(tx2gene), "entries\n\n")
  }

  # -------------------------------------------------
  # STEP 3: TXIMPORT
  # -------------------------------------------------
  cat("Step 3: Importing with tximport...\n")

  txi <- tryCatch({
    if (!level_config$tx_out) {
      tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
    } else {
      tximport(files, type = "salmon", txIn = TRUE, txOut = TRUE,
               ignoreTxVersion = FALSE, ignoreAfterBar = FALSE)
    }
  }, error = function(e) {
    cat("ERROR in tximport:", conditionMessage(e), "\n")
    NULL
  })

  if (is.null(txi)) {
    cat("Skipping level:", level_name, "\n\n")
    next
  }

  entity_type <- if (level_config$tx_out) "transcripts" else "genes"
  cat("Imported:", ncol(txi$counts), "samples,", nrow(txi$counts), entity_type, "\n\n")

  # -------------------------------------------------
  # STEP 4: SAMPLE METADATA
  # -------------------------------------------------
  cat("Step 4: Building sample metadata...\n")
  conditions <- SAMPLE_LABELS[current_ids]
  missing_labels <- is.na(conditions)
  if (any(missing_labels)) {
    cat("WARNING: No labels for:", paste(current_ids[missing_labels], collapse = ", "), "\n")
    conditions[missing_labels] <- current_ids[missing_labels]
  }
  sample_data <- data.frame(
    SampleID  = current_ids,
    Condition = unname(conditions),
    row.names = current_ids,
    stringsAsFactors = FALSE
  )
  cat("Conditions per tissue:\n")
  print(table(sample_data$Condition))
  cat("\n")

  # -------------------------------------------------
  # STEP 5: EXTRACT MATRICES
  # -------------------------------------------------
  raw_counts <- txi$counts
  tpm_matrix <- txi$abundance

  # -------------------------------------------------
  # STEP 6: SAVE FULL REFERENCE MATRIX
  # -------------------------------------------------
  cat("Step 6: Saving full reference matrices...\n")
  level_output_dir <- file.path(MATRICES_DIR, MASTER_REFERENCE, level_name)
  dir.create(level_output_dir, recursive = TRUE, showWarnings = FALSE)

  save_count_matrix(raw_counts, level_output_dir,
                    base_name    = MASTER_REFERENCE,
                    master_ref   = MASTER_REFERENCE,
                    level_suffix = level_config$output_suffix,
                    sample_labels = SAMPLE_LABELS,
                    tpm_matrix   = tpm_matrix)
  cat("\n")

  # -------------------------------------------------
  # STEP 7: PROCESS GENE GROUPS
  # -------------------------------------------------
  cat("Step 7: Processing gene groups...\n")

  gene_group_files <- list.files(GENE_GROUPS_DIR, pattern = "\\.(csv|txt|tsv)$",
                                  full.names = TRUE)

  # Filter to only configured gene groups
  gene_groups_str <- Sys.getenv("GENE_GROUPS_STR", unset = "")
  if (nzchar(gene_groups_str)) {
    enabled_groups <- trimws(strsplit(gene_groups_str, " ")[[1]])
    gene_group_files <- gene_group_files[
      tools::file_path_sans_ext(basename(gene_group_files)) %in% enabled_groups
    ]
    cat("Gene groups to process:", paste(enabled_groups, collapse = ", "), "\n")
  }

  if (length(gene_group_files) == 0) {
    cat("No matching gene group files found in", GENE_GROUPS_DIR, "\n")
  } else {
    cat("Found", length(gene_group_files), "gene group file(s)\n\n")

    CURRENT_DATASET <- Sys.getenv("CURRENT_DATASET", unset = "")
    successful_groups <- 0

    for (gene_group_file in gene_group_files) {
      gene_group_name  <- tools::file_path_sans_ext(basename(gene_group_file))
      output_folder_name <- if (nzchar(CURRENT_DATASET)) {
        paste0(gene_group_name, "_in_", CURRENT_DATASET)
      } else {
        gene_group_name
      }
      cat("  Processing:", gene_group_name, "->", output_folder_name, "\n")

      # Read gene list
      gene_list <- tryCatch({
        if (grepl("\\.csv$", gene_group_file, ignore.case = TRUE)) {
          gdf <- read.csv(gene_group_file, stringsAsFactors = FALSE, header = TRUE)
          if ("Gene_ID" %in% colnames(gdf)) gdf$Gene_ID else gdf[[1]]
        } else {
          raw_lines <- suppressWarnings(readLines(gene_group_file))
          raw_lines <- raw_lines[
            !grepl("^#|^Gene_ID", raw_lines, ignore.case = TRUE) & nzchar(raw_lines)
          ]
          raw_lines
        }
      }, error = function(e) {
        cat("    Error reading gene list:", conditionMessage(e), "\n")
        character(0)
      })
      gene_list <- trimws(gene_list[nzchar(gene_list)])

      if (length(gene_list) == 0) {
        cat("    Skipping: empty gene list\n")
        next
      }

      # Match genes - handle version suffixes (e.g. GENE.1, GENE.1.01)
      data_rownames <- rownames(raw_counts)
      base_ids <- sub("\\.[0-9]+\\.[0-9]+$", "", data_rownames)
      base_ids <- sub("\\.[0-9]+$", "", base_ids)
      base_to_full <- setNames(data_rownames, base_ids)

      genes_in_data <- character(0)
      for (gene in gene_list) {
        if (gene %in% data_rownames) {
          genes_in_data <- c(genes_in_data, gene)
        } else if (gene %in% names(base_to_full)) {
          genes_in_data <- c(genes_in_data, base_to_full[[gene]])
        }
      }
      genes_in_data <- unique(genes_in_data)

      if (length(genes_in_data) == 0) {
        cat("    No matching", entity_type, "found\n")
        next
      }
      cat("    Matched", length(genes_in_data), "/", length(gene_list), entity_type, "\n")

      # Save gene-group subset
      gene_group_dir <- file.path(level_output_dir, output_folder_name)
      dir.create(gene_group_dir, recursive = TRUE, showWarnings = FALSE)

      subset_counts <- raw_counts[genes_in_data, , drop = FALSE]
      subset_tpm    <- tpm_matrix[genes_in_data, , drop = FALSE]

      save_count_matrix(subset_counts, gene_group_dir,
                        base_name    = output_folder_name,
                        master_ref   = MASTER_REFERENCE,
                        level_suffix = level_config$output_suffix,
                        sample_labels = SAMPLE_LABELS,
                        tpm_matrix   = subset_tpm)
      successful_groups <- successful_groups + 1
    }

    cat("\n  Summary:", successful_groups, "/", length(gene_group_files),
        "gene groups processed\n")
  }

  cat("\n", paste(rep("=", 70), collapse = ""), "\n")
  cat(level_config$label, "COMPLETE\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")
}

cat(paste(rep("=", 70), collapse = ""), "\n")
cat("ALL LEVELS COMPLETE\n")
cat("Output directory:", MATRICES_DIR, "/", MASTER_REFERENCE, "\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")
