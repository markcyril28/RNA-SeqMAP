#!/usr/bin/env Rscript
# ===============================================
# BASIC HEATMAP GENERATOR - METHOD 3 (STAR + SALMON)
# ===============================================
# Generates heatmaps for gene expression data from tximport
# Input: TSV count matrix from tximport (Salmon quant)

# ===============================================
# LOAD DEPENDENCIES
# ===============================================

script_dir <- dirname(sys.frame(1)$ofile)

# Load shared configurations
config_path <- file.path(script_dir, "0_shared_config.R")
if (!file.exists(config_path)) {
  stop("Configuration file not found: ", config_path)
}
source(config_path)

# Load utility functions
utility_path <- file.path(script_dir, "2_utility_functions.R")
if (!file.exists(utility_path)) {
  stop("Utility functions file not found: ", utility_path)
}
source(utility_path)

# Load processing engine
processing_path <- file.path(script_dir, "1_processing_engine.R")
if (!file.exists(processing_path)) {
  stop("Processing engine file not found: ", processing_path)
}
source(processing_path)

cat("=== Basic Heatmap Generator (Method 3 - STAR/Salmon) ===\n")

# ===============================================
# PROCESSING FUNCTION
# ===============================================

process_single_heatmap <- function(data_matrix, output_path, title, count_type, label_type, 
                                   normalization_scheme, transpose, sort_by_expression) {
  return(generate_heatmap_violet(
    data_matrix = data_matrix,
    output_path = output_path,
    title = title,
    count_type = count_type,
    label_type = label_type,
    normalization_type = switch(normalization_scheme,
      "raw" = "Raw Counts",
      "count_type_normalized" = "Log2_Normalized",
      "tpm_normalized" = "TPM_Normalized",
      "zscore" = "Z-score_Normalized",
      "zscore_scaled_to_ten" = "Z-score_Scaled_to_Ten",
      normalization_scheme
    ),
    transpose = transpose,
    sort_by_expression = sort_by_expression
  ))
}

# ===============================================
# MAIN EXECUTION
# ===============================================

process_count_matrices(
  processing_func = process_single_heatmap,
  output_subdir = "1_heatmaps",
  output_suffix = "_heatmap.png"
)

cat("\n=== Heatmap generation complete ===\n")
