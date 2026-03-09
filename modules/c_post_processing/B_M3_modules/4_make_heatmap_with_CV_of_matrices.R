#!/usr/bin/env Rscript
# ===============================================
# HEATMAP WITH CV GENERATOR - METHOD 3 (STAR + SALMON)
# ===============================================
# Generates heatmaps with Coefficient of Variation annotations
# CV calculated on normalized counts

# ===============================================
# LOAD DEPENDENCIES
# ===============================================

script_dir <- dirname(sys.frame(1)$ofile)

source(file.path(script_dir, "0_shared_config.R"))
source(file.path(script_dir, "2_utility_functions.R"))
source(file.path(script_dir, "1_processing_engine.R"))

cat("=== Heatmap with CV Generator (Method 3 - STAR/Salmon) ===\n")

# ===============================================
# PROCESSING FUNCTION
# ===============================================

process_single_heatmap_cv <- function(data_matrix, output_path, title, count_type, label_type, 
                                      normalization_scheme, transpose, sort_by_expression) {
  # Get normalized data for CV calculation
  raw_matrix <- preprocess_for_count_type_normalized(data_matrix, count_type)
  vis_matrix <- apply_normalization(data_matrix, normalization_scheme, count_type)
  
  return(generate_heatmap_with_cv(
    data_matrix = vis_matrix,
    output_path = output_path,
    title = title,
    count_type = count_type,
    label_type = label_type,
    normalization_scheme = normalization_scheme,
    transpose = transpose,
    sort_by_expression = sort_by_expression,
    raw_data_matrix = raw_matrix
  ))
}

# ===============================================
# MAIN EXECUTION
# ===============================================

process_count_matrices(
  processing_func = process_single_heatmap_cv,
  output_subdir = "2_heatmaps_with_cv",
  output_suffix = "_heatmap_cv.png"
)

cat("\n=== Heatmap with CV generation complete ===\n")
