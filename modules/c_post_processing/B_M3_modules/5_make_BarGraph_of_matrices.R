#!/usr/bin/env Rscript
# ===============================================
# BAR GRAPH GENERATOR - METHOD 3 (STAR + SALMON)
# ===============================================
# Generates individual bar graphs for each gene

# ===============================================
# LOAD DEPENDENCIES
# ===============================================

script_dir <- dirname(sys.frame(1)$ofile)

source(file.path(script_dir, "0_shared_config.R"))
source(file.path(script_dir, "2_utility_functions.R"))
source(file.path(script_dir, "1_processing_engine.R"))

cat("=== Bar Graph Generator (Method 3 - STAR/Salmon) ===\n")

# ===============================================
# PROCESSING FUNCTION
# ===============================================

process_bargraphs_for_matrix <- function(data_matrix, output_path, title, count_type, label_type, 
                                         normalization_scheme, transpose, sort_by_expression) {
  output_dir <- dirname(output_path)
  base_name <- tools::file_path_sans_ext(basename(output_path))
  bargraph_dir <- file.path(output_dir, base_name)
  if (!dir.exists(bargraph_dir)) dir.create(bargraph_dir, recursive = TRUE)
  
  normalized_data <- apply_normalization(data_matrix, normalization_scheme, count_type)
  
  n_successful <- generate_gene_bargraphs(
    data_matrix = normalized_data,
    output_dir = bargraph_dir,
    title_prefix = title,
    count_type = count_type,
    label_type = label_type,
    normalization_scheme = normalization_scheme,
    sort_by_expression = sort_by_expression
  )
  
  cat("  Generated", n_successful, "bar graphs for:", title, "\n")
  return(n_successful > 0)
}

# ===============================================
# MAIN EXECUTION
# ===============================================

process_count_matrices(
  processing_func = process_bargraphs_for_matrix,
  output_subdir = "3_bargraphs",
  output_suffix = "_bargraphs"
)

cat("\n=== Bar graph generation complete ===\n")
