# ===============================================
# SHARED R CONFIGURATIONS FOR METHOD 3 (STAR + SALMON)
# ===============================================
# Central configuration file - edit settings here only
# Method 3: STAR Alignment + Salmon Quantification + tximport

# ===============================================
# LIBRARY DEPENDENCIES
# ===============================================

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(ggplot2)
})

# ===============================================
# SAMPLE LABELS CONFIGURATION
# ===============================================
# Maps SRR accessions to human-readable organ/tissue labels
# Used when label_type = "Organ" is selected

SAMPLE_LABELS <- c(
  # Eggplant developmental stages
  "SRR20722387" = "Fert_Ovary",
  "SRR20722232" = "Fert_Anther",
  "SRR23909869" = "Old_Leaf",
  "SRR21010466" = "Fruit_Ripening",
  "SRR3884686" = "Young_Fruit_Peel_A",
  "SRR3884687" = "Young_Fruit_Peel_B"
)

# ===============================================
# INPUT/OUTPUT PATHS
# ===============================================

# Directory structure - relative to the post-processing folder
SCRIPT_DIR <- dirname(sys.frame(1)$ofile)
BASE_DIR <- dirname(SCRIPT_DIR)

# Input directories
# tximport outputs are stored in count_matrices_from_STAR/ or alignments/
MATRIX_DIR <- file.path(BASE_DIR, "count_matrices_from_STAR")
if (!dir.exists(MATRIX_DIR)) {
  # Fallback - check for alternative locations
  alt_dirs <- c(
    file.path(BASE_DIR, "matrices"),
    file.path(BASE_DIR, "alignments", "matrices")
  )
  for (alt in alt_dirs) {
    if (dir.exists(alt)) {
      MATRIX_DIR <- alt
      break
    }
  }
}

# Gene groups input
GENE_GROUPS_DIR <- file.path(BASE_DIR, "A_GeneGroups_InputList")

# Output directories
OUTPUT_DIR <- file.path(BASE_DIR, "Figure_Outputs")
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# ===============================================
# COUNT TYPES CONFIGURATION
# ===============================================
# Method 3 uses Salmon quantification via tximport
# Available count types: counts, abundance (TPM), length-scaled

COUNT_TYPES <- c(
  "counts"     # Raw counts from tximport (for DESeq2)
  # "abundance"  # TPM from Salmon (normalized by length + library)
)

# ===============================================
# LABEL TYPES
# ===============================================
# Controls what's shown on heatmap axes

LABEL_TYPES <- c(
  "Organ"    # Use friendly organ names from SAMPLE_LABELS
  # "SRR"    # Use raw SRR accession numbers
)

# ===============================================
# NORMALIZATION SCHEMES
# ===============================================
# Method 3 uses tximport counts - can apply DESeq2-style normalization

NORM_SCHEMES <- c(
  "count_type_normalized",    # Log2(counts + 1) - appropriate for tximport counts
  "zscore_scaled_to_ten"      # Z-score scaled to 0-10 range
  # "raw",                    # Raw counts (large values, hard to visualize)
  # "zscore"                  # Standard z-score
)

# ===============================================
# DISPLAY OPTIONS
# ===============================================

TRANSPOSE <- c(
  FALSE    # Genes as rows, samples as columns (default)
  # TRUE   # Samples as rows, genes as columns
)

SORT_OPTIONS <- c(
  TRUE,    # Sort by mean expression level
  FALSE    # Keep original order (by organ/sample name)
)

# ===============================================
# VISUALIZATION SETTINGS
# ===============================================

LEGEND_POSITION <- "bottom"  # "bottom" or "right"

# ===============================================
# GENE GROUP FILTERING
# ===============================================
# Set to FALSE to skip filtering (process all genes in matrix)
# Set to TRUE to filter only genes listed in gene group files

SKIP_FILTERING <- FALSE

# ===============================================
# CONSOLIDATED OUTPUT PATHS FOR ADVANCED ANALYSES
# ===============================================

CONSOLIDATED_BASE_DIR <- OUTPUT_DIR
MASTER_REFERENCE <- "All_Smel_Genes"

# ===============================================
# HELPER FUNCTIONS
# ===============================================

# Read configuration from wrapper script files
read_config_file <- function(file_path, default_value, is_boolean = FALSE) {
  if (!file.exists(file_path)) return(default_value)
  value <- trimws(readLines(file_path, warn = FALSE))
  if (is_boolean) return(tolower(value[1]) == "true")
  return(value[nzchar(value)])
}

# Load gene groups and overwrite settings
load_runtime_config <- function() {
  gene_groups <- read_config_file(".gene_groups_temp.txt",
    default_value = c("SmelGRFs", "SmelGIFs"))
  master_reference <- read_config_file(".master_reference_temp.txt",
    default_value = MASTER_REFERENCE)
  if (length(master_reference) > 1) master_reference <- master_reference[1]
  overwrite <- read_config_file(".overwrite_temp.txt",
    default_value = TRUE, is_boolean = TRUE)
  list(gene_groups = gene_groups, master_reference = master_reference, overwrite_existing = overwrite)
}

# Create output directory safely
ensure_output_dir <- function(dir_path, clean = FALSE) {
  if (clean && dir.exists(dir_path)) unlink(dir_path, recursive = TRUE)
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
}

# Validate and read matrix
validate_and_read_matrix <- function(input_file, min_rows = 2) {
  if (!file.exists(input_file)) return(list(success = FALSE, reason = "file not found"))
  matrix_data <- read_count_matrix(input_file)
  if (is.null(matrix_data)) return(list(success = FALSE, reason = "failed to read"))
  if (nrow(matrix_data) < min_rows) return(list(success = FALSE, reason = paste0("need ≥", min_rows, " rows")))
  list(success = TRUE, data = matrix_data, n_genes = nrow(matrix_data))
}
