#!/usr/bin/env Rscript

# ===============================================
# CROSS-METHOD CONCORDANCE - CONFIGURATION
# ===============================================
# Shared configuration for all concordance analysis modules.
# Sources the main pipeline's shared_config and utility_functions,
# then sets concordance-specific parameters.

# -----------------------------------------------
# Source pipeline shared infrastructure
# -----------------------------------------------

ANALYSIS_MODULES_DIR <- Sys.getenv("ANALYSIS_MODULES_DIR", ".")

# Source pipeline shared infrastructure if available; concordance defines its own
# helpers below and reads critical values from environment variables, so these are
# optional rather than fatal.
.shared_config_path <- file.path(ANALYSIS_MODULES_DIR, "0_shared_config.R")
.utility_funcs_path <- file.path(ANALYSIS_MODULES_DIR, "1_utility_functions.R")
if (file.exists(.shared_config_path)) {
  source(.shared_config_path)
} else {
  message("[CONCORDANCE CONFIG] WARN: shared config not found: ", .shared_config_path)
  message("[CONCORDANCE CONFIG]       Using environment variables for configuration.")
}
if (file.exists(.utility_funcs_path)) {
  source(.utility_funcs_path)
} else {
  message("[CONCORDANCE CONFIG] WARN: utility functions not found: ", .utility_funcs_path)
}
rm(.shared_config_path, .utility_funcs_path)

# Source match_gene_ids if not already defined (needed by ranking stability).
# The function may be in 1_utility_functions.R (sourced above) or in the standalone
# utilities/match_gene_ids.R file. Try the standalone file first, then 1_utility_functions.R.
if (!exists("match_gene_ids", mode = "function")) {
  .match_ids_standalone <- file.path(dirname(ANALYSIS_MODULES_DIR), "utilities", "match_gene_ids.R")
  .util_funcs_fallback <- file.path(ANALYSIS_MODULES_DIR, "1_utility_functions.R")
  if (file.exists(.match_ids_standalone)) {
    source(.match_ids_standalone)
  } else if (file.exists(.util_funcs_fallback)) {
    source(.util_funcs_fallback)
  } else {
    stop("[CONCORDANCE CONFIG] match_gene_ids() not defined and neither match_gene_ids.R nor 1_utility_functions.R found.",
         "\n  Searched: ", .match_ids_standalone,
         "\n  Searched: ", .util_funcs_fallback,
         "\n  This function is required by ranking stability analysis.")
  }
  rm(.match_ids_standalone, .util_funcs_fallback)
}

# -----------------------------------------------
# Concordance-specific configuration
# -----------------------------------------------

CONCORDANCE_SCRIPT_DIR <- Sys.getenv("CONCORDANCE_SCRIPT_DIR", ".")
OUTPUT_DIR <- Sys.getenv("OUTPUT_DIR", "3_POST_PROC/cross_method_concordance")
ALIGNMENT_BASE <- Sys.getenv("ALIGNMENT_BASE", "2_ALIGNMENT_RESULTs")
POST_PROC_BASE <- Sys.getenv("POST_PROC_BASE", "3_POST_PROC")
BASE_DIR <- Sys.getenv("BASE_DIR", ".")

# Ensure MASTER_REFERENCE is set (may already be defined by 0_shared_config.R)
if (!exists("MASTER_REFERENCE") || !nzchar(MASTER_REFERENCE)) {
  MASTER_REFERENCE <- Sys.getenv("MASTER_REFERENCE", "GPE001970_genome")
}

# Ensure GENE_GROUPS_DIR is set (may already be defined by 0_shared_config.R)
if (!exists("GENE_GROUPS_DIR") || !nzchar(GENE_GROUPS_DIR)) {
  GENE_GROUPS_DIR <- Sys.getenv("GENE_GROUPS_DIR", "")
  if (!nzchar(GENE_GROUPS_DIR)) {
    GENE_GROUPS_DIR <- file.path(BASE_DIR, "inputs", "gene_groups_csv")
  }
}

# Methods to compare
METHODS_STR <- Sys.getenv("METHODS", "M1_HISAT2_RefGuided M2_HISAT2_DeNovo M3_STAR_Align M4_Salmon_Saf M5_RSEM_Bowtie2")
CONCORDANCE_METHODS <- trimws(strsplit(trimws(METHODS_STR), "\\s+")[[1]])
CONCORDANCE_METHODS <- CONCORDANCE_METHODS[nzchar(CONCORDANCE_METHODS)]

# Short labels for display
METHOD_SHORT_NAMES <- c(
  M1_HISAT2_RefGuided = "M1:HISAT2-Ref",
  M2_HISAT2_DeNovo    = "M2:HISAT2-DeNovo",
  M3_STAR_Align       = "M3:STAR-Salmon",
  M4_Salmon_Saf       = "M4:Salmon-SAF",
  M5_RSEM_Bowtie2     = "M5:Bowtie2-RSEM"
)

# Parse method-specific reference directory mapping
# Format: "M1_HISAT2_RefGuided=GPE001970_genome;M2_HISAT2_DeNovo=GPE001970_transcripts;..."
METHOD_REF_DIRS_STR <- Sys.getenv("METHOD_REF_DIRS_STR", "")
METHOD_REF_DIRS <- list()
if (nzchar(METHOD_REF_DIRS_STR)) {
  pairs <- strsplit(METHOD_REF_DIRS_STR, ";")[[1]]
  for (pair in pairs) {
    kv <- strsplit(pair, "=")[[1]]
    if (length(kv) == 2) {
      METHOD_REF_DIRS[[trimws(kv[1])]] <- trimws(kv[2])
    }
  }
}

# Gene groups for ranking stability
GENE_GROUPS_STR <- Sys.getenv("GENE_GROUPS", "SmelDMPs_v5_with_18s_and_HAP2,Selected_SmelGRF-GIF_with_two_GIF")
CONCORDANCE_GENE_GROUPS <- trimws(strsplit(GENE_GROUPS_STR, ",")[[1]])

# Figure output paths
FIGURES_DIR <- file.path(OUTPUT_DIR, "figures")
TABLES_DIR <- file.path(OUTPUT_DIR, "tables")
dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLES_DIR, recursive = TRUE, showWarnings = FALSE)

# Harmonized data intermediate file
HARMONIZED_RDS <- file.path(OUTPUT_DIR, "harmonized_tpm_matrices.rds")

# Analysis parameters
CONCORDANCE_MIN_EXPR <- 0.1       # Minimum TPM to consider a gene "expressed"
CONCORDANCE_MIN_SAMPLES <- 3      # Gene must be expressed in at least this many samples
RANKING_CHANGE_THRESHOLD <- 0.3   # Fractional rank change above this is "drastic"
CORRELATION_MIN_GENES <- 10       # Minimum nonzero genes per sample for pairwise correlation

# -----------------------------------------------
# Helper: get reference directory for a method
# -----------------------------------------------

get_method_ref_dir <- function(method) {
  if (method %in% names(METHOD_REF_DIRS)) {
    return(METHOD_REF_DIRS[[method]])
  }
  # Fallback: use MASTER_REFERENCE as-is
  return(MASTER_REFERENCE)
}

# -----------------------------------------------
# Helper: get short method name for display
# -----------------------------------------------

get_short_name <- function(method) {
  if (method %in% names(METHOD_SHORT_NAMES)) {
    return(METHOD_SHORT_NAMES[[method]])
  }
  return(method)
}

cat("[CONCORDANCE CONFIG] Methods:", paste(CONCORDANCE_METHODS, collapse = ", "), "\n")
cat("[CONCORDANCE CONFIG] Reference:", MASTER_REFERENCE, "\n")
cat("[CONCORDANCE CONFIG] Gene groups:", paste(CONCORDANCE_GENE_GROUPS, collapse = ", "), "\n")
cat("[CONCORDANCE CONFIG] Output:", OUTPUT_DIR, "\n")
