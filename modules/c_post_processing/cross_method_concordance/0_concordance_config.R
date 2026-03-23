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

# -----------------------------------------------
# Concordance-specific configuration
# -----------------------------------------------

CONCORDANCE_SCRIPT_DIR <- Sys.getenv("CONCORDANCE_SCRIPT_DIR", ".")
OUTPUT_DIR <- Sys.getenv("OUTPUT_DIR", "3_POST_PROC")
ALIGNMENT_BASE <- Sys.getenv("ALIGNMENT_BASE", "2_ALIGNMENT_RESULTs")
POST_PROC_BASE <- Sys.getenv("POST_PROC_BASE", "3_POST_PROC")
BASE_DIR <- Sys.getenv("BASE_DIR", ".")

# Ensure MASTER_REFERENCE is set (may already be defined by 0_shared_config.R).
# Default matches 0_shared_config.R ("Eggplant_V4.1") to avoid silent mismatch.
if (!exists("MASTER_REFERENCE") || !nzchar(MASTER_REFERENCE)) {
  MASTER_REFERENCE <- Sys.getenv("MASTER_REFERENCE", "Eggplant_V4.1")
}

# Ensure GENE_GROUPS_DIR is set (may already be defined by 0_shared_config.R)
if (!exists("GENE_GROUPS_DIR") || !nzchar(GENE_GROUPS_DIR)) {
  GENE_GROUPS_DIR <- Sys.getenv("GENE_GROUPS_DIR", "")
  if (!nzchar(GENE_GROUPS_DIR)) {
    GENE_GROUPS_DIR <- file.path(BASE_DIR, "inputs", "3_post_proc_inputs", "gene_groups_csv")
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
  pairs <- trimws(strsplit(METHOD_REF_DIRS_STR, ";")[[1]])
  kv_list <- strsplit(pairs, "=")
  # Vectorize: trimws each side once via vapply instead of per-iteration trimws
  for (kv in kv_list) {
    if (length(kv) == 2) {
      METHOD_REF_DIRS[[kv[1]]] <- trimws(kv[2])
    }
  }
}

# -----------------------------------------------
# Concordance mode: cross_method (default), cross_genome, cross_gene_group
# -----------------------------------------------

CONCORDANCE_MODE <- Sys.getenv("CONCORDANCE_MODE", "cross_method")

# Fixed method for cross_genome and cross_gene_group modes
FIXED_METHOD <- Sys.getenv("FIXED_METHOD", "")

# Genomes to compare in cross_genome mode (semicolon-separated)
CONCORDANCE_GENOMES_STR <- Sys.getenv("CONCORDANCE_GENOMES", "")
CONCORDANCE_GENOMES <- if (nzchar(CONCORDANCE_GENOMES_STR)) {
  trimws(strsplit(CONCORDANCE_GENOMES_STR, ";")[[1]])
} else {
  character(0)
}

# Gene groups for concordance analysis
GENE_GROUPS_STR <- Sys.getenv("GENE_GROUPS", "SmelDMPs_v5_with_18s_and_HAP2,Selected_SmelGRF-GIF_with_two_GIF")
CONCORDANCE_GENE_GROUPS <- trimws(strsplit(GENE_GROUPS_STR, ",")[[1]])

# Per-genome gene group CSV mapping for cross-genome orthology (positional row correspondence).
# Format: "genome1=csv1,csv2;genome2=csv3,csv4"
# Parsed into a named list: list(genome1 = c("csv1","csv2"), genome2 = c("csv3","csv4"))
GENOME_GENE_GROUPS_MAP_STR <- Sys.getenv("GENOME_GENE_GROUPS_MAP", "")
GENOME_GENE_GROUPS_MAP <- list()
if (nzchar(GENOME_GENE_GROUPS_MAP_STR)) {
  .ggm_pairs <- strsplit(GENOME_GENE_GROUPS_MAP_STR, ";")[[1]]
  .ggm_kv <- strsplit(.ggm_pairs, "=")
  for (.kv in .ggm_kv) {
    if (length(.kv) == 2) {
      .genome <- trimws(.kv[1])
      .csvs <- trimws(strsplit(trimws(.kv[2]), ",")[[1]])
      GENOME_GENE_GROUPS_MAP[[.genome]] <- .csvs
    }
  }
  # Guard: .genome/.csvs are only assigned when at least one entry has a valid "=" separator
  rm(list = intersect(c(".ggm_pairs", ".ggm_kv", ".kv", ".genome", ".csvs"), ls(all.names = TRUE)))
}
if (length(GENOME_GENE_GROUPS_MAP) > 0) {
  cat("[CONCORDANCE CONFIG] Genome gene groups map:\n")
  for (.g in names(GENOME_GENE_GROUPS_MAP)) {
    cat("  ", .g, ":", paste(GENOME_GENE_GROUPS_MAP[[.g]], collapse = ", "), "\n")
  }
  if (exists(".g")) rm(.g)
}

# Figure resolution (DPI) — reuse value from 0_shared_config.R if already set,
# otherwise compute from env var (300–600, default 300). Avoids redundant Sys.getenv().
if (!exists("FIGURE_DPI") || !is.integer(FIGURE_DPI)) {
  FIGURE_DPI <- as.integer(Sys.getenv("FIGURE_DPI", "300"))
  if (is.na(FIGURE_DPI) || FIGURE_DPI < 72) FIGURE_DPI <- 300L
  FIGURE_DPI <- max(300L, min(600L, FIGURE_DPI))
}

# ---------------------------------------------------------------------------
# Auto-sizing helper for ComplexHeatmap figures
# ---------------------------------------------------------------------------
# Computes png pixel dimensions and draw() padding from content metrics so
# labels, titles, dendrograms, and legends are never clipped.
#
# Returns list(width, height, padding) where width/height are in pixels at
# FIGURE_DPI, and padding is a unit(c(bottom,left,top,right), "mm") vector.
#
# Arguments:
#   row_labels    — character vector of row labels
#   col_labels    — character vector of column labels
#   col_rot       — column label rotation in degrees (default 45)
#   hm_body_cm    — heatmap body size in cm (width, height). Scalar or length-2.
#   has_dendro    — whether dendrograms are shown (adds space)
#   has_title     — whether a column_title is shown
#   legend_width_cm — approximate legend width in cm (default 4)
#   font_size     — label font size in pt (used to estimate text extents)
calc_figure_layout <- function(row_labels, col_labels,
                               col_rot = 45, hm_body_cm = c(14, 14),
                               has_dendro = TRUE, has_title = TRUE,
                               legend_width_cm = 4, font_size = 13,
                               title_lines = 1) {
  if (length(hm_body_cm) == 1) hm_body_cm <- rep(hm_body_cm, 2)

  # Estimate max label width in cm (~0.022 cm per character per pt at 300 DPI)
  char_width_cm <- font_size * 0.022
  max_row_chars <- max(nchar(row_labels), na.rm = TRUE)
  max_col_chars <- max(nchar(col_labels), na.rm = TRUE)

  row_label_cm <- max_row_chars * char_width_cm
  # Rotated column labels project both vertically and horizontally
  col_rot_rad <- col_rot * pi / 180
  col_label_height_cm <- max_col_chars * char_width_cm * abs(sin(col_rot_rad))
  col_label_width_cm  <- max_col_chars * char_width_cm * abs(cos(col_rot_rad))

  # Padding in mm: bottom (col labels), left (row labels), top (title+dendro), right (legend)
  # title_lines scales top padding: each line ~8mm at 14pt; floor at 15mm for single-line titles
  pad_bottom  <- max(20, col_label_height_cm * 10 + 15)
  pad_left    <- max(15, row_label_cm * 10 + 10)
  title_h_mm  <- if (has_title) max(15, title_lines * 8 + 5) else 0
  pad_top     <- 15 + title_h_mm + (if (has_dendro) 10 else 0)
  pad_right   <- max(20, legend_width_cm * 10 + 15)

  # Total figure size in cm
  total_w_cm <- pad_left / 10 + hm_body_cm[1] + pad_right / 10 +
                (if (has_dendro) 1.5 else 0)
  total_h_cm <- pad_top / 10 + hm_body_cm[2] + pad_bottom / 10 +
                (if (has_dendro) 1.5 else 0)

  # Convert cm to pixels at FIGURE_DPI
  cm_to_px <- FIGURE_DPI / 2.54
  width_px  <- ceiling(total_w_cm * cm_to_px)
  height_px <- ceiling(total_h_cm * cm_to_px)

  list(
    width   = width_px,
    height  = height_px,
    padding = unit(c(pad_bottom, pad_left, pad_top, pad_right), "mm")
  )
}

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
# Helper: get short display name for a comparison item
# -----------------------------------------------
# Works for methods, genomes, and gene groups depending on mode.

get_short_name <- function(item) {
  if (item %in% names(METHOD_SHORT_NAMES)) {
    return(METHOD_SHORT_NAMES[[item]])
  }
  # For genomes: strip common suffixes for brevity
  short <- sub("_genome$", "", item)
  short <- sub("_transcripts.*$", ":tx", short)
  return(short)
}

# Mode-aware heatmap title helper
get_concordance_title <- function() {
  switch(CONCORDANCE_MODE,
    cross_method          = "Cross-Method Concordance (Median Spearman Correlation)",
    cross_genome          = paste0("Cross-Genome Concordance (Median Spearman, ", get_short_name(FIXED_METHOD), ")"),
    cross_gene_group      = paste0("Cross-Gene-Group Concordance (Median Spearman, ", get_short_name(FIXED_METHOD), ")"),
    cross_equivalent_gene = paste0("Equivalent-Gene Concordance Across Genomes (", get_short_name(FIXED_METHOD), ")"),
    "Concordance (Median Spearman Correlation)"
  )
}

# -----------------------------------------------
# Shared fast CSV reader (used by gene group loading in multiple scripts)
# data.table::fread is 10-50x faster than read.csv for larger files
# -----------------------------------------------
.conc_use_dt <- requireNamespace("data.table", quietly = TRUE)

.fast_read_csv <- function(path, ...) {
  if (.conc_use_dt) {
    data.table::fread(path, data.table = FALSE, ...)
  } else {
    read.csv(path, stringsAsFactors = FALSE, header = TRUE, ...)
  }
}

cat("[CONCORDANCE CONFIG] Mode:", CONCORDANCE_MODE, "\n")
if (nzchar(FIXED_METHOD)) cat("[CONCORDANCE CONFIG] Fixed method:", FIXED_METHOD, "\n")
cat("[CONCORDANCE CONFIG] Methods:", paste(CONCORDANCE_METHODS, collapse = ", "), "\n")
cat("[CONCORDANCE CONFIG] Reference:", MASTER_REFERENCE, "\n")
if (length(CONCORDANCE_GENOMES) > 0) cat("[CONCORDANCE CONFIG] Genomes:", paste(CONCORDANCE_GENOMES, collapse = ", "), "\n")
cat("[CONCORDANCE CONFIG] Gene groups:", paste(CONCORDANCE_GENE_GROUPS, collapse = ", "), "\n")
cat("[CONCORDANCE CONFIG] Output:", OUTPUT_DIR, "\n")
