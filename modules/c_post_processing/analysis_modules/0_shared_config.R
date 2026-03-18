#!/usr/bin/env Rscript

# ===============================================
# SHARED CONFIGURATION FOR POST-PROCESSING ANALYSES
# ===============================================
# Centralized configuration to avoid duplication across scripts
# This file is sourced by all analysis modules
#
# ORGANIZATION:
#   1. Global Configuration Variables
#   2. Directory Constants
#   3. Analysis Parameters
#   4. Functions (GPU, Method Detection, Utilities)
#   5. Initialization (runs at load time)

# ===============================================
# SECTION 1: GLOBAL CONFIGURATION VARIABLES
# ===============================================
# Override via environment variables

# Method being processed (set by bash wrapper)
CURRENT_METHOD <- Sys.getenv("CURRENT_METHOD", unset = "M5_RSEM_Bowtie2")
if (!nzchar(Sys.getenv("CURRENT_METHOD", unset = ""))) {
  warning("CURRENT_METHOD env var is not set \u2014 defaulting to 'M5_RSEM_Bowtie2'. ",
          "Export CURRENT_METHOD before running analysis scripts.")
}

# Current dataset being processed (set by bash wrapper)
CURRENT_DATASET <- Sys.getenv("CURRENT_DATASET", unset = "")

# Master reference genome/transcriptome
MASTER_REFERENCE <- Sys.getenv("MASTER_REFERENCE", unset = "Eggplant_V4.1")

# GPU acceleration flag (isTRUE guards against NA from empty/malformed env var)
ENABLE_GPU <- isTRUE(as.logical(Sys.getenv("ENABLE_GPU", unset = "FALSE")))

# Thread count
THREADS <- as.integer(Sys.getenv("THREADS", unset = "8"))

# Available RAM (GB) - used to determine if memory-intensive operations are safe
AVAILABLE_RAM_GB <- as.integer(Sys.getenv("AVAILABLE_RAM_GB", unset = "24"))

# Available GPU VRAM (GB) - used for GPU memory management
GPU_VRAM_GB <- as.integer(Sys.getenv("GPU_VRAM_GB", unset = "8"))

# Global random seed for reproducibility across all stochastic operations
# (t-SNE, UMAP, GSEA permutations, DESeq2 shrinkage, WGCNA layout, etc.)
# Override via environment variable GLOBAL_RANDOM_SEED if needed.
GLOBAL_RANDOM_SEED <- as.integer(Sys.getenv("GLOBAL_RANDOM_SEED", unset = "42"))

# Memory-aware settings (with 24GB+ RAM: prioritize accuracy over memory conservation)
HIGH_MEMORY_MODE <- AVAILABLE_RAM_GB >= 16

# GPU status variables (initialized later by detect_gpu())
GPU_AVAILABLE <- FALSE
GPU_BACKEND <- "none"  # "none", "cuda", or "torch"

# ===============================================
# SECTION 2: DIRECTORY CONSTANTS
# ===============================================

# Base directories (relative to method folder)
# Different methods have different quantification output structures:
#   M1/M2 (HISAT2+StringTie): count_matrices_from_stringtie/
#   M3 (STAR): count_matrices_from_STAR/
#   M4 (Salmon): count_matrices_from_Salmon_Quant/
#   M5 (RSEM): count_matrices_from_RSEM_Quant/
# Use get_matrices_dir() for method-specific directory names
CONSOLIDATED_BASE_DIR <- "Figure_Outputs"

# Gene groups directory - use environment variable if set, otherwise derive from BASE_DIR
GENE_GROUPS_DIR <- Sys.getenv("GENE_GROUPS_DIR", unset = "")
if (GENE_GROUPS_DIR == "") {
  base_dir_fallback <- Sys.getenv("BASE_DIR", unset = "")
  if (nzchar(base_dir_fallback)) {
    GENE_GROUPS_DIR <- file.path(base_dir_fallback, "inputs", "gene_groups_csv")
  } else {
    # Last-resort: walk up three levels from ANALYSIS_MODULES_DIR to reach project root
    ANALYSIS_MODULES_DIR <- Sys.getenv("ANALYSIS_MODULES_DIR", unset = normalizePath(".", mustWork = FALSE))
    GENE_GROUPS_DIR <- file.path(dirname(dirname(dirname(ANALYSIS_MODULES_DIR))), "inputs", "gene_groups_csv")
  }
}

# SRR CSV directory for sample labels
SRR_CSV_DIR <- Sys.getenv("SRR_CSV_DIR", unset = "")
if (SRR_CSV_DIR == "") {
  base_dir_fallback <- Sys.getenv("BASE_DIR", unset = "")
  if (nzchar(base_dir_fallback)) {
    SRR_CSV_DIR <- file.path(base_dir_fallback, "inputs", "SRR_csv")
  } else {
    ANALYSIS_MODULES_DIR <- Sys.getenv("ANALYSIS_MODULES_DIR", unset = normalizePath(".", mustWork = FALSE))
    SRR_CSV_DIR <- file.path(dirname(dirname(dirname(ANALYSIS_MODULES_DIR))), "inputs", "SRR_csv")
  }
}

# Output subdirectories (include MASTER_REFERENCE as leaf for per-reference isolation)
OUTPUT_SUBDIRS <- list(
  MATRIX_CREATION = file.path("0_Matrix_Creation", MASTER_REFERENCE),
  BASIC_HEATMAP = file.path("I_Basic_Heatmap", MASTER_REFERENCE),
  CV_HEATMAP = file.path("II_Heatmap_with_CV", MASTER_REFERENCE),
  BAR_GRAPH = file.path("III_Bar_Graphs", MASTER_REFERENCE),
  WGCNA = file.path("IV_Coexpression_WGCNA", MASTER_REFERENCE),
  DEA = file.path("V_Differential_Expression", MASTER_REFERENCE),
  GSEA = file.path("VI_Gene_Set_Enrichment", MASTER_REFERENCE),
  DIM_REDUCTION = file.path("VII_PCA", MASTER_REFERENCE),
  CORRELATION = file.path("VIII_Sample_Clustering", MASTER_REFERENCE),
  TISSUE_SPEC = file.path("IX_Tissue_Specificity", MASTER_REFERENCE)
)

# Cache package availability once at load time (avoids repeated requireNamespace() per function call)
.HAS_DATATABLE <- requireNamespace("data.table", quietly = TRUE)

# ===============================================
# SECTION 3: ANALYSIS PARAMETERS
# ===============================================

# Count types - AUTO-SELECTED based on method:
#   - StringTie (M1/M2 HISAT2): tpm, fpkm, coverage
#   - Salmon (M4): tpm, NumReads (length-corrected, bias-corrected)
#   - RSEM (M5): tpm, expected_count (fpkm excluded from tximport pipeline)
# For DESeq2: use expected_count/NumReads (raw integer counts)
#   NOTE: StringTie coverage is a per-base abundance metric, NOT raw fragment counts.
#         M1 DESeq2 uses prepDE.py integer counts (gene_count_matrix.csv), not coverage.
# For visualization: TPM is preferred for cross-sample comparison
#
# NOTE: COUNT_TYPES is now dynamically set at initialization based on CURRENT_METHOD
#       Use get_count_types() function for method-specific types
COUNT_TYPES <- c()  # Initialized in SECTION 5 based on method

# Gene name display options (column names in gene_groups_csv/*.csv files)
GENE_TYPES <- c(
  #"Gene_ID",        # Column 1: Original gene identifier (e.g., SMEL4.1_01g005840)
  "Shortened_Name"   # Column 2: Human-readable short name (e.g., SmelDMP01.840)
)

# Sample label display options (column names in SRR_csv/*.csv files)
LABEL_TYPES <- c(
  #"SRR_ID",  # Column 1: SRA accession number (e.g., SRR3884686)
  "Organ"     # Column 2: Tissue/organ description (e.g., Flower_Buds)
)

# Processing levels
PROCESSING_LEVELS <- c(
  "gene_level",
  "isoform_level"
)

# Normalization schemes - AUTO-SELECTED based on count type:
#   - "raw": No transformation (only for raw counts: expected_count, NumReads, coverage)
#   - "count_type_normalized": Log2(counts + 1) transformation (universal)
#   - "deseq2_normalized": DESeq2-style median-of-ratios + log2 (only for raw counts)
#   - "zscore": Z-score normalization (after log2 transformation) (universal)
#   - "zscore_scaled_to_ten": Z-score scaled to 0-10 range (universal)
#   - "cpm": Counts Per Million (library-size normalized) + log2 (only for raw counts)
# NOTE: For DEA, use raw counts with DESeq2 internal normalization.
#       For visualization (heatmaps, PCA), use log2 or VST-transformed data.
#
# RECOMMENDATION for HISAT2/StringTie:
#   - Use TPM with "count_type_normalized" to see actual expression magnitudes
#   - Use "zscore" only when comparing relative patterns (hides absolute levels)
#   - Control genes will show variation with zscore even if expression is stable
#
# NORM_SCHEMES is now dynamically filtered per count_type using get_norm_schemes()
NORM_SCHEMES <- c(
  "count_type_normalized",  # Log2(TPM+1) - best for seeing actual expression levels
  "zscore",                 # Z-score after log2 (global) - relative pattern comparison
  #"zscore_row",            # Z-score per-gene (row-wise) - makes all genes equally visible (hides control stability)
  "zscore_scaled_to_ten"    # 0-10 scaled relative patterns - good for visual comparison
)  # All visualization schemes enabled; use get_norm_schemes() for count-type filtering

# Minimum sample/gene thresholds for various analyses
MIN_SAMPLES <- 3       # Minimum samples for correlation/clustering
MIN_GENES_DEA <- 10    # Minimum genes for differential expression
MIN_GENES_WGCNA <- 10  # Reduced to allow smaller gene groups (minimum recommended: 10)

# Sample labels (initialized later by load_sample_labels_from_csv())
SAMPLE_LABELS <- c()
SAMPLE_IDS <- c()

# ===============================================
# SECTION 4: FUNCTIONS
# ===============================================

# -----------------------------------------------
# 4.1 GPU Detection Functions
# -----------------------------------------------

# Get CUDA version - prefer conda CUDA over system
get_cuda_version <- function() {
  # First check for conda CUDA environment variable
  conda_cuda <- Sys.getenv("TORCH_CUDA_VERSION", unset = "")
  if (nzchar(conda_cuda)) {
    return(list(version = conda_cuda, source = "conda"))
  }
  
  # Check nvcc first (more accurate for installed toolkit)
  nvcc_version <- tryCatch({
    output <- system("nvcc --version 2>/dev/null", intern = TRUE)
    ver_line <- grep("release", output, value = TRUE)
    if (length(ver_line) > 0) {
      version <- gsub(".*release ([0-9]+\\.[0-9]+).*", "\\1", ver_line[1])
      return(list(version = version, source = "nvcc"))
    }
    NULL
  }, error = function(e) NULL)
  
  if (!is.null(nvcc_version)) return(nvcc_version)
  
  # Fallback: detect from nvidia-smi (driver's max supported version)
  tryCatch({
    cuda_output <- system("nvidia-smi 2>/dev/null", intern = TRUE)
    cuda_line <- grep("CUDA Version", cuda_output, value = TRUE)
    if (length(cuda_line) > 0) {
      version <- gsub(".*CUDA Version: ([0-9]+\\.[0-9]+).*", "\\1", cuda_line[1])
      return(list(version = version, source = "nvidia-smi"))
    }
    return(list(version = NULL, source = "none"))
  }, error = function(e) list(version = NULL, source = "none"))
}

# Detect GPU and compatible packages
# Note: R torch supports CUDA 11.7-11.8 and 12.1
# For other CUDA versions, GPU ops fall back to CPU
detect_gpu <- function() {
  if (!ENABLE_GPU) {
    return(list(available = FALSE, backend = "cpu", message = "GPU disabled by configuration"))
  }
  
  # Check for GPU hardware via nvidia-smi
  gpu_available <- tryCatch({
    result <- system("nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null", intern = TRUE)
    length(result) > 0 && !grepl("error|fail", result[1], ignore.case = TRUE)
  }, error = function(e) FALSE, warning = function(w) FALSE)
  
  if (!gpu_available) {
    return(list(available = FALSE, backend = "cpu", message = "No NVIDIA GPU detected"))
  }
  
  # Check CUDA version
  cuda_info <- get_cuda_version()
  cuda_version <- if (!is.null(cuda_info$version)) suppressWarnings(as.numeric(cuda_info$version)) else 0
  if (is.na(cuda_version)) cuda_version <- 0  # Guard against unparseable version strings (e.g., "12.1.1")

  # R torch supports CUDA 11.7, 11.8, and 12.1
  torch_supported <- (cuda_version >= 11.7 && cuda_version <= 11.8) || (cuda_version >= 12.1 && cuda_version < 12.2)

  if (!is.null(cuda_info$version)) {
    message("[GPU] CUDA ", cuda_info$version, " detected (", cuda_info$source, ")")
    if (!torch_supported) {
      message("[GPU] Note: R torch requires CUDA 11.7, 11.8, or 12.1, found ", cuda_info$version)
      message("[GPU] GPU matrix operations will use CPU. This does not affect pipeline results.")
    }
  }
  
  # Check for torch (preferred for matrix operations)
  if (requireNamespace("torch", quietly = TRUE)) {
    cuda_works <- tryCatch({
      torch::cuda_is_available()
    }, error = function(e) FALSE)
    
    if (cuda_works) {
      return(list(available = TRUE, backend = "torch", 
                  message = paste0("GPU available via torch (", torch::cuda_device_count(), " device(s))")))
    } else if (torch_supported) {
      message("[GPU] torch package found but CUDA backend not installed")
      message("[GPU] Run: Rscript -e 'torch::install_torch(type=\"cuda\")'")
    }
  }
  
  # Check for gpuR (alternative for matrix operations)
  if (requireNamespace("gpuR", quietly = TRUE)) {
    tryCatch({
      gpuR::detectGPUs()
      return(list(available = TRUE, backend = "gpuR", message = "GPU available via gpuR"))
    }, error = function(e) NULL)
  }
  
  return(list(available = FALSE, backend = "cpu", 
              message = "GPU detected but R GPU acceleration not available (using CPU - this is fine)"))
}

# -----------------------------------------------
# 4.2 Method Type Detection Functions
# -----------------------------------------------

# Detect quantification method type from CURRENT_METHOD
get_method_type <- function(method = CURRENT_METHOD) {
  if (grepl("HISAT2|StringTie|M1_|M2_|^M1$|^M2$", method, ignore.case = TRUE)) {
    return("stringtie")
  } else if (grepl("Salmon|M4_|^M4$", method, ignore.case = TRUE)) {
    return("salmon")
  } else if (grepl("RSEM|M5_|^M5$", method, ignore.case = TRUE)) {
    return("rsem")
  } else if (grepl("STAR|M3_|^M3$", method, ignore.case = TRUE)) {
    return("star")
  }
  return("unknown")
}

# Get method-specific quant directory
get_quant_dir <- function(method = CURRENT_METHOD) {
  method_type <- get_method_type(method)
  switch(method_type,
    "stringtie" = "stringtie_WD",
    "salmon" = "Salmon_Quant",
    "rsem" = "RSEM_Quant_WD",
    "star" = "STAR_alignment_WD",
    "quant_WD"  # default
  )
}

# Get method-appropriate count types
# For visualization: TPM is preferred (comparable across samples)
# For DESeq2: use raw counts (coverage for StringTie)
get_count_types <- function(method = CURRENT_METHOD) {
  method_type <- get_method_type(method)
  switch(method_type,
    # StringTie: TPM for cross-sample comparison, FPKM similar, coverage as abundance metric
    # NOTE: coverage is per-base depth (not raw counts); DESeq2 uses prepDE.py integers instead
    "stringtie" = c("tpm", "fpkm", "coverage"),
    # Salmon: TPM for visualization, NumReads for DESeq2
    "salmon" = c("tpm", "NumReads"),
    # RSEM: TPM for visualization, expected_count for DESeq2
    # NOTE: RSEM also produces FPKM (in .genes.results column 7), but tximport
    # only outputs counts + TPM. The bash script saves FPKM to genes.FPKM.not_cross_norm
    # for reference, but the tximport-based analysis pipeline does not use it.
    "rsem" = c("tpm", "expected_count"),
    # STAR+Salmon: TPM for visualization, NumReads (Salmon offset counts) for DESeq2
    "star" = c("tpm", "NumReads"),
    c("tpm")  # default
  )
}

# Get appropriate normalization schemes for a given count type
# Raw count types (expected_count, NumReads) support all normalizations including CPM/DESeq2
# Abundance metrics (coverage) from StringTie are NOT raw integer counts — they are per-base
#   coverage values already derived from alignment depth. CPM and DESeq2 median-of-ratios
#   assume raw fragment/read counts as input, so they are invalid for coverage.
# Pre-normalized types (tpm, fpkm) should NOT use raw/deseq2/cpm (already normalized)
get_norm_schemes <- function(count_type) {
  # True raw integer count types (from counting reads/fragments)
  raw_count_types <- c("expected_count", "numreads", "counts")
  # StringTie coverage: abundance-based metric, not raw counts
  abundance_count_types <- c("coverage")
  is_raw <- tolower(count_type) %in% raw_count_types
  is_abundance <- tolower(count_type) %in% abundance_count_types

  if (is_raw) {
    # Raw integer counts support all normalization schemes
    return(c(
      "count_type_normalized",  # Log2(x+1) - always useful
      "raw",                    # No transformation - for DESeq2 input
      "deseq2_normalized",      # Median-of-ratios + log2
      "zscore",                 # Z-score after log2
      "zscore_row",             # Z-score per-gene (row-wise)
      "zscore_scaled_to_ten",   # Z-score scaled 0-10
      "cpm"                     # Counts Per Million + log2
    ))
  } else if (is_abundance) {
    # StringTie coverage: per-base abundance metric (not raw fragment counts)
    # Log2 and z-score transformations are valid; CPM/DESeq2 are not
    # (those methods assume raw read/fragment counts for library-size correction)
    return(c(
      "count_type_normalized",  # Log2(x+1) - useful for visualization
      "raw",                    # No transformation - for inspection
      "zscore",                 # Z-score after log2 (global)
      "zscore_row",             # Z-score per-gene (row-wise)
      "zscore_scaled_to_ten"    # Z-score scaled 0-10
    ))
  } else {
    # Pre-normalized types (TPM, FPKM) - skip redundant normalizations
    # TPM/FPKM are already library-size normalized, so raw/deseq2/cpm don't make sense
    return(c(
      "count_type_normalized",  # Log2(x+1) - still useful for visualization
      "zscore",                 # Z-score after log2 (global) - useful for heatmaps
      "zscore_row",             # Z-score per-gene (row-wise) - makes all genes equally visible
      "zscore_scaled_to_ten"    # Z-score scaled 0-10 - useful for heatmaps
    ))
  }
}

# Get the raw integer count type appropriate for DESeq2 input
# DESeq2 requires raw counts, NOT pre-normalized TPM/FPKM.
# Returns NULL for methods where DEA uses a separate input path (e.g., M1 prepDE.py).
get_raw_count_type <- function(method = CURRENT_METHOD) {
  method_type <- get_method_type(method)
  switch(method_type,
    "salmon" = "NumReads",
    "rsem" = "expected_count",
    "star" = "NumReads",
    "stringtie" = NULL,  # M1 uses prepDE.py counts (handled separately); M2 lacks DEA-suitable counts
    NULL  # default
  )
}

# Check if a count_type + norm_scheme combination is valid
is_valid_norm_for_count <- function(count_type, norm_scheme) {
  valid_schemes <- get_norm_schemes(count_type)
  return(tolower(norm_scheme) %in% tolower(valid_schemes))
}

# Get method-specific matrices directory
get_matrices_dir <- function(method = CURRENT_METHOD) {
  method_type <- get_method_type(method)
  switch(method_type,
    "stringtie" = "count_matrices_from_stringtie",
    "salmon" = "count_matrices_from_Salmon_Quant",
    "rsem" = "count_matrices_from_RSEM_Quant",
    "star" = "count_matrices_from_STAR",
    "count_matrices"  # default
  )
}

# Generate combined output folder name: GeneGroup_in_Dataset
# Example: SmelDMPs_with_1_18s_rRNA_in_PRJNA328564_selected
get_output_folder_name <- function(gene_group, dataset = CURRENT_DATASET) {
  if (nzchar(dataset)) {
    return(paste0(gene_group, "_in_", dataset))
  } else {
    return(gene_group)  # fallback to gene group only if no dataset specified
  }
}

# -----------------------------------------------
# 4.3 Sample Labels Functions
# -----------------------------------------------

# Load sample labels from all CSV files in SRR_csv directory
# CSV format: SRR_ID,Organ,Notes
# If SRR_COMBINED_LIST_STR env var is set, filter to only those samples AND preserve order
load_sample_labels_from_csv <- function(srr_csv_dir = SRR_CSV_DIR) {
  labels <- c()
  if (!dir.exists(srr_csv_dir)) return(labels)
  
  csv_files <- list.files(srr_csv_dir, pattern = "\\.csv$", full.names = TRUE)
  for (csv_file in csv_files) {
    tryCatch({
      df <- if (.HAS_DATATABLE) {
        data.table::fread(csv_file, header = TRUE, data.table = FALSE)
      } else {
        read.csv(csv_file, stringsAsFactors = FALSE, header = TRUE, comment.char = "#")
      }
      if ("SRR_ID" %in% colnames(df) && "Organ" %in% colnames(df)) {
        df <- df[!is.na(df$SRR_ID) & nzchar(trimws(df$SRR_ID)), ]
        new_labels <- setNames(df$Organ, df$SRR_ID)
        labels <- c(labels, new_labels[!names(new_labels) %in% names(labels)])
      }
    }, error = function(e) NULL)
  }
  
  # Filter to only samples specified in SRR_COMBINED_LIST_STR (from bash config)
  # Format may be "SRR123 SRR456" or "SRR123:Organ1 SRR456:Organ2"
  # IMPORTANT: Preserve the order from the env var (which follows CSV file order)
  srr_list_str <- Sys.getenv("SRR_COMBINED_LIST_STR", unset = "")
  if (nzchar(srr_list_str)) {
    enabled_samples <- trimws(strsplit(srr_list_str, " ")[[1]])
    # Extract just the SRR ID (before colon if present)
    enabled_samples <- sapply(strsplit(enabled_samples, ":"), `[`, 1)
    original_count <- length(labels)
    # Reorder labels to match the order in SRR_COMBINED_LIST_STR (CSV order)
    labels <- labels[enabled_samples[enabled_samples %in% names(labels)]]
    cat("[CONFIG] Filtering samples: ", original_count, " -> ", length(labels), 
        " (from SRR_COMBINED_LIST_STR)\n", sep = "")
  }
  
  return(labels)
}

# -----------------------------------------------
# 4.4 Runtime Configuration Functions
# -----------------------------------------------

read_config_file <- function(file_path, default_value, is_boolean = FALSE) {
  if (!file.exists(file_path)) return(default_value)
  value <- trimws(readLines(file_path, warn = FALSE))
  if (length(value) == 0 || !nzchar(value[1])) return(default_value)
  if (is_boolean) {
    return(tolower(value[1]) == "true")
  } else {
    return(value[nzchar(value)])
  }
}

load_runtime_config <- function(method_modules_dir = ".") {
  gene_groups <- read_config_file(
    file.path(method_modules_dir, ".gene_groups_temp.txt"),
    default_value = c("SmelDMPs", "SmelGRF-GIF")
  )
  
  master_reference <- read_config_file(
    file.path(method_modules_dir, ".master_reference_temp.txt"),
    default_value = MASTER_REFERENCE
  )
  if (length(master_reference) > 1) master_reference <- master_reference[1]
  
  # Check environment variable first (set by bash wrapper), then fall back to temp file
  overwrite_env <- Sys.getenv("OVERWRITE_EXISTING", unset = "")
  if (nzchar(overwrite_env)) {
    overwrite <- tolower(overwrite_env) == "true"
  } else {
    overwrite <- read_config_file(
      file.path(method_modules_dir, ".overwrite_temp.txt"),
      default_value = TRUE,
      is_boolean = TRUE
    )
  }
  
  list(
    gene_groups = gene_groups,
    master_reference = master_reference,
    overwrite_existing = overwrite
  )
}

# -----------------------------------------------
# 4.5 Output and Printing Functions
# -----------------------------------------------

print_separator <- function(char = "=", width = 60) {
  cat("\n", paste(rep(char, width), collapse = ""), "\n")
}

print_config_summary <- function(title, config) {
  print_separator()
  cat(title, "\n")
  print_separator()
  cat("\nConfiguration:\n")
  cat("  * Method:", CURRENT_METHOD, "(", get_method_type(CURRENT_METHOD), ")\n")
  cat("  * Master Reference:", config$master_reference, "\n")
  cat("  * Overwrite existing:", config$overwrite_existing, "\n")
  cat("  * Gene groups:", paste(config$gene_groups, collapse = ", "), "\n")
  cat("  * Count types:", paste(COUNT_TYPES, collapse = ", "), "\n")
  cat("  * Threads:", THREADS, "\n")
  cat("  * Random seed:", GLOBAL_RANDOM_SEED, "\n")
  if (GPU_AVAILABLE) {
    cat("  * GPU:", "ENABLED (", GPU_BACKEND, ")\n", sep = "")
  } else if (ENABLE_GPU) {
    cat("  * GPU: REQUESTED but unavailable (using CPU)\n")
  } else {
    cat("  * GPU: disabled\n")
  }
  cat("\n")
}

print_summary <- function(successful, total, skipped = NULL) {
  print_separator()
  cat("SUMMARY:", successful, "/", total, "items generated")
  if (!is.null(skipped) && skipped > 0) {
    cat(" (", skipped, " skipped)\n", sep = "")
  } else {
    cat("\n")
  }
  print_separator()
}

# -----------------------------------------------
# 4.6 File and Path Utility Functions
# -----------------------------------------------

ensure_output_dir <- function(dir_path, clean = FALSE) {
  if (clean && dir.exists(dir_path)) {
    unlink(dir_path, recursive = TRUE)
  }
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
}

build_input_path <- function(gene_group, processing_level, count_type, gene_type, 
                             matrices_dir = NULL, master_ref = MASTER_REFERENCE,
                             method = CURRENT_METHOD, label_type = "Organ",
                             dataset = CURRENT_DATASET) {
  if (is.null(matrices_dir)) {
    matrices_dir <- get_matrices_dir(method)
  }
  
  method_type <- get_method_type(method)
  
  # Get combined folder name (GeneGroup_in_Dataset) if dataset is specified
  folder_name <- get_output_folder_name(gene_group, dataset)
  
  # StringTie (M1/M2) has different file naming and structure
  if (method_type == "stringtie") {
    # StringTie matrix builder only creates "geneName" files (raw Gene_IDs as row names);
    # name conversion to Shortened_Name is applied at runtime by apply_labels().
    # Always use "geneName" regardless of gene_type to match actual file names.
    stringtie_gene_type <- "geneName"
    stringtie_label_type <- if (label_type == "Organ") "Organ" else "SRR"
    # File name uses folder_name (includes dataset suffix)
    file.path(matrices_dir, folder_name,
              paste0(folder_name, "_", count_type, "_counts_", stringtie_gene_type, 
                     "_", stringtie_label_type, "_from_", master_ref, ".tsv"))
  } else {
    # Tximport path (Salmon/RSEM): count_matrices_from_*/{master_ref}/{level}/{folder_name}/
    # Both full-genome (gene_group == master_ref) and gene group subsets use the same structure;
    # folder_name = get_output_folder_name(gene_group, dataset) = "{gene_group}_in_{dataset}"
    file.path(matrices_dir, master_ref, processing_level, folder_name,
              paste0(folder_name, "_", count_type, "_", gene_type,
                     "_from_", master_ref, "_", processing_level, ".tsv"))
  }
}

build_title_base <- function(gene_group, count_type, gene_type, label_type, 
                             processing_level, norm_scheme, master_ref = MASTER_REFERENCE,
                             dataset = CURRENT_DATASET) {
  # Include dataset in title if specified
  base_name <- get_output_folder_name(gene_group, dataset)
  paste0(base_name, "_", count_type, "_", gene_type, "_",
         label_type, "_from_", master_ref, "_", processing_level, "_", norm_scheme)
}

validate_and_read_matrix <- function(input_file, min_rows = 2) {
  if (!file.exists(input_file)) {
    return(list(success = FALSE, reason = "file not found", data = NULL, n_genes = 0L))
  }
  matrix_data <- read_count_matrix(input_file)
  if (is.null(matrix_data)) {
    return(list(success = FALSE, reason = "failed to read", data = NULL, n_genes = 0L))
  }
  if (nrow(matrix_data) < min_rows) {
    return(list(success = FALSE, reason = paste0("need >=", min_rows, " rows"), data = NULL, n_genes = 0L))
  }
  list(success = TRUE, data = matrix_data, n_genes = nrow(matrix_data))
}

should_skip_existing <- function(output_path, overwrite) {
  !overwrite && file.exists(output_path)
}

# ===============================================
# SECTION 5: INITIALIZATION (runs at load time)
# ===============================================

# Initialize GPU detection — skip the (slow) nvidia-smi / nvcc probes when
# the user has explicitly disabled GPU support, saving ~0.5-1s per script load.
if (ENABLE_GPU) {
  .gpu_info <- detect_gpu()
  GPU_AVAILABLE <- .gpu_info$available
  GPU_BACKEND <- .gpu_info$backend
  if (GPU_AVAILABLE || interactive()) {
    message("[GPU] ", .gpu_info$message)
  }
} else {
  GPU_AVAILABLE <- FALSE
  GPU_BACKEND <- "none"
}

# Initialize COUNT_TYPES based on current method
COUNT_TYPES <- get_count_types(CURRENT_METHOD)
message("[CONFIG] Method: ", CURRENT_METHOD, " -> Count types: ", paste(COUNT_TYPES, collapse = ", "))

# Initialize sample labels from CSV files
SAMPLE_LABELS <- load_sample_labels_from_csv()
SAMPLE_IDS <- names(SAMPLE_LABELS)

if (length(SAMPLE_LABELS) == 0) {
  warning("No sample labels loaded from CSV files in: ", SRR_CSV_DIR,
          "\nEnsure CSV files have SRR_ID and Organ columns.")
}
