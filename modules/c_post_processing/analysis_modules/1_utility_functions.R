#!/usr/bin/env Rscript

# ===============================================
# UTILITY FUNCTIONS FOR POST-PROCESSING ANALYSES
# ===============================================
# Comprehensive utility functions for all analysis modules

# Guard: only source 0_shared_config.R if it hasn't been loaded yet (CURRENT_METHOD is its sentinel).
# Use ANALYSIS_MODULES_DIR env var (set by the bash wrapper) with a safe fallback to
# sys.frame(1)$ofile for interactive/direct-source calls.
if (!exists("CURRENT_METHOD")) {
  .utils_dir <- tryCatch(
    dirname(sys.frame(1)$ofile),
    error = function(e) Sys.getenv("ANALYSIS_MODULES_DIR", ".")
  )
  source(file.path(.utils_dir, "0_shared_config.R"))
}

# ===============================================
# GPU-ACCELERATED COMPUTATION FUNCTIONS
# ===============================================
# These functions use GPU when available, with automatic CPU fallback
# GPU VRAM limit: GPU_VRAM_GB (default 8GB) - matrices larger than this fall back to CPU

# Estimate matrix memory usage in GB (for float32)
estimate_matrix_memory_gb <- function(nrow, ncol, dtype_bytes = 4) {
  (nrow * ncol * dtype_bytes) / (1024^3)
}

# Check if matrix fits in GPU VRAM (with safety margin for intermediate results)
fits_in_gpu_vram <- function(nrow, ncol, safety_factor = 0.7) {
  # For operations like correlation, we need space for input + output + intermediates
  # Correlation of NxM matrix produces NxN output, plus intermediate centered matrix
  estimated_usage <- estimate_matrix_memory_gb(nrow, ncol) * 3  # 3x for safety
  max_allowed <- GPU_VRAM_GB * safety_factor
  return(estimated_usage <= max_allowed)
}

# GPU-accelerated correlation matrix computation
# WGCNA cor() is the bottleneck for large datasets
gpu_cor <- function(x, method = "pearson") {
  # GPU path only supports Pearson; fall back to CPU for Spearman/Kendall
  if (!GPU_AVAILABLE || nrow(x) < 100 || !fits_in_gpu_vram(nrow(x), ncol(x)) || method != "pearson") {
    if (GPU_AVAILABLE && method != "pearson") {
      message("[GPU] ", method, " correlation not supported on GPU, using CPU")
    } else if (GPU_AVAILABLE && !fits_in_gpu_vram(nrow(x), ncol(x))) {
      message("[GPU] Matrix too large for ", GPU_VRAM_GB, "GB VRAM, using CPU")
    }
    # Use CPU for small matrices, large matrices, non-Pearson, or when GPU unavailable
    return(cor(x, method = method, use = "pairwise.complete.obs"))
  }

  tryCatch({
    if (GPU_BACKEND == "torch") {
      # Use torch for GPU column-column correlation (matching R's cor())
      x_tensor <- torch::torch_tensor(as.matrix(x), device = "cuda")
      n <- x_tensor$size(0)
      if (n < 2) stop("Need at least 2 rows for correlation")
      # Center columns: mean(dim=0) averages across rows for each column
      x_centered <- x_tensor - x_tensor$mean(dim = 0, keepdim = TRUE)
      # Column-column covariance: t(X_c) %*% X_c / (nrow - 1)
      cov_matrix <- torch::torch_mm(x_centered$t(), x_centered) / (n - 1)
      # Compute standard deviations of columns (clamp to avoid division by zero for constant columns)
      std_dev <- torch::torch_sqrt(torch::torch_diag(cov_matrix))
      std_dev <- torch::torch_clamp(std_dev, min = 1e-12)
      # Compute correlation
      cor_matrix <- cov_matrix / torch::torch_outer(std_dev, std_dev)
      result <- as.matrix(cor_matrix$cpu())
      rownames(result) <- colnames(x)
      colnames(result) <- colnames(x)
      return(result)
    } else if (GPU_BACKEND == "gpuR") {
      # Use gpuR for GPU correlation (cov() computes column-column covariance)
      x_gpu <- gpuR::vclMatrix(as.matrix(x), type = "float")
      result <- as.matrix(gpuR::cov(x_gpu))
      # Convert covariance to correlation
      std_dev <- sqrt(diag(result))
      result <- result / outer(std_dev, std_dev)
      rownames(result) <- colnames(x)
      colnames(result) <- colnames(x)
      return(result)
    }
  }, error = function(e) {
    message("[GPU] Correlation failed, falling back to CPU: ", e$message)
  })
  
  # Fallback to CPU
  return(cor(x, method = method, use = "pairwise.complete.obs"))
}

# GPU-accelerated PCA
gpu_prcomp <- function(x, center = TRUE, scale. = FALSE, rank. = NULL) {
  # Check if matrix fits in GPU VRAM
  if (!GPU_AVAILABLE || nrow(x) < 50 || !fits_in_gpu_vram(nrow(x), ncol(x))) {
    return(prcomp(x, center = center, scale. = scale., rank. = rank.))
  }
  
  tryCatch({
    if (GPU_BACKEND == "torch") {
      x_mat <- as.matrix(x)
      
      # Center and scale on CPU first (small operation)
      if (center) {
        col_means <- colMeans(x_mat, na.rm = TRUE)
        x_mat <- sweep(x_mat, 2, col_means, "-")
      }
      if (scale.) {
        col_sds <- if (requireNamespace("matrixStats", quietly = TRUE)) {
          matrixStats::colSds(x_mat, na.rm = TRUE)
        } else {
          apply(x_mat, 2, sd, na.rm = TRUE)
        }
        col_sds[col_sds == 0] <- 1
        x_mat <- sweep(x_mat, 2, col_sds, "/")
      }
      
      # SVD on GPU
      x_tensor <- torch::torch_tensor(x_mat, device = "cuda", dtype = torch::torch_float32())
      svd_result <- torch::linalg_svd(x_tensor, full_matrices = FALSE)
      
      # Extract results
      n <- nrow(x_mat)
      sdev <- as.numeric(svd_result[[2]]$cpu()) / sqrt(max(1, n - 1))
      rotation <- as.matrix(svd_result[[3]]$t()$cpu())
      # Scores = X_centered %*% V; svd_result[[3]] = Vh (k×p), so Vh$t() = V (p×k)
      x_scores <- as.matrix(torch::torch_mm(x_tensor, svd_result[[3]]$t())$cpu())
      
      # Build prcomp-compatible result
      result <- list(
        sdev = sdev,
        rotation = rotation,
        x = x_scores,
        center = if (center) col_means else FALSE,
        scale = if (scale.) col_sds else FALSE
      )
      class(result) <- "prcomp"
      rownames(result$x) <- rownames(x)
      colnames(result$x) <- paste0("PC", seq_len(ncol(result$x)))
      colnames(result$rotation) <- paste0("PC", seq_len(ncol(result$rotation)))
      rownames(result$rotation) <- colnames(x)
      
      return(result)
    }
  }, error = function(e) {
    message("[GPU] PCA failed, falling back to CPU: ", e$message)
  })
  
  # Fallback to CPU
  return(prcomp(x, center = center, scale. = scale., rank. = rank.))
}

# GPU-accelerated matrix multiplication (for large TOM calculations in WGCNA)
gpu_matmult <- function(A, B) {
  # Check if matrices fit in GPU VRAM
  if (!GPU_AVAILABLE || nrow(A) < 100 || !fits_in_gpu_vram(nrow(A), ncol(B))) {
    return(A %*% B)
  }
  
  tryCatch({
    if (GPU_BACKEND == "torch") {
      A_tensor <- torch::torch_tensor(as.matrix(A), device = "cuda", dtype = torch::torch_float32())
      B_tensor <- torch::torch_tensor(as.matrix(B), device = "cuda", dtype = torch::torch_float32())
      result <- as.matrix(torch::torch_mm(A_tensor, B_tensor)$cpu())
      rownames(result) <- rownames(A)
      colnames(result) <- colnames(B)
      return(result)
    } else if (GPU_BACKEND == "gpuR") {
      A_gpu <- gpuR::vclMatrix(as.matrix(A), type = "float")
      B_gpu <- gpuR::vclMatrix(as.matrix(B), type = "float")
      result <- as.matrix(A_gpu %*% B_gpu)
      rownames(result) <- rownames(A)
      colnames(result) <- colnames(B)
      return(result)
    }
  }, error = function(e) {
    message("[GPU] Matrix multiplication failed, falling back to CPU: ", e$message)
  })
  
  return(A %*% B)
}

# GPU-accelerated Euclidean distance matrix (for clustering)
gpu_dist <- function(x, method = "euclidean") {
  # Only GPU-accelerate Euclidean distance; fallback for other methods
  if (!GPU_AVAILABLE || method != "euclidean" || nrow(x) < 50 || !fits_in_gpu_vram(nrow(x), ncol(x))) {
    return(dist(x, method = method))
  }
  
  tryCatch({
    if (GPU_BACKEND == "torch") {
      x_tensor <- torch::torch_tensor(as.matrix(x), device = "cuda", dtype = torch::torch_float32())
      # Compute pairwise squared distances: ||a-b||^2 = ||a||^2 + ||b||^2 - 2*a.b
      sq_norms <- torch::torch_sum(x_tensor^2, dim = 2, keepdim = TRUE)
      distances_sq <- sq_norms + sq_norms$t() - 2 * torch::torch_mm(x_tensor, x_tensor$t())
      # Clamp to avoid negative values from numerical precision issues
      distances_sq <- torch::torch_clamp(distances_sq, min = 0)
      distances <- torch::torch_sqrt(distances_sq)
      result <- as.matrix(distances$cpu())
      rownames(result) <- rownames(x)
      colnames(result) <- rownames(x)
      return(as.dist(result))
    }
  }, error = function(e) {
    message("[GPU] Distance calculation failed, falling back to CPU: ", e$message)
  })
  
  return(dist(x, method = method))
}

# ===============================================
# FILE I/O FUNCTIONS
# ===============================================

read_count_matrix <- function(file_path) {
  tryCatch({
    # .rds cache: readRDS is ~10-50x faster than TSV parsing for repeated reads.
    # Each analysis module (heatmap, PCA, DEA, etc.) re-reads the same matrix files;
    # the .rds cache eliminates redundant parsing after the first read.
    rds_path <- paste0(file_path, ".rds")
    if (file.exists(rds_path) &&
        file.mtime(rds_path) >= file.mtime(file_path)) {
      return(readRDS(rds_path))
    }

    # data.table::fread() is 10-50x faster than read.table() for large matrices
    data <- if (.HAS_DATATABLE) {
      data.table::fread(file_path, header = TRUE, sep = "\t",
                        na.strings = c("", "NA", "null"),
                        data.table = FALSE)
    } else {
      read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                 check.names = FALSE, na.strings = c("", "NA", "null"))
    }
    # Normalize duplicate column names: fread uses "_1" suffixes but downstream
    # code (apply_labels, replicate averaging) expects ".1" from make.unique()
    if (any(duplicated(colnames(data)))) {
      colnames(data) <- make.unique(colnames(data))
    }
    if (any(duplicated(data[, 1]))) {
      data[, 1] <- make.unique(as.character(data[, 1]), sep = "_")
    }
    rownames(data) <- data[, 1]
    data <- data[, -1, drop = FALSE]
    # Vectorized type coercion: identify non-numeric columns in one pass, convert in bulk
    non_num <- which(!vapply(data, is.numeric, logical(1)))
    if (length(non_num) > 0) {
      data[non_num] <- lapply(data[non_num], function(x) suppressWarnings(as.numeric(as.character(x))))
    }
    data_matrix <- as.matrix(data)
    data_matrix[is.na(data_matrix)] <- 0

    # Save .rds cache for subsequent reads by other analysis modules
    tryCatch(saveRDS(data_matrix, rds_path), error = function(e) NULL)

    return(data_matrix)
  }, error = function(e) {
    cat("Error reading file:", file_path, "-", e$message, "\n")
    return(NULL)
  })
}

save_matrix_data <- function(data_matrix, output_path, metadata = NULL) {
  tryCatch({
    base_path <- tools::file_path_sans_ext(output_path)
    matrix_df <- data.frame(
      Gene_ID = rownames(data_matrix),
      data_matrix,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    if (.HAS_DATATABLE) {
      data.table::fwrite(matrix_df, file = paste0(base_path, ".tsv"),
                         sep = "\t", quote = FALSE)
    } else {
      write.table(matrix_df, file = paste0(base_path, ".tsv"),
                  sep = "\t", quote = FALSE, row.names = FALSE)
    }
    return(TRUE)
  }, error = function(e) {
    cat("Error saving matrix:", e$message, "\n")
    return(FALSE)
  })
}

# Read gene list from file (supports CSV and plain text formats)
# CSV format expected: Gene_ID,Shortened_Name,... (uses Gene_ID column)
# Plain text format: one gene ID per line
read_gene_list_from_file <- function(gene_list_file) {
  if (!file.exists(gene_list_file)) return(NULL)
  
  file_ext <- tolower(tools::file_ext(gene_list_file))
  
  tryCatch({
    if (file_ext == "csv") {
      gene_df <- read.csv(gene_list_file, stringsAsFactors = FALSE, header = TRUE)
      if ("Gene_ID" %in% colnames(gene_df)) {
        gene_list <- gene_df$Gene_ID
      } else {
        gene_list <- gene_df[[1]]  # Fallback to first column
      }
    } else {
      gene_list <- suppressWarnings(readLines(gene_list_file))
      gene_list <- gene_list[!grepl("^#|^Gene_ID", gene_list, ignore.case = TRUE) & nzchar(gene_list)]
    }
    return(trimws(gene_list))
  }, error = function(e) {
    cat("Error reading gene list:", e$message, "\n")
    return(NULL)
  })
}

truncate_labels <- function(labels, max_length = 25) {
  labels <- as.character(labels)
  labels[is.na(labels)] <- ""
  too_long <- nchar(labels) > max_length
  labels[too_long] <- paste0(substr(labels[too_long], 1, max_length - 3), "...")
  labels
}

# ===============================================
# GENE ID MATCHING
# ===============================================
# Shared logic for matching a gene list to matrix row names, handling
# version suffixes common in eggplant IDs (e.g., SMEL4.1_06g023900.1.01).
# Used by: filter_by_gene_group(), tximport_salmon_to_matrices.R,
#           tximport_star_to_matrices.R

match_gene_ids <- function(gene_list, data_rownames) {
  # Precompute base IDs by stripping version suffixes (.X.XX then .X)
  base_ids <- sub("\\.[0-9]+\\.[0-9]+$", "", data_rownames)
  base_ids <- sub("\\.[0-9]+$", "", base_ids)

  # Build lookup: base_id -> data_rownames indices (vectorized, O(n) via split)
  # Use environment as hash map for O(1) lookups
  base_to_rows <- new.env(hash = TRUE, parent = emptyenv(), size = length(data_rownames))
  # split() + seq_along avoids O(n²) c() concatenation that occurs with incremental appends
  idx_groups <- split(seq_along(data_rownames), base_ids)
  for (nm in names(idx_groups)) {
    base_to_rows[[nm]] <- idx_groups[[nm]]
  }
  # Resolve indices to row names (deferred to avoid repeated string concatenation)
  .resolve_rows <- function(key) {
    idx <- base_to_rows[[key]]
    if (is.null(idx)) return(NULL)
    data_rownames[idx]
  }
  rowname_set <- new.env(hash = TRUE, parent = emptyenv(), size = length(data_rownames))
  for (rn in data_rownames) rowname_set[[rn]] <- TRUE

  # Vectorized: exact matches first
  exact_mask <- gene_list %in% data_rownames
  matched_list <- list(gene_list[exact_mask])

  # Non-exact: try forward match (gene_list ID is base -> find suffixed data rows)
  # Pre-allocate list to avoid O(n²) c() concatenation
  non_exact <- gene_list[!exact_mask]
  if (length(non_exact) > 0) {
    ne_results <- vector("list", length(non_exact))
    for (i in seq_along(non_exact)) {
      gene <- non_exact[i]
      hits <- .resolve_rows(gene)
      if (!is.null(hits)) {
        ne_results[[i]] <- hits
      } else {
        # Reverse: strip suffix from gene_list ID to match base-level row IDs
        gene_base <- sub("\\.[0-9]+$", "", gene)
        if (gene_base != gene) {
          hits2 <- .resolve_rows(gene_base)
          if (!is.null(hits2)) {
            ne_results[[i]] <- hits2
          } else if (!is.null(rowname_set[[gene_base]])) {
            ne_results[[i]] <- gene_base
          }
        }
      }
    }
    matched_list <- c(matched_list, ne_results)
  }
  unique(unlist(matched_list, use.names = FALSE))
}

# ===============================================
# NORMALIZATION FUNCTIONS
# ===============================================

preprocess_for_raw <- function(data_matrix) {
  if (is.null(data_matrix) || nrow(data_matrix) == 0) return(NULL)
  data_matrix[is.na(data_matrix) | is.infinite(data_matrix)] <- 0
  return(data_matrix)
}

preprocess_for_cpm <- function(data_matrix, count_type = "expected_count") {
  # Counts Per Million (CPM): library-size normalization
  # Formula: CPM = (count / total_library_count) * 1e6
  # Note: CPM accounts for sequencing depth differences between samples
  # Log2(CPM + 1) applied for variance stabilization (pseudocount of 1)
  # Suitable for visualization; for DEA, use DESeq2's internal normalization
  if (is.null(data_matrix) || nrow(data_matrix) == 0) return(NULL)
  # CPM requires raw read/fragment counts — not valid for coverage (per-base depth)
  if (tolower(count_type) == "coverage") {
    warning("CPM normalization is invalid for coverage data (per-base abundance, not raw counts). ",
            "Falling back to log2(x+1).", call. = FALSE)
    return(preprocess_for_count_type_normalized(data_matrix, count_type))
  }
  lib_sizes <- colSums(data_matrix, na.rm = TRUE)
  zero_libs <- lib_sizes == 0
  if (any(zero_libs)) {
    warning("CPM: ", sum(zero_libs), " sample(s) have zero total counts: ",
            paste(head(names(lib_sizes)[zero_libs], 5), collapse = ", "),
            " -- setting to 1 to avoid division by zero.", call. = FALSE)
    lib_sizes[zero_libs] <- 1
  }
  data_cpm <- sweep(data_matrix, 2, lib_sizes/1e6, FUN = "/")
  data_cpm[is.na(data_cpm) | is.infinite(data_cpm)] <- 0
  return(log2(data_cpm + 1))  # +1 pseudocount before log
}

preprocess_for_count_type_normalized <- function(data_matrix, count_type) {
  if (is.null(data_matrix) || nrow(data_matrix) == 0) return(NULL)
  data_processed <- log2(data_matrix + 1)
  data_processed[is.na(data_processed) | is.infinite(data_processed)] <- 0
  return(data_processed)
}

preprocess_for_zscore <- function(data_matrix, count_type, .log2_cache = NULL) {
  # Z-score normalization: GLOBAL standardization (across entire matrix)
  # All values scaled using global mean and sd
  # This preserves relative gene stability - housekeeping genes will show consistent values
  # while variable genes will show high/low extremes
  if (is.null(data_matrix) || nrow(data_matrix) == 0) return(NULL)
  data_norm <- if (!is.null(.log2_cache)) .log2_cache else preprocess_for_count_type_normalized(data_matrix, count_type)
  if (nrow(data_norm) > 1 && ncol(data_norm) > 1) {
    # Global z-score (preserves relative stability across genes)
    # Single as.matrix() conversion (was called twice — saves one full matrix copy)
    .mat_vals <- as.matrix(data_norm)
    global_mean <- mean(.mat_vals, na.rm = TRUE)
    global_sd <- sd(.mat_vals, na.rm = TRUE)
    if (global_sd == 0 || !is.finite(global_sd)) global_sd <- 1
    data_zscore <- (data_norm - global_mean) / global_sd
    data_zscore[is.na(data_zscore) | is.infinite(data_zscore)] <- 0
    return(data_zscore)
  }
  return(data_norm)
}

preprocess_for_zscore_row <- function(data_matrix, count_type, .log2_cache = NULL) {
  # Z-score normalization: ROW-WISE (per-gene) standardization
  # Each gene is scaled to its own mean and sd across samples
  # This makes all genes equally visible - good for pattern comparison
  # but hides absolute expression level differences between genes
  if (is.null(data_matrix) || nrow(data_matrix) == 0) return(NULL)
  data_norm <- if (!is.null(.log2_cache)) .log2_cache else preprocess_for_count_type_normalized(data_matrix, count_type)
  if (nrow(data_norm) > 1 && ncol(data_norm) > 1) {
    # Row-wise z-score (each gene normalized independently)
    # Vectorized SD: avoid apply() loop using matrix arithmetic
    row_means <- rowMeans(data_norm, na.rm = TRUE)
    n_c <- ncol(data_norm)
    row_sds <- sqrt(rowSums((data_norm - row_means)^2, na.rm = TRUE) / (n_c - 1))
    row_sds[row_sds == 0 | !is.finite(row_sds)] <- 1  # Avoid division by zero
    data_zscore <- (data_norm - row_means) / row_sds
    data_zscore[is.na(data_zscore) | is.infinite(data_zscore)] <- 0
    return(data_zscore)
  }
  return(data_norm)
}

preprocess_for_zscore_scaled_to_ten <- function(data_matrix, count_type, .log2_cache = NULL) {
  # Z-score normalization scaled to 0-10 range
  # First applies global z-score (like preprocess_for_zscore), then rescales to 0-10
  # This preserves z-score patterns but with an intuitive 0-10 scale
  if (is.null(data_matrix) || nrow(data_matrix) == 0) return(NULL)

  # First apply global z-score normalization (pass cache to avoid redundant log2)
  data_zscore <- preprocess_for_zscore(data_matrix, count_type, .log2_cache)
  if (is.null(data_zscore)) return(NULL)
  
  # Then scale z-scores to 0-10 range
  min_val <- min(data_zscore, na.rm = TRUE)
  max_val <- max(data_zscore, na.rm = TRUE)
  val_range <- max_val - min_val
  
  if (val_range > 0 && is.finite(val_range)) {
    data_scaled <- ((data_zscore - min_val) / val_range) * 10
  } else {
    data_scaled <- matrix(5, nrow = nrow(data_zscore), ncol = ncol(data_zscore))
    rownames(data_scaled) <- rownames(data_zscore)
    colnames(data_scaled) <- colnames(data_zscore)
  }
  data_scaled[is.na(data_scaled) | is.infinite(data_scaled)] <- 0
  return(data_scaled)
}

preprocess_for_deseq2_normalized <- function(data_matrix, count_type = "expected_count") {
  # DESeq2-style median-of-ratios normalization
  # Computes size factors based on geometric mean of each gene across samples
  if (is.null(data_matrix) || nrow(data_matrix) == 0) return(NULL)
  # Median-of-ratios assumes raw fragment counts — not valid for coverage (per-base depth)
  if (tolower(count_type) == "coverage") {
    warning("DESeq2 median-of-ratios normalization is invalid for coverage data (per-base abundance). ",
            "Falling back to log2(x+1).", call. = FALSE)
    return(preprocess_for_count_type_normalized(data_matrix, count_type))
  }
  
  # Replace zeros/NAs with small value for geometric mean calculation
  data_clean <- data_matrix
  data_clean[data_clean == 0 | is.na(data_clean)] <- 0.5
  
  # Calculate geometric mean per gene (row)
  log_data <- log(data_clean)
  geo_means <- exp(rowMeans(log_data, na.rm = TRUE))
  
  # Remove genes with zero geometric mean
  valid_genes <- geo_means > 0 & is.finite(geo_means)
  if (sum(valid_genes) == 0) {
    # Fallback to simple log2 if no valid genes
    return(preprocess_for_count_type_normalized(data_matrix, count_type))
  }
  
  # Calculate size factors (median of ratios for each sample)
  ratios <- sweep(data_clean[valid_genes, , drop = FALSE], 1, geo_means[valid_genes], FUN = "/")
  size_factors <- if (requireNamespace("matrixStats", quietly = TRUE)) {
    matrixStats::colMedians(ratios, na.rm = TRUE)
  } else {
    apply(ratios, 2, median, na.rm = TRUE)
  }
  size_factors[size_factors == 0 | !is.finite(size_factors)] <- 1
  
  # Normalize counts by size factors
  data_normalized <- sweep(data_matrix, 2, size_factors, FUN = "/")
  
  # Apply log2 transformation
  data_normalized <- log2(data_normalized + 1)
  data_normalized[is.na(data_normalized) | is.infinite(data_normalized)] <- 0
  
  return(data_normalized)
}

apply_normalization <- function(data_matrix, normalization_scheme, count_type, .log2_cache = NULL) {
  # .log2_cache: optional pre-computed log2(data_matrix + 1) to avoid redundant computation
  # when calling multiple normalization schemes on the same raw data.
  switch(normalization_scheme,
    "raw" = preprocess_for_raw(data_matrix),
    "count_type_normalized" = {
      if (!is.null(.log2_cache)) .log2_cache
      else preprocess_for_count_type_normalized(data_matrix, count_type)
    },
    "deseq2_normalized" = preprocess_for_deseq2_normalized(data_matrix, count_type),
    "zscore" = preprocess_for_zscore(data_matrix, count_type, .log2_cache),
    "zscore_row" = preprocess_for_zscore_row(data_matrix, count_type, .log2_cache),
    "zscore_scaled_to_ten" = preprocess_for_zscore_scaled_to_ten(data_matrix, count_type, .log2_cache),
    "cpm" = preprocess_for_cpm(data_matrix, count_type),
    stop("Unknown normalization scheme: ", normalization_scheme)
  )
}

# ===============================================
# HEATMAP LEGEND HELPERS
# ===============================================

get_legend_layout <- function(position = "bottom") {
  if (position == "bottom") {
    list(direction = "horizontal", height = grid::unit(2, "cm"), width = grid::unit(10, "cm"),
         grid_height = grid::unit(0.8, "cm"), grid_width = grid::unit(2.5, "cm"))
  } else {
    list(direction = "vertical", height = grid::unit(10, "cm"), width = grid::unit(2, "cm"),
         grid_height = grid::unit(2.5, "cm"), grid_width = grid::unit(0.8, "cm"))
  }
}

get_legend_title <- function(normalization_type, count_type = NULL) {
  switch(normalization_type,
    "raw" = if (!is.null(count_type)) paste0("Raw ", toupper(count_type)) else "Raw Counts",
    "count_type_normalized" = "Log2 Expression",
    "deseq2_normalized" = "DESeq2 Norm. Expression",
    "zscore" = "Z-score (Global)",
    "zscore_row" = "Z-score (Per-Gene)",
    "zscore_scaled_to_ten" = "Z-Score [0-10]",
    "cpm" = "Log2(CPM+1)",
    "Expression"
  )
}

get_norm_display_name <- function(norm_scheme) {
  switch(norm_scheme,
    "raw" = "Raw",
    "count_type_normalized" = "Count-Type_Normalized",
    "deseq2_normalized" = "DESeq2_Normalized",
    "zscore" = "Z-score_Global",
    "zscore_row" = "Z-score_Per-Gene",
    "zscore_scaled_to_ten" = "Z-score_Scaled_to_Ten",
    "cpm" = "CPM_Normalized",
    norm_scheme
  )
}

# ===============================================
# SAMPLE LABEL CONVERSION
# ===============================================

convert_to_organ_labels <- function(counts_matrix) {
  current_cols <- colnames(counts_matrix)
  # Short-circuit: if columns are already organ labels (not SRR IDs), skip conversion
  if (all(current_cols %in% SAMPLE_LABELS) || !any(current_cols %in% SAMPLE_IDS)) {
    return(counts_matrix)
  }
  colnames_organ <- SAMPLE_LABELS[current_cols]
  colnames_organ[is.na(colnames_organ)] <- current_cols[is.na(colnames_organ)]
  result <- counts_matrix
  colnames(result) <- colnames_organ
  result
}

# ===============================================
# GENE NAME CONVERSION
# ===============================================

# Cache for gene_groups_csv directory listing (avoids repeated list.files() calls)
.gene_groups_csv_cache <- new.env(hash = TRUE, parent = emptyenv())
# Cache for loaded gene name mappings (avoids re-reading CSV on every convert_to_shortened_names call)
.gene_name_mapping_cache <- new.env(hash = TRUE, parent = emptyenv())

.find_gene_group_csv <- function(gene_group, gene_groups_dir) {
  # Check top-level first (fast path)
  csv_file <- file.path(gene_groups_dir, paste0(gene_group, ".csv"))
  if (file.exists(csv_file)) return(csv_file)

  # Build/reuse cached directory listing (one list.files call per gene_groups_dir)
  cache_key <- gene_groups_dir
  if (is.null(.gene_groups_csv_cache[[cache_key]])) {
    all_files <- list.files(gene_groups_dir, pattern = "\\.csv$",
                            recursive = TRUE, full.names = TRUE)
    .fnames <- tools::file_path_sans_ext(basename(all_files))
    .first <- !duplicated(.fnames)
    file_map <- setNames(all_files[.first], .fnames[.first])
    .gene_groups_csv_cache[[cache_key]] <- file_map
  }
  file_map <- .gene_groups_csv_cache[[cache_key]]
  if (gene_group %in% names(file_map)) return(file_map[[gene_group]])
  return(csv_file)  # return original (non-existent) path for downstream file.exists check
}

# Load gene name mapping from gene_groups CSV files or reference gene_info.csv
# Supported CSV formats:
#   - Gene_ID,Shortened_Name,...
#   - Gene,Shortened_Name,...
#   - Gene_ID,Name,... (for reference gene_info.csv files)
load_gene_name_mapping <- function(gene_group, gene_groups_dir = GENE_GROUPS_DIR) {
  # Return cached mapping if available (avoids re-reading CSV per call)
  cache_key <- paste0(gene_group, "|", gene_groups_dir)
  cached <- .gene_name_mapping_cache[[cache_key]]
  if (!is.null(cached)) return(cached)

  # Use cached directory listing for fast CSV lookup
  csv_file <- .find_gene_group_csv(gene_group, gene_groups_dir)

  # If not found, check if it's a reference with a gene_info.csv
  if (!file.exists(csv_file)) {
    # Check for gene_info.csv companion file in INPUT_FASTAs
    input_fastas_dir <- Sys.getenv("INPUT_FASTAS_DIR", "")
    if (input_fastas_dir == "") {
      # Try to find it relative to workspace
      potential_paths <- c(
        file.path(dirname(dirname(dirname(gene_groups_dir))), "inputs"),
        file.path(dirname(dirname(gene_groups_dir)), "inputs"),
        "../../../../inputs"
      )
      for (path in potential_paths) {
        if (dir.exists(path)) {
          input_fastas_dir <- normalizePath(path, mustWork = FALSE)
          break
        }
      }
    }

    if (input_fastas_dir != "") {
      # Look for gene_info.csv matching the gene_group name
      gene_info_file <- file.path(input_fastas_dir, "mapping", paste0(gene_group, ".gene_info.csv"))
      if (file.exists(gene_info_file)) {
        csv_file <- gene_info_file
      }
    }
  }

  if (!file.exists(csv_file)) return(NULL)

  tryCatch({
    # data.table::fread is 5-10x faster than read.csv for larger gene group files
    df <- if (.HAS_DATATABLE) {
      data.table::fread(csv_file, header = TRUE, data.table = FALSE)
    } else {
      read.csv(csv_file, stringsAsFactors = FALSE, header = TRUE)
    }
    # Support both "Gene_ID" and "Gene" as the ID column
    gene_col <- if ("Gene_ID" %in% colnames(df)) "Gene_ID" else if ("Gene" %in% colnames(df)) "Gene" else NULL
    # Support both "Shortened_Name" and "Name" as the display name column
    name_col <- if ("Shortened_Name" %in% colnames(df)) "Shortened_Name" else if ("Name" %in% colnames(df)) "Name" else NULL

    if (!is.null(gene_col) && !is.null(name_col)) {
      result <- setNames(trimws(df[[name_col]]), trimws(df[[gene_col]]))
      .gene_name_mapping_cache[[cache_key]] <- result
      return(result)
    }
    return(NULL)
  }, error = function(e) NULL)
}

# Convert Gene_ID to Shortened_Name in row names
# Handles both gene-level IDs (SMEL4.1_06g023900) and 
# transcript-level IDs (SMEL4.1_06g023900.1.01)
convert_to_shortened_names <- function(counts_matrix, gene_group) {
  mapping <- load_gene_name_mapping(gene_group)
  if (is.null(mapping)) return(counts_matrix)
  
  current_rownames <- rownames(counts_matrix)

  # Vectorized multi-level suffix stripping: compute all variants at once,
  # then cascade matches (exact → strip .X.XX → strip .X → reverse lookup)
  base_double <- sub("\\.[0-9]+\\.[0-9]+$", "", current_rownames)
  base_single <- sub("\\.[0-9]+$", "", base_double)

  # Cascade: first match wins (avoids repeated subset+reassign passes)
  new_rownames <- mapping[current_rownames]
  na_mask <- is.na(new_rownames)
  if (any(na_mask)) {
    hits2 <- mapping[base_double[na_mask]]
    new_rownames[na_mask] <- hits2
    na_mask <- is.na(new_rownames)
  }
  if (any(na_mask)) {
    hits3 <- mapping[base_single[na_mask]]
    new_rownames[na_mask] <- hits3
    na_mask <- is.na(new_rownames)
  }

  # Reverse lookup: row ID is shorter than mapping key (e.g., gene-level "SMEL5_06g022750"
  # when mapping has transcript-level "SMEL5_06g022750.1" as key).
  # Build a stripped-key mapping and try matching.
  if (any(na_mask)) {
    mapping_keys <- names(mapping)
    stripped_keys <- sub("\\.[0-9]+$", "", mapping_keys)
    # Only use entries where stripping actually changed the key (avoids false matches)
    changed <- stripped_keys != mapping_keys
    if (any(changed)) {
      reverse_mapping <- setNames(mapping[changed], stripped_keys[changed])
      # Remove duplicates (keep first occurrence)
      reverse_mapping <- reverse_mapping[!duplicated(names(reverse_mapping))]
      reverse_hits <- reverse_mapping[current_rownames[na_mask]]
      new_rownames[na_mask] <- reverse_hits
    }
  }

  # Keep original name if still no mapping found
  na_mask <- is.na(new_rownames)
  new_rownames[na_mask] <- current_rownames[na_mask]
  
  result <- counts_matrix
  rownames(result) <- new_rownames
  result
}

# Apply label transformations based on gene_type and label_type
# Also reorders columns to match SAMPLE_IDS order (from CSV file)
# NOTE: Does NOT filter out columns - only reorders columns that match SAMPLE_IDS
#       Columns not in SAMPLE_IDS are appended at the end in their original order
apply_labels <- function(counts_matrix, gene_group, gene_type, label_type) {
  result <- counts_matrix
  
  # Reorder columns to match SAMPLE_IDS order (preserves CSV file order)
  # Important: Keep ALL columns - just reorder those that match
  if (length(SAMPLE_IDS) > 0 && ncol(result) > 0) {
    current_cols <- colnames(result)
    
    # Check if columns are organ labels (matching SAMPLE_LABELS values)
    # Handle R's make.unique suffixes (.1, .2, etc.) from read.table on duplicate organ names
    base_cols <- sub("\\.[0-9]+$", "", current_cols)
    if (any(base_cols %in% SAMPLE_LABELS)) {
      # Columns are organ labels - reorder using SAMPLE_LABELS values order
      ordered_organs <- SAMPLE_LABELS[SAMPLE_IDS]
      ordered_organs <- ordered_organs[!is.na(ordered_organs)]
      # Match against base organ names (stripped of .1/.2 suffixes) so that
      # duplicates like "Flower_Buds.1" correctly match "Flower_Buds".
      # Vectorized approach: for each ordered_organ, find its first unused column match.
      # When organs repeat, we need sequential consumption (first match, then second, etc.)
      # Group column indices by base organ name for O(1) lookup per organ.
      col_groups <- split(seq_along(base_cols), base_cols)
      # Track consumption position within each group
      group_pos <- integer(length(col_groups))
      names(group_pos) <- names(col_groups)
      matched_indices <- integer(length(ordered_organs))
      n_matched <- 0L
      for (organ in ordered_organs) {
        if (!is.null(col_groups[[organ]])) {
          pos <- group_pos[[organ]] + 1L
          indices <- col_groups[[organ]]
          if (pos <= length(indices)) {
            n_matched <- n_matched + 1L
            matched_indices[n_matched] <- indices[pos]
            group_pos[[organ]] <- pos
          }
        }
      }
      matched_indices <- matched_indices[seq_len(n_matched)]
      used <- logical(ncol(result))
      used[matched_indices] <- TRUE
      # Columns not matched by the expected order (keep at end)
      unmatched_indices <- which(!used)
      ordered_indices <- c(matched_indices, unmatched_indices)
      if (length(ordered_indices) == ncol(result) &&
          !identical(ordered_indices, seq_len(ncol(result)))) {
        result <- result[, ordered_indices, drop = FALSE]
      }
    } else if (any(current_cols %in% SAMPLE_IDS)) {
      # Columns are SRR IDs - reorder using SAMPLE_IDS order directly
      matched_cols <- SAMPLE_IDS[SAMPLE_IDS %in% current_cols]
      unmatched_cols <- current_cols[!current_cols %in% matched_cols]
      ordered_cols <- c(matched_cols, unmatched_cols)
      if (length(ordered_cols) > 0 && !identical(ordered_cols, current_cols)) {
        result <- result[, ordered_cols, drop = FALSE]
      }
    }
    # If columns don't match either pattern, keep original order (no filtering)
  }
  
  # Apply row label transformation (gene names)
  if (gene_type == "Shortened_Name") {
    result <- convert_to_shortened_names(result, gene_group)
  }
  # gene_type == "Gene_ID" means keep original row names
  
  # Apply column label transformation (sample names)
  if (label_type == "Organ") {
    result <- convert_to_organ_labels(result)
  }
  # label_type == "SRR_ID" means keep original column names
  
  return(result)
}

# ===============================================
# TISSUE GROUP MAPPING (DRY - used by PCA, Correlation, DEA)
# ===============================================
# Maps organ names from CSV to biological tissue groups
# Update this function when CSV organ names change

map_tissue_to_group <- function(tissue_name) {
  if (grepl("Root|Stem|Leaf|Leaves|Senescent", tissue_name, ignore.case = TRUE)) {
    return("Vegetative")
  }
  if (grepl("Flower|Bud|Pistil|Stamen", tissue_name, ignore.case = TRUE)) {
    return("Reproductive")
  }
  if (grepl("Fruit|peduncle", tissue_name, ignore.case = TRUE)) {
    return("Fruit")
  }
  if (grepl("Radicle|Cotyledon", tissue_name, ignore.case = TRUE)) {
    return("Seedling")
  }
  return("Other")
}

# Vectorized version using grepl on entire vector (avoids per-element sapply overhead)
get_tissue_groups <- function(tissue_names) {
  result <- rep("Other", length(tissue_names))
  # Later matches override earlier ones, so check in reverse priority
  mask_seedling     <- grepl("Radicle|Cotyledon", tissue_names, ignore.case = TRUE)
  mask_fruit        <- grepl("Fruit|peduncle", tissue_names, ignore.case = TRUE)
  mask_reproductive <- grepl("Flower|Bud|Pistil|Stamen", tissue_names, ignore.case = TRUE)
  mask_vegetative   <- grepl("Root|Stem|Leaf|Leaves|Senescent", tissue_names, ignore.case = TRUE)
  result[mask_seedling]     <- "Seedling"
  result[mask_fruit]        <- "Fruit"
  result[mask_reproductive] <- "Reproductive"
  result[mask_vegetative]   <- "Vegetative"
  result
}

# ===============================================
# COLOR FUNCTIONS
# ===============================================

get_violet_color_scale <- function(n_breaks = 100) {
  # Light lavender to deep purple gradient (matching image palette)
  #colorRampPalette(c("#E8D5F0", "#D4B5E3", "#C095D6", "#AC75C9", "#9855BC", "#8435AF", "#6F1FA2", "#5A0F8F", "#45007C"))(n_breaks)
  colorRampPalette(c("#dab3ddff", "#d9afe0ff", "#c57fd1ff", "#ac44beff", "#8E24AA", "#6A1B9A", "#4A148C", "#2F1B69"))(n_breaks)
}

get_blue_red_color_scale <- function(n_breaks = 100) {
  colorRampPalette(c("#2166AC", "#67A9CF", "#F7F7F7", "#EF8A62", "#B2182B"))(n_breaks)
}

get_cv_color_scale <- function(n_breaks = 100) {
  # Deep violet to light lavender for CV: low CV = deep (stable), high CV = pale (variable)
  # colorRamp2(seq(min,max,...), colors) maps min→first color, max→last color
  colorRampPalette(c("#4A148C", "#7B1FA2", "#9C27B0", "#AB47BC",
                     "#BA68C8", "#CE93D8", "#E1BEE7", "#F3E5F5"))(n_breaks)
}
