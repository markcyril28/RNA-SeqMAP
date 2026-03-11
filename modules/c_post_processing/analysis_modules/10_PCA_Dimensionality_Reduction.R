#!/usr/bin/env Rscript

# ===============================================
# PCA & DIMENSIONALITY REDUCTION MODULE
# ===============================================
# PCA, t-SNE, UMAP for sample visualization

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(RColorBrewer)
  library(Rtsne)
  library(umap)
  library(factoextra)
})

SCRIPT_DIR <- Sys.getenv("ANALYSIS_MODULES_DIR", ".")
source(file.path(SCRIPT_DIR, "0_shared_config.R"))
source(file.path(SCRIPT_DIR, "1_utility_functions.R"))

# ===============================================
# CONFIGURATION
# ===============================================

DIM_REDUCTION_OUT_DIR <- file.path(CONSOLIDATED_BASE_DIR, OUTPUT_SUBDIRS$DIM_REDUCTION)

# PCA parameters (MAXIMUM ACCURACY - using all 32GB RAM)
# Use 0 or Inf to include ALL genes (no subsampling for maximum accuracy)
N_TOP_VARIABLE_GENES <- 0          # 0 = use all genes, no filtering
CENTER_DATA <- TRUE
SCALE_DATA <- TRUE
N_PCS_TO_PLOT <- 10                # More principal components for thorough analysis

# t-SNE parameters (MAXIMUM ACCURACY - exact algorithm)
# Perplexity must be < (n_samples - 1) / 3 to avoid errors
# For n=10 samples: max perplexity = 3; for n=30: max perplexity = 9
# Dynamically adjusted in run_tsne_analysis()
TSNE_PERPLEXITY <- 15              # Higher perplexity for global structure
TSNE_MAX_ITER <- 5000              # Maximum iterations for full convergence
TSNE_SEED <- GLOBAL_RANDOM_SEED     # Reproducibility (from shared config)
MIN_SAMPLES_TSNE <- 4              # t-SNE needs >= 4 samples
TSNE_THETA <- 0.0                  # Exact t-SNE (0 = exact, no approximation)

# UMAP parameters (MAXIMUM ACCURACY)
UMAP_N_NEIGHBORS <- 15             # More neighbors for global structure
UMAP_MIN_DIST <- 0.05              # Tighter clustering for better separation
UMAP_METRIC <- "cosine"            # Cosine similarity (best for expression data)
UMAP_SEED <- GLOBAL_RANDOM_SEED
UMAP_N_EPOCHS <- 1000              # Maximum epochs for full convergence

# Figure toggles
GENERATE_DIM_FIGURES <- list(
  pca_scree = TRUE,
  pca_scatter = TRUE,
  pca_loadings = TRUE,
  tsne_scatter = TRUE,
  umap_scatter = TRUE,
  combined_panel = TRUE
)

# ===============================================
# DATA PREPROCESSING
# ===============================================

select_variable_genes <- function(data_matrix, n_top = N_TOP_VARIABLE_GENES) {
  if (n_top <= 0 || n_top >= nrow(data_matrix)) return(data_matrix)
  
  gene_vars <- apply(data_matrix, 1, var, na.rm = TRUE)
  gene_vars[is.na(gene_vars)] <- 0
  
  top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:min(n_top, length(gene_vars))]
  return(data_matrix[top_genes, , drop = FALSE])
}

prepare_sample_metadata <- function(sample_ids) {
  tissues <- SAMPLE_LABELS[sample_ids]
  tissues[is.na(tissues)] <- sample_ids[is.na(tissues)]
  
  # Strip R's make.unique suffixes (.1, .2, etc.) from organ labels for display
  tissues <- sub("\\.[0-9]+$", "", tissues)
  
  # Use centralized tissue group mapping from 1_utility_functions.R
  tissue_groups <- get_tissue_groups(tissues)
  
  data.frame(
    Sample = sample_ids,
    Tissue = tissues,
    TissueGroup = tissue_groups,
    stringsAsFactors = FALSE
  )
}

# ===============================================
# PCA FUNCTIONS
# ===============================================

run_pca_analysis <- function(data_matrix, metadata, output_dir, gene_group) {
  cat("    Running PCA")
  if (GPU_AVAILABLE) cat(" (GPU accelerated)")
  cat("\n")

  # PCA needs at least 3 samples for meaningful 2D projection
  if (ncol(data_matrix) < 3) {
    cat("    Too few samples for PCA (need >= 3, have", ncol(data_matrix), ")\n")
    return(NULL)
  }

  data_t <- t(data_matrix)

  col_vars <- apply(data_t, 2, var, na.rm = TRUE)
  data_t <- data_t[, col_vars > 0, drop = FALSE]

  if (ncol(data_t) < 3) {
    cat("    Too few variable genes for PCA\n")
    return(NULL)
  }
  
  # Use GPU-accelerated PCA when available
  pca_result <- gpu_prcomp(data_t, center = CENTER_DATA, scale. = SCALE_DATA)
  var_explained <- summary(pca_result)$importance[2, ] * 100
  
  if (GENERATE_DIM_FIGURES$pca_scree) {
    png(file.path(output_dir, paste0(gene_group, "_PCA_scree.png")),
        width = 800, height = 600, res = 100)
    print(fviz_eig(pca_result, addlabels = TRUE) +
      ggtitle(paste0(gene_group, " - Variance Explained")) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold")))
    dev.off()
  }
  
  pca_df <- as.data.frame(pca_result$x)
  pca_df$Sample <- rownames(pca_df)
  pca_df <- merge(pca_df, metadata, by = "Sample")
  
  write.table(pca_df, file.path(output_dir, paste0(gene_group, "_PCA_coords.tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  if (GENERATE_DIM_FIGURES$pca_scatter) {
    p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = TissueGroup)) +
      geom_point(size = 4) +
      geom_text_repel(aes(label = Tissue), size = 3) +
      labs(title = paste0(gene_group, " - PCA"),
           x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
           y = paste0("PC2 (", round(var_explained[2], 1), "%)")) +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    
    ggsave(file.path(output_dir, paste0(gene_group, "_PCA.png")),
           p, width = 10, height = 8, dpi = 150)
  }
  
  return(list(pca = pca_result, df = pca_df, var = var_explained))
}

# ===============================================
# t-SNE FUNCTIONS
# ===============================================

run_tsne_analysis <- function(data_matrix, metadata, output_dir, gene_group) {
  if (!GENERATE_DIM_FIGURES$tsne_scatter) return(NULL)
  
  cat("    Running t-SNE (exact algorithm for accuracy)\n")
  
  data_t <- t(data_matrix)
  
  # Dynamic perplexity: must be < (n - 1) / 3
  n_samples <- nrow(data_t)
  max_perplexity <- floor((n_samples - 1) / 3)
  effective_perplexity <- min(TSNE_PERPLEXITY, max_perplexity)
  
  if (n_samples < MIN_SAMPLES_TSNE || effective_perplexity < 1) {
    cat("    Too few samples for t-SNE (need >=", MIN_SAMPLES_TSNE, ")\n")
    return(NULL)
  }
  
  set.seed(TSNE_SEED)
  # Use exact t-SNE (theta=0) for maximum accuracy when sample count is manageable
  theta_val <- if (exists("TSNE_THETA")) TSNE_THETA else 0.0
  
  tsne_result <- Rtsne(data_t, dims = 2, 
                       perplexity = effective_perplexity,
                       max_iter = TSNE_MAX_ITER, 
                       theta = theta_val,          # 0 = exact, higher = faster approximation
                       check_duplicates = FALSE,
                       pca = TRUE,                 # PCA initialization for stability
                       normalize = TRUE)           # Normalize input for better results
  
  tsne_df <- data.frame(
    Sample = rownames(data_t),
    tSNE1 = tsne_result$Y[, 1],
    tSNE2 = tsne_result$Y[, 2]
  )
  tsne_df <- merge(tsne_df, metadata, by = "Sample")
  
  p <- ggplot(tsne_df, aes(x = tSNE1, y = tSNE2, color = TissueGroup)) +
    geom_point(size = 4) +
    geom_text_repel(aes(label = Tissue), size = 3) +
    labs(title = paste0(gene_group, " - t-SNE")) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  ggsave(file.path(output_dir, paste0(gene_group, "_tSNE.png")),
         p, width = 10, height = 8, dpi = 150)
  
  return(tsne_df)
}

# ===============================================
# UMAP FUNCTIONS
# ===============================================

run_umap_analysis <- function(data_matrix, metadata, output_dir, gene_group) {
  if (!GENERATE_DIM_FIGURES$umap_scatter) return(NULL)
  
  cat("    Running UMAP\n")
  
  data_t <- t(data_matrix)
  
  # UMAP needs at least 4 samples; n_neighbors is clamped below via min()
  if (nrow(data_t) < 4) {
    cat("    Too few samples for UMAP (need >= 4, have", nrow(data_t), ")\n")
    return(NULL)
  }
  
  set.seed(UMAP_SEED)
  umap_config <- umap.defaults
  umap_config$n_neighbors <- min(UMAP_N_NEIGHBORS, nrow(data_t) - 1)
  umap_config$min_dist <- UMAP_MIN_DIST
  umap_config$metric <- UMAP_METRIC
  umap_config$n_epochs <- if (exists("UMAP_N_EPOCHS")) UMAP_N_EPOCHS else 200
  
  umap_result <- umap(data_t, config = umap_config)
  
  umap_df <- data.frame(
    Sample = rownames(data_t),
    UMAP1 = umap_result$layout[, 1],
    UMAP2 = umap_result$layout[, 2]
  )
  umap_df <- merge(umap_df, metadata, by = "Sample")
  
  p <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = TissueGroup)) +
    geom_point(size = 4) +
    geom_text_repel(aes(label = Tissue), size = 3) +
    labs(title = paste0(gene_group, " - UMAP")) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  ggsave(file.path(output_dir, paste0(gene_group, "_UMAP.png")),
         p, width = 10, height = 8, dpi = 150)
  
  return(umap_df)
}

# ===============================================
# MAIN FUNCTION
# ===============================================

run_dimensionality_reduction <- function(config = NULL, matrices_dir = NULL) {
  # Get method base directory from environment for config file loading
  method_base_dir <- Sys.getenv("METHOD_BASE_DIR", unset = ".")
  if (is.null(config)) config <- load_runtime_config(method_base_dir)
  if (is.null(matrices_dir)) {
    matrices_dir <- file.path(method_base_dir, get_matrices_dir(CURRENT_METHOD))
  }
  ensure_output_dir(DIM_REDUCTION_OUT_DIR)
  
  print_config_summary("PCA & DIMENSIONALITY REDUCTION", config)
  
  # Log GPU status for PCA computation
  if (GPU_AVAILABLE) {
    cat("GPU acceleration: ENABLED for PCA (", GPU_VRAM_GB, "GB VRAM)\n", sep = "")
  }
  
  successful <- 0
  total <- 0
  
  for (gene_group in config$gene_groups) {
    cat("Processing:", gene_group, "\n")
    total <- total + 1
    
    output_folder_name <- get_output_folder_name(gene_group, CURRENT_DATASET)
    output_dir <- file.path(DIM_REDUCTION_OUT_DIR, output_folder_name)
    ensure_output_dir(output_dir)
    
    input_file <- build_input_path(gene_group, PROCESSING_LEVELS[1],
                                   COUNT_TYPES[1], GENE_TYPES[1],
                                   matrices_dir, config$master_reference)
    
    validation <- validate_and_read_matrix(input_file, 5)
    if (!validation$success) {
      cat("  Skipped:", validation$reason, "\n")
      next
    }
    
    data_matrix <- validation$data
    if (N_TOP_VARIABLE_GENES > 0 && nrow(data_matrix) > N_TOP_VARIABLE_GENES) {
      data_matrix <- select_variable_genes(data_matrix, N_TOP_VARIABLE_GENES)
    }
    
    data_matrix <- apply_normalization(data_matrix, NORM_SCHEMES[1], COUNT_TYPES[1])
    
    metadata <- prepare_sample_metadata(colnames(data_matrix))
    
    tryCatch(run_pca_analysis(data_matrix, metadata, output_dir, gene_group),
             error = function(e) cat("  PCA failed:", e$message, "\n"))
    tryCatch(run_tsne_analysis(data_matrix, metadata, output_dir, gene_group),
             error = function(e) cat("  t-SNE failed:", e$message, "\n"))
    tryCatch(run_umap_analysis(data_matrix, metadata, output_dir, gene_group),
             error = function(e) cat("  UMAP failed:", e$message, "\n"))
    
    successful <- successful + 1
    cat("  Complete\n")
  }
  
  print_summary(successful, total)
}

if (!interactive() && identical(environment(), globalenv())) {
  run_dimensionality_reduction()
}
