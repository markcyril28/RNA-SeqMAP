#!/usr/bin/env Rscript

# ===============================================
# SAMPLE CORRELATION & CLUSTERING MODULE
# ===============================================
# Sample-sample correlation and hierarchical clustering

suppressPackageStartupMessages({
  library(pheatmap)
  library(RColorBrewer)
  library(ggplot2)
})

SCRIPT_DIR <- Sys.getenv("ANALYSIS_MODULES_DIR", ".")
source(file.path(SCRIPT_DIR, "0_shared_config.R"))
source(file.path(SCRIPT_DIR, "1_utility_functions.R"))

# ===============================================
# CONFIGURATION
# ===============================================

CORRELATION_OUT_DIR <- file.path(CONSOLIDATED_BASE_DIR, OUTPUT_SUBDIRS$CORRELATION)

# Correlation parameters
# Spearman is preferred for RNA-seq: robust to outliers, handles non-linear relationships
# Pearson assumes normality and linear relationships (use on log-transformed data)
# For RNA-seq QC, correlation > 0.9 between replicates indicates good quality
CORRELATION_METHOD <- "spearman"  # "pearson", "spearman", "kendall"
CLUSTERING_METHOD <- "complete"   # "complete", "ward.D2", "average"
SHOW_NUMBERS <- TRUE
NUMBER_SIZE <- 8
# MIN_SAMPLES defined in 0_shared_config.R (requires >= 3 for correlation)

# Figure toggles
GENERATE_COR_FIGURES <- list(
  correlation_heatmap = TRUE,
  dendrogram = TRUE,
  scatterplot_matrix = FALSE
)

# ===============================================
# CORRELATION FUNCTIONS
# ===============================================

calculate_sample_correlation <- function(data_matrix, method = CORRELATION_METHOD) {
  # Use GPU-accelerated correlation when available (for Pearson)
  # Falls back to CPU for Spearman/Kendall or when GPU unavailable
  if (method == "pearson" && GPU_AVAILABLE) {
    cor_matrix <- gpu_cor(data_matrix, method = method)
  } else {
    cor_matrix <- cor(data_matrix, use = "pairwise.complete.obs", method = method)
  }
  cor_matrix[is.na(cor_matrix)] <- 0
  return(cor_matrix)
}

create_correlation_heatmap <- function(cor_matrix, output_path, title, 
                                        annotation_df = NULL, show_numbers = SHOW_NUMBERS) {
  tryCatch({
    # Color palette
    colors <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100)
    
    # Annotation colors
    ann_colors <- NULL
    if (!is.null(annotation_df)) {
      tissue_groups <- unique(annotation_df$TissueGroup)
      n_groups <- length(tissue_groups)
      group_colors <- if (n_groups <= 8) {
        brewer.pal(max(3, n_groups), "Set2")[1:n_groups]
      } else {
        colorRampPalette(brewer.pal(8, "Set2"))(n_groups)
      }
      names(group_colors) <- tissue_groups
      ann_colors <- list(TissueGroup = group_colors)
    }
    
    png(output_path, width = 1000, height = 900, res = 100)
    # Strip R's make.unique suffixes (.1, .2) from organ labels for display
    clean_labels <- sub("\\.[0-9]+$", "", rownames(cor_matrix))
    pheatmap(
      cor_matrix,
      main = title,
      color = colors,
      breaks = seq(-1, 1, length.out = 101),
      clustering_method = CLUSTERING_METHOD,
      display_numbers = show_numbers,
      number_format = "%.2f",
      fontsize_number = NUMBER_SIZE,
      labels_row = clean_labels,
      labels_col = clean_labels,
      annotation_row = annotation_df,
      annotation_col = annotation_df,
      annotation_colors = ann_colors,
      border_color = NA
    )
    dev.off()
    
    cat("    Generated:", basename(output_path), "\n")
    return(TRUE)
  }, error = function(e) {
    cat("    Error:", e$message, "\n")
    return(FALSE)
  })
}

create_dendrogram <- function(cor_matrix, output_path, title) {
  if (!GENERATE_COR_FIGURES$dendrogram) return(NULL)
  
  tryCatch({
    # Convert correlation to distance (1 - cor)
    # Use GPU-accelerated distance if working with raw data
    dist_matrix <- as.dist(1 - cor_matrix)
    hc <- hclust(dist_matrix, method = CLUSTERING_METHOD)
    # Strip R's make.unique suffixes (.1, .2) from organ labels for display
    hc$labels <- sub("\\.[0-9]+$", "", hc$labels)
    
    png(output_path, width = 1000, height = 600, res = 100)
    plot(hc, main = title, xlab = "", sub = "", hang = -1)
    dev.off()
    
    return(hc)
  }, error = function(e) {
    cat("    Dendrogram error:", e$message, "\n")
    return(NULL)
  })
}

# ===============================================
# MAIN FUNCTION
# ===============================================

run_sample_correlation <- function(config = NULL, matrices_dir = NULL) {
  # Get method base directory from environment for config file loading
  method_base_dir <- Sys.getenv("METHOD_BASE_DIR", unset = ".")
  if (is.null(config)) config <- load_runtime_config(method_base_dir)
  if (is.null(matrices_dir)) {
    matrices_dir <- file.path(method_base_dir, get_matrices_dir(CURRENT_METHOD))
  }
  ensure_output_dir(CORRELATION_OUT_DIR)
  
  print_config_summary("SAMPLE CORRELATION & CLUSTERING", config)
  
  # Log GPU status
  if (GPU_AVAILABLE && CORRELATION_METHOD == "pearson") {
    cat("GPU acceleration: ENABLED for Pearson correlation (", GPU_VRAM_GB, "GB VRAM)\n", sep = "")
  } else if (CORRELATION_METHOD != "pearson") {
    cat("GPU acceleration: N/A (", CORRELATION_METHOD, " requires CPU)\n", sep = "")
  }
  
  successful <- 0
  total <- 0
  
  for (gene_group in config$gene_groups) {
    cat("Processing:", gene_group, "\n")
    total <- total + 1
    
    output_folder_name <- get_output_folder_name(gene_group, CURRENT_DATASET)
    output_dir <- file.path(CORRELATION_OUT_DIR, output_folder_name)
    ensure_output_dir(output_dir)
    
    input_file <- build_input_path(gene_group, PROCESSING_LEVELS[1],
                                   COUNT_TYPES[1], GENE_TYPES[1],
                                   matrices_dir, config$master_reference)
    
    validation <- validate_and_read_matrix(input_file, 5)
    if (!validation$success) {
      cat("  Skipped:", validation$reason, "\n")
      next
    }
    
    data_matrix <- apply_normalization(validation$data, NORM_SCHEMES[1], COUNT_TYPES[1])
    
    # Calculate correlation
    cor_matrix <- calculate_sample_correlation(data_matrix, CORRELATION_METHOD)
    
    # Save correlation matrix
    cor_df <- as.data.frame(cor_matrix)
    cor_df <- cbind(Sample = rownames(cor_df), cor_df)
    write.table(cor_df, file.path(output_dir, paste0(gene_group, "_correlation_matrix.tsv")),
                sep = "\t", row.names = FALSE, quote = FALSE)
    
    # Prepare annotation
    sample_ids <- colnames(data_matrix)
    tissues <- SAMPLE_LABELS[sample_ids]
    tissues[is.na(tissues)] <- sample_ids[is.na(tissues)]
    
    # Use centralized tissue group mapping from 1_utility_functions.R
    tissue_groups <- get_tissue_groups(tissues)
    
    annotation_df <- data.frame(
      TissueGroup = tissue_groups,
      row.names = colnames(cor_matrix)
    )
    
    # Create heatmap
    if (GENERATE_COR_FIGURES$correlation_heatmap) {
      create_correlation_heatmap(
        cor_matrix,
        file.path(output_dir, paste0(gene_group, "_correlation_heatmap.png")),
        paste0(gene_group, " - Sample Correlation (", CORRELATION_METHOD, ")"),
        annotation_df
      )
    }
    
    # Create dendrogram
    create_dendrogram(
      cor_matrix,
      file.path(output_dir, paste0(gene_group, "_dendrogram.png")),
      paste0(gene_group, " - Sample Clustering")
    )
    
    successful <- successful + 1
    cat("  Complete\n")
  }
  
  print_summary(successful, total)
}

if (!interactive() && identical(environment(), globalenv())) {
  run_sample_correlation()
}
