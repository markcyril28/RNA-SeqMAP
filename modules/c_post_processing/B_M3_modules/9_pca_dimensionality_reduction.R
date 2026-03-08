#!/usr/bin/env Rscript

# ===============================================
# PCA & DIMENSIONALITY REDUCTION FOR METHOD 3 (STAR + SALMON)
# ===============================================
# Performs PCA, t-SNE, and UMAP for sample visualization
# Detects batch effects and outliers

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(RColorBrewer)
  library(Rtsne)
  library(umap)
  library(factoextra)
})

script_dir <- dirname(sys.frame(1)$ofile)
source(file.path(script_dir, "0_shared_config.R"))
source(file.path(script_dir, "2_utility_functions.R"))

# ===============================================
# CONFIGURATION
# ===============================================

DIM_REDUCTION_SUBDIR <- "VII_Dimensionality_Reduction"
DIM_REDUCTION_OUT_DIR <- file.path(CONSOLIDATED_BASE_DIR, DIM_REDUCTION_SUBDIR)

# PCA parameters
N_TOP_VARIABLE_GENES <- 1000       # Top variable genes for PCA (0 = use all)
CENTER_DATA <- TRUE                # Center data before PCA
SCALE_DATA <- TRUE                 # Scale data before PCA
N_PCS_TO_PLOT <- 4                 # Number of PCs to plot pairwise

# t-SNE parameters
TSNE_PERPLEXITY <- 5               # Perplexity (5-50, lower for small datasets)
TSNE_MAX_ITER <- 1000              # Max iterations
TSNE_SEED <- 42                    # Random seed for reproducibility

# UMAP parameters
UMAP_N_NEIGHBORS <- 5              # Number of neighbors (2-100)
UMAP_MIN_DIST <- 0.3               # Minimum distance (0.0-0.99)
UMAP_METRIC <- "euclidean"         # Distance metric
UMAP_SEED <- 42                    # Random seed

# Load runtime configuration
config <- load_runtime_config()
ensure_output_dir(DIM_REDUCTION_OUT_DIR)

# FIGURE GENERATION CONTROL
GENERATE_DIM_FIGURES <- list(
  pca_scree = TRUE,                  # Scree plot (variance explained)
  pca_scatter = TRUE,                # PC1 vs PC2, PC3 vs PC4
  pca_loadings = TRUE,               # Top gene loadings per PC
  pca_biplot = FALSE,                # Biplot with sample and gene arrows
  tsne_scatter = TRUE,               # t-SNE 2D plot
  umap_scatter = TRUE,               # UMAP 2D plot
  combined_panel = TRUE              # Combined PCA, t-SNE, UMAP panel
)

# ===============================================
# DATA PREPROCESSING
# ===============================================

# Select top variable genes
select_variable_genes <- function(data_matrix, n_top = N_TOP_VARIABLE_GENES) {
  if (n_top <= 0 || n_top >= nrow(data_matrix)) {
    return(data_matrix)
  }
  
  # Calculate variance per gene
  gene_vars <- apply(data_matrix, 1, var, na.rm = TRUE)
  gene_vars[is.na(gene_vars)] <- 0
  
  # Select top variable
  top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:min(n_top, length(gene_vars))]
  
  return(data_matrix[top_genes, , drop = FALSE])
}

# Prepare sample metadata
prepare_sample_metadata <- function(sample_ids) {
  tissues <- SAMPLE_LABELS[sample_ids]
  tissues[is.na(tissues)] <- sample_ids[is.na(tissues)]
  
  # Extract tissue group
  tissue_groups <- sapply(tissues, function(t) {
    if (grepl("Root|Stem|Leaf|Leaves", t, ignore.case = TRUE)) return("Vegetative")
    if (grepl("Flower|Bud|Pistil", t, ignore.case = TRUE)) return("Reproductive")
    if (grepl("Fruit", t, ignore.case = TRUE)) return("Fruit")
    if (grepl("Radicle|Cotyledon|Seed", t, ignore.case = TRUE)) return("Seedling")
    return("Other")
  })
  
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
  cat("    Running PCA\n")
  
  # Transpose: samples as rows, genes as columns
  data_t <- t(data_matrix)
  
  # Remove constant columns
  col_vars <- apply(data_t, 2, var, na.rm = TRUE)
  data_t <- data_t[, col_vars > 0, drop = FALSE]
  
  if (ncol(data_t) < 3) {
    cat("    Warning: Too few variable genes for PCA\n")
    return(NULL)
  }
  
  # Run PCA
  pca_result <- prcomp(data_t, center = CENTER_DATA, scale. = SCALE_DATA)
  
  # Variance explained
  var_explained <- summary(pca_result)$importance[2, ] * 100
  
  # Scree plot
  if (GENERATE_DIM_FIGURES$pca_scree) {
    png(file.path(output_dir, paste0(gene_group, "_PCA_scree.png")), 
        width = 800, height = 600, res = 100)
    fviz_eig(pca_result, addlabels = TRUE, ylim = c(0, max(var_explained) * 1.1)) +
      ggtitle(paste0(gene_group, " - PCA Variance Explained")) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    dev.off()
  }
  
  # Create PCA data frame
  pca_df <- as.data.frame(pca_result$x)
  pca_df$Sample <- rownames(pca_df)
  pca_df <- merge(pca_df, metadata, by = "Sample")
  
  # Export PCA coordinates
  write.table(pca_df, 
              file = file.path(output_dir, paste0(gene_group, "_PCA_coordinates.tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  # PCA scatter plots
  if (GENERATE_DIM_FIGURES$pca_scatter) {
    # PC1 vs PC2
    p1 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = TissueGroup, shape = TissueGroup)) +
      geom_point(size = 4, alpha = 0.8) +
      geom_text_repel(aes(label = Tissue), size = 2.5, max.overlaps = 20) +
      scale_color_brewer(palette = "Set1") +
      labs(title = paste0(gene_group, " - PCA"),
           x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
           y = paste0("PC2 (", round(var_explained[2], 1), "%)")) +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"),
            legend.position = "right")
    
    ggsave(file.path(output_dir, paste0(gene_group, "_PCA_PC1_PC2.png")),
           p1, width = 12, height = 10, dpi = 150)
    
    # PC3 vs PC4 if available
    if (ncol(pca_result$x) >= 4) {
      p2 <- ggplot(pca_df, aes(x = PC3, y = PC4, color = TissueGroup, shape = TissueGroup)) +
        geom_point(size = 4, alpha = 0.8) +
        geom_text_repel(aes(label = Tissue), size = 2.5, max.overlaps = 20) +
        scale_color_brewer(palette = "Set1") +
        labs(title = paste0(gene_group, " - PCA"),
             x = paste0("PC3 (", round(var_explained[3], 1), "%)"),
             y = paste0("PC4 (", round(var_explained[4], 1), "%)")) +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))
      
      ggsave(file.path(output_dir, paste0(gene_group, "_PCA_PC3_PC4.png")),
             p2, width = 12, height = 10, dpi = 150)
    }
  }
  
  # Gene loadings
  if (GENERATE_DIM_FIGURES$pca_loadings) {
    loadings <- as.data.frame(pca_result$rotation[, 1:min(4, ncol(pca_result$rotation))])
    loadings$Gene <- rownames(loadings)
    
    # Top loadings per PC
    for (pc in paste0("PC", 1:min(2, ncol(pca_result$rotation)))) {
      top_pos <- head(loadings[order(-loadings[[pc]]), ], 10)
      top_neg <- head(loadings[order(loadings[[pc]]), ], 10)
      top_loadings <- rbind(top_pos, top_neg)
      top_loadings <- top_loadings[order(top_loadings[[pc]]), ]
      
      p <- ggplot(top_loadings, aes(x = reorder(Gene, get(pc)), y = get(pc))) +
        geom_bar(stat = "identity", fill = ifelse(top_loadings[[pc]] > 0, "#B2182B", "#2166AC")) +
        coord_flip() +
        labs(title = paste0(gene_group, " - Top ", pc, " Loadings"),
             x = "Gene", y = paste0(pc, " Loading")) +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))
      
      ggsave(file.path(output_dir, paste0(gene_group, "_", pc, "_loadings.png")),
             p, width = 10, height = 8, dpi = 150)
    }
    
    write.table(loadings,
                file = file.path(output_dir, paste0(gene_group, "_PCA_loadings.tsv")),
                sep = "\t", row.names = FALSE, quote = FALSE)
  }
  
  return(list(pca = pca_result, df = pca_df, var = var_explained))
}

# ===============================================
# t-SNE FUNCTIONS
# ===============================================

run_tsne_analysis <- function(data_matrix, metadata, output_dir, gene_group) {
  if (!GENERATE_DIM_FIGURES$tsne_scatter) return(NULL)
  
  cat("    Running t-SNE\n")
  
  # Transpose: samples as rows
  data_t <- t(data_matrix)
  
  # Adjust perplexity for small datasets
  perp <- min(TSNE_PERPLEXITY, floor((nrow(data_t) - 1) / 3))
  if (perp < 2) {
    cat("    Skipped t-SNE: too few samples\n")
    return(NULL)
  }
  
  # Run t-SNE
  set.seed(TSNE_SEED)
  tsne_result <- tryCatch({
    Rtsne(data_t, dims = 2, perplexity = perp, 
          max_iter = TSNE_MAX_ITER, check_duplicates = FALSE, pca = TRUE)
  }, error = function(e) {
    cat("    t-SNE failed:", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(tsne_result)) return(NULL)
  
  # Create data frame
  tsne_df <- data.frame(
    tSNE1 = tsne_result$Y[, 1],
    tSNE2 = tsne_result$Y[, 2],
    Sample = rownames(data_t)
  )
  tsne_df <- merge(tsne_df, metadata, by = "Sample")
  
  # Plot
  p <- ggplot(tsne_df, aes(x = tSNE1, y = tSNE2, color = TissueGroup, shape = TissueGroup)) +
    geom_point(size = 4, alpha = 0.8) +
    geom_text_repel(aes(label = Tissue), size = 2.5, max.overlaps = 20) +
    scale_color_brewer(palette = "Set1") +
    labs(title = paste0(gene_group, " - t-SNE (perplexity=", perp, ")"),
         x = "t-SNE 1", y = "t-SNE 2") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  ggsave(file.path(output_dir, paste0(gene_group, "_tSNE.png")),
         p, width = 12, height = 10, dpi = 150)
  
  write.table(tsne_df,
              file = file.path(output_dir, paste0(gene_group, "_tSNE_coordinates.tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  return(tsne_df)
}

# ===============================================
# UMAP FUNCTIONS
# ===============================================

run_umap_analysis <- function(data_matrix, metadata, output_dir, gene_group) {
  if (!GENERATE_DIM_FIGURES$umap_scatter) return(NULL)
  
  cat("    Running UMAP\n")
  
  # Transpose: samples as rows
  data_t <- t(data_matrix)
  
  # Adjust n_neighbors for small datasets
  n_neighbors <- min(UMAP_N_NEIGHBORS, nrow(data_t) - 1)
  if (n_neighbors < 2) {
    cat("    Skipped UMAP: too few samples\n")
    return(NULL)
  }
  
  # Run UMAP
  set.seed(UMAP_SEED)
  umap_config <- umap.defaults
  umap_config$n_neighbors <- n_neighbors
  umap_config$min_dist <- UMAP_MIN_DIST
  umap_config$metric <- UMAP_METRIC
  
  umap_result <- tryCatch({
    umap(data_t, config = umap_config)
  }, error = function(e) {
    cat("    UMAP failed:", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(umap_result)) return(NULL)
  
  # Create data frame
  umap_df <- data.frame(
    UMAP1 = umap_result$layout[, 1],
    UMAP2 = umap_result$layout[, 2],
    Sample = rownames(data_t)
  )
  umap_df <- merge(umap_df, metadata, by = "Sample")
  
  # Plot
  p <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = TissueGroup, shape = TissueGroup)) +
    geom_point(size = 4, alpha = 0.8) +
    geom_text_repel(aes(label = Tissue), size = 2.5, max.overlaps = 20) +
    scale_color_brewer(palette = "Set1") +
    labs(title = paste0(gene_group, " - UMAP"),
         x = "UMAP 1", y = "UMAP 2") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  ggsave(file.path(output_dir, paste0(gene_group, "_UMAP.png")),
         p, width = 12, height = 10, dpi = 150)
  
  write.table(umap_df,
              file = file.path(output_dir, paste0(gene_group, "_UMAP_coordinates.tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  return(umap_df)
}

# ===============================================
# MAIN EXECUTION
# ===============================================

run_dimensionality_reduction <- function() {
  print_separator()
  cat("DIMENSIONALITY REDUCTION ANALYSIS\n")
  print_separator()
  
  for (gene_group in config$gene_groups) {
    cat("\nProcessing gene group:", gene_group, "\n")
    
    group_out_dir <- file.path(DIM_REDUCTION_OUT_DIR, gene_group)
    ensure_output_dir(group_out_dir)
    
    # Load count matrix
    input_file <- build_input_path(gene_group, "gene_level", "expected_count", "geneName")
    
    result <- validate_and_read_matrix(input_file, min_rows = 5)
    if (!result$success) {
      cat("  Skipped:", result$reason, "\n")
      next
    }
    
    count_matrix <- result$data
    cat("  Loaded matrix:", nrow(count_matrix), "genes x", ncol(count_matrix), "samples\n")
    
    # Normalize (log2 CPM)
    data_norm <- preprocess_for_cpm(count_matrix)
    
    # Select variable genes
    data_var <- select_variable_genes(data_norm, N_TOP_VARIABLE_GENES)
    cat("  Using", nrow(data_var), "variable genes\n")
    
    # Prepare metadata
    metadata <- prepare_sample_metadata(colnames(data_var))
    
    # Run analyses
    pca_result <- run_pca_analysis(data_var, metadata, group_out_dir, gene_group)
    tsne_result <- run_tsne_analysis(data_var, metadata, group_out_dir, gene_group)
    umap_result <- run_umap_analysis(data_var, metadata, group_out_dir, gene_group)
    
    cat("  Completed dimensionality reduction\n")
  }
  
  print_separator()
  cat("Dimensionality reduction complete. Output:", DIM_REDUCTION_OUT_DIR, "\n")
  print_separator()
}

# Run if executed directly
if (!interactive()) {
  run_dimensionality_reduction()
}
