#!/usr/bin/env Rscript

# ===============================================
# SAMPLE CORRELATION & CLUSTERING FOR METHOD 3 (STAR + SALMON)
# ===============================================
# Computes sample-sample correlations and hierarchical clustering
# Identifies outliers and batch effects

suppressPackageStartupMessages({
  library(pheatmap)
  library(ggplot2)
  library(RColorBrewer)
  library(dendextend)
  library(corrplot)
})

script_dir <- dirname(sys.frame(1)$ofile)
source(file.path(script_dir, "0_shared_config.R"))
source(file.path(script_dir, "2_utility_functions.R"))

# ===============================================
# CONFIGURATION
# ===============================================

CORRELATION_SUBDIR <- "VIII_Sample_Correlation"
CORRELATION_OUT_DIR <- file.path(CONSOLIDATED_BASE_DIR, CORRELATION_SUBDIR)

# Correlation parameters
CORRELATION_METHOD <- "spearman"   # "pearson" or "spearman"
CLUSTERING_METHOD <- "complete"    # "complete", "average", "ward.D2"
DISTANCE_METHOD <- "euclidean"     # "euclidean", "manhattan", "correlation"

# Outlier detection
OUTLIER_SD_THRESHOLD <- 2.5        # Samples with mean correlation < mean - threshold*SD
MIN_CORRELATION <- 0.7             # Minimum acceptable correlation

# Bootstrap resampling for cluster stability
BOOTSTRAP_ITERATIONS <- 100        # 0 to skip bootstrap analysis

# Load runtime configuration
config <- load_runtime_config()
ensure_output_dir(CORRELATION_OUT_DIR)

# FIGURE GENERATION CONTROL
GENERATE_CORR_FIGURES <- list(
  correlation_heatmap = TRUE,        # Sample-sample correlation heatmap
  correlation_matrix_plot = TRUE,    # corrplot style visualization
  dendrogram = TRUE,                 # Hierarchical clustering dendrogram
  distance_heatmap = TRUE,           # Sample distance heatmap
  outlier_barplot = TRUE,            # Mean correlation per sample
  silhouette_plot = FALSE            # Cluster silhouette analysis
)

# ===============================================
# CORRELATION FUNCTIONS
# ===============================================

# Compute sample correlation matrix
compute_sample_correlation <- function(data_matrix, method = CORRELATION_METHOD) {
  cor_matrix <- cor(data_matrix, method = method, use = "pairwise.complete.obs")
  cor_matrix[is.na(cor_matrix)] <- 0
  return(cor_matrix)
}

# Compute sample distance matrix
compute_sample_distance <- function(data_matrix, method = DISTANCE_METHOD) {
  if (method == "correlation") {
    cor_matrix <- cor(data_matrix, method = CORRELATION_METHOD, use = "pairwise.complete.obs")
    dist_matrix <- as.dist(1 - cor_matrix)
  } else {
    dist_matrix <- dist(t(data_matrix), method = method)
  }
  return(dist_matrix)
}

# Create correlation heatmap
create_correlation_heatmap <- function(cor_matrix, metadata, output_dir, gene_group) {
  if (!GENERATE_CORR_FIGURES$correlation_heatmap) return(NULL)
  
  # Annotation
  annotation_row <- data.frame(
    Tissue = metadata$TissueGroup,
    row.names = metadata$Sample
  )
  
  # Color palette
  col_palette <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)
  
  # Tissue group colors
  tissue_colors <- brewer.pal(length(unique(metadata$TissueGroup)), "Set1")
  names(tissue_colors) <- unique(metadata$TissueGroup)
  annotation_colors <- list(Tissue = tissue_colors)
  
  png(file.path(output_dir, paste0(gene_group, "_correlation_heatmap.png")),
      width = 1000, height = 1000, res = 100)
  
  pheatmap(cor_matrix,
           clustering_method = CLUSTERING_METHOD,
           annotation_row = annotation_row,
           annotation_col = annotation_row,
           annotation_colors = annotation_colors,
           color = col_palette,
           main = paste0(gene_group, " - Sample Correlation (", CORRELATION_METHOD, ")"),
           fontsize = 8,
           display_numbers = nrow(cor_matrix) <= 15,
           number_format = "%.2f",
           number_color = "black")
  
  dev.off()
  
  # Export correlation matrix
  cor_df <- as.data.frame(cor_matrix)
  cor_df$Sample <- rownames(cor_df)
  cor_df <- cor_df[, c("Sample", setdiff(names(cor_df), "Sample"))]
  write.table(cor_df,
              file = file.path(output_dir, paste0(gene_group, "_correlation_matrix.tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)
}

# Create corrplot style visualization
create_corrplot <- function(cor_matrix, output_dir, gene_group) {
  if (!GENERATE_CORR_FIGURES$correlation_matrix_plot) return(NULL)
  
  png(file.path(output_dir, paste0(gene_group, "_corrplot.png")),
      width = 1000, height = 1000, res = 100)
  
  corrplot(cor_matrix, method = "color", type = "upper",
           order = "hclust", hclust.method = CLUSTERING_METHOD,
           tl.col = "black", tl.srt = 45, tl.cex = 0.7,
           addCoef.col = if (nrow(cor_matrix) <= 15) "black" else NULL,
           number.cex = 0.6,
           col = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
           title = paste0(gene_group, " - Sample Correlation"),
           mar = c(0, 0, 2, 0))
  
  dev.off()
}

# Create hierarchical clustering dendrogram
create_dendrogram <- function(data_matrix, metadata, output_dir, gene_group) {
  if (!GENERATE_CORR_FIGURES$dendrogram) return(NULL)
  
  # Compute distance and cluster
  dist_matrix <- compute_sample_distance(data_matrix)
  hc <- hclust(dist_matrix, method = CLUSTERING_METHOD)
  
  # Create dendrogram with colored labels
  dend <- as.dendrogram(hc)
  
  # Color by tissue group
  tissue_groups <- metadata$TissueGroup[match(labels(dend), metadata$Sample)]
  unique_groups <- unique(tissue_groups)
  group_colors <- brewer.pal(max(3, length(unique_groups)), "Set1")[1:length(unique_groups)]
  names(group_colors) <- unique_groups
  label_colors <- group_colors[tissue_groups]
  
  dend <- color_labels(dend, col = label_colors)
  dend <- set(dend, "labels_cex", 0.7)
  
  png(file.path(output_dir, paste0(gene_group, "_dendrogram.png")),
      width = 1200, height = 600, res = 100)
  
  par(mar = c(10, 4, 4, 2))
  plot(dend, main = paste0(gene_group, " - Hierarchical Clustering"))
  
  # Add legend
  legend("topright", legend = names(group_colors), 
         fill = group_colors, title = "Tissue Group", cex = 0.8)
  
  dev.off()
  
  return(hc)
}

# Create distance heatmap
create_distance_heatmap <- function(data_matrix, metadata, output_dir, gene_group) {
  if (!GENERATE_CORR_FIGURES$distance_heatmap) return(NULL)
  
  dist_matrix <- compute_sample_distance(data_matrix)
  dist_matrix_full <- as.matrix(dist_matrix)
  
  annotation_row <- data.frame(
    Tissue = metadata$TissueGroup,
    row.names = metadata$Sample
  )
  
  png(file.path(output_dir, paste0(gene_group, "_distance_heatmap.png")),
      width = 1000, height = 1000, res = 100)
  
  pheatmap(dist_matrix_full,
           clustering_method = CLUSTERING_METHOD,
           annotation_row = annotation_row,
           annotation_col = annotation_row,
           color = colorRampPalette(c("white", "#2166AC"))(100),
           main = paste0(gene_group, " - Sample Distance (", DISTANCE_METHOD, ")"),
           fontsize = 8)
  
  dev.off()
}

# Detect outliers based on mean correlation
detect_outliers <- function(cor_matrix, output_dir, gene_group) {
  # Mean correlation per sample (excluding self)
  diag(cor_matrix) <- NA
  mean_cors <- rowMeans(cor_matrix, na.rm = TRUE)
  
  # Outlier detection
  global_mean <- mean(mean_cors, na.rm = TRUE)
  global_sd <- sd(mean_cors, na.rm = TRUE)
  outlier_threshold <- global_mean - OUTLIER_SD_THRESHOLD * global_sd
  
  outlier_df <- data.frame(
    Sample = names(mean_cors),
    Mean_Correlation = mean_cors,
    Is_Outlier = mean_cors < outlier_threshold | mean_cors < MIN_CORRELATION,
    stringsAsFactors = FALSE
  )
  outlier_df <- outlier_df[order(outlier_df$Mean_Correlation), ]
  
  # Export
  write.table(outlier_df,
              file = file.path(output_dir, paste0(gene_group, "_sample_correlations.tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Bar plot
  if (GENERATE_CORR_FIGURES$outlier_barplot) {
    p <- ggplot(outlier_df, aes(x = reorder(Sample, Mean_Correlation), 
                                 y = Mean_Correlation,
                                 fill = Is_Outlier)) +
      geom_bar(stat = "identity") +
      geom_hline(yintercept = outlier_threshold, linetype = "dashed", color = "red") +
      geom_hline(yintercept = MIN_CORRELATION, linetype = "dotted", color = "orange") +
      scale_fill_manual(values = c("FALSE" = "#2166AC", "TRUE" = "#B2182B")) +
      coord_flip() +
      labs(title = paste0(gene_group, " - Mean Sample Correlation"),
           x = "Sample", y = paste0("Mean ", CORRELATION_METHOD, " Correlation"),
           fill = "Potential Outlier") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    
    ggsave(file.path(output_dir, paste0(gene_group, "_mean_correlation_barplot.png")),
           p, width = 10, height = 8, dpi = 150)
  }
  
  # Report outliers
  n_outliers <- sum(outlier_df$Is_Outlier)
  if (n_outliers > 0) {
    cat("    Potential outliers:", paste(outlier_df$Sample[outlier_df$Is_Outlier], collapse = ", "), "\n")
  }
  
  return(outlier_df)
}

# ===============================================
# CLUSTER STABILITY ANALYSIS
# ===============================================

assess_cluster_stability <- function(data_matrix, n_bootstrap = BOOTSTRAP_ITERATIONS,
                                      output_dir, gene_group) {
  if (n_bootstrap <= 0) return(NULL)
  
  cat("    Running bootstrap cluster stability (", n_bootstrap, " iterations)\n")
  
  n_genes <- nrow(data_matrix)
  original_hc <- hclust(compute_sample_distance(data_matrix), method = CLUSTERING_METHOD)
  original_order <- original_hc$order
  
  # Bootstrap
  cophenetic_cors <- numeric(n_bootstrap)
  
  for (i in 1:n_bootstrap) {
    # Resample genes with replacement
    boot_idx <- sample(1:n_genes, n_genes, replace = TRUE)
    boot_data <- data_matrix[boot_idx, ]
    
    boot_hc <- hclust(compute_sample_distance(boot_data), method = CLUSTERING_METHOD)
    
    # Cophenetic correlation
    cophenetic_cors[i] <- cor(cophenetic(original_hc), cophenetic(boot_hc), method = "spearman")
  }
  
  stability_stats <- data.frame(
    Mean_Cophenetic_Cor = mean(cophenetic_cors, na.rm = TRUE),
    SD_Cophenetic_Cor = sd(cophenetic_cors, na.rm = TRUE),
    CI_Lower = quantile(cophenetic_cors, 0.025, na.rm = TRUE),
    CI_Upper = quantile(cophenetic_cors, 0.975, na.rm = TRUE)
  )
  
  write.table(stability_stats,
              file = file.path(output_dir, paste0(gene_group, "_cluster_stability.tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  cat("    Cluster stability (cophenetic cor):", round(stability_stats$Mean_Cophenetic_Cor, 3), "\n")
  
  return(stability_stats)
}

# ===============================================
# MAIN EXECUTION
# ===============================================

run_sample_correlation_analysis <- function() {
  print_separator()
  cat("SAMPLE CORRELATION & CLUSTERING ANALYSIS\n")
  print_separator()
  
  for (gene_group in config$gene_groups) {
    cat("\nProcessing gene group:", gene_group, "\n")
    
    group_out_dir <- file.path(CORRELATION_OUT_DIR, gene_group)
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
    
    # Normalize
    data_norm <- preprocess_for_cpm(count_matrix)
    
    # Prepare metadata
    metadata <- data.frame(
      Sample = colnames(data_norm),
      Tissue = SAMPLE_LABELS[colnames(data_norm)],
      stringsAsFactors = FALSE
    )
    metadata$Tissue[is.na(metadata$Tissue)] <- metadata$Sample[is.na(metadata$Tissue)]
    metadata$TissueGroup <- sapply(metadata$Tissue, function(t) {
      if (grepl("Root|Stem|Leaf|Leaves", t, ignore.case = TRUE)) return("Vegetative")
      if (grepl("Flower|Bud|Pistil", t, ignore.case = TRUE)) return("Reproductive")
      if (grepl("Fruit", t, ignore.case = TRUE)) return("Fruit")
      if (grepl("Radicle|Cotyledon|Seed", t, ignore.case = TRUE)) return("Seedling")
      return("Other")
    })
    
    # Compute correlation
    cor_matrix <- compute_sample_correlation(data_norm)
    
    # Generate outputs
    create_correlation_heatmap(cor_matrix, metadata, group_out_dir, gene_group)
    create_corrplot(cor_matrix, group_out_dir, gene_group)
    create_dendrogram(data_norm, metadata, group_out_dir, gene_group)
    create_distance_heatmap(data_norm, metadata, group_out_dir, gene_group)
    detect_outliers(cor_matrix, group_out_dir, gene_group)
    assess_cluster_stability(data_norm, BOOTSTRAP_ITERATIONS, group_out_dir, gene_group)
    
    cat("  Completed correlation analysis\n")
  }
  
  print_separator()
  cat("Correlation analysis complete. Output:", CORRELATION_OUT_DIR, "\n")
  print_separator()
}

# Run if executed directly
if (!interactive()) {
  run_sample_correlation_analysis()
}
