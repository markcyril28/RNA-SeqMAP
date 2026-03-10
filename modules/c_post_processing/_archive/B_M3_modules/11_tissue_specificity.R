#!/usr/bin/env Rscript

# ===============================================
# TISSUE SPECIFICITY ANALYSIS FOR METHOD 3 (STAR + SALMON)
# ===============================================
# Calculates tissue specificity scores (Tau, entropy, etc.)
# Identifies tissue-enriched and housekeeping genes

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
  library(dplyr)
  library(tidyr)
})

script_dir <- dirname(sys.frame(1)$ofile)
source(file.path(script_dir, "0_shared_config.R"))
source(file.path(script_dir, "2_utility_functions.R"))

# ===============================================
# CONFIGURATION
# ===============================================

TISSUE_SPEC_SUBDIR <- "IX_Tissue_Specificity"
TISSUE_SPEC_OUT_DIR <- file.path(CONSOLIDATED_BASE_DIR, TISSUE_SPEC_SUBDIR)

# Tau index thresholds
TAU_TISSUE_SPECIFIC <- 0.85        # Tau >= this = tissue-specific
TAU_HOUSEKEEPING <- 0.25           # Tau <= this = housekeeping/ubiquitous
TAU_MODERATE <- 0.5                # Intermediate specificity threshold

# Expression thresholds
MIN_EXPRESSION <- 1                # Minimum expression (CPM) to consider expressed
ENRICHMENT_FOLD <- 2               # Fold enrichment to call tissue-enriched

# Number of top genes to report per tissue
N_TOP_GENES_PER_TISSUE <- 20

# Load runtime configuration
config <- load_runtime_config()
ensure_output_dir(TISSUE_SPEC_OUT_DIR)

# FIGURE GENERATION CONTROL
GENERATE_TISSUE_SPEC_FIGURES <- list(
  tau_histogram = TRUE,              # Distribution of Tau values
  tau_vs_expression = TRUE,          # Tau vs mean expression scatter
  tissue_enriched_heatmap = TRUE,    # Heatmap of top tissue-enriched genes
  specificity_barplot = TRUE,        # Genes per tissue category
  expression_profile = TRUE          # Expression profiles of top specific genes
)

# ===============================================
# TISSUE SPECIFICITY METRICS
# ===============================================

# Calculate Tau index (Yanai et al. 2005)
# Tau ranges from 0 (ubiquitous) to 1 (tissue-specific)
calculate_tau <- function(expression_vector) {
  n <- length(expression_vector)
  if (n < 2) return(NA)
  
  # Normalize to max expression
  max_expr <- max(expression_vector, na.rm = TRUE)
  if (max_expr == 0 || is.na(max_expr)) return(NA)
  
  x <- expression_vector / max_expr
  x[is.na(x)] <- 0
  
  # Tau formula
  tau <- sum(1 - x) / (n - 1)
  return(tau)
}

# Calculate Shannon entropy-based specificity
# Lower entropy = more specific
calculate_entropy_specificity <- function(expression_vector) {
  total <- sum(expression_vector, na.rm = TRUE)
  if (total == 0) return(NA)
  
  p <- expression_vector / total
  p[p == 0] <- NA  # Avoid log(0)
  
  H <- -sum(p * log2(p), na.rm = TRUE)
  H_max <- log2(length(expression_vector))
  
  # Normalize to 0-1 (1 = specific, 0 = ubiquitous)
  specificity <- 1 - (H / H_max)
  return(specificity)
}

# Calculate Gini coefficient
calculate_gini <- function(expression_vector) {
  x <- sort(expression_vector[!is.na(expression_vector)])
  n <- length(x)
  if (n == 0 || sum(x) == 0) return(NA)
  
  gini <- (2 * sum(seq_along(x) * x) / (n * sum(x))) - ((n + 1) / n)
  return(gini)
}

# Calculate all specificity metrics
calculate_all_specificity_metrics <- function(data_matrix) {
  cat("  Calculating specificity metrics\n")
  
  metrics_df <- data.frame(
    Gene = rownames(data_matrix),
    Mean_Expression = rowMeans(data_matrix, na.rm = TRUE),
    Max_Expression = apply(data_matrix, 1, max, na.rm = TRUE),
    N_Tissues_Expressed = apply(data_matrix, 1, function(x) sum(x >= MIN_EXPRESSION, na.rm = TRUE)),
    Tau = apply(data_matrix, 1, calculate_tau),
    Entropy_Specificity = apply(data_matrix, 1, calculate_entropy_specificity),
    Gini = apply(data_matrix, 1, calculate_gini),
    stringsAsFactors = FALSE
  )
  
  # Identify max tissue
  metrics_df$Max_Tissue <- colnames(data_matrix)[apply(data_matrix, 1, which.max)]
  
  # Classify genes
  metrics_df$Specificity_Class <- "Moderate"
  metrics_df$Specificity_Class[metrics_df$Tau >= TAU_TISSUE_SPECIFIC] <- "Tissue_Specific"
  metrics_df$Specificity_Class[metrics_df$Tau <= TAU_HOUSEKEEPING] <- "Housekeeping"
  
  return(metrics_df)
}

# ===============================================
# TISSUE ENRICHMENT ANALYSIS
# ===============================================

# Identify tissue-enriched genes (expression in tissue X >> other tissues)
identify_tissue_enriched_genes <- function(data_matrix, fold_threshold = ENRICHMENT_FOLD) {
  n_tissues <- ncol(data_matrix)
  enriched_list <- list()
  
  for (tissue in colnames(data_matrix)) {
    tissue_expr <- data_matrix[, tissue]
    other_mean <- rowMeans(data_matrix[, colnames(data_matrix) != tissue, drop = FALSE], na.rm = TRUE)
    other_mean[other_mean == 0] <- 0.01  # Avoid division by zero
    
    fold_enrichment <- tissue_expr / other_mean
    
    # Filter: expressed in tissue AND enriched over others
    enriched_idx <- which(tissue_expr >= MIN_EXPRESSION & fold_enrichment >= fold_threshold)
    
    if (length(enriched_idx) > 0) {
      enriched_df <- data.frame(
        Gene = rownames(data_matrix)[enriched_idx],
        Tissue = tissue,
        Expression = tissue_expr[enriched_idx],
        Fold_Enrichment = fold_enrichment[enriched_idx],
        stringsAsFactors = FALSE
      )
      enriched_df <- enriched_df[order(-enriched_df$Fold_Enrichment), ]
      enriched_list[[tissue]] <- enriched_df
    }
  }
  
  all_enriched <- do.call(rbind, enriched_list)
  if (!is.null(all_enriched)) {
    rownames(all_enriched) <- NULL
  }
  
  return(all_enriched)
}

# ===============================================
# VISUALIZATION FUNCTIONS
# ===============================================

create_tau_histogram <- function(metrics_df, output_dir, gene_group) {
  if (!GENERATE_TISSUE_SPEC_FIGURES$tau_histogram) return(NULL)
  
  p <- ggplot(metrics_df, aes(x = Tau, fill = Specificity_Class)) +
    geom_histogram(bins = 50, alpha = 0.8, color = "white") +
    geom_vline(xintercept = c(TAU_HOUSEKEEPING, TAU_TISSUE_SPECIFIC), 
               linetype = "dashed", color = c("blue", "red")) +
    scale_fill_manual(values = c("Housekeeping" = "#2166AC", 
                                  "Moderate" = "grey60", 
                                  "Tissue_Specific" = "#B2182B")) +
    labs(title = paste0(gene_group, " - Tau Index Distribution"),
         x = "Tau (Tissue Specificity Index)",
         y = "Number of Genes",
         fill = "Classification") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  ggsave(file.path(output_dir, paste0(gene_group, "_tau_histogram.png")),
         p, width = 10, height = 6, dpi = 150)
}

create_tau_vs_expression <- function(metrics_df, output_dir, gene_group) {
  if (!GENERATE_TISSUE_SPEC_FIGURES$tau_vs_expression) return(NULL)
  
  plot_data <- metrics_df[metrics_df$Mean_Expression > 0, ]
  
  p <- ggplot(plot_data, aes(x = log10(Mean_Expression + 1), y = Tau, 
                              color = Specificity_Class)) +
    geom_point(alpha = 0.5, size = 1.5) +
    geom_hline(yintercept = c(TAU_HOUSEKEEPING, TAU_TISSUE_SPECIFIC), 
               linetype = "dashed", color = "grey40") +
    scale_color_manual(values = c("Housekeeping" = "#2166AC", 
                                   "Moderate" = "grey60", 
                                   "Tissue_Specific" = "#B2182B")) +
    labs(title = paste0(gene_group, " - Tau vs Expression"),
         x = "Log10(Mean Expression + 1)",
         y = "Tau (Tissue Specificity)",
         color = "Classification") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  ggsave(file.path(output_dir, paste0(gene_group, "_tau_vs_expression.png")),
         p, width = 10, height = 8, dpi = 150)
}

create_tissue_enriched_heatmap <- function(data_matrix, enriched_df, output_dir, 
                                            gene_group, n_per_tissue = 5) {
  if (!GENERATE_TISSUE_SPEC_FIGURES$tissue_enriched_heatmap) return(NULL)
  if (is.null(enriched_df) || nrow(enriched_df) == 0) return(NULL)
  
  # Select top genes per tissue
  top_genes <- enriched_df %>%
    group_by(Tissue) %>%
    slice_max(order_by = Fold_Enrichment, n = n_per_tissue) %>%
    pull(Gene) %>%
    unique()
  
  if (length(top_genes) < 2) return(NULL)
  
  # Subset and normalize
  heatmap_data <- data_matrix[top_genes, , drop = FALSE]
  heatmap_data <- t(scale(t(log2(heatmap_data + 1))))
  heatmap_data[is.na(heatmap_data)] <- 0
  
  png(file.path(output_dir, paste0(gene_group, "_tissue_enriched_heatmap.png")),
      width = 1000, height = 1200, res = 100)
  
  pheatmap(heatmap_data,
           cluster_rows = TRUE, cluster_cols = TRUE,
           show_rownames = nrow(heatmap_data) <= 50,
           main = paste0(gene_group, " - Top Tissue-Enriched Genes"),
           color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
           fontsize_row = 6)
  
  dev.off()
}

create_specificity_barplot <- function(metrics_df, output_dir, gene_group) {
  if (!GENERATE_TISSUE_SPEC_FIGURES$specificity_barplot) return(NULL)
  
  # Count genes per tissue for tissue-specific genes
  tissue_counts <- metrics_df %>%
    filter(Specificity_Class == "Tissue_Specific") %>%
    count(Max_Tissue, name = "N_Genes") %>%
    arrange(desc(N_Genes))
  
  if (nrow(tissue_counts) == 0) return(NULL)
  
  p <- ggplot(tissue_counts, aes(x = reorder(Max_Tissue, N_Genes), y = N_Genes)) +
    geom_bar(stat = "identity", fill = "#B2182B", alpha = 0.8) +
    coord_flip() +
    labs(title = paste0(gene_group, " - Tissue-Specific Genes per Tissue"),
         x = "Tissue", y = "Number of Tissue-Specific Genes") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  ggsave(file.path(output_dir, paste0(gene_group, "_tissue_specific_counts.png")),
         p, width = 10, height = 8, dpi = 150)
}

create_expression_profiles <- function(data_matrix, metrics_df, output_dir, 
                                        gene_group, n_genes = 6) {
  if (!GENERATE_TISSUE_SPEC_FIGURES$expression_profile) return(NULL)
  
  # Select top tissue-specific genes
  top_specific <- metrics_df %>%
    filter(Specificity_Class == "Tissue_Specific") %>%
    arrange(desc(Tau)) %>%
    head(n_genes) %>%
    pull(Gene)
  
  if (length(top_specific) == 0) return(NULL)
  
  # Prepare long format data
  plot_data <- data_matrix[top_specific, , drop = FALSE]
  plot_data <- as.data.frame(plot_data)
  plot_data$Gene <- rownames(plot_data)
  plot_data_long <- pivot_longer(plot_data, cols = -Gene, 
                                  names_to = "Tissue", values_to = "Expression")
  
  p <- ggplot(plot_data_long, aes(x = Tissue, y = Expression, fill = Tissue)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    facet_wrap(~ Gene, scales = "free_y", ncol = 2) +
    labs(title = paste0(gene_group, " - Top Tissue-Specific Gene Profiles"),
         x = "", y = "Expression") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
          legend.position = "none")
  
  ggsave(file.path(output_dir, paste0(gene_group, "_expression_profiles.png")),
         p, width = 12, height = 10, dpi = 150)
}

# ===============================================
# MAIN EXECUTION
# ===============================================

run_tissue_specificity_analysis <- function() {
  print_separator()
  cat("TISSUE SPECIFICITY ANALYSIS\n")
  print_separator()
  
  for (gene_group in config$gene_groups) {
    cat("\nProcessing gene group:", gene_group, "\n")
    
    group_out_dir <- file.path(TISSUE_SPEC_OUT_DIR, gene_group)
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
    
    # Convert to CPM for analysis
    lib_sizes <- colSums(count_matrix, na.rm = TRUE)
    lib_sizes[lib_sizes == 0] <- 1
    cpm_matrix <- sweep(count_matrix, 2, lib_sizes / 1e6, FUN = "/")
    
    # Replace SRR IDs with tissue names if possible
    tissue_names <- SAMPLE_LABELS[colnames(cpm_matrix)]
    tissue_names[is.na(tissue_names)] <- colnames(cpm_matrix)[is.na(tissue_names)]
    colnames(cpm_matrix) <- tissue_names
    
    # Calculate specificity metrics
    metrics_df <- calculate_all_specificity_metrics(cpm_matrix)
    
    # Summary
    n_specific <- sum(metrics_df$Specificity_Class == "Tissue_Specific", na.rm = TRUE)
    n_housekeeping <- sum(metrics_df$Specificity_Class == "Housekeeping", na.rm = TRUE)
    cat("  Tissue-specific genes (Tau >=", TAU_TISSUE_SPECIFIC, "):", n_specific, "\n")
    cat("  Housekeeping genes (Tau <=", TAU_HOUSEKEEPING, "):", n_housekeeping, "\n")
    
    # Export metrics
    metrics_df <- metrics_df[order(-metrics_df$Tau), ]
    write.table(metrics_df,
                file = file.path(group_out_dir, paste0(gene_group, "_tissue_specificity_metrics.tsv")),
                sep = "\t", row.names = FALSE, quote = FALSE)
    
    # Identify tissue-enriched genes
    enriched_df <- identify_tissue_enriched_genes(cpm_matrix)
    if (!is.null(enriched_df)) {
      write.table(enriched_df,
                  file = file.path(group_out_dir, paste0(gene_group, "_tissue_enriched_genes.tsv")),
                  sep = "\t", row.names = FALSE, quote = FALSE)
      cat("  Tissue-enriched genes identified:", nrow(enriched_df), "\n")
    }
    
    # Generate figures
    create_tau_histogram(metrics_df, group_out_dir, gene_group)
    create_tau_vs_expression(metrics_df, group_out_dir, gene_group)
    create_tissue_enriched_heatmap(cpm_matrix, enriched_df, group_out_dir, gene_group)
    create_specificity_barplot(metrics_df, group_out_dir, gene_group)
    create_expression_profiles(cpm_matrix, metrics_df, group_out_dir, gene_group)
    
    cat("  Completed tissue specificity analysis\n")
  }
  
  print_separator()
  cat("Tissue specificity analysis complete. Output:", TISSUE_SPEC_OUT_DIR, "\n")
  print_separator()
}

# Run if executed directly
if (!interactive()) {
  run_tissue_specificity_analysis()
}
