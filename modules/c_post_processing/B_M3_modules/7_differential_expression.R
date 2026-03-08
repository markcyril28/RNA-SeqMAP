#!/usr/bin/env Rscript

# ===============================================
# DIFFERENTIAL EXPRESSION ANALYSIS FOR METHOD 3 (STAR + SALMON)
# ===============================================
# Performs DESeq2-based differential expression analysis
# Generates volcano plots, MA plots, and significant gene tables

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
})

script_dir <- dirname(sys.frame(1)$ofile)
source(file.path(script_dir, "0_shared_config.R"))
source(file.path(script_dir, "2_utility_functions.R"))

# ===============================================
# CONFIGURATION
# ===============================================

DEA_SUBDIR <- "V_Differential_Expression"
DEA_OUT_DIR <- file.path(CONSOLIDATED_BASE_DIR, DEA_SUBDIR)

# DESeq2 parameters
PADJ_THRESHOLD <- 0.05             # Adjusted p-value cutoff for significance
LFC_THRESHOLD <- 1.0               # Log2 fold change threshold (1 = 2-fold change)
SHRINKAGE_TYPE <- "apeglm"         # LFC shrinkage method: "apeglm", "ashr", or "normal"
MIN_COUNT_FILTER <- 10             # Minimum total counts per gene across samples
INDEPENDENT_FILTERING <- TRUE      # Use DESeq2 independent filtering

# Comparison configuration - define contrasts as list of (condition1, condition2)
# Will compare condition1 vs condition2 (positive LFC = higher in condition1)
TISSUE_GROUPS <- list(
  "Vegetative" = c("Root", "Stem", "Mature_Leaves", "Young_Leaves"),
  "Reproductive" = c("Flower_Buds", "Flowers", "Opened_Buds", "Pistils"),
  "Fruit" = c("Young_Fruits", "Mature_Fruits"),
  "Seedling" = c("Radicles", "Cotyledons")
)

# Load runtime configuration
config <- load_runtime_config()
ensure_output_dir(DEA_OUT_DIR)

# FIGURE GENERATION CONTROL
GENERATE_DEA_FIGURES <- list(
  volcano_plot = TRUE,               # Volcano plot (LFC vs -log10 p-value)
  ma_plot = TRUE,                    # MA plot (mean expression vs LFC)
  top_genes_heatmap = TRUE,          # Heatmap of top DE genes
  pvalue_histogram = FALSE           # P-value distribution diagnostic
)

# ===============================================
# DEA CORE FUNCTIONS
# ===============================================

# Prepare DESeq2 dataset from count matrix
prepare_deseq_dataset <- function(count_matrix, sample_info) {
  # Ensure counts are integers
  count_matrix <- round(count_matrix)
  storage.mode(count_matrix) <- "integer"
  
  # Filter low-count genes
  keep <- rowSums(count_matrix) >= MIN_COUNT_FILTER
  count_matrix <- count_matrix[keep, , drop = FALSE]
  
  if (nrow(count_matrix) < 10) {
    cat("    Warning: Too few genes after filtering (", nrow(count_matrix), ")\n")
    return(NULL)
  }
  
  # Create DESeq2 dataset
  dds <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = sample_info,
    design = ~ condition
  )
  
  return(dds)
}

# Run DESeq2 analysis
run_deseq2 <- function(dds, contrast_name, output_dir) {
  # Run differential expression
  dds <- DESeq(dds, quiet = TRUE)
  
  # Get results with shrinkage if possible
  res <- tryCatch({
    if (SHRINKAGE_TYPE == "apeglm") {
      # apeglm requires coefficient name
      coef_name <- resultsNames(dds)[2]
      lfcShrink(dds, coef = coef_name, type = "apeglm", quiet = TRUE)
    } else {
      results(dds, alpha = PADJ_THRESHOLD)
    }
  }, error = function(e) {
    cat("    LFC shrinkage failed, using standard results\n")
    results(dds, alpha = PADJ_THRESHOLD)
  })
  
  # Convert to data frame and add gene names
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  res_df <- res_df[, c("gene", setdiff(names(res_df), "gene"))]
  
  # Add significance classification
  res_df$significance <- "NS"
  res_df$significance[res_df$padj < PADJ_THRESHOLD & res_df$log2FoldChange > LFC_THRESHOLD] <- "Up"
  res_df$significance[res_df$padj < PADJ_THRESHOLD & res_df$log2FoldChange < -LFC_THRESHOLD] <- "Down"
  res_df$significance <- factor(res_df$significance, levels = c("Down", "NS", "Up"))
  
  # Sort by adjusted p-value
  res_df <- res_df[order(res_df$padj), ]
  
  return(list(dds = dds, results = res_df))
}

# Create volcano plot
create_volcano_plot <- function(res_df, contrast_name, output_dir) {
  if (!GENERATE_DEA_FIGURES$volcano_plot) return(NULL)
  
  # Prepare data
  plot_data <- res_df[!is.na(res_df$padj), ]
  plot_data$neglog10p <- -log10(plot_data$padj)
  plot_data$neglog10p[is.infinite(plot_data$neglog10p)] <- max(plot_data$neglog10p[is.finite(plot_data$neglog10p)]) + 10
  
  # Top genes to label
  n_label <- min(10, sum(plot_data$significance != "NS"))
  top_genes <- head(plot_data[plot_data$significance != "NS", ], n_label)
  
  # Count significant genes
  n_up <- sum(plot_data$significance == "Up", na.rm = TRUE)
  n_down <- sum(plot_data$significance == "Down", na.rm = TRUE)
  
  # Create plot
  p <- ggplot(plot_data, aes(x = log2FoldChange, y = neglog10p, color = significance)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("Down" = "#2166AC", "NS" = "grey60", "Up" = "#B2182B"),
                       labels = c(paste0("Down (", n_down, ")"), 
                                  "NS", 
                                  paste0("Up (", n_up, ")"))) +
    geom_vline(xintercept = c(-LFC_THRESHOLD, LFC_THRESHOLD), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(PADJ_THRESHOLD), linetype = "dashed", color = "grey40") +
    labs(title = paste0("Volcano Plot: ", contrast_name),
         x = "Log2 Fold Change",
         y = "-Log10 Adjusted P-value",
         color = "Significance") +
    theme_minimal() +
    theme(legend.position = "bottom",
          plot.title = element_text(hjust = 0.5, face = "bold"))
  
  if (nrow(top_genes) > 0) {
    p <- p + geom_text_repel(data = top_genes, aes(label = gene),
                              size = 2.5, max.overlaps = 15, segment.size = 0.2)
  }
  
  ggsave(file.path(output_dir, paste0(contrast_name, "_volcano.png")),
         p, width = 10, height = 8, dpi = 150)
  
  return(p)
}

# Create MA plot
create_ma_plot <- function(res_df, contrast_name, output_dir) {
  if (!GENERATE_DEA_FIGURES$ma_plot) return(NULL)
  
  plot_data <- res_df[!is.na(res_df$baseMean) & !is.na(res_df$log2FoldChange), ]
  plot_data$log_baseMean <- log10(plot_data$baseMean + 1)
  
  p <- ggplot(plot_data, aes(x = log_baseMean, y = log2FoldChange, color = significance)) +
    geom_point(alpha = 0.5, size = 1.2) +
    scale_color_manual(values = c("Down" = "#2166AC", "NS" = "grey60", "Up" = "#B2182B")) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black") +
    geom_hline(yintercept = c(-LFC_THRESHOLD, LFC_THRESHOLD), linetype = "dashed", color = "grey40") +
    labs(title = paste0("MA Plot: ", contrast_name),
         x = "Log10 Mean Expression",
         y = "Log2 Fold Change",
         color = "Significance") +
    theme_minimal() +
    theme(legend.position = "bottom",
          plot.title = element_text(hjust = 0.5, face = "bold"))
  
  ggsave(file.path(output_dir, paste0(contrast_name, "_MA.png")),
         p, width = 10, height = 8, dpi = 150)
  
  return(p)
}

# Create heatmap of top DE genes
create_de_heatmap <- function(dds, res_df, contrast_name, output_dir, n_top = 50) {
  if (!GENERATE_DEA_FIGURES$top_genes_heatmap) return(NULL)
  
  # Get normalized counts
  vst_data <- tryCatch({
    assay(vst(dds, blind = FALSE))
  }, error = function(e) {
    log2(counts(dds, normalized = TRUE) + 1)
  })
  
  # Select top DE genes
  sig_genes <- res_df[res_df$significance != "NS" & !is.na(res_df$padj), ]
  top_genes <- head(sig_genes$gene, n_top)
  
  if (length(top_genes) < 2) {
    cat("    Too few significant genes for heatmap\n")
    return(NULL)
  }
  
  # Subset and scale
  heatmap_data <- vst_data[top_genes, , drop = FALSE]
  heatmap_data <- t(scale(t(heatmap_data)))
  
  # Annotation
  annotation_col <- data.frame(
    Condition = colData(dds)$condition,
    row.names = colnames(heatmap_data)
  )
  
  png(file.path(output_dir, paste0(contrast_name, "_top_DE_heatmap.png")),
      width = 1000, height = 1200, res = 100)
  
  pheatmap(heatmap_data,
           annotation_col = annotation_col,
           cluster_rows = TRUE, cluster_cols = TRUE,
           show_rownames = nrow(heatmap_data) <= 50,
           main = paste0("Top ", min(n_top, nrow(heatmap_data)), " DE Genes: ", contrast_name),
           fontsize_row = 6)
  
  dev.off()
}

# Export significant genes
export_significant_genes <- function(res_df, contrast_name, output_dir) {
  # Full results
  write.table(res_df, 
              file = file.path(output_dir, paste0(contrast_name, "_all_results.tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Significant only
  sig_genes <- res_df[res_df$significance != "NS" & !is.na(res_df$padj), ]
  write.table(sig_genes,
              file = file.path(output_dir, paste0(contrast_name, "_significant_genes.tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Summary
  summary_stats <- data.frame(
    Comparison = contrast_name,
    Total_Genes = nrow(res_df),
    Significant_Total = sum(res_df$significance != "NS", na.rm = TRUE),
    Upregulated = sum(res_df$significance == "Up", na.rm = TRUE),
    Downregulated = sum(res_df$significance == "Down", na.rm = TRUE),
    Padj_Threshold = PADJ_THRESHOLD,
    LFC_Threshold = LFC_THRESHOLD
  )
  
  return(summary_stats)
}

# ===============================================
# PAIRWISE COMPARISON FUNCTION
# ===============================================

run_pairwise_comparison <- function(count_matrix, condition1_samples, condition2_samples,
                                     condition1_name, condition2_name, output_dir) {
  contrast_name <- paste0(condition1_name, "_vs_", condition2_name)
  cat("  Processing:", contrast_name, "\n")
  
  # Subset samples
  all_samples <- c(condition1_samples, condition2_samples)
  available_samples <- intersect(all_samples, colnames(count_matrix))
  
  if (length(available_samples) < 4) {
    cat("    Skipped: insufficient samples (", length(available_samples), ")\n")
    return(NULL)
  }
  
  sub_matrix <- count_matrix[, available_samples, drop = FALSE]
  
  # Create sample info
  sample_info <- data.frame(
    row.names = available_samples,
    condition = factor(ifelse(available_samples %in% condition1_samples, 
                               condition1_name, condition2_name),
                        levels = c(condition2_name, condition1_name))
  )
  
  # Prepare DESeq dataset
  dds <- prepare_deseq_dataset(sub_matrix, sample_info)
  if (is.null(dds)) return(NULL)
  
  # Run analysis
  result <- run_deseq2(dds, contrast_name, output_dir)
  
  # Generate outputs
  create_volcano_plot(result$results, contrast_name, output_dir)
  create_ma_plot(result$results, contrast_name, output_dir)
  create_de_heatmap(result$dds, result$results, contrast_name, output_dir)
  summary_stats <- export_significant_genes(result$results, contrast_name, output_dir)
  
  cat("    Significant genes: ", sum(result$results$significance != "NS", na.rm = TRUE), "\n")
  
  return(summary_stats)
}

# ===============================================
# MAIN EXECUTION
# ===============================================

run_differential_expression_analysis <- function() {
  print_separator()
  cat("DIFFERENTIAL EXPRESSION ANALYSIS (DESeq2)\n")
  print_separator()
  
  all_summaries <- data.frame()
  
  for (gene_group in config$gene_groups) {
    cat("\nProcessing gene group:", gene_group, "\n")
    
    group_out_dir <- file.path(DEA_OUT_DIR, gene_group)
    ensure_output_dir(group_out_dir)
    
    # Load count matrix (gene-level, expected_count)
    input_file <- build_input_path(gene_group, "gene_level", "expected_count", "geneName")
    
    result <- validate_and_read_matrix(input_file, min_rows = 10)
    if (!result$success) {
      cat("  Skipped:", result$reason, "\n")
      next
    }
    
    count_matrix <- result$data
    cat("  Loaded matrix:", nrow(count_matrix), "genes x", ncol(count_matrix), "samples\n")
    
    # Map sample IDs to tissue names for grouping
    sample_tissues <- SAMPLE_LABELS[colnames(count_matrix)]
    names(sample_tissues) <- colnames(count_matrix)
    
    # Run tissue group comparisons
    group_names <- names(TISSUE_GROUPS)
    for (i in 1:(length(group_names) - 1)) {
      for (j in (i + 1):length(group_names)) {
        g1 <- group_names[i]
        g2 <- group_names[j]
        
        # Find samples belonging to each group
        g1_samples <- names(sample_tissues)[sample_tissues %in% TISSUE_GROUPS[[g1]]]
        g2_samples <- names(sample_tissues)[sample_tissues %in% TISSUE_GROUPS[[g2]]]
        
        if (length(g1_samples) >= 2 && length(g2_samples) >= 2) {
          summary <- run_pairwise_comparison(count_matrix, g1_samples, g2_samples, g1, g2, group_out_dir)
          if (!is.null(summary)) {
            summary$GeneGroup <- gene_group
            all_summaries <- rbind(all_summaries, summary)
          }
        }
      }
    }
  }
  
  # Save combined summary
  if (nrow(all_summaries) > 0) {
    write.table(all_summaries,
                file = file.path(DEA_OUT_DIR, "DEA_summary_all_comparisons.tsv"),
                sep = "\t", row.names = FALSE, quote = FALSE)
  }
  
  print_separator()
  cat("DEA complete. Output:", DEA_OUT_DIR, "\n")
  print_separator()
}

# Run if executed directly
if (!interactive()) {
  run_differential_expression_analysis()
}
