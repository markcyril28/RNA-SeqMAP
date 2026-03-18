#!/usr/bin/env Rscript

# ===============================================
# TISSUE SPECIFICITY MODULE
# ===============================================
# Calculates tissue specificity indices (Tau, Z-score)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(pheatmap)
  library(RColorBrewer)
})

SCRIPT_DIR <- Sys.getenv("ANALYSIS_MODULES_DIR", ".")
source(file.path(SCRIPT_DIR, "0_shared_config.R"))
source(file.path(SCRIPT_DIR, "1_utility_functions.R"))

# ===============================================
# CONFIGURATION
# ===============================================

TISSUE_SPEC_OUT_DIR <- file.path(CONSOLIDATED_BASE_DIR, OUTPUT_SUBDIRS$TISSUE_SPEC)

# Tau threshold for tissue-specific genes
TAU_THRESHOLD <- 0.8
MIN_EXPRESSION <- 1  # Minimum expression to consider

# Figure toggles
GENERATE_TISSUE_FIGURES <- list(
  tau_distribution = TRUE,
  specificity_heatmap = TRUE,
  top_specific_genes = TRUE,
  tissue_barplot = TRUE
)

# ===============================================
# TISSUE SPECIFICITY CALCULATIONS
# ===============================================

# Calculate Tau index (0 = ubiquitous, 1 = tissue-specific)
# Reference: Yanai et al., 2005 - Genome-wide midrange transcription profiles
# Tau = sum(1 - x_i/x_max) / (n - 1) where x_i is expression in tissue i
calculate_tau <- function(expression_vector) {
  x <- as.numeric(expression_vector)
  x[x < 0 | is.na(x)] <- 0  # Handle negative and NA values
  x_max <- max(x, na.rm = TRUE)
  
  if (x_max == 0 || is.na(x_max) || !is.finite(x_max)) return(NA)
  
  # Normalize to max expression
  x_norm <- x / x_max
  n <- length(x)
  
  # Need at least 2 tissues for meaningful Tau
  if (n < 2) return(NA)
  
  tau <- sum(1 - x_norm, na.rm = TRUE) / (n - 1)
  
  # Tau should be between 0 and 1
  tau <- max(0, min(1, tau))
  
  return(tau)
}

# Calculate tissue specificity for all genes
calculate_tissue_specificity <- function(data_matrix) {
  # Vectorized Tau calculation: avoid 3 separate apply() loops
  data_clean <- data_matrix
  data_clean[data_clean < 0 | is.na(data_clean)] <- 0
  # apply(max) avoids the O(n*m) data.frame copy that do.call(pmax, as.data.frame()) creates
  row_max <- apply(data_clean, 1, max, na.rm = TRUE)
  n_tissues <- ncol(data_matrix)

  # Tau = sum(1 - x_i/x_max) / (n - 1), vectorized across all genes
  valid_rows <- row_max > 0 & is.finite(row_max)
  tau_values <- rep(NA_real_, nrow(data_matrix))
  if (any(valid_rows)) {
    data_norm <- data_clean[valid_rows, , drop = FALSE] / row_max[valid_rows]
    tau_values[valid_rows] <- rowSums(1 - data_norm, na.rm = TRUE) / (n_tissues - 1)
    tau_values[valid_rows] <- pmax(0, pmin(1, tau_values[valid_rows]))
  }

  # Vectorized max tissue (column index of max value per row)
  max_col_idx <- max.col(data_clean, ties.method = "first")
  max_tissue <- colnames(data_matrix)[max_col_idx]
  max_tissue[!valid_rows] <- NA_character_

  # Vectorized max expression
  max_expr <- row_max
  max_expr[!valid_rows] <- NA_real_
  
  result <- data.frame(
    Gene = rownames(data_matrix),
    Tau = tau_values,
    MaxTissue = max_tissue,
    MaxExpression = max_expr,
    IsSpecific = !is.na(tau_values) & tau_values >= TAU_THRESHOLD,
    stringsAsFactors = FALSE
  )
  
  result <- result[order(-result$Tau), ]
  return(result)
}

# ===============================================
# VISUALIZATION FUNCTIONS
# ===============================================

create_tau_distribution <- function(specificity_df, output_path, title) {
  if (!GENERATE_TISSUE_FIGURES$tau_distribution) return(NULL)
  
  p <- ggplot(specificity_df, aes(x = Tau)) +
    geom_histogram(bins = 50, fill = "#6C3483", color = "white", alpha = 0.8) +
    geom_vline(xintercept = TAU_THRESHOLD, linetype = "dashed", color = "red", linewidth = 1) +
    annotate("text", x = TAU_THRESHOLD + 0.05, y = Inf, label = paste0("Tau >= ", TAU_THRESHOLD),
             hjust = 0, vjust = 2, color = "red") +
    labs(title = title, x = "Tau (Tissue Specificity Index)", 
         y = "Number of Genes",
         subtitle = paste0("N = ", nrow(specificity_df), " genes")) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  ggsave(output_path, p, width = 10, height = 6, dpi = 150)
  return(p)
}

create_tissue_specificity_heatmap <- function(data_matrix, specificity_df, output_path, 
                                               title, n_top = 50) {
  if (!GENERATE_TISSUE_FIGURES$specificity_heatmap) return(NULL)
  
  # Filter to specific genes
  specific_genes <- specificity_df$Gene[specificity_df$IsSpecific]
  specific_genes <- head(specific_genes, n_top)
  
  if (length(specific_genes) == 0) {
    cat("    No tissue-specific genes found\n")
    return(NULL)
  }
  
  hm_data <- data_matrix[specific_genes, , drop = FALSE]
  
  # Row-scale for visualization
  hm_scaled <- t(scale(t(hm_data)))
  hm_scaled[is.na(hm_scaled)] <- 0
  
  png(output_path, width = 1000, height = 800, res = 100)
  on.exit(try(dev.off(), silent = TRUE), add = TRUE)
  # Strip R's make.unique suffixes (.1, .2) from organ labels for display
  clean_col_labels <- sub("\\.[0-9]+$", "", colnames(hm_scaled))
  pheatmap(
    hm_scaled,
    main = title,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = nrow(hm_scaled) <= 30,
    labels_col = clean_col_labels,
    color = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100),
    border_color = NA
  )
  dev.off()
  on.exit(NULL)
  
  return(TRUE)
}

create_tissue_gene_counts <- function(specificity_df, output_path, title) {
  if (!GENERATE_TISSUE_FIGURES$tissue_barplot) return(NULL)
  
  specific_df <- specificity_df[specificity_df$IsSpecific, ]
  if (nrow(specific_df) == 0) return(NULL)
  
  tissue_counts <- specific_df %>%
    group_by(MaxTissue) %>%
    summarise(Count = n(), .groups = "drop") %>%
    arrange(desc(Count))
  
  p <- ggplot(tissue_counts, aes(x = reorder(MaxTissue, Count), y = Count, fill = MaxTissue)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(title = title, x = "Tissue", y = "Number of Specific Genes") +
    theme_minimal() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold")) +
    scale_fill_manual(values = colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(nrow(tissue_counts)))
  
  ggsave(output_path, p, width = 10, height = 6, dpi = 150)
  return(p)
}

# ===============================================
# MAIN FUNCTION
# ===============================================

run_tissue_specificity <- function(config = NULL, matrices_dir = NULL) {
  # Get method base directory from environment for config file loading
  method_base_dir <- Sys.getenv("METHOD_BASE_DIR", unset = ".")
  if (is.null(config)) config <- load_runtime_config(method_base_dir)
  if (is.null(matrices_dir)) {
    matrices_dir <- file.path(method_base_dir, get_matrices_dir(CURRENT_METHOD))
  }
  ensure_output_dir(TISSUE_SPEC_OUT_DIR)
  
  print_config_summary("TISSUE SPECIFICITY ANALYSIS", config)
  
  successful <- 0
  total <- 0
  
  for (gene_group in config$gene_groups) {
    cat("Processing:", gene_group, "\n")
    total <- total + 1
    
    output_folder_name <- get_output_folder_name(gene_group, CURRENT_DATASET)
    output_dir <- file.path(TISSUE_SPEC_OUT_DIR, output_folder_name)
    ensure_output_dir(output_dir)
    
    input_file <- build_input_path(gene_group, PROCESSING_LEVELS[1],
                                   COUNT_TYPES[1], GENE_TYPES[1],
                                   matrices_dir, config$master_reference)
    
    validation <- validate_and_read_matrix(input_file, 5)
    if (!validation$success) {
      cat("  Skipped:", validation$reason, "\n")
      next
    }
    
    # Tau index (Yanai et al., 2005) requires LINEAR-scale expression values.
    # Tau = sum(1 - x_i/x_max) / (n-1) — log transformation compresses dynamic range
    # and artificially reduces Tau, making genes appear less tissue-specific.
    # Use raw TPM/FPKM for Tau; log-transform only for visualization (heatmaps).
    #
    # Average biological replicates within each tissue before Tau calculation.
    # Tau measures specificity across DISTINCT tissues, not individual samples.
    # Without averaging, replicates inflate the denominator (n-1) while contributing
    # near-identical expression values, systematically underestimating Tau.
    # Convert SRR IDs to organ labels so that:
    #   1. Biological replicates are correctly grouped by organ for averaging
    #   2. MaxTissue reports organ names (e.g., "Flower_Buds"), not SRR IDs
    #   3. Heatmaps display organ labels
    organ_data <- convert_to_organ_labels(validation$data)
    tissue_data <- organ_data
    col_names <- colnames(tissue_data)
    # Strip R's make.unique suffixes (.1, .2, etc.) to recover base tissue names
    base_tissues <- sub("\\.[0-9]+$", "", col_names)
    unique_tissues <- unique(base_tissues)
    if (length(unique_tissues) < length(col_names)) {
      # Replicates present — average expression within each tissue
      averaged_matrix <- sapply(unique_tissues, function(tissue) {
        tissue_cols <- which(base_tissues == tissue)
        if (length(tissue_cols) == 1) {
          tissue_data[, tissue_cols]
        } else {
          rowMeans(tissue_data[, tissue_cols, drop = FALSE], na.rm = TRUE)
        }
      })
      # sapply returns a vector (not a matrix) when tissue_data has a single row;
      # force back to matrix so downstream apply()/rownames() calls work correctly.
      if (!is.matrix(averaged_matrix)) {
        averaged_matrix <- matrix(averaged_matrix, nrow = nrow(tissue_data),
                                  dimnames = list(rownames(tissue_data), unique_tissues))
      }
      rownames(averaged_matrix) <- rownames(tissue_data)
      tissue_data <- averaged_matrix
    }
    specificity <- calculate_tissue_specificity(tissue_data)

    # Log-transformed data for heatmap visualization (organ-labeled, per-sample)
    data_matrix <- apply_normalization(organ_data, NORM_SCHEMES[1], COUNT_TYPES[1])
    
    # Save results
    write.table(specificity, file.path(output_dir, paste0(gene_group, "_tissue_specificity.tsv")),
                sep = "\t", row.names = FALSE, quote = FALSE)
    
    n_specific <- sum(specificity$IsSpecific, na.rm = TRUE)
    cat("  Found", n_specific, "tissue-specific genes (Tau >=", TAU_THRESHOLD, ")\n")
    
    # Generate figures
    create_tau_distribution(
      specificity,
      file.path(output_dir, paste0(gene_group, "_tau_distribution.png")),
      paste0(gene_group, " - Tau Distribution")
    )
    
    create_tissue_specificity_heatmap(
      data_matrix, specificity,
      file.path(output_dir, paste0(gene_group, "_specific_genes_heatmap.png")),
      paste0(gene_group, " - Tissue-Specific Genes")
    )
    
    create_tissue_gene_counts(
      specificity,
      file.path(output_dir, paste0(gene_group, "_tissue_gene_counts.png")),
      paste0(gene_group, " - Genes per Tissue")
    )
    
    successful <- successful + 1
    cat("  Complete\n")
  }
  
  print_summary(successful, total)
}

if (!interactive() && identical(environment(), globalenv())) {
  run_tissue_specificity()
}
