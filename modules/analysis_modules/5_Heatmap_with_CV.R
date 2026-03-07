#!/usr/bin/env Rscript

# ===============================================
# HEATMAP WITH CV (COEFFICIENT OF VARIATION) MODULE
# ===============================================
# Generates heatmaps with CV annotations per gene

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(grid)
})

SCRIPT_DIR <- Sys.getenv("ANALYSIS_MODULES_DIR", ".")
source(file.path(SCRIPT_DIR, "0_shared_config.R"))
source(file.path(SCRIPT_DIR, "1_utility_functions.R"))
source(file.path(SCRIPT_DIR, "2_processing_engine.R"))

# ===============================================
# CONFIGURATION
# ===============================================

LEGEND_POSITION <- "bottom"
CV_HEATMAP_OUT_DIR <- file.path(CONSOLIDATED_BASE_DIR, OUTPUT_SUBDIRS$CV_HEATMAP)

# Export raw values alongside heatmap images
# When TRUE, saves the data matrix (with CV column) as a TSV file with the same name as the PNG
EXPORT_RAW_VALUES <- TRUE

# ===============================================
# CV CALCULATION
# ===============================================

# Calculate Coefficient of Variation: CV = (SD / mean) * 100
# NOTE: CV is most meaningful on raw or linear-scale data, not log-transformed.
#       For log data, consider using SD directly as a variability measure.
calculate_cv <- function(data_matrix, is_log_scale = FALSE) {
  apply(data_matrix, 1, function(row) {
    # For log-transformed data, return SD as variability measure
    if (is_log_scale) {
      row_sd <- sd(row, na.rm = TRUE)
      return(if (is.finite(row_sd)) row_sd else NA)
    }
    # For linear data, calculate true CV
    row_mean <- mean(row, na.rm = TRUE)
    row_sd <- sd(row, na.rm = TRUE)
    if (row_mean > 0 && is.finite(row_mean) && is.finite(row_sd)) {
      return(row_sd / row_mean * 100)
    }
    return(NA)
  })
}

# ===============================================
# CV HEATMAP GENERATION
# ===============================================

generate_heatmap_with_cv <- function(data_matrix, output_path, title,
                                      count_type, label_type, normalization_type,
                                      transpose = FALSE, sort_by_expression = FALSE) {
  tryCatch({
    # Determine if data is log-scale based on normalization type
    is_log_scale <- normalization_type %in% c("count_type_normalized", "cpm", 
                                               "deseq2_normalized", "Count-Type_Normalized")
    
    # Calculate CV before any transformation (use raw scale CV if possible)
    cv_values <- calculate_cv(data_matrix, is_log_scale = is_log_scale)
    
    if (transpose) {
      data_matrix <- t(data_matrix)
    }
    
    # Sort by mean expression when sort_by_expression is TRUE
    # Goal: Put high expression genes closer to the organ/sample labels
    if (sort_by_expression) {
      if (transpose) {
        # Organs_as_Rows: genes are columns, sort columns so high expression is LEFT (near row labels)
        col_means <- colMeans(data_matrix, na.rm = TRUE)
        data_matrix <- data_matrix[, order(col_means, decreasing = TRUE), drop = FALSE]
        # Also reorder CV values to match (CV is per-gene, now in columns)
        cv_values <- cv_values[order(col_means, decreasing = TRUE)]
      } else {
        # Genes_as_Rows: genes are rows, sort rows so high expression is TOP (near column labels)
        row_means <- rowMeans(data_matrix, na.rm = TRUE)
        data_matrix <- data_matrix[order(row_means, decreasing = TRUE), , drop = FALSE]
        # Also reorder CV values to match
        cv_values <- cv_values[order(row_means, decreasing = TRUE)]
      }
    }
    
    # Color scale with quantile-based range for better visibility
    # Use 2nd to 98th percentile to avoid extreme values dominating the scale
    # This makes low-expression genes more visible while keeping patterns accurate
    data_values <- as.vector(data_matrix)
    data_values <- data_values[is.finite(data_values)]
    
    if (length(data_values) < 2) {
      cat("      Warning: Insufficient data for heatmap\n")
      return(FALSE)
    }
    
    # Use quantiles for color scale bounds (more robust than min/max)
    color_min <- quantile(data_values, 0.02, na.rm = TRUE)
    color_max <- quantile(data_values, 0.98, na.rm = TRUE)
    
    # Handle case where quantiles are identical (no variation)
    if (color_min == color_max) {
      # Fall back to full range
      color_min <- min(data_values)
      color_max <- max(data_values)
      if (color_min == color_max) {
        cat("      Warning: No data variation, skipping heatmap\n")
        return(FALSE)
      }
    }
    
    # Configure legend breaks based on normalization type
    # Use case-insensitive matching and check for key patterns
    is_zscore_scaled <- grepl("zscore.*scaled.*ten|z-score.*scaled.*ten", normalization_type, ignore.case = TRUE)
    
    if (is_zscore_scaled) {
      # For zscore_scaled_to_ten: data is already scaled to 0-10
      # Use fixed 0-10 color scale and legend with 5 parts (0, 2.5, 5, 7.5, 10)
      legend_breaks <- c(0, 2.5, 5, 7.5, 10)
      legend_labels <- as.character(legend_breaks)
      color_min <- 0
      color_max <- 10
    } else {
      # Divide range into 5 equal parts for other normalizations
      legend_breaks <- seq(color_min, color_max, length.out = 5)
      legend_labels <- sprintf("%.1f", legend_breaks)
    }
    
    color_fun <- colorRamp2(
      seq(color_min, color_max, length.out = 100),
      get_violet_color_scale(100)
    )
    
    # CV color scale using blue-red diverging palette from utility functions
    cv_range <- range(cv_values, na.rm = TRUE)
    # Handle case where all CV values are identical
    if (cv_range[1] == cv_range[2] || any(is.na(cv_range))) {
      cv_range <- c(0, 100)  # Default CV range
    }
    cv_color_fun <- colorRamp2(
      seq(cv_range[1], cv_range[2], length.out = 100),
      get_blue_red_color_scale(100)
    )
    
    # Legend layout
    legend_layout <- get_legend_layout(LEGEND_POSITION)
    legend_title <- get_legend_title(normalization_type, count_type)
    
    # Calculate dimensions for square cells with auto-sizing
    n_rows <- nrow(data_matrix)
    n_cols <- ncol(data_matrix)
    cell_size <- unit(12, "mm")  # Square cell size
    
    # Auto-calculate image dimensions based on heatmap size
    # Add extra margins for CV annotation column
    margin_width <- 550   # Space for row names, CV annotation, and legend
    margin_height <- 550  # Space for column names, title, legend, and top/bottom padding
    img_width <- max(900, n_cols * 60 + margin_width)
    img_height <- max(800, n_rows * 60 + margin_height)
    
    # Main heatmap without dendrograms, with square cells and visible borders
    ht <- Heatmap(
      data_matrix,
      name = legend_title,
      col = color_fun,
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      show_row_dend = FALSE,
      show_column_dend = FALSE,
      show_row_names = nrow(data_matrix) <= 50,
      show_column_names = TRUE,
      row_names_side = "left",
      row_names_gp = gpar(fontsize = 11),
      column_names_gp = gpar(fontsize = 11),
      column_names_rot = 45,
      column_title = title,
      column_title_gp = gpar(fontsize = 14, fontface = "bold"),
      width = n_cols * cell_size,
      height = n_rows * cell_size,
      rect_gp = gpar(col = "white", lwd = 0.25),
      heatmap_legend_param = list(
        direction = legend_layout$direction,
        legend_height = legend_layout$height,
        legend_width = legend_layout$width,
        title_gp = gpar(fontsize = 11),
        labels_gp = gpar(fontsize = 10),
        at = legend_breaks,
        labels = legend_labels
      )
    )
    
    # CV annotation (only if not transposed - CV is per-gene and genes are rows)
    if (!transpose && length(cv_values) == nrow(data_matrix)) {
      # Calculate CV legend breaks
      cv_legend_breaks <- seq(cv_range[1], cv_range[2], length.out = 5)
      cv_legend_labels <- sprintf("%.0f", cv_legend_breaks)
      
      cv_anno <- rowAnnotation(
        CV = cv_values,
        col = list(CV = cv_color_fun),
        annotation_legend_param = list(
          CV = list(
            title = if (is_log_scale) "SD (%)" else "CV (%)",
            direction = legend_layout$direction,
            legend_width = legend_layout$width,
            at = cv_legend_breaks,
            labels = cv_legend_labels
          )
        )
      )
      ht <- ht + cv_anno
    }
    
    # Save with auto-adjusted dimensions
    png(output_path, width = img_width, height = img_height, res = 150)
    draw(ht, heatmap_legend_side = LEGEND_POSITION)
    dev.off()
    
    # Export raw values with CV as TSV alongside the PNG
    if (exists("EXPORT_RAW_VALUES") && EXPORT_RAW_VALUES) {
      tsv_path <- sub("\\.png$", "_values.tsv", output_path)
      # Convert to data frame with row names and CV column
      if (transpose) {
        # When transposed, genes are columns - add CV as a row
        export_df <- data.frame(GeneID = rownames(data_matrix), data_matrix, check.names = FALSE)
      } else {
        # When not transposed, genes are rows - add CV as a column
        export_df <- data.frame(
          GeneID = rownames(data_matrix),
          CV = cv_values,
          data_matrix,
          check.names = FALSE
        )
      }
      write.table(export_df, tsv_path, sep = "\t", row.names = FALSE, quote = FALSE)
      cat("      Exported values:", basename(tsv_path), "\n")
    }
    
    cat("      Generated:", basename(output_path), "\n")
    return(TRUE)
  }, error = function(e) {
    cat("      Error:", e$message, "\n")
    return(FALSE)
  })
}

# ===============================================
# PROCESSING CALLBACK
# ===============================================

process_cv_heatmap <- function(gene_group, gene_group_output_dir, processing_level,
                               count_type, gene_type, label_type, norm_scheme,
                               raw_data_matrix, normalized_data, overwrite, extra_options) {
  
  local_total <- 0
  local_successful <- 0
  local_skipped <- 0
  
  norm_display <- get_norm_display_name(norm_scheme)
  
  for (orient in get_orientation_options()) {
    for (sorting in get_sorting_options()) {
      
      version_dir <- file.path(
        gene_group_output_dir, processing_level, count_type, gene_type,
        norm_scheme, orient$orient_name, sorting$sort_name
      )
      ensure_output_dir(version_dir)
      
      title_base <- build_title_base(gene_group, count_type, gene_type,
                                     label_type, processing_level, norm_scheme)
      output_path <- file.path(version_dir, paste0(title_base, "_with_CV.png"))
      
      local_total <- local_total + 1
      
      if (should_skip_existing(output_path, overwrite)) {
        local_skipped <- local_skipped + 1
        next
      }
      
      success <- generate_heatmap_with_cv(
        data_matrix = normalized_data,
        output_path = output_path,
        title = paste0(title_base, " (with CV)"),
        count_type = count_type,
        label_type = label_type,
        normalization_type = norm_display,
        transpose = orient$transpose,
        sort_by_expression = sorting$sort
      )
      
      if (success) local_successful <- local_successful + 1
    }
  }
  
  list(total = local_total, successful = local_successful, skipped = local_skipped)
}

# ===============================================
# MAIN
# ===============================================

run_cv_heatmap <- function(config = NULL, matrices_dir = NULL) {
  if (is.null(config)) config <- load_runtime_config()
  ensure_output_dir(CV_HEATMAP_OUT_DIR)
  
  print_config_summary("HEATMAP WITH CV GENERATION", config)
  
  results <- process_all_combinations(
    config = config,
    output_base_dir = CV_HEATMAP_OUT_DIR,
    processing_callback = process_cv_heatmap,
    matrices_dir = matrices_dir
  )
  
  print_summary(results$successful, results$total, results$skipped)
}

if (!interactive() && identical(environment(), globalenv())) {
  run_cv_heatmap()
}
