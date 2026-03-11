#!/usr/bin/env Rscript

# ===============================================
# BASIC HEATMAP GENERATION MODULE
# ===============================================
# Generates standard heatmaps from expression matrices

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(grid)
})

# Pre-initialize fontconfig to suppress "using without calling FcInit()" warning
# that fires on the first graphics device creation in a new session.
invisible(suppressWarnings({
  tmp <- tempfile(fileext = ".png")
  png(tmp, width = 1, height = 1); dev.off()
  file.remove(tmp)
}))

SCRIPT_DIR <- Sys.getenv("ANALYSIS_MODULES_DIR", ".")
source(file.path(SCRIPT_DIR, "0_shared_config.R"))
source(file.path(SCRIPT_DIR, "1_utility_functions.R"))
source(file.path(SCRIPT_DIR, "2_processing_engine.R"))

# ===============================================
# CONFIGURATION
# ===============================================

LEGEND_POSITION <- "bottom"
HEATMAP_OUT_DIR <- file.path(CONSOLIDATED_BASE_DIR, OUTPUT_SUBDIRS$BASIC_HEATMAP)

# Export raw values alongside heatmap images
# When TRUE, saves the data matrix as a TSV file with the same name as the PNG
EXPORT_RAW_VALUES <- TRUE

# ===============================================
# HEATMAP GENERATION FUNCTION
# ===============================================

generate_heatmap_violet <- function(data_matrix, output_path, title,
                                     count_type, label_type, normalization_type,
                                     norm_scheme = NULL,
                                     transpose = FALSE, sort_by_expression = FALSE) {
  tryCatch({
    # Prepare data
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
      } else {
        # Genes_as_Rows: genes are rows, sort rows so high expression is TOP (near column labels)
        row_means <- rowMeans(data_matrix, na.rm = TRUE)
        data_matrix <- data_matrix[order(row_means, decreasing = TRUE), , drop = FALSE]
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
    # All schemes use quantile-based color bounds so visual intensity is consistent.
    # zscore_scaled_to_ten is a linear rescaling of zscore to [0,10]; using the same
    # quantile approach ensures identical color patterns between the two.
    # For zscore_scaled_to_ten: show a complete 0-10 legend (increment of 2) so the
    # reader sees the full intuitive scale, even though colors are quantile-mapped.
    is_zscore_scaled <- grepl("zscore.*scaled.*ten|z-score.*scaled.*ten",
                              normalization_type, ignore.case = TRUE)
    if (is_zscore_scaled) {
      legend_breaks <- seq(0, 10, by = 2)
      legend_labels <- as.character(legend_breaks)
    } else {
      legend_breaks <- seq(color_min, color_max, length.out = 5)
      legend_labels <- sprintf("%.1f", legend_breaks)
    }
    
    color_fun <- colorRamp2(
      seq(color_min, color_max, length.out = 100),
      get_violet_color_scale(100)
    )
    
    # Legend
    legend_layout <- get_legend_layout(LEGEND_POSITION)
    # Use internal norm_scheme for get_legend_title() (which switches on internal names);
    # normalization_type is the display name used for file naming and the is_zscore_scaled check.
    legend_title <- get_legend_title(if (!is.null(norm_scheme)) norm_scheme else normalization_type, count_type)
    
    # Calculate dimensions for square cells with auto-sizing
    n_rows <- nrow(data_matrix)
    n_cols <- ncol(data_matrix)
    cell_size <- unit(12, "mm")  # Square cell size
    
    # Auto-calculate image dimensions based on heatmap size
    # Add margins for labels, title, and legend
    margin_width <- 500   # Space for row names and legend
    margin_height <- 550  # Space for column names, title, legend, and top/bottom padding
    img_width <- max(800, n_cols * 60 + margin_width)
    img_height <- max(800, n_rows * 60 + margin_height)
    
    # Clean organ label suffixes (.1, .2) added by R's make.unique on duplicate names
    # Only strip from the axis that carries organ/tissue labels
    if (label_type == "Organ") {
      if (transpose) {
        clean_row_labels <- sub("\\.[0-9]+$", "", rownames(data_matrix))
        clean_col_labels <- colnames(data_matrix)
      } else {
        clean_row_labels <- rownames(data_matrix)
        clean_col_labels <- sub("\\.[0-9]+$", "", colnames(data_matrix))
      }
    } else {
      clean_row_labels <- rownames(data_matrix)
      clean_col_labels <- colnames(data_matrix)
    }

    # Create heatmap without dendrograms, with square cells and visible borders
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
      row_labels = clean_row_labels,
      column_labels = clean_col_labels,
      column_names_side = if (transpose) "top" else "bottom",
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
    
    # Save with auto-adjusted dimensions
    png(output_path, width = img_width, height = img_height, res = 150)
    draw(ht, heatmap_legend_side = LEGEND_POSITION)
    dev.off()
    
    # Export raw values as TSV alongside the PNG
    if (exists("EXPORT_RAW_VALUES") && EXPORT_RAW_VALUES) {
      tsv_path <- sub("\\.png$", "_values.tsv", output_path)
      # Convert to data frame with row names as first column
      row_id_label <- if (transpose) {
        if (label_type == "Organ") "OrganID" else "SampleID"
      } else {
        "GeneID"
      }
      export_df <- data.frame(V1 = rownames(data_matrix), data_matrix, check.names = FALSE)
      names(export_df)[1] <- row_id_label
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

process_basic_heatmap <- function(gene_group, gene_group_output_dir, processing_level,
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
      output_path <- file.path(version_dir, paste0(title_base, ".png"))
      
      local_total <- local_total + 1
      
      if (should_skip_existing(output_path, overwrite)) {
        cat("      Skipping (exists):", basename(output_path), "\n")
        local_skipped <- local_skipped + 1
        next
      }
      
      success <- generate_heatmap_violet(
        data_matrix = normalized_data,
        output_path = output_path,
        title = title_base,
        count_type = count_type,
        label_type = label_type,
        normalization_type = norm_display,
        norm_scheme = norm_scheme,
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

run_basic_heatmap <- function(config = NULL, matrices_dir = NULL) {
  if (is.null(config)) config <- load_runtime_config()
  ensure_output_dir(HEATMAP_OUT_DIR)
  
  print_config_summary("BASIC HEATMAP GENERATION", config)
  
  results <- process_all_combinations(
    config = config,
    output_base_dir = HEATMAP_OUT_DIR,
    processing_callback = process_basic_heatmap,
    matrices_dir = matrices_dir
  )
  
  print_summary(results$successful, results$total, results$skipped)
  cat("Output directory:", HEATMAP_OUT_DIR, "\n")
}

# Run if executed directly
if (!interactive() && identical(environment(), globalenv())) {
  run_basic_heatmap()
}
