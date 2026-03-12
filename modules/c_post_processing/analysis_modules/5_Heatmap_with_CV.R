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
# margin: 1 = row-wise (per gene), 2 = column-wise (per sample)
calculate_cv <- function(data_matrix, is_log_scale = FALSE, margin = 1) {
  # Vectorized CV: avoid apply() loop over rows/columns
  if (margin == 1) {
    # Row-wise (per gene)
    if (is_log_scale) {
      rm <- rowMeans(data_matrix, na.rm = TRUE)
      n_c <- ncol(data_matrix)
      row_sds <- sqrt(rowSums((data_matrix - rm)^2, na.rm = TRUE) / (n_c - 1))
      row_sds[!is.finite(row_sds)] <- NA
      return(row_sds)
    }
    rm <- rowMeans(data_matrix, na.rm = TRUE)
    n_c <- ncol(data_matrix)
    row_sds <- sqrt(rowSums((data_matrix - rm)^2, na.rm = TRUE) / (n_c - 1))
    cv <- ifelse(rm > 0 & is.finite(rm) & is.finite(row_sds), row_sds / rm * 100, NA)
    return(cv)
  } else {
    # Column-wise (per sample)
    if (is_log_scale) {
      cm <- colMeans(data_matrix, na.rm = TRUE)
      n_r <- nrow(data_matrix)
      col_sds <- sqrt(colSums((data_matrix - rep(cm, each = n_r))^2, na.rm = TRUE) / (n_r - 1))
      col_sds[!is.finite(col_sds)] <- NA
      return(col_sds)
    }
    cm <- colMeans(data_matrix, na.rm = TRUE)
    n_r <- nrow(data_matrix)
    col_sds <- sqrt(colSums((data_matrix - rep(cm, each = n_r))^2, na.rm = TRUE) / (n_r - 1))
    cv <- ifelse(cm > 0 & is.finite(cm) & is.finite(col_sds), col_sds / cm * 100, NA)
    return(cv)
  }
}

# ===============================================
# CV HEATMAP GENERATION
# ===============================================

generate_heatmap_with_cv <- function(data_matrix, output_path, title,
                                      count_type, label_type, normalization_type,
                                      norm_scheme = NULL,
                                      transpose = FALSE, sort_by_expression = FALSE,
                                      raw_data_matrix = NULL) {
  tryCatch({
    # CV is most meaningful on raw linear-scale data (TPM/FPKM/coverage).
    # When raw_data_matrix is provided (from the processing callback), use it so that
    # z-score and other centered normalizations don't produce all-NA CV values.
    # Fallback: use data_matrix directly (covers legacy / direct calls).
    cv_source <- if (!is.null(raw_data_matrix)) raw_data_matrix else data_matrix
    # Per-gene CV (row-wise on original orientation) and per-sample CV (column-wise)
    gene_cv <- calculate_cv(cv_source, is_log_scale = FALSE, margin = 1)
    sample_cv <- calculate_cv(cv_source, is_log_scale = FALSE, margin = 2)

    if (transpose) {
      data_matrix <- t(data_matrix)
      # After transpose: rows = samples, columns = genes
      # row_cv = per-sample CV, col_cv = per-gene CV
      row_cv <- sample_cv
      col_cv <- gene_cv
      row_cv_label <- "Organ CV"
      col_cv_label <- "Gene CV"
    } else {
      # Default: rows = genes, columns = samples
      row_cv <- gene_cv
      col_cv <- sample_cv
      row_cv_label <- "Gene CV"
      col_cv_label <- "Organ CV"
    }

    # Sort by mean expression when sort_by_expression is TRUE
    # Goal: Put high expression genes closer to the organ/sample labels
    if (sort_by_expression) {
      if (transpose) {
        # Organs_as_Rows: genes are columns, sort columns so high expression is LEFT (near row labels)
        col_means <- colMeans(data_matrix, na.rm = TRUE)
        sort_order <- order(col_means, decreasing = TRUE)
        data_matrix <- data_matrix[, sort_order, drop = FALSE]
        col_cv <- col_cv[sort_order]
      } else {
        # Genes_as_Rows: genes are rows, sort rows so high expression is TOP (near column labels)
        row_means <- rowMeans(data_matrix, na.rm = TRUE)
        sort_order <- order(row_means, decreasing = TRUE)
        data_matrix <- data_matrix[sort_order, , drop = FALSE]
        row_cv <- row_cv[sort_order]
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
    
    # CV color scales: violet gradient matching heatmap (low CV = deep, high CV = pale)
    # Row CV color scale
    row_cv_range <- range(row_cv, na.rm = TRUE)
    if (row_cv_range[1] == row_cv_range[2] || any(is.na(row_cv_range))) {
      row_cv_range <- c(0, 100)
    }
    row_cv_color_fun <- colorRamp2(
      seq(row_cv_range[1], row_cv_range[2], length.out = 100),
      get_cv_color_scale(100)
    )
    # Column CV color scale
    col_cv_range <- range(col_cv, na.rm = TRUE)
    if (col_cv_range[1] == col_cv_range[2] || any(is.na(col_cv_range))) {
      col_cv_range <- c(0, 100)
    }
    col_cv_color_fun <- colorRamp2(
      seq(col_cv_range[1], col_cv_range[2], length.out = 100),
      get_cv_color_scale(100)
    )
    
    # Legend layout
    legend_layout <- get_legend_layout(LEGEND_POSITION)
    # Use internal norm_scheme for get_legend_title() (which switches on internal names);
    # normalization_type is the display name used for file naming and log-scale detection.
    legend_title <- get_legend_title(if (!is.null(norm_scheme)) norm_scheme else normalization_type, count_type)
    
    # Calculate dimensions for square cells with auto-sizing
    n_rows <- nrow(data_matrix)
    n_cols <- ncol(data_matrix)
    cell_size <- unit(12, "mm")  # Square cell size
    
    # Auto-calculate image dimensions based on heatmap size
    # Add extra margins for CV annotation columns/rows (both row CV and column CV)
    margin_width <- 650   # Space for row names, row CV annotation, and legend
    margin_height <- 650  # Space for column names, column CV annotation, title, legend
    img_width <- max(900, n_cols * 60 + margin_width)
    img_height <- max(800, n_rows * 60 + margin_height)
    
    # Build column CV annotation BEFORE the Heatmap call so it can be passed
    # as top_annotation / bottom_annotation (avoids %v% mixing with row annotations)
    top_cv_anno <- NULL
    bottom_cv_anno <- NULL
    if (length(col_cv) == ncol(data_matrix)) {
      col_cv_text <- sprintf("%.1f", ifelse(is.na(col_cv), 0, col_cv))

      cv_col_anno <- HeatmapAnnotation(
        `CV` = unname(col_cv),
        `CV%` = anno_text(col_cv_text, gp = gpar(fontsize = 9), rot = 45,
                          location = 0.5, just = "center"),
        col = list(`CV` = col_cv_color_fun),
        show_legend = FALSE,
        annotation_label = c(`CV` = col_cv_label, `CV%` = ""),
        annotation_name_side = "left",
        annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),
        gap = unit(2, "mm")
      )

      # Place column CV on the opposite side from the column labels
      if (transpose) {
        bottom_cv_anno <- cv_col_anno
      } else {
        top_cv_anno <- cv_col_anno
      }
    }

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
      top_annotation = top_cv_anno,
      bottom_annotation = bottom_cv_anno,
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

    # Row CV annotation: color strip + numeric text (right side of heatmap)
    if (length(row_cv) == nrow(data_matrix)) {
      row_cv_text <- sprintf("%.1f", ifelse(is.na(row_cv), 0, row_cv))

      cv_row_anno <- rowAnnotation(
        `CV` = unname(row_cv),
        `CV%` = anno_text(row_cv_text, gp = gpar(fontsize = 9), location = 0.5, just = "center"),
        col = list(`CV` = row_cv_color_fun),
        show_legend = FALSE,
        annotation_label = c(`CV` = row_cv_label, `CV%` = ""),
        annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),
        gap = unit(2, "mm")
      )
      ht <- ht + cv_row_anno
    }
    
    # Save with auto-adjusted dimensions
    png(output_path, width = img_width, height = img_height, res = 150)
    draw(ht, heatmap_legend_side = LEGEND_POSITION)
    dev.off()
    
    # Export raw values with CV as TSV alongside the PNG
    if (exists("EXPORT_RAW_VALUES") && EXPORT_RAW_VALUES) {
      tryCatch({
        tsv_path <- sub("\\.png$", "_values.tsv", output_path)
        if (transpose) {
          # Organs as rows, genes as columns — row CV = sample CV
          row_id_label <- if (label_type == "Organ") "OrganID" else "SampleID"
          export_df <- data.frame(
            row_id_label = rownames(data_matrix),
            Sample_CV = row_cv,
            data_matrix,
            check.names = FALSE
          )
          names(export_df)[1] <- row_id_label
          # Append a summary row with per-gene (column) CV
          if (length(col_cv) == ncol(data_matrix)) {
            col_cv_row <- c("Gene_CV", NA, sprintf("%.1f", col_cv))
            export_df <- rbind(export_df, setNames(as.list(col_cv_row), names(export_df)))
          }
        } else {
          # Genes as rows, samples as columns — row CV = gene CV
          export_df <- data.frame(
            GeneID = rownames(data_matrix),
            Gene_CV = row_cv,
            data_matrix,
            check.names = FALSE
          )
          # Append a summary row with per-sample (column) CV
          if (length(col_cv) == ncol(data_matrix)) {
            col_cv_row <- c("Sample_CV", NA, sprintf("%.1f", col_cv))
            export_df <- rbind(export_df, setNames(as.list(col_cv_row), names(export_df)))
          }
        }
        write.table(export_df, tsv_path, sep = "\t", row.names = FALSE, quote = FALSE)
        cat("      Exported values:", basename(tsv_path), "\n")
      }, error = function(e) {
        cat("      Warning: TSV export failed:", e$message, "\n")
      })
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
  
  for (orient in get_orientation_options(gene_group)) {
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
        norm_scheme = norm_scheme,
        transpose = orient$transpose,
        sort_by_expression = sorting$sort,
        raw_data_matrix = raw_data_matrix
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
  # Get method base directory from environment for config file loading
  method_base_dir <- Sys.getenv("METHOD_BASE_DIR", unset = ".")
  if (is.null(config)) config <- load_runtime_config(method_base_dir)
  if (is.null(matrices_dir)) {
    matrices_dir <- file.path(method_base_dir, get_matrices_dir(CURRENT_METHOD))
  }
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
