# ===============================================
# UTILITY FUNCTIONS FOR METHOD 3 - STAR + SALMON QUANTIFICATION
# ===============================================
# Comprehensive utility functions (shared with M4/M5)
# M3 uses: STAR alignment + Salmon quantification + tximport
# Supports: basic heatmaps, CV heatmaps, and bar graphs

# ===============================================
# FILE I/O AND DATA PREPARATION
# ===============================================

read_count_matrix <- function(file_path) {
  tryCatch({
    data <- read.table(file_path, header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE,
                      na.strings = c("", " ", "NA", "null"))

    if (any(duplicated(data[, 1]))) {
      data[, 1] <- make.unique(as.character(data[, 1]), sep = "_")
    }

    rownames(data) <- data[, 1]
    data <- data[, -1, drop = FALSE]

    for (i in seq_len(ncol(data))) {
      if (!is.numeric(data[, i])) {
        data[, i] <- suppressWarnings(as.numeric(as.character(data[, i])))
      }
    }

    data_matrix <- as.matrix(data)
    data_matrix[is.na(data_matrix)] <- 0

    if (!is.numeric(data_matrix)) {
      data_matrix <- suppressWarnings(apply(data_matrix, 2, as.numeric))
      data_matrix[is.na(data_matrix)] <- 0
    }

    return(data_matrix)
  }, error = function(e) {
    cat("Error reading file:", file_path, "\n")
    cat("Error message:", e$message, "\n")
    return(NULL)
  })
}

save_matrix_data <- function(data_matrix, output_path, metadata = NULL) {
  tryCatch({
    base_path <- tools::file_path_sans_ext(output_path)

    if (!is.null(metadata)) {
      meta_parts <- c()
      if (!is.null(metadata$count_type)) {
        meta_parts <- c(meta_parts, toupper(metadata$count_type))
      }
      if (!is.null(metadata$gene_type)) {
        meta_parts <- c(meta_parts, metadata$gene_type)
      }
      if (!is.null(metadata$label_type)) {
        meta_parts <- c(meta_parts, paste0(metadata$label_type, "_labels"))
      }
      if (!is.null(metadata$normalization)) {
        meta_parts <- c(meta_parts, metadata$normalization)
      }
      if (!is.null(metadata$sorting)) {
        sort_text <- ifelse(metadata$sorting, "sorted_by_expression", "sorted_by_organ")
        meta_parts <- c(meta_parts, sort_text)
      }

      if (length(meta_parts) > 0) {
        metadata_suffix <- paste0("_", paste(meta_parts, collapse = "_"))
        matrix_path <- paste0(base_path, metadata_suffix, ".tsv")
      } else {
        matrix_path <- paste0(base_path, ".tsv")
      }
    } else {
      matrix_path <- paste0(base_path, ".tsv")
    }

    matrix_df <- data.frame(
      Gene_ID = rownames(data_matrix),
      data_matrix,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    write.table(matrix_df, file = matrix_path, sep = "\t",
                row.names = FALSE, col.names = TRUE, quote = FALSE)

    return(TRUE)
  }, error = function(e) {
    cat("Error saving matrix:", e$message, "\n")
    return(FALSE)
  })
}

truncate_labels <- function(labels, max_length = 25) {
  if (!is.vector(labels) && !is.null(labels)) {
    labels <- as.character(labels)
  }

  if (is.null(labels) || length(labels) == 0) {
    return(character(0))
  }

  labels <- as.character(labels)

  sapply(labels, function(x) {
    if (is.na(x) || is.null(x)) {
      return("")
    }
    x <- as.character(x)
    if (nchar(x) > max_length) {
      paste0(substr(x, 1, max_length - 3), "...")
    } else {
      x
    }
  }, USE.NAMES = FALSE)
}

# ===============================================
# NORMALIZATION FUNCTIONS
# ===============================================

preprocess_for_raw <- function(data_matrix) {
  if (is.null(data_matrix) || nrow(data_matrix) == 0) return(NULL)
  data_processed <- data_matrix
  data_processed[is.na(data_processed) | is.infinite(data_processed)] <- 0
  return(data_processed)
}

preprocess_for_raw_normalized <- function(data_matrix) {
  if (is.null(data_matrix) || nrow(data_matrix) == 0) return(NULL)
  lib_sizes <- colSums(data_matrix, na.rm = TRUE)
  lib_sizes[lib_sizes == 0] <- 1
  data_normalized <- sweep(data_matrix, 2, lib_sizes/1e6, FUN = "/")
  data_normalized[is.na(data_normalized) | is.infinite(data_normalized)] <- 0
  return(data_normalized)
}

preprocess_for_cpm <- function(data_matrix, count_type = "expected_count") {
  if (is.null(data_matrix) || nrow(data_matrix) == 0) return(NULL)
  lib_sizes <- colSums(data_matrix, na.rm = TRUE)
  lib_sizes[lib_sizes == 0] <- 1
  data_cpm <- sweep(data_matrix, 2, lib_sizes/1e6, FUN = "/")
  data_cpm[is.na(data_cpm) | is.infinite(data_cpm)] <- 0
  data_processed <- log2(data_cpm + 1)
  data_processed[is.na(data_processed) | is.infinite(data_processed)] <- 0
  return(data_processed)
}

preprocess_for_count_type_normalized <- function(data_matrix, count_type) {
  if (is.null(data_matrix) || nrow(data_matrix) == 0) return(NULL)
  data_processed <- log2(data_matrix + 1)
  data_processed[is.na(data_processed) | is.infinite(data_processed)] <- 0
  return(data_processed)
}

preprocess_for_zscore <- function(data_matrix, count_type) {
  if (is.null(data_matrix) || nrow(data_matrix) == 0) return(NULL)
  data_norm <- preprocess_for_count_type_normalized(data_matrix, count_type)
  if (nrow(data_norm) > 1 && ncol(data_norm) > 1) {
    global_mean <- mean(data_norm, na.rm = TRUE)
    global_sd <- sd(as.vector(data_norm), na.rm = TRUE)

    if (global_sd > 0 && is.finite(global_sd)) {
      data_zscore <- (data_norm - global_mean) / global_sd
    } else {
      data_zscore <- data_norm - global_mean
    }

    data_zscore[is.na(data_zscore) | is.infinite(data_zscore)] <- 0
    return(data_zscore)
  } else {
    return(data_norm)
  }
}

preprocess_for_zscore_scaled_to_ten <- function(data_matrix, count_type) {
  if (is.null(data_matrix) || nrow(data_matrix) == 0) return(NULL)
  data_processed <- preprocess_for_count_type_normalized(data_matrix, count_type)
  if (nrow(data_processed) > 1 && ncol(data_processed) > 1) {
    global_mean <- mean(data_processed, na.rm = TRUE)
    global_sd <- sd(as.vector(data_processed), na.rm = TRUE)

    if (global_sd > 0 && is.finite(global_sd)) {
      data_processed <- (data_processed - global_mean) / global_sd
    } else {
      data_processed <- data_processed - global_mean
    }

    data_processed[is.na(data_processed) | is.infinite(data_processed)] <- 0
  }

  min_val <- min(data_processed, na.rm = TRUE)
  max_val <- max(data_processed, na.rm = TRUE)
  val_range <- max_val - min_val

  if (val_range > 0 && is.finite(val_range)) {
    data_scaled <- ((data_processed - min_val) / val_range) * 10
  } else {
    data_scaled <- matrix(0, nrow = nrow(data_processed), ncol = ncol(data_processed))
  }

  data_scaled[is.na(data_scaled) | is.infinite(data_scaled)] <- 0
  return(data_scaled)
}

apply_normalization <- function(data_matrix, normalization_scheme, count_type) {
  switch(normalization_scheme,
    "raw" = preprocess_for_raw(data_matrix),
    "raw_normalized" = preprocess_for_raw_normalized(data_matrix),
    "count_type_normalized" = preprocess_for_count_type_normalized(data_matrix, count_type),
    "zscore" = preprocess_for_zscore(data_matrix, count_type),
    "zscore_scaled_to_ten" = preprocess_for_zscore_scaled_to_ten(data_matrix, count_type),
    "cpm" = preprocess_for_cpm(data_matrix, count_type),
    stop("Unknown normalization scheme: ", normalization_scheme)
  )
}

# ===============================================
# BASIC HEATMAP GENERATION
# ===============================================

generate_heatmap_violet <- function(data_matrix, output_path, title, count_type, label_type, normalization_type = "raw", transpose = FALSE, sort_by_expression = FALSE) {
  if (is.null(data_matrix) || nrow(data_matrix) == 0 || nrow(data_matrix) < 2) return(FALSE)
  if (any(is.na(data_matrix)) || any(is.infinite(data_matrix))) {
    data_matrix[is.na(data_matrix) | is.infinite(data_matrix)] <- 0
  }
  if (length(unique(as.vector(data_matrix))) == 1) return(FALSE)

  col_labels <- colnames(data_matrix)
  if (label_type == "Organ" && exists("SAMPLE_LABELS")) {
    col_labels <- ifelse(colnames(data_matrix) %in% names(SAMPLE_LABELS),
                         SAMPLE_LABELS[colnames(data_matrix)],
                         colnames(data_matrix))
    colnames(data_matrix) <- col_labels
  }

  gene_labels <- rownames(data_matrix)

  if (sort_by_expression) {
    if (transpose) {
      row_means <- rowMeans(data_matrix, na.rm = TRUE)
      row_order <- order(row_means)
      data_matrix <- data_matrix[row_order, , drop = FALSE]
      gene_labels <- gene_labels[row_order]
    } else {
      col_means <- colMeans(data_matrix, na.rm = TRUE)
      col_order <- order(col_means)
      data_matrix <- data_matrix[, col_order, drop = FALSE]
      col_labels <- col_labels[col_order]
      colnames(data_matrix) <- col_labels
    }
  }

  if (transpose) {
    data_matrix <- t(data_matrix)
    row_labels_display <- col_labels
    col_labels_display <- gene_labels
  } else {
    row_labels_display <- gene_labels
    col_labels_display <- col_labels
  }

  rownames(data_matrix) <- row_labels_display
  colnames(data_matrix) <- col_labels_display

  save_matrix_data(data_matrix, output_path)

  # Choose palette: special case for Z-score scaled to ten (0..10) vs others
  if (normalization_type == "Z-score_Scaled_to_Ten") {
    # For Z-score scaled to ten: use fixed breaks in [0,10] with sequential palette
    eggplant_colors <- colorRamp2(
      seq(0, 10, length = 8),
      c("#dab3ddff","#d9afe0ff","#c57fd1ff","#ac44beff",
        "#8E24AA","#6A1B9A","#4A148C","#2F1B69")
    )
  } else {
    # For other schemes, map colors across observed data range
    eggplant_colors <- colorRamp2(
      seq(min(data_matrix), max(data_matrix), length = 8),
      c("#FFFFFF","#F3E5F5","#CE93D8","#AB47BC",
        "#8E24AA","#6A1B9A","#4A148C","#2F1B69")
    )
  }

  # Calculate evenly spaced breaks for legend
  data_range <- range(data_matrix, na.rm = TRUE)
  legend_breaks <- if (normalization_type == "Z-score_Scaled_to_Ten") {
    seq(0, 10, by = 2)
  } else {
    pretty(data_range, n = 5)
  }

  legend_title <- switch(normalization_type,
    "raw" = paste0("Raw ", toupper(count_type)),
    "Count-Type_Normalized" = "Log2 Expression",
    "Z-score_Normalized" = "Z-score",
    "Z-score_Scaled_to_Ten" = "Z-Score Scaled to [0-10]",
    "Expression"
  )

  if (LEGEND_POSITION == "bottom") {
    legend_direction <- "horizontal"
    legend_height <- unit(2, "cm")
    legend_width <- unit(10, "cm")
    grid_height <- unit(0.8, "cm")
    grid_width <- unit(2.5, "cm")
  } else {  # "right" (default)
    legend_direction <- "vertical"
    legend_height <- unit(10, "cm")
    legend_width <- unit(2, "cm")
    grid_height <- unit(2.5, "cm")
    grid_width <- unit(0.8, "cm")
  }

  tryCatch({
    # Truncate row names for display if many rows
    max_length <- if (nrow(data_matrix) > 50) 20 else 30

    if (is.null(rownames(data_matrix)) || length(rownames(data_matrix)) == 0) {
      rownames(data_matrix) <- paste0("Row_", seq_len(nrow(data_matrix)))
    }

    row_labels <- truncate_labels(rownames(data_matrix), max_length = max_length)
    rownames(data_matrix) <- row_labels

    if (transpose) {
      show_row_names <- TRUE
      show_col_names <- ifelse(ncol(data_matrix) > 50, FALSE, TRUE)
      row_fontsize <- 20
      col_fontsize <- if (ncol(data_matrix) > 50) 7 else 20
      row_side <- "left"
      col_side <- "top"
      col_rot <- 45
      row_label <- "Organs"
      col_label <- "Genes"
    } else {
      show_row_names <- ifelse(nrow(data_matrix) > 50, FALSE, TRUE)
      show_col_names <- TRUE
      row_fontsize <- if (nrow(data_matrix) > 50) 7 else 20
      col_fontsize <- 20
      row_side <- "left"
      col_side <- "top"
      col_rot <- 45
      row_label <- "Genes"
      col_label <- "Organs"
    }

    clean_title <- title
    clean_title <- gsub("_geneName.*", "", clean_title)
    clean_title <- gsub("_geneID.*", "", clean_title)
    clean_title <- gsub("_SRR.*", "", clean_title)
    clean_title <- gsub("_Organ.*", "", clean_title)
    clean_title <- gsub("_", " ", clean_title)

    sorting_text <- ifelse(sort_by_expression, "Sorted_by_Expression", "Sorted_by_Organ")
    full_title <- paste(clean_title, "-", normalization_type, "-", sorting_text)
    if (normalization_type == "raw") {
      full_title <- paste(clean_title, "-", sorting_text)
    }

    ht <- Heatmap(
      data_matrix,
      name = legend_title,
      col = eggplant_colors,
      show_row_names = show_row_names,
      row_names_side = row_side,
      row_names_gp = gpar(fontsize = row_fontsize),
      row_names_max_width = unit(10, "cm"),
      row_title = row_label,
      row_title_gp = gpar(fontsize = 14, fontface = "bold"),
      row_title_side = "left",
      cluster_rows = FALSE,
      show_column_names = show_col_names,
      column_names_side = col_side,
      column_names_gp = gpar(fontsize = col_fontsize),
      column_names_rot = col_rot,
      column_title_side = "top",
      cluster_columns = FALSE,
      rect_gp = gpar(col = "white", lwd = if (max(nrow(data_matrix), ncol(data_matrix)) > 50) 0.3 else 0.5),
      heatmap_legend_param = list(
        title = legend_title,
        legend_direction = legend_direction,
        legend_height = legend_height,
        legend_width = legend_width,
        title_gp = gpar(fontsize = 12, fontface = "bold"),
        title_position = "topcenter",
        labels_gp = gpar(fontsize = 12),
        grid_height = grid_height,
        grid_width = grid_width,
        at = legend_breaks,
        just = c("left", "top")
      ),
      column_title = paste0(full_title, "\n", col_label),
      column_title_gp = gpar(fontsize = 14, fontface = "bold")
    )

    n_rows <- nrow(data_matrix)
    n_cols <- ncol(data_matrix)
    base_cell_size <- 1
    base_width <- n_cols * base_cell_size
    base_height <- n_rows * base_cell_size
    padding <- 0.75
    final_width <- max(8, min(20, base_width + 2 * padding))
    final_height <- max(6, min(16, base_height + 2 * padding))

    png(output_path, width = final_width, height = final_height, units = "in", res = 1080)
    draw(ht,
         heatmap_legend_side = LEGEND_POSITION,
         annotation_legend_side = LEGEND_POSITION,
         merge_legends = TRUE,
         padding = unit(c(padding, padding, padding, padding), "inches"))
    dev.off()

    return(TRUE)
  }, error = function(e) {
    cat("Error generating heatmap for:", title, "\n")
    cat("Error details:", e$message, "\n")
    return(FALSE)
  })
}

# ===============================================
# HEATMAP WITH CV GENERATION
# ===============================================

generate_heatmap_with_cv <- function(data_matrix, output_path, title, count_type, label_type, normalization_scheme = "count_type_normalized", sort_by_expression = FALSE, raw_data_matrix = NULL) {
  TEXT_SIZE <- 14

  if (is.null(data_matrix) || nrow(data_matrix) == 0) return(FALSE)
  if (nrow(data_matrix) < 2) return(FALSE)
  if (any(is.na(data_matrix)) || any(is.infinite(data_matrix))) {
    data_matrix[is.na(data_matrix) | is.infinite(data_matrix)] <- 0
  }
  if (length(unique(as.vector(data_matrix))) == 1) return(FALSE)

  col_labels <- colnames(data_matrix)
  if (label_type == "SRR") {
    col_labels <- colnames(data_matrix)
  } else if (label_type == "Organ" && exists("SAMPLE_LABELS")) {
    col_labels <- ifelse(colnames(data_matrix) %in% names(SAMPLE_LABELS),
                         SAMPLE_LABELS[colnames(data_matrix)],
                         colnames(data_matrix))
  }

  colnames(data_matrix) <- col_labels

  if (!is.null(raw_data_matrix)) {
    colnames(raw_data_matrix) <- col_labels
  }

  gene_labels <- rownames(data_matrix)

  if (sort_by_expression) {
    col_means <- colMeans(data_matrix, na.rm = TRUE)
    col_order <- order(col_means)
    data_matrix <- data_matrix[, col_order, drop = FALSE]
    col_labels <- col_labels[col_order]
    colnames(data_matrix) <- col_labels

    if (!is.null(raw_data_matrix)) {
      raw_data_matrix <- raw_data_matrix[, col_order, drop = FALSE]
      colnames(raw_data_matrix) <- col_labels
    }
  }

  if (is.null(raw_data_matrix)) {
    stop("ERROR: Cannot calculate CV without raw data.\n",
         "       raw_data_matrix parameter is required for all normalizations.\n",
         "       CV must be calculated on untransformed counts for valid interpretation.\n",
         "       Current normalization scheme: ", normalization_scheme)
  }

  cv_values <- apply(raw_data_matrix, 1, function(x) {
    mean_val <- mean(x, na.rm = TRUE)
    sd_val <- sd(x, na.rm = TRUE)
    if (mean_val == 0 || is.na(mean_val) || is.na(sd_val)) {
      return(NA_real_)
    }
    return(sd_val / mean_val)
  })
  cv_labels <- gene_labels

  data_scaled <- data_matrix
  data_scaled[is.na(data_scaled)] <- 0

  row_labels_display <- gene_labels
  col_labels_display <- col_labels

  rownames(data_scaled) <- row_labels_display
  colnames(data_scaled) <- col_labels_display

  gene_type <- if (grepl("geneName", title, ignore.case = TRUE)) {
    "geneName"
  } else if (grepl("geneID", title, ignore.case = TRUE)) {
    "geneID"
  } else {
    "genes"
  }

  save_matrix_data(data_scaled, output_path, metadata = list(
    count_type = count_type,
    gene_type = gene_type,
    label_type = label_type,
    normalization = normalization_scheme,
    sorting = sort_by_expression
  ))

  cv_df <- data.frame(
    Label = cv_labels,
    CV = round(cv_values, 3),
    stringsAsFactors = FALSE
  )

  max_length <- if (nrow(data_scaled) > 50) 15 else 20

  if (is.null(rownames(data_scaled)) || length(rownames(data_scaled)) == 0) {
    rownames(data_scaled) <- paste0("Row_", seq_len(nrow(data_scaled)))
  }

  row_labels <- truncate_labels(rownames(data_scaled), max_length = max_length)
  rownames(data_scaled) <- row_labels

  # Choose color mapping depending on normalization scheme
  if (normalization_scheme == "zscore") {
    # For pure z-score, use symmetric diverging palette around zero
    max_abs_zscore <- max(abs(data_scaled), na.rm = TRUE)
    if (is.finite(max_abs_zscore) && max_abs_zscore > 0) {
      zscore_limit <- max(2, ceiling(max_abs_zscore))
      heatmap_colors <- colorRamp2(
        seq(-zscore_limit, zscore_limit, length = 10),
        c("#FFFFFF", "#FFEB3B", "#CDDC39", "#8BC34A", "#4CAF50", 
          "#009688", "#00BCD4", "#3F51B5", "#673AB7", "#4A148C")
      )
    } else {
      heatmap_colors <- colorRamp2(c(0, 1), c("#FFFFFF", "#4A148C"))
    }
  } else if (normalization_scheme == "zscore_scaled_to_ten") {
    # For Z-score scaled to [0,10], use fixed range with sequential palette
    heatmap_colors <- colorRamp2(
      seq(0, 10, length = 10),
      c("#FFFFFF", "#FFEB3B", "#CDDC39", "#8BC34A", "#4CAF50", 
        "#009688", "#00BCD4", "#3F51B5", "#673AB7", "#4A148C")
    )
  } else {
    # For other normalizations, use data range with sequential palette
    data_range <- range(data_scaled, na.rm = TRUE, finite = TRUE)
    if (length(data_range) == 2 && data_range[1] != data_range[2]) {
      heatmap_colors <- colorRamp2(
        seq(data_range[1], data_range[2], length = 10),
        c("#FFFFFF", "#FFEB3B", "#CDDC39", "#8BC34A", "#4CAF50", 
          "#009688", "#00BCD4", "#3F51B5", "#673AB7", "#4A148C")
      )
    } else {
      heatmap_colors <- colorRamp2(c(0, 1), c("#FFFFFF", "#4A148C"))
    }
  }

  # Calculate legend ticks
  legend_breaks_cv <- if (normalization_scheme == "zscore_scaled_to_ten") {
    seq(0, 10, by = 2)
  } else if (normalization_scheme == "zscore") {
    max_abs_zscore <- max(abs(data_scaled), na.rm = TRUE)
    if (is.finite(max_abs_zscore) && max_abs_zscore > 0) {
      zscore_limit <- max(2, ceiling(max_abs_zscore))
      pretty(seq(-zscore_limit, zscore_limit, length.out = 5), n = 5)
    } else {
      NULL
    }
  } else {
    data_range <- range(data_scaled, na.rm = TRUE)
    pretty(data_range, n = 5)
  }

  tryCatch({
    cv_annotation <- rowAnnotation(
      CV = anno_text(
        sprintf("%.3f", cv_values),
        gp = gpar(fontsize = TEXT_SIZE, col = "black"),
        just = "left",
        width = unit(1.8, "cm")
      ),
      annotation_name_gp = gpar(fontsize = TEXT_SIZE, fontface = "bold"),
      annotation_name_side = "top",
      annotation_name_rot = 0,
      gap = unit(0.2, "cm"),
      simple_anno_size = unit(1.8, "cm"),
      annotation_label = "Coefficient of Variation",
      show_annotation_name = TRUE
    )

    legend_title <- switch(normalization_scheme,
      "raw" = "Raw Counts (Z-scored for viz)",
      "raw_normalized" = "CPM (Z-scored for viz)",
      "count_type_normalized" = "Log2 Expression (Z-scored for viz)",
      "zscore" = "Z-score",
      "zscore_scaled_to_ten" = "ZScore Scaled to [0-10]",
      "cpm" = "Log2(CPM+1) (Z-scored for viz)",
      "Expression"
    )

    # Configure legend placement parameters
    if (exists("LEGEND_POSITION") && LEGEND_POSITION == "bottom") {
      legend_direction <- "horizontal"
      legend_height <- unit(1.5, "cm")
      legend_width <- unit(6, "cm")
      grid_height <- unit(0.6, "cm")
      grid_width <- unit(1.2, "cm")
    } else {  # "top" or "right"
      legend_direction <- "vertical"
      legend_height <- unit(8, "cm")
      legend_width <- unit(2, "cm")
      grid_height <- unit(1.8, "cm")
      grid_width <- unit(0.8, "cm")
    }

    show_row_names <- TRUE
    show_col_names <- TRUE
    row_fontsize <- TEXT_SIZE * 0.90
    col_fontsize <- TEXT_SIZE * 0.90
    col_rot <- 45
    row_label <- "Genes"
    col_label <- "Organs"

    clean_title <- title
    clean_title <- gsub("_geneName.*", "", clean_title)
    clean_title <- gsub("_geneID.*", "", clean_title)
    clean_title <- gsub("_SRR.*", "", clean_title)
    clean_title <- gsub("_Organ.*", "", clean_title)
    clean_title <- gsub("_", " ", clean_title)

    count_type_text <- toupper(count_type)
    normalization_text <- normalization_scheme
    sorting_text <- ifelse(sort_by_expression, "sorted by expression", "sorted by organ")

    full_title <- paste0(clean_title, " (", count_type_text, " | ",
                         gene_type, " | ", label_type, " labels | ",
                         normalization_text, " | ", sorting_text, ")")

    # Construct heatmap object with improved settings
    ht <- Heatmap(
      data_scaled,
      name = legend_title,
      col = heatmap_colors,
      
      # Row (gene) settings
      show_row_names = show_row_names,
      row_names_side = "left",
      row_names_gp = gpar(fontsize = row_fontsize),
      row_names_max_width = unit(6, "cm"),
      row_title = row_label,
      row_title_gp = gpar(fontsize = TEXT_SIZE, fontface = "bold"),
      row_title_side = "left",
      cluster_rows = FALSE,
      
      # Column (sample) settings
      show_column_names = show_col_names,
      column_names_side = "bottom",
      column_names_gp = gpar(fontsize = col_fontsize),
      column_names_rot = col_rot,
      column_title_side = "top",
      cluster_columns = FALSE,
      
      # Heatmap body settings
      rect_gp = gpar(col = "white", lwd = 0.3),
      
      # Legend settings (dynamic based on LEGEND_POSITION)
      heatmap_legend_param = list(
        title = legend_title,
        legend_direction = legend_direction,
        legend_height = legend_height,
        legend_width = legend_width,
        title_gp = gpar(fontsize = TEXT_SIZE * 0.86, fontface = "bold"),
        labels_gp = gpar(fontsize = TEXT_SIZE * 0.71),
        grid_height = grid_height,
        grid_width = grid_width,
        at = legend_breaks_cv
      ),
      
      # Title
      column_title = full_title,
      column_title_gp = gpar(fontsize = TEXT_SIZE, fontface = "bold"),
      
      # CV annotation on the right
      right_annotation = cv_annotation
    )

    # Compute dynamic image size and render
    n_rows <- nrow(data_scaled)
    n_cols <- ncol(data_scaled)
    base_cell_size <- 0.5
    min_width <- 8
    min_height <- 6
    max_width <- 24
    max_height <- 20
    label_width <- 4
    label_height <- 3
    legend_space <- if (exists("LEGEND_POSITION") && LEGEND_POSITION == "bottom") 2 else 3
    final_width <- max(min_width, min(max_width, n_cols * base_cell_size + label_width + legend_space))
    final_height <- max(min_height, min(max_height, n_rows * base_cell_size + label_height + legend_space))

    png(output_path, width = final_width, height = final_height, units = "in", res = 300)
    draw(ht, 
         heatmap_legend_side = "top", 
         annotation_legend_side = "top",
         merge_legends = TRUE,
         gap = unit(0.5, "cm"),
         padding = unit(c(1.5, 0.75, 0.75, 0.75), "inches"))

    # Add axis label at the bottom for samples/organs
    grid.text(col_label, x = 0.5, y = 0.17, just = "center",
              gp = gpar(fontsize = TEXT_SIZE, fontface = "bold"))

    dev.off()

    base_path <- tools::file_path_sans_ext(output_path)
    sort_text <- ifelse(sort_by_expression, "sorted_by_expression", "sorted_by_organ")
    cv_output_path <- paste0(base_path, "_CV_data_",
                             toupper(count_type), "_",
                             gene_type, "_",
                             label_type, "_labels_",
                             normalization_scheme, "_",
                             sort_text, ".tsv")
    write.table(cv_df, file = cv_output_path, sep = "\t",
                row.names = FALSE, col.names = TRUE, quote = FALSE)

    return(TRUE)
  }, error = function(e) {
    cat("Error generating heatmap with CV for:", title, "\n")
    cat("Error details:", e$message, "\n")
    return(FALSE)
  })
}

# ===============================================
# BAR GRAPH GENERATION
# ===============================================

generate_gene_bargraphs <- function(data_matrix, output_dir, title_prefix, count_type, label_type, normalization_scheme = "count_type_normalized", sort_by_expression = FALSE, overwrite = TRUE) {
  TEXT_SIZE <- 14

  if (is.null(data_matrix) || nrow(data_matrix) == 0) return(0)
  if (ncol(data_matrix) == 0) return(0)

  col_labels <- colnames(data_matrix)
  if (label_type == "SRR") {
    col_labels <- colnames(data_matrix)
  } else if (label_type == "Organ" && exists("SAMPLE_LABELS")) {
    col_labels <- ifelse(colnames(data_matrix) %in% names(SAMPLE_LABELS),
                         SAMPLE_LABELS[colnames(data_matrix)],
                         colnames(data_matrix))
  }

  successful_plots <- 0

  for (i in seq_len(nrow(data_matrix))) {
    gene_name <- rownames(data_matrix)[i]
    gene_values <- suppressWarnings(as.numeric(data_matrix[i, ]))

    if (all(gene_values == 0)) next

    plot_data <- data.frame(
      Stage = col_labels,
      Expression = gene_values,
      stringsAsFactors = FALSE
    )

    if (sort_by_expression) {
      plot_data <- plot_data[order(plot_data$Expression), ]
      plot_data$Stage <- factor(plot_data$Stage, levels = plot_data$Stage)
    } else {
      plot_data$Stage <- factor(plot_data$Stage, levels = col_labels)
    }

    y_label <- switch(normalization_scheme,
      "raw" = "Raw counts",
      "raw_normalized" = "CPM (not log-transformed)",
      "count_type_normalized" = "Log2(Expression + 1)",
      "zscore" = "Z-score",
      "zscore_scaled_to_ten" = "Z-Score Scaled to [0-10]",
      "cpm" = "Log2(CPM + 1)",
      "Normalized Expression"
    )

    # Color ramp for fill aesthetic - consistent with heatmaps
    custom_colors <- c("#eeecd6ff", "#FFEB3B", "#CDDC39", "#8BC34A", "#4CAF50", 
                      "#009688", "#00BCD4", "#3F51B5", "#673AB7", "#4A148C")

    # Build the ggplot bar graph
    p <- ggplot(plot_data, aes(x = .data$Stage, y = .data$Expression, fill = .data$Expression)) +
      geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
      labs(title = gene_name, x = "", y = y_label) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = TEXT_SIZE, hjust = 0.5, margin = margin(b = 10)),
        axis.text.x = element_text(angle = 45, hjust = 1, size = TEXT_SIZE * 0.90),
        axis.text.y = element_text(size = TEXT_SIZE * 0.90),
        axis.title.y = element_text(size = TEXT_SIZE * 0.90, margin = margin(r = 10)),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        legend.background = element_rect(fill = "white", color = NA),
        plot.margin = margin(10, 10, 10, 10)
      ) +
      scale_fill_gradientn(colors = custom_colors, name = "Expression") +
      scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

    clean_gene_name <- gsub("[^A-Za-z0-9_-]", "_", gene_name)
    output_file <- file.path(output_dir, paste0(clean_gene_name, "_", normalization_scheme, "_bargraph.png"))

    # Check if file exists and skip if not overwriting
    if (!overwrite && file.exists(output_file)) {
      successful_plots <- successful_plots + 1  # Count as successful
      next
    }

    tryCatch({
      ggsave(output_file, plot = p, width = 8, height = 6, units = "in", dpi = 300, bg = "white")
      successful_plots <- successful_plots + 1
    }, error = function(e) {
      cat("Error saving bar graph for gene:", gene_name, "\n")
      cat("Error details:", e$message, "\n")
    })
  }

  return(successful_plots)
}
