#!/usr/bin/env Rscript

# ===============================================
# BAR GRAPH GENERATION MODULE
# ===============================================
# Generates bar graphs showing expression per gene/sample

suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  library(RColorBrewer)
})

SCRIPT_DIR <- Sys.getenv("ANALYSIS_MODULES_DIR", ".")
source(file.path(SCRIPT_DIR, "0_shared_config.R"))
source(file.path(SCRIPT_DIR, "1_utility_functions.R"))
source(file.path(SCRIPT_DIR, "2_processing_engine.R"))

# ===============================================
# CONFIGURATION
# ===============================================

BAR_GRAPH_OUT_DIR <- file.path(CONSOLIDATED_BASE_DIR, OUTPUT_SUBDIRS$BAR_GRAPH)

# Bar graph style options
BAR_FILL_PALETTE <- "Set2"
SHOW_ERROR_BARS <- TRUE
# Error bar type: "SD" for standard deviation, "SE" for standard error
ERROR_BAR_TYPE <- "SD"
FACET_WRAP_COLS <- 4

# ===============================================
# BAR GRAPH GENERATION
# ===============================================

generate_bar_graph <- function(data_matrix, output_path, title, 
                                y_label = "Expression", facet_by_gene = TRUE) {
  tryCatch({
    # Convert to long format
    df <- as.data.frame(data_matrix)
    df$Gene <- rownames(df)
    
    df_long <- tidyr::pivot_longer(
      df, 
      cols = -Gene, 
      names_to = "Sample", 
      values_to = "Expression"
    )
    
    # Color palette
    n_samples <- length(unique(df_long$Sample))
    if (n_samples <= 8) {
      colors <- brewer.pal(max(3, n_samples), BAR_FILL_PALETTE)
    } else {
      colors <- colorRampPalette(brewer.pal(8, BAR_FILL_PALETTE))(n_samples)
    }
    
    # Create plot
    if (facet_by_gene && nrow(data_matrix) <= 20) {
      p <- ggplot(df_long, aes(x = Sample, y = Expression, fill = Sample)) +
        geom_bar(stat = "identity", position = "dodge") +
        facet_wrap(~ Gene, scales = "free_y", ncol = FACET_WRAP_COLS) +
        scale_fill_manual(values = colors) +
        labs(title = title, x = "Sample", y = y_label) +
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          strip.text = element_text(face = "bold"),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold")
        )
    } else {
      # Aggregate by sample for large gene sets
      df_summary <- df_long %>%
        group_by(Sample) %>%
        summarise(
          Mean = mean(Expression, na.rm = TRUE),
          SD = sd(Expression, na.rm = TRUE),
          .groups = "drop"
        )
      
      p <- ggplot(df_summary, aes(x = Sample, y = Mean, fill = Sample)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values = colors) +
        labs(title = paste0(title, " (Mean across ", nrow(data_matrix), " genes)"),
             x = "Sample", y = paste0("Mean ", y_label)) +
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold")
        )
      
      if (SHOW_ERROR_BARS) {
        p <- p + geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), 
                               width = 0.2, alpha = 0.7)
      }
    }
    
    # Dynamic sizing
    n_genes <- nrow(data_matrix)
    height <- if (facet_by_gene && n_genes <= 20) max(8, ceiling(n_genes / FACET_WRAP_COLS) * 2.5) else 6
    width <- max(10, n_samples * 0.8)
    
    ggsave(output_path, p, width = width, height = height, dpi = 150)
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

process_bar_graph <- function(gene_group, gene_group_output_dir, processing_level,
                              count_type, gene_type, label_type, norm_scheme,
                              raw_data_matrix, normalized_data, overwrite, extra_options) {
  
  local_total <- 0
  local_successful <- 0
  local_skipped <- 0
  
  version_dir <- file.path(gene_group_output_dir, processing_level, count_type, 
                           gene_type, norm_scheme)
  ensure_output_dir(version_dir)
  
  title_base <- build_title_base(gene_group, count_type, gene_type,
                                 label_type, processing_level, norm_scheme)
  output_path <- file.path(version_dir, paste0(title_base, "_barplot.png"))
  
  local_total <- 1
  
  if (should_skip_existing(output_path, overwrite)) {
    local_skipped <- 1
    return(list(total = local_total, successful = local_successful, skipped = local_skipped))
  }
  
  success <- generate_bar_graph(
    data_matrix = normalized_data,
    output_path = output_path,
    title = gene_group,
    y_label = get_legend_title(norm_scheme, count_type)
  )
  
  if (success) local_successful <- 1
  
  list(total = local_total, successful = local_successful, skipped = local_skipped)
}

# ===============================================
# MAIN
# ===============================================

run_bar_graph <- function(config = NULL, matrices_dir = NULL) {
  # Get method base directory from environment for config file loading
  method_base_dir <- Sys.getenv("METHOD_BASE_DIR", unset = ".")
  if (is.null(config)) config <- load_runtime_config(method_base_dir)
  if (is.null(matrices_dir)) {
    matrices_dir <- file.path(method_base_dir, get_matrices_dir(CURRENT_METHOD))
  }
  ensure_output_dir(BAR_GRAPH_OUT_DIR)
  
  print_config_summary("BAR GRAPH GENERATION", config)
  
  results <- process_all_combinations(
    config = config,
    output_base_dir = BAR_GRAPH_OUT_DIR,
    processing_callback = process_bar_graph,
    matrices_dir = matrices_dir
  )
  
  print_summary(results$successful, results$total, results$skipped)
}

if (!interactive() && identical(environment(), globalenv())) {
  run_bar_graph()
}
