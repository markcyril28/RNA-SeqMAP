#!/usr/bin/env Rscript

# ===============================================
# PROCESSING ENGINE FOR ALL ANALYSIS MODULES
# ===============================================
# Generic processing loop to avoid repetition across scripts

# ===============================================
# GENERIC PROCESSING FUNCTION
# ===============================================

process_all_combinations <- function(
  config,
  output_base_dir,
  processing_callback,
  min_rows = 2,
  extra_options = NULL,
  matrices_dir = NULL  # Auto-detect based on method if NULL
) {
  counters <- list(total = 0, successful = 0, skipped = 0)
  
  # Auto-detect matrices directory based on current method
  if (is.null(matrices_dir)) {
    matrices_dir <- get_matrices_dir(CURRENT_METHOD)
  }
  
  # Get method-appropriate count types (filter COUNT_TYPES to valid ones for method)
  method_type <- get_method_type(CURRENT_METHOD)
  valid_count_types <- get_count_types(CURRENT_METHOD)
  
  # Map user-configured count types to method-equivalent types
  # When a user requests a count type that doesn't exist for a method,
  # map it to the closest available alternative.
  count_type_mapping <- list(
    "stringtie" = list(
      "expected_count" = "coverage",  # fallback only — coverage is per-base depth, NOT raw fragment counts
      "tpm" = "tpm",
      "fpkm" = "fpkm",
      "coverage" = "coverage"
    ),
    "salmon" = list(
      "expected_count" = "NumReads",
      "tpm" = "tpm",
      "NumReads" = "NumReads"
    ),
    "rsem" = list(
      "expected_count" = "expected_count",
      "tpm" = "tpm",
      "fpkm" = "fpkm"
    ),
    "star" = list(
      "expected_count" = "NumReads",
      "tpm" = "tpm",
      "NumReads" = "NumReads"
    )
  )
  
  # Get mapping for current method (or use identity mapping)
  method_mapping <- count_type_mapping[[method_type]]
  if (is.null(method_mapping)) {
    method_mapping <- setNames(as.list(valid_count_types), valid_count_types)
  }
  
  # Map configured count types to method-appropriate ones
  active_count_types <- unique(unlist(lapply(COUNT_TYPES, function(ct) {
    if (!is.null(method_mapping[[ct]])) {
      method_mapping[[ct]]
    } else if (ct %in% valid_count_types) {
      ct
    } else {
      NULL
    }
  })))
  
  if (length(active_count_types) == 0) {
    cat("  Warning: No valid count types for method", CURRENT_METHOD, 
        "- configured:", paste(COUNT_TYPES, collapse = ","), "\n")
  }
  
  # For StringTie, processing_level is not used in path, so use a placeholder
  active_processing_levels <- if (method_type == "stringtie") c("gene_level") else PROCESSING_LEVELS
  
  for (gene_group in config$gene_groups) {
    # Get combined output folder name (GeneGroup_in_Dataset)
    output_folder_name <- get_output_folder_name(gene_group, CURRENT_DATASET)
    cat("Processing gene group:", gene_group, "\n")
    if (nzchar(CURRENT_DATASET)) {
      cat("  Dataset:", CURRENT_DATASET, "-> Output folder:", output_folder_name, "\n")
    }
    
    gene_group_output_dir <- file.path(output_base_dir, output_folder_name)
    ensure_output_dir(gene_group_output_dir, clean = FALSE)
    
    for (processing_level in active_processing_levels) {
      for (count_type in active_count_types) {
        for (gene_type in GENE_TYPES) {
          for (label_type in LABEL_TYPES) {
            
            input_file <- build_input_path(gene_group, processing_level, count_type, 
                                           gene_type, matrices_dir, config$master_reference,
                                           CURRENT_METHOD, label_type)
            
            validation <- validate_and_read_matrix(input_file, min_rows)
            
            if (!validation$success) {
              cat("  Skipping (", validation$reason, "):", processing_level, "|", 
                  basename(input_file), "\n", sep = "")
              next
            }
            
            cat("  Processing:", processing_level, "|", count_type, "|", gene_type, "|", 
                label_type, "- Found", validation$n_genes, "genes\n")
            
            # Sequential processing for normalization schemes
            # Filter configured NORM_SCHEMES to only those valid for this count type
            # (e.g., "raw"/"cpm"/"deseq2_normalized" are invalid for pre-normalized TPM/FPKM)
            active_norm_schemes <- NORM_SCHEMES[sapply(NORM_SCHEMES, function(ns) is_valid_norm_for_count(count_type, ns))]
            # Pre-compute log2 once for all norm schemes that need it (avoids redundant computation)
            .log2_cache <- NULL
            needs_log2 <- any(active_norm_schemes %in% c("count_type_normalized", "zscore", "zscore_row", "zscore_scaled_to_ten"))
            if (needs_log2) {
              .log2_cache <- preprocess_for_count_type_normalized(validation$data, count_type)
            }
            for (norm_scheme in active_norm_schemes) {
              data_normalized <- apply_normalization(validation$data, norm_scheme, count_type, .log2_cache)
              
              if (is.null(data_normalized)) {
                cat("    Failed normalization:", norm_scheme, "\n")
                next
              }
              
              # Apply gene/sample label transformations based on gene_type and label_type
              raw_labeled <- apply_labels(validation$data, gene_group, gene_type, label_type)
              normalized_labeled <- apply_labels(data_normalized, gene_group, gene_type, label_type)
              
              callback_result <- processing_callback(
                gene_group = gene_group,
                gene_group_output_dir = gene_group_output_dir,
                processing_level = processing_level,
                count_type = count_type,
                gene_type = gene_type,
                label_type = label_type,
                norm_scheme = norm_scheme,
                raw_data_matrix = raw_labeled,
                normalized_data = normalized_labeled,
                overwrite = config$overwrite_existing,
                extra_options = extra_options
              )
              
              counters$total <- counters$total + callback_result$total
              counters$successful <- counters$successful + callback_result$successful
              counters$skipped <- counters$skipped + callback_result$skipped
            }
          }
        }
      }
    }
  }
  
  return(counters)
}

# ===============================================
# ORIENTATION AND SORTING OPTIONS
# ===============================================

get_sorting_options <- function() {
  list(
    list(sort = FALSE, sort_name = "Original_Order"),
    list(sort = TRUE, sort_name = "Sorted_by_Expression")
  )
}

get_orientation_options <- function(gene_group = NULL) {
  # SmelDMP gene groups: organs-as-rows layout only
  if (!is.null(gene_group) && grepl("SmelDMP", gene_group, ignore.case = TRUE)) {
    return(list(
      list(transpose = TRUE, orient_name = "Organs_as_Rows")
    ))
  }
  list(
    list(transpose = FALSE, orient_name = "Genes_as_Rows"),
    list(transpose = TRUE, orient_name = "Organs_as_Rows")
  )
}
