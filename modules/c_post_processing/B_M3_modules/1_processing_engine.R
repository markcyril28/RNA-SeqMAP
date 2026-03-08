# ===============================================
# GENERIC PROCESSING ENGINE - METHOD 3 (STAR + SALMON)
# ===============================================
# Loops through all parameter combinations and applies
# the specified processing function to each count matrix

# ===============================================
# MAIN PROCESSING LOOP
# ===============================================

process_count_matrices <- function(processing_func, output_subdir, output_suffix) {
  
  # Get configuration from shared config
  count_types <- if (exists("COUNT_TYPES")) COUNT_TYPES else c("counts")
  label_types <- if (exists("LABEL_TYPES")) LABEL_TYPES else c("SRR")
  norm_schemes <- if (exists("NORM_SCHEMES")) NORM_SCHEMES else c("zscore_scaled_to_ten")
  transpose_opts <- if (exists("TRANSPOSE")) TRANSPOSE else c(FALSE)
  sort_opts <- if (exists("SORT_OPTIONS")) SORT_OPTIONS else c(TRUE, FALSE)
  
  # Get gene group files (if gene-specific filtering is needed)
  gene_groups_files <- if (exists("GENE_GROUPS_DIR") && dir.exists(GENE_GROUPS_DIR)) {
    list.files(GENE_GROUPS_DIR, pattern = "\\.(csv|txt|tsv)$", full.names = TRUE)
  } else {
    character(0)
  }
  
  # Find input matrices - tximport produces TSV files
  matrix_files <- list.files(MATRIX_DIR, pattern = "\\.(tsv|csv)$", full.names = TRUE)
  
  if (length(matrix_files) == 0) {
    cat("No matrix files found in:", MATRIX_DIR, "\n")
    cat("Run the main M3 pipeline first to generate tximport outputs.\n")
    return(invisible(NULL))
  }
  
  cat("Found", length(matrix_files), "matrix file(s)\n")
  cat("Gene groups:", length(gene_groups_files), "\n")
  
  total_processed <- 0
  total_success <- 0
  
  for (matrix_file in matrix_files) {
    matrix_name <- tools::file_path_sans_ext(basename(matrix_file))
    cat("\n--- Processing:", matrix_name, "---\n")
    
    # Read raw data
    raw_data <- read_count_matrix(matrix_file)
    if (is.null(raw_data)) {
      cat("  Failed to read matrix\n")
      next
    }
    
    # Determine gene groups to process
    if (length(gene_groups_files) > 0 && !isTRUE(SKIP_FILTERING)) {
      groups_to_process <- gene_groups_files
    } else {
      # Process all genes (create a pseudo gene group)
      groups_to_process <- list(list(name = "AllGenes", genes = rownames(raw_data)))
    }
    
    for (gene_group_item in groups_to_process) {
      # Handle both file paths and list objects
      if (is.character(gene_group_item)) {
        gene_group_name <- tools::file_path_sans_ext(basename(gene_group_item))
        # Support CSV format with Gene_ID column
        file_ext <- tolower(tools::file_ext(gene_group_item))
        if (file_ext == "csv") {
          gene_df <- tryCatch(read.csv(gene_group_item, stringsAsFactors = FALSE, header = TRUE), error = function(e) NULL)
          if (!is.null(gene_df) && "Gene_ID" %in% colnames(gene_df)) {
            gene_list <- trimws(gene_df$Gene_ID)
          } else if (!is.null(gene_df)) {
            gene_list <- trimws(gene_df[[1]])
          } else {
            gene_list <- character(0)
          }
        } else {
          gene_list <- readLines(gene_group_item, warn = FALSE)
          gene_list <- gene_list[!grepl("^#|^Gene_ID", gene_list, ignore.case = TRUE)]
          gene_list <- trimws(gene_list)
        }
        gene_list <- gene_list[gene_list != ""]
      } else {
        gene_group_name <- gene_group_item$name
        gene_list <- gene_group_item$genes
      }
      
      # Filter matrix to gene group
      matching_genes <- intersect(gene_list, rownames(raw_data))
      if (length(matching_genes) == 0) {
        cat("  No matching genes for:", gene_group_name, "\n")
        next
      }
      
      filtered_data <- raw_data[matching_genes, , drop = FALSE]
      cat("  Gene group:", gene_group_name, "(", length(matching_genes), "genes )\n")
      
      # Loop through all parameter combinations
      for (count_type in count_types) {
        for (label_type in label_types) {
          for (norm_scheme in norm_schemes) {
            for (transpose in transpose_opts) {
              for (sort_by_expr in sort_opts) {
                
                total_processed <- total_processed + 1
                
                # Apply normalization
                normalized_data <- apply_normalization(filtered_data, norm_scheme, count_type)
                if (is.null(normalized_data)) next
                
                # Build output path
                sort_label <- if (sort_by_expr) "sorted_expr" else "sorted_organ"
                trans_label <- if (transpose) "transposed" else "normal"
                
                output_name <- paste(gene_group_name, count_type, label_type, 
                                    norm_scheme, trans_label, sort_label, sep = "_")
                
                output_subdir_full <- file.path(OUTPUT_DIR, output_subdir, gene_group_name)
                if (!dir.exists(output_subdir_full)) dir.create(output_subdir_full, recursive = TRUE)
                
                output_path <- file.path(output_subdir_full, paste0(output_name, output_suffix))
                
                # Call processing function
                result <- tryCatch({
                  processing_func(
                    data_matrix = normalized_data,
                    output_path = output_path,
                    title = paste(gene_group_name, "-", count_type),
                    count_type = count_type,
                    label_type = label_type,
                    normalization_scheme = norm_scheme,
                    transpose = transpose,
                    sort_by_expression = sort_by_expr
                  )
                }, error = function(e) {
                  cat("    Error:", e$message, "\n")
                  FALSE
                })
                
                if (isTRUE(result)) total_success <- total_success + 1
              }
            }
          }
        }
      }
    }
  }
  
  cat("\n=== Processing Summary ===\n")
  cat("Total combinations:", total_processed, "\n")
  cat("Successful:", total_success, "\n")
  cat("Failed:", total_processed - total_success, "\n")
}
