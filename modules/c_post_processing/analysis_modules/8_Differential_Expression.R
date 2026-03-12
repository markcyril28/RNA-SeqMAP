#!/usr/bin/env Rscript

# ===============================================
# DIFFERENTIAL EXPRESSION ANALYSIS MODULE
# ===============================================
# DESeq2-based differential expression analysis
#
# INPUT SOURCE (differs from visualization modules):
#   This module uses RAW INTEGER COUNTS, not the TPM/FPKM/Coverage matrices
#   consumed by heatmap, PCA, and other visualization modules.
#   - M1 HISAT2 RefGuided: prepDE.py integer counts (gene_count_matrix.csv)
#     staged by prepde_matrix_linker.sh into .../deseq2_input/
#   - M3/M4/M5 (Salmon/RSEM): Prefer DESeqDataSetFromTximport() using the saved
#     tximport RDS object, which preserves transcript-length offsets (Soneson et al.
#     2015). Falls back to DESeqDataSetFromMatrix() with rounded counts if RDS
#     is not available.
#   DESeq2 performs its own internal normalization (median-of-ratios) on these
#   raw counts — no external normalization is applied before DESeq2.

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
})

SCRIPT_DIR <- Sys.getenv("ANALYSIS_MODULES_DIR", ".")
source(file.path(SCRIPT_DIR, "0_shared_config.R"))
source(file.path(SCRIPT_DIR, "1_utility_functions.R"))

# ===============================================
# CONFIGURATION
# ===============================================

DEA_OUT_DIR <- file.path(CONSOLIDATED_BASE_DIR, OUTPUT_SUBDIRS$DEA)

# DESeq2 parameters
# Note: padj threshold of 0.05 and LFC of 1.0 are standard stringent cutoffs
# apeglm shrinkage (Zhu et al. 2019) provides better LFC estimates than normal
PADJ_THRESHOLD <- 0.05
LFC_THRESHOLD <- 1.0
SHRINKAGE_TYPE <- "apeglm"  # Options: "apeglm" (recommended), "ashr", "normal"
MIN_COUNT_FILTER <- 10      # Min reads in >= 2 samples (pre-filtering)
INDEPENDENT_FILTERING <- TRUE  # DESeq2's automatic low-count filtering

# Tissue groups for comparisons
# IMPORTANT: Values must match 'Organ' column in SRR_csv/*.csv files exactly
TISSUE_GROUPS <- list(
  "Vegetative" = c("Roots", "Stems", "Leaves", "Senescent_leaves"),
  "Reproductive" = c("Buds_0.7cm", "Opened_Buds", "Flowers", "Pistils"),
  "Fruit" = c("Fruits_1cm", "Fruits_Stage_1", "Fruits_6cm", "Fruits_Skin_Stage_2", 
              "Fruits_Flesh_Stage_2", "Fruits_Calyx_Stage_2", "Fruits_Skin_Stage_3",
              "Fruits_Flesh_Stage_3", "Fruits_peduncle"),
  "Seedling" = c("Radicles", "Cotyledons")
)

# Figure toggles
GENERATE_DEA_FIGURES <- list(
  volcano_plot = TRUE,
  ma_plot = TRUE,
  top_genes_heatmap = TRUE,
  pvalue_histogram = FALSE
)

# ===============================================
# DEA CORE FUNCTIONS
# ===============================================

prepare_deseq_dataset <- function(count_matrix, sample_info) {
  # DESeq2 requires raw integer counts - round expected_count from Salmon/RSEM
  # NOTE: For M4 Salmon / M5 RSEM, rounding expected counts and using
  # DESeqDataSetFromMatrix loses tximport's gene-length offset corrections.
  # The gold-standard approach is DESeqDataSetFromTximport(txi, ...).
  # This simplified path still gives valid results but may be slightly less
  # accurate for genes with large transcript-length variation across samples.
  count_matrix <- round(count_matrix)
  storage.mode(count_matrix) <- "integer"
  
  # Validate minimum samples per condition (need >= 2 for variance estimation)
  condition_counts <- table(sample_info$condition)
  if (any(condition_counts < 2)) {
    cat("  Warning: Some conditions have < 2 replicates (unreliable statistics)\n")
  }
  if (ncol(count_matrix) < 4) {
    cat("  Too few samples for reliable DEA (need >= 4 total)\n")
    return(NULL)
  }
  
  # Remove genes with zero counts across all samples (uninformative)
  count_matrix <- count_matrix[rowSums(count_matrix) > 0, , drop = FALSE]
  
  # Apply minimum count filter: require MIN_COUNT_FILTER in at least 2 samples
  # This pre-filtering reduces multiple testing burden and improves power
  keep <- rowSums(count_matrix >= MIN_COUNT_FILTER) >= 2
  count_matrix <- count_matrix[keep, , drop = FALSE]
  
  if (nrow(count_matrix) < MIN_GENES_DEA) {
    cat("  Too few genes after filtering (need >=", MIN_GENES_DEA, ")\n")
    return(NULL)
  }
  
  dds <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = sample_info,
    design = ~ condition
  )

  return(dds)
}

# Create DESeqDataSet from tximport object (preserves transcript-length offsets).
# For M3/M4/M5: uses average transcript length correction from Salmon/RSEM via
# DESeqDataSetFromTximport() — the DESeq2-recommended approach (Love et al. 2014,
# Soneson et al. 2015). This accounts for differential isoform usage across samples
# that changes the effective gene length, improving accuracy over rounded counts.
prepare_deseq_from_tximport <- function(txi, sample_ids, sample_info, gene_ids = NULL) {
  # Subset tximport to requested samples
  txi_sub <- txi
  txi_sub$counts    <- txi$counts[, sample_ids, drop = FALSE]
  txi_sub$abundance <- txi$abundance[, sample_ids, drop = FALSE]
  txi_sub$length    <- txi$length[, sample_ids, drop = FALSE]

  # Optionally filter to gene group
  if (!is.null(gene_ids)) {
    matched <- match_gene_ids(gene_ids, rownames(txi_sub$counts))
    if (length(matched) == 0) {
      cat("  No gene IDs matched in tximport object\n")
      return(NULL)
    }
    txi_sub$counts    <- txi_sub$counts[matched, , drop = FALSE]
    txi_sub$abundance <- txi_sub$abundance[matched, , drop = FALSE]
    txi_sub$length    <- txi_sub$length[matched, , drop = FALSE]
  }

  # Pre-filter low-count genes (same criteria as prepare_deseq_dataset)
  raw_counts <- round(txi_sub$counts)
  keep <- rowSums(raw_counts) > 0 & rowSums(raw_counts >= MIN_COUNT_FILTER) >= 2
  txi_sub$counts    <- txi_sub$counts[keep, , drop = FALSE]
  txi_sub$abundance <- txi_sub$abundance[keep, , drop = FALSE]
  txi_sub$length    <- txi_sub$length[keep, , drop = FALSE]

  if (nrow(txi_sub$counts) < MIN_GENES_DEA) {
    cat("  Too few genes after filtering (need >=", MIN_GENES_DEA, ")\n")
    return(NULL)
  }

  condition_counts <- table(sample_info$condition)
  if (any(condition_counts < 2)) {
    cat("  Warning: Some conditions have < 2 replicates (unreliable statistics)\n")
  }
  if (ncol(txi_sub$counts) < 4) {
    cat("  Too few samples for reliable DEA (need >= 4 total)\n")
    return(NULL)
  }

  dds <- DESeqDataSetFromTximport(
    txi = txi_sub,
    colData = sample_info,
    design = ~ condition
  )

  cat("    Using tximport offsets (transcript-length corrected)\n")
  return(dds)
}

run_deseq2 <- function(dds, contrast_name, output_dir) {
  set.seed(GLOBAL_RANDOM_SEED)  # Reproducibility for apeglm shrinkage estimation
  dds <- DESeq(dds, quiet = TRUE)
  
  res <- tryCatch({
    if (SHRINKAGE_TYPE == "apeglm") {
      # Find the contrast coefficient — skip the intercept (first element)
      coef_names <- resultsNames(dds)
      contrast_coefs <- coef_names[grepl("^condition_", coef_names)]
      coef_name <- if (length(contrast_coefs) > 0) contrast_coefs[1] else if (length(coef_names) > 1) coef_names[2] else NULL
      if (is.null(coef_name)) {
        cat("    Warning: no condition coefficient found for apeglm shrinkage, using unshrunken results\n")
        results(dds, alpha = PADJ_THRESHOLD)
      } else {
        lfcShrink(dds, coef = coef_name, type = "apeglm", quiet = TRUE)
      }
    } else {
      results(dds, alpha = PADJ_THRESHOLD)
    }
  }, error = function(e) {
    cat("    Note: apeglm shrinkage failed, using unshrunken results:", e$message, "\n")
    results(dds, alpha = PADJ_THRESHOLD)
  })
  
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  res_df <- res_df[, c("gene", setdiff(names(res_df), "gene"))]
  
  res_df$significance <- "NS"
  res_df$significance[res_df$padj < PADJ_THRESHOLD & res_df$log2FoldChange > LFC_THRESHOLD] <- "Up"
  res_df$significance[res_df$padj < PADJ_THRESHOLD & res_df$log2FoldChange < -LFC_THRESHOLD] <- "Down"
  res_df$significance <- factor(res_df$significance, levels = c("Down", "NS", "Up"))
  
  res_df <- res_df[order(res_df$padj), ]
  
  return(list(dds = dds, results = res_df))
}

create_volcano_plot <- function(res_df, contrast_name, output_dir) {
  if (!GENERATE_DEA_FIGURES$volcano_plot) return(NULL)
  
  plot_data <- res_df[!is.na(res_df$padj), ]
  plot_data$neglog10p <- -log10(plot_data$padj)
  finite_vals <- plot_data$neglog10p[is.finite(plot_data$neglog10p)]
  inf_replacement <- if (length(finite_vals) > 0) max(finite_vals) + 10 else 300
  plot_data$neglog10p[is.infinite(plot_data$neglog10p)] <- inf_replacement
  
  n_label <- min(10, sum(plot_data$significance != "NS"))
  top_genes <- head(plot_data[plot_data$significance != "NS", ], n_label)
  
  n_up <- sum(plot_data$significance == "Up", na.rm = TRUE)
  n_down <- sum(plot_data$significance == "Down", na.rm = TRUE)
  
  p <- ggplot(plot_data, aes(x = log2FoldChange, y = neglog10p, color = significance)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("Down" = "#2166AC", "NS" = "grey60", "Up" = "#B2182B"),
                       labels = c(paste0("Down (", n_down, ")"), "NS", paste0("Up (", n_up, ")"))) +
    geom_vline(xintercept = c(-LFC_THRESHOLD, LFC_THRESHOLD), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(PADJ_THRESHOLD), linetype = "dashed", color = "grey40") +
    labs(title = paste0("Volcano: ", contrast_name),
         x = "Log2 Fold Change", y = "-Log10 Adj P-value") +
    theme_minimal() +
    theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5, face = "bold"))
  
  if (nrow(top_genes) > 0) {
    p <- p + geom_text_repel(data = top_genes, aes(label = gene),
                              size = 2.5, max.overlaps = 15)
  }
  
  ggsave(file.path(output_dir, paste0(contrast_name, "_volcano.png")),
         p, width = 10, height = 8, dpi = 150)
  
  return(p)
}

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
         x = "Log10 Mean Expression", y = "Log2 Fold Change") +
    theme_minimal() +
    theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5, face = "bold"))
  
  ggsave(file.path(output_dir, paste0(contrast_name, "_MA.png")),
         p, width = 10, height = 8, dpi = 150)
  
  return(p)
}

create_de_heatmap <- function(dds, res_df, contrast_name, output_dir, n_top = 50) {
  if (!GENERATE_DEA_FIGURES$top_genes_heatmap) return(NULL)
  
  vst_data <- tryCatch({
    assay(vst(dds, blind = FALSE))
  }, error = function(e) {
    log2(counts(dds, normalized = TRUE) + 1)
  })
  
  sig_genes <- res_df[res_df$significance != "NS", "gene"]
  if (length(sig_genes) == 0) return(NULL)
  
  top_genes <- head(sig_genes, n_top)
  hm_data <- vst_data[top_genes, , drop = FALSE]
  
  png(file.path(output_dir, paste0(contrast_name, "_top_DE_heatmap.png")),
      width = 1000, height = 800, res = 100)
  on.exit(try(dev.off(), silent = TRUE), add = TRUE)
  pheatmap(hm_data, scale = "row", cluster_rows = TRUE, cluster_cols = TRUE,
           show_rownames = nrow(hm_data) <= 30,
           main = paste0("Top DE Genes: ", contrast_name))
  dev.off()
  on.exit(NULL)
}

# ===============================================
# M1-SPECIFIC INPUT FUNCTIONS
# ===============================================

# For M1 HISAT2 RefGuided: load the whole-genome integer count matrix produced by prepDE.py
# (staged by prepde_matrix_linker.sh), then filter down to the requested gene group.
load_m1_gene_group_counts <- function(gene_group, matrices_dir, master_ref, cached_gene_ids = NULL,
                                      cached_full_matrix = NULL) {
  # Use cached full matrix if available (avoids re-reading per gene group)
  if (!is.null(cached_full_matrix)) {
    full_matrix <- cached_full_matrix
  } else {
    # Full genome matrix staged by prepde_matrix_linker.sh
    deseq2_csv <- file.path(matrices_dir, master_ref, "deseq2_input", "gene_count_matrix.csv")
    if (!file.exists(deseq2_csv)) {
      return(list(success = FALSE, reason = paste("M1 count matrix not found:", deseq2_csv)))
    }

    full_matrix <- tryCatch(
      as.matrix(data.table::fread(deseq2_csv, header = TRUE, data.table = FALSE,
                                  check.names = FALSE)),
      error = function(e) NULL
    )
    if (!is.null(full_matrix) && ncol(full_matrix) > 0) {
      rownames(full_matrix) <- full_matrix[, 1]
      full_matrix <- full_matrix[, -1, drop = FALSE]
      storage.mode(full_matrix) <- "numeric"
    }
    if (is.null(full_matrix) || nrow(full_matrix) == 0) {
      return(list(success = FALSE, reason = "failed to read M1 count matrix"))
    }
  }

  # Filter to gene group — use pre-cached gene IDs to avoid repeated list.files()
  gene_ids <- cached_gene_ids
  if (is.null(gene_ids)) {
    gene_group_csv <- file.path(GENE_GROUPS_DIR, paste0(gene_group, ".csv"))
    if (!file.exists(gene_group_csv)) {
      hits <- list.files(GENE_GROUPS_DIR, pattern = paste0("^", gene_group, "\\.csv$"),
                         recursive = TRUE, full.names = TRUE)
      if (length(hits) > 0) gene_group_csv <- hits[1]
    }
    if (file.exists(gene_group_csv)) {
      gdf <- tryCatch(data.table::fread(gene_group_csv, header = TRUE, data.table = FALSE),
                     error = function(e) NULL)
      if (!is.null(gdf) && nrow(gdf) > 0) {
        gene_ids <- if ("Gene_ID" %in% colnames(gdf)) trimws(gdf$Gene_ID) else trimws(gdf[[1]])
      }
    }
  }
  if (!is.null(gene_ids) && length(gene_ids) > 0) {
      gene_ids <- gene_ids[nzchar(gene_ids)]
      matched <- match_gene_ids(gene_ids, rownames(full_matrix))
      if (length(matched) == 0) {
        return(list(success = FALSE, reason = paste("no gene IDs matched in M1 matrix for", gene_group)))
      }
      full_matrix <- full_matrix[matched, , drop = FALSE]
      cat("  M1: filtered to", nrow(full_matrix), "genes for", gene_group, "\n")
  }

  if (nrow(full_matrix) < MIN_GENES_DEA) {
    return(list(success = FALSE,
                reason = paste0("too few genes after filtering (", nrow(full_matrix),
                                " < ", MIN_GENES_DEA, ")")))
  }
  list(success = TRUE, data = full_matrix, n_genes = nrow(full_matrix))
}

# ===============================================
# MAIN DEA FUNCTION
# ===============================================

run_differential_expression <- function(config = NULL, matrices_dir = NULL) {
  # M2 HISAT2 De Novo does not produce raw integer counts required by DESeq2.
  # StringTie abundance outputs (TPM/FPKM/coverage) are pre-normalized metrics.
  if (grepl("M2_HISAT2_DeNovo", CURRENT_METHOD, ignore.case = TRUE)) {
    cat("Differential expression is not supported for M2 HISAT2 De Novo.\n")
    cat("M2 produces only TPM/FPKM/coverage (not raw integer counts).\n")
    cat("Use M1 (prepDE.py counts), M3/M4 (Salmon NumReads), or M5 (RSEM expected_count).\n")
    return(invisible(NULL))
  }

  # Get method base directory from environment for config file loading
  method_base_dir <- Sys.getenv("METHOD_BASE_DIR", unset = ".")
  if (is.null(config)) config <- load_runtime_config(method_base_dir)
  ensure_output_dir(DEA_OUT_DIR)

  if (is.null(matrices_dir)) {
    matrices_dir <- file.path(method_base_dir, get_matrices_dir(CURRENT_METHOD))
  }

  print_config_summary("DIFFERENTIAL EXPRESSION ANALYSIS", config)

  # M1 RefGuided: use prepDE.py integer count matrix from deseq2_input/
  is_m1 <- grepl("M1_HISAT2_RefGuided", CURRENT_METHOD, ignore.case = TRUE)

  # For M3/M4/M5: load tximport RDS if available (saved by Matrix_Creation scripts).
  # DESeqDataSetFromTximport preserves transcript-length offsets for more accurate DEA.
  method_type <- get_method_type(CURRENT_METHOD)
  txi_rds_path <- file.path(matrices_dir, config$master_reference,
                             "gene_level", "tximport_gene_level.rds")
  txi_obj <- NULL
  if (method_type %in% c("salmon", "rsem", "star") && file.exists(txi_rds_path)) {
    txi_obj <- tryCatch(readRDS(txi_rds_path), error = function(e) {
      cat("  Warning: Failed to load tximport RDS:", e$message, "\n")
      NULL
    })
    if (!is.null(txi_obj)) {
      cat("  Loaded tximport RDS (transcript-length offsets available)\n")
    }
  } else if (method_type %in% c("salmon", "rsem", "star")) {
    warning("tximport RDS not found at: ", txi_rds_path,
            " — using rounded counts (re-run Matrix_Creation to enable tximport offsets)",
            call. = FALSE)
    cat("  WARNING: tximport RDS not found, falling back to rounded counts\n")
    cat("  Re-run Matrix_Creation to enable transcript-length offset corrections\n")
  }

  successful <- 0
  total <- 0

  if (length(config$gene_groups) == 0) {
    cat("ERROR: No gene groups configured. Check .gene_groups_temp.txt\n")
    return(invisible(NULL))
  }

  # Pre-cache gene group CSV gene IDs (avoids re-reading per contrast pair)
  # Single recursive list.files() call instead of one per gene group
  .gene_group_ids_cache <- list()
  .all_csvs <- list.files(GENE_GROUPS_DIR, pattern = "\\.csv$",
                          recursive = TRUE, full.names = TRUE)
  .csv_lookup <- setNames(.all_csvs, tools::file_path_sans_ext(basename(.all_csvs)))
  for (.gg in config$gene_groups) {
    .gg_csv <- file.path(GENE_GROUPS_DIR, paste0(.gg, ".csv"))
    if (!file.exists(.gg_csv) && .gg %in% names(.csv_lookup)) {
      .gg_csv <- .csv_lookup[[.gg]]
    }
    if (file.exists(.gg_csv)) {
      .gdf <- tryCatch(data.table::fread(.gg_csv, header = TRUE),
                       error = function(e) NULL)
      if (!is.null(.gdf) && nrow(.gdf) > 0) {
        .gene_group_ids_cache[[.gg]] <- if ("Gene_ID" %in% colnames(.gdf)) trimws(.gdf$Gene_ID) else trimws(.gdf[[1]])
      }
    }
  }

  # Pre-load M1 genome-wide matrix once (avoids re-reading per gene group)
  .m1_full_matrix_cache <- NULL
  if (is_m1) {
    .m1_csv <- file.path(matrices_dir, config$master_reference, "deseq2_input", "gene_count_matrix.csv")
    if (file.exists(.m1_csv)) {
      .m1_full_matrix_cache <- tryCatch({
        .m <- as.matrix(data.table::fread(.m1_csv, header = TRUE, data.table = FALSE, check.names = FALSE))
        rownames(.m) <- .m[, 1]
        .m <- .m[, -1, drop = FALSE]
        storage.mode(.m) <- "numeric"
        .m
      }, error = function(e) NULL)
      if (!is.null(.m1_full_matrix_cache))
        cat("  M1: pre-loaded genome matrix (", nrow(.m1_full_matrix_cache), " genes)\n")
    }
  }

  for (gene_group in config$gene_groups) {
    cat("Processing:", gene_group, "\n")

    output_folder_name <- get_output_folder_name(gene_group, CURRENT_DATASET)
    output_dir <- file.path(DEA_OUT_DIR, output_folder_name)
    ensure_output_dir(output_dir)

    if (is_m1) {
      validation <- load_m1_gene_group_counts(gene_group, matrices_dir, config$master_reference,
                                                cached_gene_ids = .gene_group_ids_cache[[gene_group]],
                                                cached_full_matrix = .m1_full_matrix_cache)
    } else {
      # DESeq2 requires raw integer counts (NumReads/expected_count), NOT TPM.
      # get_raw_count_type() returns the appropriate raw count type for each method.
      raw_ct <- get_raw_count_type(CURRENT_METHOD)
      if (is.null(raw_ct)) {
        cat("  Skipped: no raw count type defined for", CURRENT_METHOD, "\n")
        next
      }
      input_file <- build_input_path(gene_group, PROCESSING_LEVELS[1],
                                     raw_ct, "Gene_ID",
                                     matrices_dir, config$master_reference)
      validation <- validate_and_read_matrix(input_file, MIN_GENES_DEA)
    }

    if (!validation$success) {
      cat("  Skipped:", validation$reason, "\n")
      next
    }

    # Create sample info with tissue groups
    sample_ids <- colnames(validation$data)
    sample_tissues <- SAMPLE_LABELS[sample_ids]

    # Warn about unmapped samples (NA labels silently drop from all tissue groups)
    na_mask <- is.na(sample_tissues)
    if (any(na_mask)) {
      cat("  Warning:", sum(na_mask), "sample(s) have no tissue label in SAMPLE_LABELS:",
          paste(sample_ids[na_mask], collapse = ", "), "\n")
      cat("  These samples will be excluded from all pairwise comparisons.\n")
    }

    # Run pairwise comparisons between tissue groups
    group_names <- names(TISSUE_GROUPS)
    if (length(group_names) < 2) {
      cat("  Skipped: fewer than 2 tissue groups defined\n")
      next
    }

    # Pre-cache sample membership per tissue group (avoids repeated %in% per pair)
    samples_by_group <- lapply(TISSUE_GROUPS, function(tissues) {
      sample_ids[sample_tissues %in% tissues]
    })

    for (i in 1:(length(group_names) - 1)) {
      for (j in (i + 1):length(group_names)) {
        total <- total + 1

        group1 <- group_names[i]
        group2 <- group_names[j]

        # Use pre-cached lookups
        samples_g1 <- samples_by_group[[group1]]
        samples_g2 <- samples_by_group[[group2]]
        
        if (length(samples_g1) < 2 || length(samples_g2) < 2) {
          cat("  Skipped", group1, "vs", group2, ": insufficient samples (",
              length(samples_g1), "vs", length(samples_g2), ")\n")
          next
        }

        # Prepare subset — verify samples exist in data columns
        all_samples <- c(samples_g1, samples_g2)
        missing_in_data <- all_samples[!all_samples %in% colnames(validation$data)]
        if (length(missing_in_data) > 0) {
          cat("  Warning:", length(missing_in_data), "sample(s) not found in count matrix:",
              paste(head(missing_in_data, 5), collapse = ", "), "\n")
          all_samples <- all_samples[all_samples %in% colnames(validation$data)]
          samples_g1 <- samples_g1[samples_g1 %in% all_samples]
          samples_g2 <- samples_g2[samples_g2 %in% all_samples]
          if (length(samples_g1) < 2 || length(samples_g2) < 2) {
            cat("  Skipped", group1, "vs", group2, ": insufficient samples after filtering\n")
            next
          }
        }
        count_subset <- validation$data[, all_samples, drop = FALSE]
        
        # Set group2 as reference level so the DESeq2 coefficient is
        # condition_{group1}_vs_{group2}, matching contrast_name direction.
        # Without explicit levels, R uses alphabetical order, making the
        # contrast direction inconsistent with the naming.
        sample_info <- data.frame(
          row.names = all_samples,
          condition = factor(c(rep(group1, length(samples_g1)),
                               rep(group2, length(samples_g2))),
                             levels = c(group2, group1))
        )
        
        contrast_name <- paste0(gene_group, "_", group1, "_vs_", group2)
        cat("  Running:", contrast_name, "\n")

        # Prefer tximport path for M3/M4/M5 (preserves length offsets)
        dds <- NULL
        if (!is.null(txi_obj)) {
          # Use cached gene IDs (loaded once per gene_group, not per contrast)
          gene_ids_for_filter <- .gene_group_ids_cache[[gene_group]]
          dds <- tryCatch(
            prepare_deseq_from_tximport(txi_obj, all_samples, sample_info,
                                        gene_ids = gene_ids_for_filter),
            error = function(e) {
              cat("    tximport path failed, falling back to matrix:", e$message, "\n")
              NULL
            }
          )
        }
        # Fallback to rounded count matrix
        if (is.null(dds)) {
          dds <- prepare_deseq_dataset(count_subset, sample_info)
        }
        if (is.null(dds)) next
        
        result <- run_deseq2(dds, contrast_name, output_dir)
        
        # Save results
        write.table(result$results,
                    file.path(output_dir, paste0(contrast_name, "_results.tsv")),
                    sep = "\t", row.names = FALSE, quote = FALSE)
        
        # Create plots
        create_volcano_plot(result$results, contrast_name, output_dir)
        create_ma_plot(result$results, contrast_name, output_dir)
        create_de_heatmap(result$dds, result$results, contrast_name, output_dir)
        
        n_sig <- sum(result$results$significance != "NS", na.rm = TRUE)
        cat("    Found", n_sig, "significant genes\n")
        
        successful <- successful + 1
      }
    }
  }
  
  print_summary(successful, total)
}

if (!interactive() && identical(environment(), globalenv())) {
  run_differential_expression()
}
