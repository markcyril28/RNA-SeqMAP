#!/usr/bin/env Rscript

# ===============================================
# GENE SET ENRICHMENT ANALYSIS MODULE
# ===============================================
# GO/KEGG enrichment using clusterProfiler

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(enrichplot)
  library(ggplot2)
  library(DOSE)
})

SCRIPT_DIR <- Sys.getenv("ANALYSIS_MODULES_DIR", ".")
source(file.path(SCRIPT_DIR, "0_shared_config.R"))
source(file.path(SCRIPT_DIR, "1_utility_functions.R"))

# ===============================================
# CONFIGURATION
# ===============================================

GSEA_OUT_DIR <- file.path(CONSOLIDATED_BASE_DIR, OUTPUT_SUBDIRS$GSEA)

# Enrichment parameters
PVALUE_CUTOFF <- 0.05
QVALUE_CUTOFF <- 0.2
MIN_GENE_SET_SIZE <- 5
MAX_GENE_SET_SIZE <- 500
N_TOP_TERMS <- 20

# GSEA-specific
GSEA_NPERM <- 1000
# Score type: "std" for standard, "pos" for only positive (upregulated), "neg" for only negative
GSEA_SCORE_TYPE <- "std"  

# ID mapping (if gene names don't match gene set database)
# Set to TRUE if using gene symbols from annotation vs gene IDs
NEED_ID_CONVERSION <- FALSE

# Custom gene set file (optional)
CUSTOM_GENE_SETS_FILE <- NULL

# Figure toggles
GENERATE_GSEA_FIGURES <- list(
  dotplot = TRUE,
  barplot = TRUE,
  cnetplot = FALSE,
  emapplot = TRUE
)

# ===============================================
# GENE SET LOADING
# ===============================================

load_custom_gene_sets <- function(file_path) {
  if (is.null(file_path) || !file.exists(file_path)) return(NULL)
  
  gs_data <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  if (!all(c("term_id", "gene") %in% colnames(gs_data))) {
    cat("  Warning: Custom gene set needs 'term_id' and 'gene' columns\n")
    return(NULL)
  }
  
  term2gene <- gs_data[, c("term_id", "gene")]
  term2name <- NULL
  if ("term_name" %in% colnames(gs_data)) {
    term2name <- unique(gs_data[, c("term_id", "term_name")])
  }
  
  return(list(term2gene = term2gene, term2name = term2name))
}

create_wgcna_gene_sets <- function(wgcna_dir) {
  module_files <- list.files(wgcna_dir, pattern = "_module_assignments.tsv$",
                              full.names = TRUE, recursive = TRUE)
  
  if (length(module_files) == 0) return(NULL)
  
  all_modules <- data.frame()
  for (f in module_files) {
    mod_data <- read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    if (all(c("Gene", "Module") %in% colnames(mod_data))) {
      mod_data$term_id <- paste0("WGCNA_", mod_data$Module)
      all_modules <- rbind(all_modules, mod_data[, c("term_id", "Gene")])
    }
  }
  
  if (nrow(all_modules) == 0) return(NULL)
  
  colnames(all_modules) <- c("term_id", "gene")
  term2name <- data.frame(
    term_id = unique(all_modules$term_id),
    term_name = gsub("WGCNA_", "Module: ", unique(all_modules$term_id))
  )
  
  return(list(term2gene = all_modules, term2name = term2name))
}

# ===============================================
# ENRICHMENT FUNCTIONS
# ===============================================

run_ora_enrichment <- function(gene_list, universe, gene_sets, term2name,
                                analysis_name, output_dir) {
  cat("    Running ORA for", length(gene_list), "genes\n")
  
  if (length(gene_list) < 3) {
    cat("    Skipped: too few genes\n")
    return(NULL)
  }
  
  enrich_result <- enricher(
    gene = gene_list,
    universe = universe,
    TERM2GENE = gene_sets,
    TERM2NAME = term2name,
    pvalueCutoff = PVALUE_CUTOFF,
    qvalueCutoff = QVALUE_CUTOFF,
    minGSSize = MIN_GENE_SET_SIZE,
    maxGSSize = MAX_GENE_SET_SIZE
  )
  
  if (is.null(enrich_result) || nrow(enrich_result@result) == 0) {
    cat("    No significant enrichments\n")
    return(NULL)
  }
  
  result_df <- as.data.frame(enrich_result)
  write.table(result_df,
              file.path(output_dir, paste0(analysis_name, "_ORA_results.tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  generate_enrichment_plots(enrich_result, analysis_name, output_dir)
  
  cat("    Found", nrow(result_df), "enriched terms\n")
  return(enrich_result)
}

run_gsea_analysis <- function(ranked_genes, gene_sets, term2name,
                               analysis_name, output_dir) {
  cat("    Running GSEA\n")
  
  if (length(ranked_genes) < 10) {
    cat("    Skipped: too few genes\n")
    return(NULL)
  }
  
  set.seed(GLOBAL_RANDOM_SEED)  # Reproducibility for permutation test
  gsea_result <- GSEA(
    ranked_genes,
    TERM2GENE = gene_sets,
    TERM2NAME = term2name,
    pvalueCutoff = PVALUE_CUTOFF,
    minGSSize = MIN_GENE_SET_SIZE,
    maxGSSize = MAX_GENE_SET_SIZE,
    nPermSimple = GSEA_NPERM,
    scoreType = GSEA_SCORE_TYPE
  )
  
  if (is.null(gsea_result) || nrow(gsea_result@result) == 0) {
    cat("    No significant enrichments\n")
    return(NULL)
  }
  
  result_df <- as.data.frame(gsea_result)
  write.table(result_df,
              file.path(output_dir, paste0(analysis_name, "_GSEA_results.tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  generate_enrichment_plots(gsea_result, analysis_name, output_dir, is_gsea = TRUE)
  
  cat("    Found", nrow(result_df), "enriched terms\n")
  return(gsea_result)
}

generate_enrichment_plots <- function(enrich_result, analysis_name, output_dir, is_gsea = FALSE) {
  n_terms <- min(N_TOP_TERMS, nrow(enrich_result@result))
  
  if (GENERATE_GSEA_FIGURES$dotplot) {
    p <- dotplot(enrich_result, showCategory = n_terms) +
      ggtitle(paste0(analysis_name, " - Top Enriched Terms")) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    ggsave(file.path(output_dir, paste0(analysis_name, "_dotplot.png")),
           p, width = 10, height = 8, dpi = 150)
  }
  
  if (GENERATE_GSEA_FIGURES$barplot) {
    p <- barplot(enrich_result, showCategory = n_terms) +
      ggtitle(paste0(analysis_name, " - Enrichment Barplot")) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    ggsave(file.path(output_dir, paste0(analysis_name, "_barplot.png")),
           p, width = 10, height = 8, dpi = 150)
  }
  
  if (GENERATE_GSEA_FIGURES$emapplot && n_terms >= 2) {
    tryCatch({
      enrich_result2 <- pairwise_termsim(enrich_result)
      p <- emapplot(enrich_result2, showCategory = n_terms) +
        ggtitle(paste0(analysis_name, " - Enrichment Map"))
      ggsave(file.path(output_dir, paste0(analysis_name, "_emapplot.png")),
             p, width = 12, height = 10, dpi = 150)
    }, error = function(e) {
      cat("    Emapplot failed:", e$message, "\n")
    })
  }
}

# ===============================================
# MAIN GSEA FUNCTION
# ===============================================

run_gene_set_enrichment <- function(config = NULL, matrices_dir = NULL) {
  # Get method base directory from environment for config file loading
  method_base_dir <- Sys.getenv("METHOD_BASE_DIR", unset = ".")
  if (is.null(config)) config <- load_runtime_config(method_base_dir)

  # Set up matrices_dir based on method_base_dir if not provided
  if (is.null(matrices_dir)) {
    matrices_dir <- file.path(method_base_dir, get_matrices_dir(CURRENT_METHOD))
  }

  ensure_output_dir(GSEA_OUT_DIR)
  
  print_config_summary("GENE SET ENRICHMENT ANALYSIS", config)
  
  # Load gene sets
  gene_sets <- NULL
  term2name <- NULL
  
  if (!is.null(CUSTOM_GENE_SETS_FILE)) {
    custom <- load_custom_gene_sets(CUSTOM_GENE_SETS_FILE)
    if (!is.null(custom)) {
      gene_sets <- custom$term2gene
      term2name <- custom$term2name
    }
  }
  
  # Try WGCNA modules as gene sets
  wgcna_dir <- file.path(CONSOLIDATED_BASE_DIR, OUTPUT_SUBDIRS$WGCNA)
  if (is.null(gene_sets) && dir.exists(wgcna_dir)) {
    wgcna_sets <- create_wgcna_gene_sets(wgcna_dir)
    if (!is.null(wgcna_sets)) {
      gene_sets <- wgcna_sets$term2gene
      term2name <- wgcna_sets$term2name
      cat("Using WGCNA modules as gene sets\n")
    }
  }
  
  if (is.null(gene_sets)) {
    cat("No gene sets available. Provide custom gene sets or run WGCNA first.\n")
    return(NULL)
  }

  # Build universe from FULL TRANSCRIPTOME (all measured genes), not gene set database.
  # Using the gene set database as universe inflates significance and produces false positives.
  # IMPORTANT: Always use "Gene_ID" gene type so that universe IDs match gene set IDs
  # (gene sets from WGCNA modules use raw Gene_IDs, not Shortened_Names).
  full_file <- build_input_path(config$master_reference, PROCESSING_LEVELS[1],
                                COUNT_TYPES[1], "Gene_ID",
                                matrices_dir, config$master_reference)
  universe <- NULL
  if (file.exists(full_file)) {
    full_data <- read_count_matrix(full_file)
    universe <- rownames(full_data)
    cat("Universe: ", length(universe), " measured genes (from full transcriptome matrix)\n")
  } else {
    cat("  Warning: Full transcriptome matrix not found at:", full_file, "\n")
    cat("  Falling back to gene set database genes as universe (may inflate significance)\n")
    universe <- unique(gene_sets$gene)
  }

  successful <- 0
  total <- 0

  for (gene_group in config$gene_groups) {
    cat("Processing:", gene_group, "\n")
    total <- total + 1

    output_folder_name <- get_output_folder_name(gene_group, CURRENT_DATASET)
    output_dir <- file.path(GSEA_OUT_DIR, output_folder_name)
    ensure_output_dir(output_dir)

    # IMPORTANT: Always use "Gene_ID" gene type for query genes to match
    # universe and gene set IDs (which are always in Gene_ID format).
    input_file <- build_input_path(gene_group, PROCESSING_LEVELS[1],
                                   COUNT_TYPES[1], "Gene_ID",
                                   matrices_dir, config$master_reference)

    validation <- validate_and_read_matrix(input_file, 5)
    if (!validation$success) {
      cat("  Skipped:", validation$reason, "\n")
      next
    }

    # Use gene IDs as the list for ORA (must match universe and gene set namespace)
    query_genes <- rownames(validation$data)

    result <- run_ora_enrichment(query_genes, universe, gene_sets, term2name,
                                  gene_group, output_dir)
    
    if (!is.null(result)) successful <- successful + 1
  }
  
  print_summary(successful, total)
}

if (!interactive() && identical(environment(), globalenv())) {
  run_gene_set_enrichment()
}
