#!/usr/bin/env Rscript

# ===============================================
# GENE SET ENRICHMENT ANALYSIS FOR METHOD 3 (STAR + SALMON)
# ===============================================
# Performs GO/KEGG enrichment analysis using clusterProfiler
# Supports ORA (over-representation) and GSEA (pre-ranked)

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(enrichplot)
  library(ggplot2)
  library(DOSE)
})

script_dir <- dirname(sys.frame(1)$ofile)
source(file.path(script_dir, "0_shared_config.R"))
source(file.path(script_dir, "2_utility_functions.R"))

# ===============================================
# CONFIGURATION
# ===============================================

GSEA_SUBDIR <- "VI_Gene_Set_Enrichment"
GSEA_OUT_DIR <- file.path(CONSOLIDATED_BASE_DIR, GSEA_SUBDIR)

# Enrichment parameters
PVALUE_CUTOFF <- 0.05              # P-value cutoff for enrichment
QVALUE_CUTOFF <- 0.2               # Q-value (FDR) cutoff
MIN_GENE_SET_SIZE <- 5             # Minimum genes in a gene set
MAX_GENE_SET_SIZE <- 500           # Maximum genes in a gene set
N_TOP_TERMS <- 20                  # Number of top terms to display

# GSEA-specific parameters
GSEA_NPERM <- 1000                 # Number of permutations for GSEA
GSEA_SCORE_TYPE <- "std"           # Score type: "std", "pos", "neg"

# Custom gene set file (optional) - tab-separated: term_id, term_name, gene
CUSTOM_GENE_SETS_FILE <- NULL      # e.g., "0_INPUT_FASTAs/docs/custom_gene_sets.tsv"

# Load runtime configuration
config <- load_runtime_config()
ensure_output_dir(GSEA_OUT_DIR)

# FIGURE GENERATION CONTROL
GENERATE_GSEA_FIGURES <- list(
  dotplot = TRUE,                    # Dot plot of enriched terms
  barplot = TRUE,                    # Bar plot of enriched terms
  cnetplot = FALSE,                  # Gene-concept network (can be slow)
  upsetplot = FALSE,                 # UpSet plot of term overlaps
  heatplot = FALSE,                  # Heatmap of genes in terms
  emapplot = TRUE                    # Enrichment map (term similarity network)
)

# ===============================================
# GENE SET LOADING FUNCTIONS
# ===============================================

# Load custom gene sets from file
load_custom_gene_sets <- function(file_path) {
  if (is.null(file_path) || !file.exists(file_path)) {
    return(NULL)
  }
  
  gs_data <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # Expected columns: term_id, term_name, gene
  if (!all(c("term_id", "gene") %in% colnames(gs_data))) {
    cat("  Warning: Custom gene set file must have 'term_id' and 'gene' columns\n")
    return(NULL)
  }
  
  # Convert to TERM2GENE format
  term2gene <- gs_data[, c("term_id", "gene")]
  
  # Create TERM2NAME if term_name exists
  term2name <- NULL
  if ("term_name" %in% colnames(gs_data)) {
    term2name <- unique(gs_data[, c("term_id", "term_name")])
  }
  
  return(list(term2gene = term2gene, term2name = term2name))
}

# Create module-based gene sets from WGCNA results
create_wgcna_gene_sets <- function(wgcna_dir) {
  module_file <- list.files(wgcna_dir, pattern = "_gene_module_assignments.tsv$", 
                            full.names = TRUE, recursive = TRUE)
  
  if (length(module_file) == 0) return(NULL)
  
  all_modules <- data.frame()
  for (f in module_file) {
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
    term_name = gsub("WGCNA_", "WGCNA Module: ", unique(all_modules$term_id))
  )
  
  return(list(term2gene = all_modules, term2name = term2name))
}

# ===============================================
# OVER-REPRESENTATION ANALYSIS (ORA)
# ===============================================

run_ora_enrichment <- function(gene_list, universe, gene_sets, term2name = NULL,
                                analysis_name, output_dir) {
  cat("    Running ORA for", length(gene_list), "genes\n")
  
  if (length(gene_list) < 3) {
    cat("    Skipped: too few genes\n")
    return(NULL)
  }
  
  # Run enrichment
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
    cat("    No significant enrichments found\n")
    return(NULL)
  }
  
  # Export results
  result_df <- as.data.frame(enrich_result)
  write.table(result_df,
              file = file.path(output_dir, paste0(analysis_name, "_ORA_results.tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Generate plots
  generate_enrichment_plots(enrich_result, analysis_name, output_dir, "ORA")
  
  cat("    Found", nrow(result_df), "enriched terms\n")
  return(enrich_result)
}

# ===============================================
# GENE SET ENRICHMENT ANALYSIS (GSEA)
# ===============================================

run_gsea_enrichment <- function(ranked_genes, gene_sets, term2name = NULL,
                                 analysis_name, output_dir) {
  cat("    Running GSEA on", length(ranked_genes), "ranked genes\n")
  
  # Sort by rank (descending)
  ranked_genes <- sort(ranked_genes, decreasing = TRUE)
  
  # Run GSEA
  gsea_result <- tryCatch({
    GSEA(
      geneList = ranked_genes,
      TERM2GENE = gene_sets,
      TERM2NAME = term2name,
      pvalueCutoff = PVALUE_CUTOFF,
      minGSSize = MIN_GENE_SET_SIZE,
      maxGSSize = MAX_GENE_SET_SIZE,
      nPermSimple = GSEA_NPERM,
      scoreType = GSEA_SCORE_TYPE,
      verbose = FALSE
    )
  }, error = function(e) {
    cat("    GSEA failed:", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(gsea_result) || nrow(gsea_result@result) == 0) {
    cat("    No significant GSEA results\n")
    return(NULL)
  }
  
  # Export results
  result_df <- as.data.frame(gsea_result)
  write.table(result_df,
              file = file.path(output_dir, paste0(analysis_name, "_GSEA_results.tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Generate plots
  generate_enrichment_plots(gsea_result, analysis_name, output_dir, "GSEA")
  
  # GSEA-specific running score plot for top terms
  if (nrow(gsea_result@result) > 0) {
    top_term <- gsea_result@result$ID[1]
    tryCatch({
      p <- gseaplot2(gsea_result, geneSetID = top_term, title = top_term)
      ggsave(file.path(output_dir, paste0(analysis_name, "_GSEA_running_score.png")),
             p, width = 10, height = 8, dpi = 150)
    }, error = function(e) NULL)
  }
  
  cat("    Found", nrow(result_df), "enriched terms\n")
  return(gsea_result)
}

# ===============================================
# VISUALIZATION FUNCTIONS
# ===============================================

generate_enrichment_plots <- function(enrich_result, analysis_name, output_dir, method) {
  n_terms <- min(N_TOP_TERMS, nrow(enrich_result@result))
  if (n_terms < 1) return()
  
  # Dot plot
  if (GENERATE_GSEA_FIGURES$dotplot) {
    tryCatch({
      p <- dotplot(enrich_result, showCategory = n_terms) +
        ggtitle(paste0(method, ": ", analysis_name)) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))
      ggsave(file.path(output_dir, paste0(analysis_name, "_", method, "_dotplot.png")),
             p, width = 10, height = 8, dpi = 150)
    }, error = function(e) NULL)
  }
  
  # Bar plot
  if (GENERATE_GSEA_FIGURES$barplot) {
    tryCatch({
      p <- barplot(enrich_result, showCategory = n_terms) +
        ggtitle(paste0(method, ": ", analysis_name)) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))
      ggsave(file.path(output_dir, paste0(analysis_name, "_", method, "_barplot.png")),
             p, width = 10, height = 8, dpi = 150)
    }, error = function(e) NULL)
  }
  
  # Enrichment map
  if (GENERATE_GSEA_FIGURES$emapplot && nrow(enrich_result@result) >= 3) {
    tryCatch({
      enrich_result_pw <- pairwise_termsim(enrich_result)
      p <- emapplot(enrich_result_pw, showCategory = min(30, n_terms)) +
        ggtitle(paste0(method, " Term Network: ", analysis_name))
      ggsave(file.path(output_dir, paste0(analysis_name, "_", method, "_emapplot.png")),
             p, width = 12, height = 10, dpi = 150)
    }, error = function(e) NULL)
  }
  
  # Gene-concept network
  if (GENERATE_GSEA_FIGURES$cnetplot && nrow(enrich_result@result) >= 2) {
    tryCatch({
      p <- cnetplot(enrich_result, showCategory = min(5, n_terms), 
                    categorySize = "pvalue", foldChange = NULL) +
        ggtitle(paste0(method, " Gene Network: ", analysis_name))
      ggsave(file.path(output_dir, paste0(analysis_name, "_", method, "_cnetplot.png")),
             p, width = 14, height = 12, dpi = 150)
    }, error = function(e) NULL)
  }
}

# ===============================================
# INTEGRATION WITH DEA RESULTS
# ===============================================

run_enrichment_from_dea <- function(dea_results_file, gene_sets, term2name,
                                     analysis_name, output_dir) {
  dea_data <- read.table(dea_results_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # Get significant genes
  if ("significance" %in% colnames(dea_data)) {
    sig_up <- dea_data$gene[dea_data$significance == "Up"]
    sig_down <- dea_data$gene[dea_data$significance == "Down"]
    all_sig <- c(sig_up, sig_down)
  } else {
    all_sig <- dea_data$gene[dea_data$padj < 0.05 & abs(dea_data$log2FoldChange) > 1]
  }
  
  universe <- dea_data$gene
  
  # Run ORA for all significant genes
  if (length(all_sig) >= 3) {
    run_ora_enrichment(all_sig, universe, gene_sets, term2name,
                       paste0(analysis_name, "_all_DE"), output_dir)
  }
  
  # Run separate ORA for up/down regulated
  if (exists("sig_up") && length(sig_up) >= 3) {
    run_ora_enrichment(sig_up, universe, gene_sets, term2name,
                       paste0(analysis_name, "_upregulated"), output_dir)
  }
  if (exists("sig_down") && length(sig_down) >= 3) {
    run_ora_enrichment(sig_down, universe, gene_sets, term2name,
                       paste0(analysis_name, "_downregulated"), output_dir)
  }
  
  # Run GSEA with ranked list (by stat or log2FC * -log10(pvalue))
  if ("stat" %in% colnames(dea_data)) {
    ranked <- setNames(dea_data$stat, dea_data$gene)
  } else if (all(c("log2FoldChange", "pvalue") %in% colnames(dea_data))) {
    dea_data$pvalue[dea_data$pvalue == 0] <- min(dea_data$pvalue[dea_data$pvalue > 0]) / 10
    ranked <- setNames(dea_data$log2FoldChange * -log10(dea_data$pvalue), dea_data$gene)
  } else {
    return()
  }
  
  ranked <- ranked[!is.na(ranked) & is.finite(ranked)]
  run_gsea_enrichment(ranked, gene_sets, term2name, analysis_name, output_dir)
}

# ===============================================
# MAIN EXECUTION
# ===============================================

run_gene_set_enrichment_analysis <- function() {
  print_separator()
  cat("GENE SET ENRICHMENT ANALYSIS\n")
  print_separator()
  
  # Load gene sets
  custom_gs <- load_custom_gene_sets(CUSTOM_GENE_SETS_FILE)
  wgcna_gs <- create_wgcna_gene_sets(file.path(CONSOLIDATED_BASE_DIR, "III_Coexpression_WGCNA"))
  
  # Combine gene sets
  all_gene_sets <- NULL
  all_term2name <- NULL
  
  if (!is.null(custom_gs)) {
    all_gene_sets <- custom_gs$term2gene
    all_term2name <- custom_gs$term2name
    cat("Loaded custom gene sets:", nrow(unique(custom_gs$term2gene[, 1, drop = FALSE])), "terms\n")
  }
  
  if (!is.null(wgcna_gs)) {
    if (is.null(all_gene_sets)) {
      all_gene_sets <- wgcna_gs$term2gene
      all_term2name <- wgcna_gs$term2name
    } else {
      all_gene_sets <- rbind(all_gene_sets, wgcna_gs$term2gene)
      all_term2name <- rbind(all_term2name, wgcna_gs$term2name)
    }
    cat("Loaded WGCNA module gene sets:", length(unique(wgcna_gs$term2gene$term_id)), "modules\n")
  }
  
  if (is.null(all_gene_sets)) {
    cat("\nNo gene sets available. Please provide:\n")
    cat("  1. Custom gene sets file (CUSTOM_GENE_SETS_FILE)\n")
    cat("  2. Or run WGCNA first to create module-based gene sets\n")
    return()
  }
  
  # Check for DEA results and run enrichment
  dea_dir <- file.path(CONSOLIDATED_BASE_DIR, "V_Differential_Expression")
  if (dir.exists(dea_dir)) {
    dea_files <- list.files(dea_dir, pattern = "_all_results.tsv$", 
                            recursive = TRUE, full.names = TRUE)
    
    cat("\nProcessing", length(dea_files), "DEA result files\n")
    
    for (dea_file in dea_files) {
      analysis_name <- gsub("_all_results.tsv$", "", basename(dea_file))
      group_dir <- dirname(dea_file)
      out_dir <- file.path(group_dir, "enrichment")
      ensure_output_dir(out_dir)
      
      cat("\n  Processing:", analysis_name, "\n")
      run_enrichment_from_dea(dea_file, all_gene_sets, all_term2name, analysis_name, out_dir)
    }
  } else {
    cat("\nNo DEA results found. Run differential expression analysis first.\n")
  }
  
  print_separator()
  cat("GSEA complete. Output:", GSEA_OUT_DIR, "\n")
  print_separator()
}

# Run if executed directly
if (!interactive()) {
  run_gene_set_enrichment_analysis()
}
