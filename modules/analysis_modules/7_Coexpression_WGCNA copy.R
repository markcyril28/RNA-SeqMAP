#!/usr/bin/env Rscript

# ===============================================
# WGCNA COEXPRESSION ANALYSIS MODULE
# ===============================================
# Weighted gene co-expression network analysis
# Runs WGCNA on FULL TRANSCRIPTOME with query genes highlighted

suppressPackageStartupMessages({
  library(WGCNA)
  library(dynamicTreeCut)
  library(fastcluster)
  library(ggplot2)
  library(RColorBrewer)
  library(igraph)
  # Interactive HTML output libraries
  library(networkD3)
  library(htmlwidgets)
  library(visNetwork)
  library(jsonlite)  # For embedding data in HTML
})

allowWGCNAThreads()

SCRIPT_DIR <- Sys.getenv("ANALYSIS_MODULES_DIR", ".")
source(file.path(SCRIPT_DIR, "0_shared_config.R"))
source(file.path(SCRIPT_DIR, "1_utility_functions.R"))

# ===============================================
# OUTPUT FORMAT TOGGLES (SET THESE!)
# ===============================================
# Toggle which outputs to generate:
OUTPUT_PNG <- TRUE                 # Static PNG images
OUTPUT_HTML <- TRUE                # Interactive HTML (zoomable, hoverable, with parameter toggles)
OUTPUT_HTML_SIMPLE <- TRUE         # Simple HTML (fast loading, minimal UI, single network view)
OUTPUT_RAW <- TRUE                 # Raw data files (TSV, RDS) for downstream analysis

# ===============================================
# CONFIGURATION
# ===============================================

WGCNA_OUT_DIR <- file.path(CONSOLIDATED_BASE_DIR, OUTPUT_SUBDIRS$WGCNA)

# WGCNA parameters (STRINGENT settings for high-confidence results)
SOFT_POWER_RANGE <- 1:30           # Extended range for optimal fit
MIN_MODULE_SIZE_DEFAULT <- 50      # STRINGENT: Larger modules only
MIN_MODULE_SIZE_SMALL <- 10        # STRINGENT: Min 10 genes per module
MERGE_CUT_HEIGHT <- 0.15           # STRINGENT: More conservative module merging
NETWORK_TYPE <- "signed"           # "signed" preserves biological interpretation
DEEP_SPLIT <- 2                    # STRINGENT: Moderate sensitivity (range 0-4)
PAM_STAGE <- TRUE                  # PAM for more accurate module assignment

# Default coexpression parameters (used for PNG output)
N_HUB_GENES <- 5                   # Top hub genes per module
N_COEXPRESSED_GENES <- 50          # Top co-expressed genes per query
N_NETWORK_GENES <- 30              # Genes in network visualization
COR_THRESHOLD <- 0.90              # Default correlation threshold
TOP_VAR_GENES <- 5000              # Top N variable genes from transcriptome

# ===============================================
# PARAMETER COMBINATIONS FOR INTERACTIVE HTML
# ===============================================
# These are precalculated and embedded in HTML for user toggling
PARAM_COMBINATIONS <- list(
  # Correlation thresholds to precalculate
  cor_thresholds = c(0.70, 0.80, 0.90, 0.95),
  
  # Network sizes (number of top correlated genes per query gene)
  network_sizes = c(10, 30, 50, 100),
  
  # Query-only network thresholds
  query_network_thresholds = c(0.50, 0.70, 0.85)
)

# Figure toggles
GENERATE_WGCNA_FIGURES <- list(
  soft_threshold_plot = TRUE,
  module_dendrogram = TRUE,
  eigengene_adjacency = TRUE,
  correlation_network = TRUE,
  query_gene_focus = TRUE
)

# Legacy compatibility (kept for backward compatibility)
OUTPUT_FORMATS <- list(
  png = OUTPUT_PNG,
  html = OUTPUT_HTML
)

# ===============================================
# WGCNA CORE FUNCTIONS
# ===============================================

pick_soft_threshold <- function(data_matrix, output_dir, gene_group) {
  powers <- SOFT_POWER_RANGE
  sft <- pickSoftThreshold(data_matrix, powerVector = powers, 
                           networkType = NETWORK_TYPE, verbose = 0)
  
  if (GENERATE_WGCNA_FIGURES$soft_threshold_plot) {
    png(file.path(output_dir, paste0(gene_group, "_soft_threshold.png")), 
        width = 1000, height = 500, res = 100)
    par(mfrow = c(1, 2))
    fit_index <- -sign(sft$fitIndices[,3]) * sft$fitIndices[,2]
    plot(sft$fitIndices[,1], fit_index,
         xlab = "Soft Threshold", ylab = "R^2",
         type = "n", main = "Scale Independence")
    text(sft$fitIndices[,1], fit_index, labels = powers, col = "red")
    abline(h = 0.80, col = "red")
    plot(sft$fitIndices[,1], sft$fitIndices[,5],
         xlab = "Soft Threshold", ylab = "Mean Connectivity",
         type = "n", main = "Mean Connectivity")
    text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, col = "red")
    dev.off()
  }
  
  fit_index <- -sign(sft$fitIndices[,3]) * sft$fitIndices[,2]
  power_selected <- which(fit_index > 0.80)[1]
  if (is.na(power_selected)) power_selected <- which.max(fit_index)
  
  return(powers[power_selected])
}

build_network_and_detect_modules <- function(data_matrix, soft_power, output_dir, 
                                              gene_group, query_genes = NULL) {
  # Dynamic module size based on gene count
  n_genes <- ncol(data_matrix)  # Genes are columns after transpose for WGCNA
  min_module_size <- if (n_genes < 100) MIN_MODULE_SIZE_SMALL else MIN_MODULE_SIZE_DEFAULT
  
  # Use PAM and deepSplit for more accurate module detection
  deep_split_val <- if (exists("DEEP_SPLIT")) DEEP_SPLIT else 2
  pam_stage_val <- if (exists("PAM_STAGE")) PAM_STAGE else TRUE
  
  net <- blockwiseModules(
    data_matrix, 
    power = soft_power,
    networkType = NETWORK_TYPE,
    TOMType = "signed",
    minModuleSize = min_module_size,
    reassignThreshold = 0,
    mergeCutHeight = MERGE_CUT_HEIGHT,
    deepSplit = deep_split_val,        # Sensitivity for module detection (0-4)
    pamStage = pam_stage_val,          # PAM for accurate gene assignment
    pamRespectsDendro = TRUE,          # PAM respects dendrogram structure
    numericLabels = TRUE,
    saveTOMs = TRUE,                   # Save TOM for downstream analysis
    saveTOMFileBase = file.path(output_dir, paste0(gene_group, "_TOM")),
    verbose = 1                        # Show progress
  )
  
  module_colors <- labels2colors(net$colors)
  
  if (GENERATE_WGCNA_FIGURES$module_dendrogram) {
    # Create query gene indicator bar (red = query gene, white = other)
    all_genes <- colnames(data_matrix)
    query_indicator <- rep("white", length(all_genes))
    if (!is.null(query_genes)) {
      matched_query <- which(all_genes %in% query_genes)
      query_indicator[matched_query] <- "red"
    }
    
    # Stack module colors and query indicator
    color_matrix <- cbind(module_colors[net$blockGenes[[1]]], 
                          query_indicator[net$blockGenes[[1]]])
    colnames(color_matrix) <- c("Module", "Query Genes")
    
    # ===== PNG OUTPUT: MODULE DENDROGRAM =====
    if (OUTPUT_PNG) {
      png(file.path(output_dir, paste0(gene_group, "_module_dendrogram.png")), 
          width = 1600, height = 800, res = 100)
      
      plotDendroAndColors(net$dendrograms[[1]], color_matrix,
                          dendroLabels = FALSE, hang = 0.03,
                          addGuide = TRUE, guideHang = 0.05,
                          main = paste0(gene_group, " - Module Dendrogram (Query genes in RED)"))
      dev.off()
    }
    
    # ===== HTML OUTPUT: INTERACTIVE MODULE SUMMARY =====
    if (OUTPUT_HTML) {
      # Create module summary table as interactive HTML
      module_summary <- data.frame(
        Gene = all_genes,
        Module = module_colors,
        Is_Query = query_indicator == "red",
        stringsAsFactors = FALSE
      )
      
      # Count genes per module
      module_counts <- as.data.frame(table(module_colors))
      colnames(module_counts) <- c("Module", "Gene_Count")
      module_counts$Has_Query <- sapply(module_counts$Module, function(m) {
        sum(module_summary$Module == m & module_summary$Is_Query) > 0
      })
      
      # Create HTML table with module info
      html_content <- paste0(
        "<!DOCTYPE html>\n<html>\n<head>\n",
        "<title>", gene_group, " - Module Summary</title>\n",
        "<style>\n",
        "body { font-family: Arial, sans-serif; margin: 20px; }\n",
        "h1 { color: #333; }\n",
        "table { border-collapse: collapse; width: 100%; margin: 20px 0; }\n",
        "th, td { border: 1px solid #ddd; padding: 12px; text-align: left; }\n",
        "th { background-color: #4CAF50; color: white; }\n",
        "tr:nth-child(even) { background-color: #f2f2f2; }\n",
        "tr:hover { background-color: #ddd; }\n",
        ".query { background-color: #ffcccc; font-weight: bold; }\n",
        ".module-tag { padding: 5px 10px; border-radius: 5px; color: white; }\n",
        "</style>\n",
        "</head>\n<body>\n",
        "<h1>", gene_group, " - WGCNA Module Summary</h1>\n",
        "<h2>Module Counts</h2>\n",
        "<table>\n<tr><th>Module</th><th>Gene Count</th><th>Contains Query Genes</th></tr>\n"
      )
      
      for (i in 1:nrow(module_counts)) {
        html_content <- paste0(html_content,
          "<tr><td><span class='module-tag' style='background-color:", module_counts$Module[i], ";'>",
          module_counts$Module[i], "</span></td><td>", module_counts$Gene_Count[i], 
          "</td><td>", ifelse(module_counts$Has_Query[i], "✓ YES", "No"), "</td></tr>\n"
        )
      }
      
      html_content <- paste0(html_content, "</table>\n",
        "<h2>Query Genes</h2>\n",
        "<table>\n<tr><th>Gene</th><th>Module</th></tr>\n"
      )
      
      query_subset <- module_summary[module_summary$Is_Query, ]
      for (i in 1:nrow(query_subset)) {
        html_content <- paste0(html_content,
          "<tr class='query'><td>", query_subset$Gene[i], 
          "</td><td><span class='module-tag' style='background-color:", query_subset$Module[i], ";'>",
          query_subset$Module[i], "</span></td></tr>\n"
        )
      }
      
      html_content <- paste0(html_content, "</table>\n</body>\n</html>")
      
      writeLines(html_content, file.path(output_dir, paste0(gene_group, "_module_summary.html")))
      cat("  Saved HTML module summary:", paste0(gene_group, "_module_summary.html"), "\n")
    }
  }
  
  return(list(net = net, module_colors = module_colors))
}

calculate_module_eigengenes <- function(data_matrix, module_colors, output_dir, gene_group) {
  ME_list <- moduleEigengenes(data_matrix, colors = module_colors)
  MEs <- ME_list$eigengenes
  
  if (ncol(MEs) >= 2 && GENERATE_WGCNA_FIGURES$eigengene_adjacency) {
    # ===== PNG OUTPUT: EIGENGENE ADJACENCY =====
    if (OUTPUT_PNG) {
      png(file.path(output_dir, paste0(gene_group, "_eigengene_adjacency.png")), 
          width = 800, height = 800, res = 100)
      plotEigengeneNetworks(MEs, "", marDendro = c(0, 4, 1, 2), marHeatmap = c(3, 4, 1, 2),
                            plotDendrograms = TRUE, xLabelsAngle = 90)
      dev.off()
    }
    
    # ===== HTML OUTPUT: INTERACTIVE EIGENGENE HEATMAP =====
    if (OUTPUT_HTML) {
      # Compute eigengene correlation
      ME_cor <- cor(MEs, use = "pairwise.complete.obs")
      
      # Create heatmaply interactive heatmap if available, otherwise use D3
      tryCatch({
        if (requireNamespace("heatmaply", quietly = TRUE)) {
          hm <- heatmaply::heatmaply(
            ME_cor,
            main = paste0(gene_group, " - Module Eigengene Correlations"),
            xlab = "Module",
            ylab = "Module",
            colors = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
            dendrogram = "both",
            margins = c(100, 100, 50, 50)
          )
          htmlwidgets::saveWidget(hm, 
                                  file.path(output_dir, paste0(gene_group, "_eigengene_heatmap.html")),
                                  selfcontained = TRUE)
          cat("  Saved interactive eigengene heatmap HTML\n")
        } else {
          # Fallback: simple HTML table
          html_content <- paste0(
            "<!DOCTYPE html>\n<html>\n<head>\n",
            "<title>", gene_group, " - Eigengene Correlations</title>\n",
            "<style>\n",
            "body { font-family: Arial, sans-serif; margin: 20px; }\n",
            "table { border-collapse: collapse; margin: 20px auto; }\n",
            "th, td { border: 1px solid #ddd; padding: 8px; text-align: center; min-width: 60px; }\n",
            "th { background-color: #4CAF50; color: white; }\n",
            "</style>\n</head>\n<body>\n",
            "<h1>", gene_group, " - Module Eigengene Correlations</h1>\n",
            "<table>\n<tr><th></th>"
          )
          for (col in colnames(ME_cor)) {
            html_content <- paste0(html_content, "<th>", col, "</th>")
          }
          html_content <- paste0(html_content, "</tr>\n")
          for (i in 1:nrow(ME_cor)) {
            html_content <- paste0(html_content, "<tr><th>", rownames(ME_cor)[i], "</th>")
            for (j in 1:ncol(ME_cor)) {
              val <- round(ME_cor[i, j], 2)
              # Color based on correlation
              if (val > 0.5) bg <- "#f4a582"
              else if (val < -0.5) bg <- "#92c5de"
              else bg <- "#f7f7f7"
              html_content <- paste0(html_content, 
                "<td style='background-color:", bg, ";'>", val, "</td>")
            }
            html_content <- paste0(html_content, "</tr>\n")
          }
          html_content <- paste0(html_content, "</table>\n</body>\n</html>")
          writeLines(html_content, file.path(output_dir, paste0(gene_group, "_eigengene_correlations.html")))
        }
      }, error = function(e) {
        cat("  Note: Could not create interactive eigengene heatmap:", e$message, "\n")
      })
    }
  }
  
  return(MEs)
}

identify_hub_genes <- function(gene_info, n_top = N_HUB_GENES) {
  modules <- unique(gene_info$Module)
  modules <- modules[modules != "grey"]
  
  hub_genes <- data.frame()
  for (mod in modules) {
    mod_genes <- gene_info[gene_info$Module == mod, ]
    kME_col <- paste0("kME", mod)
    if (kME_col %in% colnames(mod_genes)) {
      mod_genes <- mod_genes[order(-mod_genes[[kME_col]]), ]
      top_n <- min(n_top, nrow(mod_genes))
      hub_genes <- rbind(hub_genes, mod_genes[1:top_n, ])
    }
  }
  
  return(hub_genes)
}

create_correlation_network <- function(data_matrix, query_genes, output_dir, gene_group) {
  if (!GENERATE_WGCNA_FIGURES$correlation_network) return(NULL)
  
  # Compute correlation matrix (use GPU if available)
  cor_matrix <- gpu_cor(t(data_matrix))
  
  # Filter to query genes and their top correlated genes
  matched <- query_genes[query_genes %in% rownames(cor_matrix)]
  if (length(matched) == 0) {
    cat("  Warning: No query genes found in correlation matrix\n")
    return(NULL)
  }
  
  cat("  Query genes in network:", length(matched), "\n")
  
  # ===== PRECALCULATE ALL PARAMETER COMBINATIONS =====
  cat("  Precalculating parameter combinations for interactive HTML...\n")
  
  # Store all precalculated networks
  precalc_networks <- list()
  
  for (net_size in PARAM_COMBINATIONS$network_sizes) {
    # Get top correlated genes for this network size
    all_genes <- matched
    for (qg in matched) {
      cors <- cor_matrix[qg, ]
      cors <- cors[!names(cors) %in% matched]
      top_cors <- names(sort(abs(cors), decreasing = TRUE)[1:min(net_size, length(cors))])
      all_genes <- unique(c(all_genes, top_cors))
    }
    
    # Subset correlation matrix
    sub_cor <- cor_matrix[all_genes, all_genes]
    
    for (cor_thresh in PARAM_COMBINATIONS$cor_thresholds) {
      # Format threshold consistently: always 2 decimal places
      thresh_str <- sprintf("%.2f", cor_thresh)
      key <- paste0("size_", net_size, "_cor_", gsub("\\.", "", thresh_str))
      
      adj <- abs(sub_cor)
      adj[adj < cor_thresh] <- 0
      diag(adj) <- 0
      
      # Store node and edge data
      g <- graph_from_adjacency_matrix(adj, mode = "undirected", weighted = TRUE)
      if (vcount(g) == 0) next
      
      V(g)$is_query <- V(g)$name %in% matched
      
      nodes_df <- data.frame(
        id = V(g)$name,
        label = V(g)$name,
        group = ifelse(V(g)$is_query, "Query", "Correlated"),
        isQuery = V(g)$is_query,
        stringsAsFactors = FALSE
      )
      
      edge_list <- as_data_frame(g, what = "edges")
      if (nrow(edge_list) > 0) {
        edges_df <- data.frame(
          from = edge_list$from,
          to = edge_list$to,
          weight = round(edge_list$weight, 4),
          stringsAsFactors = FALSE
        )
      } else {
        edges_df <- data.frame(from = character(), to = character(), 
                               weight = numeric(), stringsAsFactors = FALSE)
      }
      
      precalc_networks[[key]] <- list(
        nodes = nodes_df,
        edges = edges_df,
        net_size = net_size,
        cor_threshold = cor_thresh,
        n_nodes = nrow(nodes_df),
        n_edges = nrow(edges_df),
        n_query = sum(nodes_df$isQuery)
      )
    }
  }
  
  # Precalculate query-only networks
  query_networks <- list()
  query_sub_cor <- cor_matrix[matched, matched]
  
  for (q_thresh in PARAM_COMBINATIONS$query_network_thresholds) {
    thresh_str <- sprintf("%.2f", q_thresh)
    key <- paste0("query_cor_", gsub("\\.", "", thresh_str))
    
    query_adj <- abs(query_sub_cor)
    query_adj[query_adj < q_thresh] <- 0
    diag(query_adj) <- 0
    
    g_q <- graph_from_adjacency_matrix(query_adj, mode = "undirected", weighted = TRUE)
    if (vcount(g_q) == 0) next
    
    nodes_q <- data.frame(
      id = V(g_q)$name,
      label = V(g_q)$name,
      stringsAsFactors = FALSE
    )
    
    edge_list_q <- as_data_frame(g_q, what = "edges")
    if (nrow(edge_list_q) > 0) {
      edges_q <- data.frame(
        from = edge_list_q$from,
        to = edge_list_q$to,
        weight = round(edge_list_q$weight, 4),
        stringsAsFactors = FALSE
      )
    } else {
      edges_q <- data.frame(from = character(), to = character(),
                            weight = numeric(), stringsAsFactors = FALSE)
    }
    
    query_networks[[key]] <- list(
      nodes = nodes_q,
      edges = edges_q,
      cor_threshold = q_thresh,
      n_nodes = nrow(nodes_q),
      n_edges = nrow(edges_q)
    )
  }
  
  cat("    Precalculated", length(precalc_networks), "full networks and", 
      length(query_networks), "query-only networks\n")
  cat("    Full network keys:", paste(names(precalc_networks), collapse = ", "), "\n")
  cat("    Query network keys:", paste(names(query_networks), collapse = ", "), "\n")
  
  # ===== PNG OUTPUT: Use default parameters =====
  if (OUTPUT_PNG) {
    # Use default settings for PNG
    all_genes <- matched
    for (qg in matched) {
      cors <- cor_matrix[qg, ]
      cors <- cors[!names(cors) %in% matched]
      top_cors <- names(sort(abs(cors), decreasing = TRUE)[1:min(N_NETWORK_GENES, length(cors))])
      all_genes <- unique(c(all_genes, top_cors))
    }
    
    sub_cor <- cor_matrix[all_genes, all_genes]
    adj <- abs(sub_cor)
    adj[adj < COR_THRESHOLD] <- 0
    diag(adj) <- 0
    
    g <- graph_from_adjacency_matrix(adj, mode = "undirected", weighted = TRUE)
    V(g)$is_query <- V(g)$name %in% matched
    
    colors <- ifelse(V(g)$is_query, "#B2182B", "#4393C3")
    sizes <- ifelse(V(g)$is_query, 18, 6)
    label_sizes <- ifelse(V(g)$is_query, 1.2, 0.5)
    label_colors <- ifelse(V(g)$is_query, "black", "gray40")
    frame_widths <- ifelse(V(g)$is_query, 3, 1)
    
    png(file.path(output_dir, paste0(gene_group, "_correlation_network.png")),
        width = 1600, height = 1400, res = 120)
    
    set.seed(42)
    layout <- layout_with_fr(g)
    
    plot(g, 
         layout = layout,
         vertex.color = colors,
         vertex.size = sizes,
         vertex.label.cex = label_sizes,
         vertex.label.color = label_colors,
         vertex.label.font = ifelse(V(g)$is_query, 2, 1),
         vertex.frame.width = frame_widths,
         vertex.frame.color = ifelse(V(g)$is_query, "black", NA),
         edge.width = E(g)$weight * 3,
         edge.color = adjustcolor("gray50", alpha.f = 0.5),
         main = paste0(gene_group, " - Correlation Network\n(Query genes in RED/BOLD, r>=", COR_THRESHOLD, ")"))
    
    legend("bottomleft", 
           legend = c(paste0("Query genes (n=", sum(V(g)$is_query), ")"), 
                      paste0("Correlated genes (n=", sum(!V(g)$is_query), ")")),
           col = c("#B2182B", "#4393C3"), 
           pch = 19, pt.cex = c(2, 1.2), bty = "n", cex = 0.9)
    
    dev.off()
    
    # Query-only network PNG
    if (GENERATE_WGCNA_FIGURES$query_gene_focus && length(matched) >= 2) {
      query_adj <- abs(query_sub_cor)
      query_adj[query_adj < 0.7] <- 0
      diag(query_adj) <- 0
      
      g_query <- graph_from_adjacency_matrix(query_adj, mode = "undirected", weighted = TRUE)
      
      png(file.path(output_dir, paste0(gene_group, "_query_genes_network.png")),
          width = 1200, height = 1000, res = 120)
      
      set.seed(42)
      layout_q <- layout_with_fr(g_query)
      
      plot(g_query,
           layout = layout_q,
           vertex.color = "#B2182B",
           vertex.size = 25,
           vertex.label.cex = 1.0,
           vertex.label.color = "black",
           vertex.label.font = 2,
           vertex.frame.width = 3,
           vertex.frame.color = "black",
           edge.width = E(g_query)$weight * 5,
           edge.color = adjustcolor("darkred", alpha.f = 0.6),
           edge.label = round(E(g_query)$weight, 2),
           edge.label.cex = 0.7,
           main = paste0(gene_group, " - Query Genes Co-expression\n(Edge weights = correlation, r>=0.70)"))
      
      dev.off()
    }
  }
  
  # ===== HTML OUTPUT: INTERACTIVE WITH PARAMETER TOGGLES =====
  if (OUTPUT_HTML && length(precalc_networks) > 0) {
    cat("  Building interactive HTML with parameter toggles...\n")
    
    # Convert to JSON for embedding
    networks_json <- toJSON(precalc_networks, auto_unbox = TRUE)
    query_networks_json <- toJSON(query_networks, auto_unbox = TRUE)
    
    html_content <- paste0('<!DOCTYPE html>
<html>
<head>
  <title>', gene_group, ' - Interactive Correlation Network</title>
  <script src="https://unpkg.com/vis-network/standalone/umd/vis-network.min.js"></script>
  <style>
    * { box-sizing: border-box; margin: 0; padding: 0; }
    body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif; background: #f5f5f5; }
    .header { background: linear-gradient(135deg, #1a365d 0%, #2c5282 100%); color: white; padding: 20px; text-align: center; }
    .header h1 { margin-bottom: 5px; font-size: 1.8em; }
    .header .subtitle { opacity: 0.9; font-size: 0.95em; }
    .container { display: flex; height: calc(100vh - 100px); }
    .sidebar { width: 320px; background: white; padding: 20px; overflow-y: auto; box-shadow: 2px 0 10px rgba(0,0,0,0.1); }
    .main { flex: 1; display: flex; flex-direction: column; }
    #network { flex: 1; background: white; margin: 15px; border-radius: 8px; box-shadow: 0 2px 10px rgba(0,0,0,0.1); }
    .control-group { margin-bottom: 20px; padding-bottom: 15px; border-bottom: 1px solid #e2e8f0; }
    .control-group:last-child { border-bottom: none; }
    .control-group h3 { color: #2d3748; margin-bottom: 10px; font-size: 0.9em; text-transform: uppercase; letter-spacing: 0.5px; }
    .control-group label { display: block; margin: 8px 0; cursor: pointer; padding: 8px 12px; border-radius: 6px; transition: background 0.2s; }
    .control-group label:hover { background: #f7fafc; }
    .control-group input[type="radio"] { margin-right: 10px; }
    .control-group input[type="checkbox"] { margin-right: 10px; }
    select { width: 100%; padding: 10px; border: 1px solid #e2e8f0; border-radius: 6px; font-size: 14px; background: white; cursor: pointer; }
    .stats { background: #edf2f7; padding: 15px; border-radius: 8px; margin-top: 15px; }
    .stats h4 { margin-bottom: 10px; color: #2d3748; font-size: 0.85em; text-transform: uppercase; }
    .stat-row { display: flex; justify-content: space-between; margin: 5px 0; font-size: 0.9em; }
    .stat-label { color: #718096; }
    .stat-value { font-weight: 600; color: #2d3748; }
    .legend { margin-top: 20px; }
    .legend-item { display: flex; align-items: center; margin: 8px 0; font-size: 0.9em; }
    .legend-dot { width: 16px; height: 16px; border-radius: 50%; margin-right: 10px; }
    .query-dot { background: #B2182B; border: 2px solid #7f1d1d; }
    .corr-dot { background: #4393C3; border: 2px solid #2563eb; }
    .tab-buttons { display: flex; margin-bottom: 15px; }
    .tab-btn { flex: 1; padding: 10px; border: none; background: #e2e8f0; cursor: pointer; font-size: 0.85em; transition: all 0.2s; }
    .tab-btn:first-child { border-radius: 6px 0 0 6px; }
    .tab-btn:last-child { border-radius: 0 6px 6px 0; }
    .tab-btn.active { background: #3182ce; color: white; }
    .help-text { font-size: 0.8em; color: #718096; margin-top: 5px; }
  </style>
</head>
<body>
  <div class="header">
    <h1>', gene_group, '</h1>
    <div class="subtitle">Interactive Co-expression Network | Query genes highlighted in RED</div>
  </div>
  
  <div class="container">
    <div class="sidebar">
      <div class="tab-buttons">
        <button class="tab-btn active" onclick="switchTab(\'full\')">Full Network</button>
        <button class="tab-btn" onclick="switchTab(\'query\')">Query Only</button>
      </div>
      
      <div id="full-controls">
        <div class="control-group">
          <h3>Correlation Threshold</h3>
          ', paste(sapply(PARAM_COMBINATIONS$cor_thresholds, function(ct) {
            checked <- if(ct == COR_THRESHOLD) "checked" else ""
            paste0('<label><input type="radio" name="corThresh" value="', ct, '" ', checked, ' onchange="updateNetwork()"> r >= ', ct, '</label>')
          }), collapse = "\n          "), '
          <p class="help-text">Higher = fewer but stronger connections</p>
        </div>
        
        <div class="control-group">
          <h3>Network Size</h3>
          ', paste(sapply(PARAM_COMBINATIONS$network_sizes, function(ns) {
            checked <- if(ns == N_NETWORK_GENES) "checked" else ""
            paste0('<label><input type="radio" name="netSize" value="', ns, '" ', checked, ' onchange="updateNetwork()"> Top ', ns, ' genes/query</label>')
          }), collapse = "\n          "), '
          <p class="help-text">Genes most correlated with query genes</p>
        </div>
      </div>
      
      <div id="query-controls" style="display:none;">
        <div class="control-group">
          <h3>Query Network Threshold</h3>
          ', paste(sapply(PARAM_COMBINATIONS$query_network_thresholds, function(qt) {
            checked <- if(qt == 0.70) "checked" else ""
            paste0('<label><input type="radio" name="queryThresh" value="', qt, '" ', checked, ' onchange="updateQueryNetwork()"> r >= ', qt, '</label>')
          }), collapse = "\n          "), '
        </div>
      </div>
      
      <div class="control-group">
        <h3>Display Options</h3>
        <label><input type="checkbox" id="showLabels" checked onchange="toggleLabels()"> Show gene labels</label>
        <label><input type="checkbox" id="physics" checked onchange="togglePhysics()"> Enable physics</label>
      </div>
      
      <div class="stats" id="stats">
        <h4>Network Statistics</h4>
        <div class="stat-row"><span class="stat-label">Total nodes:</span><span class="stat-value" id="statNodes">-</span></div>
        <div class="stat-row"><span class="stat-label">Query genes:</span><span class="stat-value" id="statQuery">-</span></div>
        <div class="stat-row"><span class="stat-label">Edges:</span><span class="stat-value" id="statEdges">-</span></div>
        <div class="stat-row"><span class="stat-label">Cor threshold:</span><span class="stat-value" id="statThresh">-</span></div>
      </div>
      
      <div class="legend">
        <h4 style="margin-bottom:10px; font-size:0.85em; color:#2d3748; text-transform:uppercase;">Legend</h4>
        <div class="legend-item"><div class="legend-dot query-dot"></div>Query Gene (larger, bold)</div>
        <div class="legend-item"><div class="legend-dot corr-dot"></div>Correlated Gene</div>
      </div>
    </div>
    
    <div class="main">
      <div id="network"><div style="display:flex;align-items:center;justify-content:center;height:100%;color:#666;">Loading network...</div></div>
    </div>
  </div>
  
  <script>
    // Precalculated network data
    let networks = {};
    let queryNetworks = {};
    
    try {
      networks = ', networks_json, ';
      queryNetworks = ', query_networks_json, ';
      console.log("Networks loaded:", Object.keys(networks).length, "full networks,", Object.keys(queryNetworks).length, "query networks");
    } catch (e) {
      console.error("Error parsing network data:", e);
    }
    
    let network = null;
    let currentTab = "full";
    
    function getSelectedValue(name) {
      const radio = document.querySelector(\'input[name="\' + name + \'"]:checked\');
      return radio ? radio.value : null;
    }
    
    function buildKey(netSize, corThresh) {
      // Format to 2 decimal places to match R key format, then remove dot
      const formatted = parseFloat(corThresh).toFixed(2);
      return "size_" + netSize + "_cor_" + formatted.replace(".", "");
    }
    
    function buildQueryKey(corThresh) {
      const formatted = parseFloat(corThresh).toFixed(2);
      return "query_cor_" + formatted.replace(".", "");
    }
    
    function createNetwork(nodes, edges, stats) {
      const container = document.getElementById("network");
      
      if (!nodes || nodes.length === 0) {
        container.innerHTML = "<div style=\\"display:flex;align-items:center;justify-content:center;height:100%;color:#666;\\">No nodes to display for this configuration</div>";
        return;
      }
      
      const visNodes = nodes.map(n => ({
        id: n.id,
        label: n.label,
        color: n.isQuery ? "#B2182B" : "#4393C3",
        size: n.isQuery ? 30 : 15,
        font: { size: n.isQuery ? 16 : 11, bold: n.isQuery, color: n.isQuery ? "#1a1a1a" : "#4a4a4a" },
        borderWidth: n.isQuery ? 3 : 1,
        borderWidthSelected: 4,
        title: "<b>" + n.label + "</b><br>" + (n.isQuery ? "Query Gene" : "Correlated Gene")
      }));
      
      const visEdges = edges.map(e => ({
        from: e.from,
        to: e.to,
        width: Math.max(1, e.weight * 5),
        color: { color: "rgba(100,100,100,0.4)", highlight: "#666" },
        title: "Correlation: " + e.weight.toFixed(3)
      }));
      
      const data = { nodes: new vis.DataSet(visNodes), edges: new vis.DataSet(visEdges) };
      
      const options = {
        layout: {
          improvedLayout: false,  // Disable for large networks
          randomSeed: 42
        },
        physics: {
          enabled: document.getElementById("physics").checked,
          solver: "forceAtlas2Based",
          forceAtlas2Based: { 
            gravitationalConstant: -50, 
            centralGravity: 0.01,
            springLength: 100,
            springConstant: 0.08,
            damping: 0.4
          },
          stabilization: { 
            enabled: true,
            iterations: 200,
            updateInterval: 25
          }
        },
        interaction: {
          hover: true,
          tooltipDelay: 100,
          navigationButtons: true,
          keyboard: true,
          zoomView: true,
          dragView: true
        },
        nodes: { shape: "dot" },
        edges: { smooth: { type: "continuous" } }
      };
      
      if (network) network.destroy();
      network = new vis.Network(container, data, options);
      
      // Fit network to view after stabilization
      network.once("stabilizationIterationsDone", function() {
        network.fit({ animation: { duration: 500 } });
      });
      
      // Update stats
      document.getElementById("statNodes").textContent = stats.n_nodes || nodes.length;
      document.getElementById("statQuery").textContent = stats.n_query || nodes.filter(n => n.isQuery).length;
      document.getElementById("statEdges").textContent = stats.n_edges || edges.length;
      document.getElementById("statThresh").textContent = "r >= " + stats.cor_threshold;
      
      // Apply label visibility
      toggleLabels();
    }
    
    function updateNetwork() {
      const netSize = getSelectedValue("netSize");
      const corThresh = getSelectedValue("corThresh");
      const key = buildKey(netSize, corThresh);
      
      console.log("Looking for network:", key);
      console.log("Available networks:", Object.keys(networks));
      
      if (networks[key]) {
        const net = networks[key];
        createNetwork(net.nodes, net.edges, net);
      } else {
        // Try to find any available network as fallback
        const availableKeys = Object.keys(networks);
        console.warn("Network not found:", key, "Available:", availableKeys);
        if (availableKeys.length > 0) {
          const fallbackKey = availableKeys[0];
          console.log("Using fallback network:", fallbackKey);
          const net = networks[fallbackKey];
          createNetwork(net.nodes, net.edges, net);
        } else {
          document.getElementById("network").innerHTML = "<div style=\\"display:flex;align-items:center;justify-content:center;height:100%;color:#666;\\">No network data available</div>";
        }
      }
    }
    
    function updateQueryNetwork() {
      const corThresh = getSelectedValue("queryThresh");
      const key = buildQueryKey(corThresh);
      
      console.log("Looking for query network:", key);
      console.log("Available query networks:", Object.keys(queryNetworks));
      
      if (queryNetworks[key]) {
        const net = queryNetworks[key];
        // Add isQuery=true for all nodes in query network
        const nodes = net.nodes.map(n => ({...n, isQuery: true}));
        createNetwork(nodes, net.edges, {...net, n_query: net.n_nodes});
      } else {
        const availableKeys = Object.keys(queryNetworks);
        if (availableKeys.length > 0) {
          const fallbackKey = availableKeys[0];
          const net = queryNetworks[fallbackKey];
          const nodes = net.nodes.map(n => ({...n, isQuery: true}));
          createNetwork(nodes, net.edges, {...net, n_query: net.n_nodes});
        }
      }
    }
    
    function switchTab(tab) {
      currentTab = tab;
      document.querySelectorAll(".tab-btn").forEach((btn, i) => {
        btn.classList.toggle("active", (i === 0 && tab === "full") || (i === 1 && tab === "query"));
      });
      document.getElementById("full-controls").style.display = tab === "full" ? "block" : "none";
      document.getElementById("query-controls").style.display = tab === "query" ? "block" : "none";
      
      if (tab === "full") updateNetwork();
      else updateQueryNetwork();
    }
    
    function toggleLabels() {
      if (!network) return;
      const show = document.getElementById("showLabels").checked;
      network.setOptions({ nodes: { font: { size: show ? 14 : 0 } } });
    }
    
    function togglePhysics() {
      if (!network) return;
      network.setOptions({ physics: { enabled: document.getElementById("physics").checked } });
    }
    
    // Initialize on page load
    document.addEventListener("DOMContentLoaded", function() {
      console.log("Page loaded, initializing network...");
      console.log("Networks object:", networks);
      console.log("Network keys:", Object.keys(networks));
      
      // Check if vis library loaded
      if (typeof vis === "undefined") {
        document.getElementById("network").innerHTML = "<div style=\\"display:flex;align-items:center;justify-content:center;height:100%;color:#c00;padding:20px;text-align:center;\\"><div><b>vis-network library failed to load</b><br><br>Please ensure you have internet access.<br>The network visualization requires the vis-network library from unpkg.com</div></div>";
        return;
      }
      
      if (Object.keys(networks).length === 0) {
        document.getElementById("network").innerHTML = "<div style=\\"display:flex;align-items:center;justify-content:center;height:100%;color:#c00;padding:20px;text-align:center;\\"><div><b>No network data found</b><br><br>This may happen if:<br>- Correlation thresholds are too stringent<br>- No edges pass the threshold<br><br>Try running with lower COR_THRESHOLD</div></div>";
        return;
      }
      
      try {
        updateNetwork();
      } catch (e) {
        console.error("Error updating network:", e);
        document.getElementById("network").innerHTML = "<div style=\\"display:flex;align-items:center;justify-content:center;height:100%;color:#c00;padding:20px;\\">Error: " + e.message + "</div>";
      }
    });
  </script>
</body>
</html>')
    
    writeLines(html_content, file.path(output_dir, paste0(gene_group, "_network_interactive.html")))
    cat("  Saved interactive HTML with toggles:", paste0(gene_group, "_network_interactive.html"), "\n")
  }
  
  # ===== SIMPLE HTML OUTPUT: Fast loading, minimal UI =====
  if (OUTPUT_HTML_SIMPLE && length(precalc_networks) > 0) {
    cat("  Building simple HTML network view...\n")
    
    # Use default network (middle settings for good balance)
    default_key <- paste0("size_", N_NETWORK_GENES, "_cor_", gsub("\\.", "", sprintf("%.2f", COR_THRESHOLD)))
    if (!default_key %in% names(precalc_networks)) {
      default_key <- names(precalc_networks)[1]
    }
    
    simple_net <- precalc_networks[[default_key]]
    simple_nodes_json <- toJSON(simple_net$nodes, auto_unbox = TRUE)
    simple_edges_json <- toJSON(simple_net$edges, auto_unbox = TRUE)
    
    simple_html <- paste0('<!DOCTYPE html>
<html>
<head>
  <title>', gene_group, ' - Network (Simple View)</title>
  <script src="https://unpkg.com/vis-network/standalone/umd/vis-network.min.js"></script>
  <style>
    * { margin: 0; padding: 0; box-sizing: border-box; }
    body { font-family: system-ui, -apple-system, sans-serif; background: #f8f9fa; }
    .header { background: #2c5282; color: white; padding: 15px 20px; display: flex; justify-content: space-between; align-items: center; }
    .header h1 { font-size: 1.3em; font-weight: 500; }
    .header .info { font-size: 0.85em; opacity: 0.9; }
    #network { width: 100%; height: calc(100vh - 60px); background: white; }
    .legend { position: fixed; bottom: 20px; left: 20px; background: white; padding: 12px 16px; border-radius: 8px; box-shadow: 0 2px 8px rgba(0,0,0,0.15); font-size: 0.85em; }
    .legend-item { display: flex; align-items: center; margin: 4px 0; }
    .legend-dot { width: 14px; height: 14px; border-radius: 50%; margin-right: 8px; }
    .query { background: #B2182B; }
    .corr { background: #4393C3; }
  </style>
</head>
<body>
  <div class="header">
    <h1>', gene_group, ' - Co-expression Network</h1>
    <div class="info">', simple_net$n_nodes, ' genes | ', simple_net$n_edges, ' edges | r \u2265 ', simple_net$cor_threshold, '</div>
  </div>
  <div id="network"></div>
  <div class="legend">
    <div class="legend-item"><div class="legend-dot query"></div>Query genes (', simple_net$n_query, ')</div>
    <div class="legend-item"><div class="legend-dot corr"></div>Correlated genes (', simple_net$n_nodes - simple_net$n_query, ')</div>
  </div>
  <script>
    const nodes = ', simple_nodes_json, ';
    const edges = ', simple_edges_json, ';
    
    const visNodes = nodes.map(n => ({
      id: n.id, label: n.label,
      color: n.isQuery ? "#B2182B" : "#4393C3",
      size: n.isQuery ? 25 : 12,
      font: { size: n.isQuery ? 14 : 9, bold: n.isQuery, color: "#333" },
      borderWidth: n.isQuery ? 2 : 1,
      title: n.label + (n.isQuery ? " (Query)" : "")
    }));
    
    const visEdges = edges.map(e => ({
      from: e.from, to: e.to,
      width: Math.max(0.5, e.weight * 3),
      color: { color: "rgba(120,120,120,0.3)" },
      smooth: false
    }));
    
    const container = document.getElementById("network");
    const data = { nodes: new vis.DataSet(visNodes), edges: new vis.DataSet(visEdges) };
    const options = {
      layout: { improvedLayout: false, randomSeed: 42 },
      physics: {
        enabled: true,
        solver: "barnesHut",
        barnesHut: { gravitationalConstant: -3000, centralGravity: 0.5, springLength: 80, springConstant: 0.05 },
        stabilization: { iterations: 40, fit: true },
        maxVelocity: 50
      },
      interaction: { hover: true, tooltipDelay: 50, zoomView: true, dragView: true },
      nodes: { shape: "dot" },
      edges: { smooth: false }
    };
    
    const network = new vis.Network(container, data, options);
    network.once("stabilizationIterationsDone", function() {
      network.setOptions({ physics: { enabled: false } });
      network.fit();
      network.moveTo({ scale: 0.6 });
    });
  </script>
</body>
</html>')
    
    writeLines(simple_html, file.path(output_dir, paste0(gene_group, "_network_simple.html")))
    cat("  Saved simple HTML:", paste0(gene_group, "_network_simple.html"), "\n")
  }
  
  return(list(precalc_networks = precalc_networks, query_networks = query_networks))
}

# ===============================================
# MAIN WGCNA FUNCTION
# ===============================================

run_wgcna <- function(config = NULL, matrices_dir = NULL) {
  # Get method base directory from environment for config file loading
  method_base_dir <- Sys.getenv("METHOD_BASE_DIR", unset = ".")
  if (is.null(config)) config <- load_runtime_config(method_base_dir)
  
  # Set up matrices_dir based on method_base_dir if not provided
  if (is.null(matrices_dir)) {
    matrices_dir <- file.path(method_base_dir, get_matrices_dir(CURRENT_METHOD))
  }
  cat("  Matrices directory:", matrices_dir, "\n")
  
  ensure_output_dir(WGCNA_OUT_DIR)
  
  print_config_summary("WGCNA COEXPRESSION ANALYSIS", config)
  cat("Mode: Full transcriptome analysis with query genes highlighted\n")
  cat("Top variable genes used:", TOP_VAR_GENES, "\n")
  
  # Log GPU status for correlation network computation
  if (GPU_AVAILABLE) {
    cat("GPU acceleration: ENABLED for correlation matrices (", GPU_VRAM_GB, "GB VRAM)\n", sep = "")
  }
  
  successful <- 0
  total <- 0
  
  for (gene_group in config$gene_groups) {
    cat("\n", paste(rep("=", 60), collapse = ""), "\n")
    cat("Processing:", gene_group, "\n")
    cat(paste(rep("=", 60), collapse = ""), "\n")
    total <- total + 1
    
    output_dir <- file.path(WGCNA_OUT_DIR, gene_group)
    ensure_output_dir(output_dir)
    
    # ===== STEP 1: Read the QUERY GENES (from gene group matrix) =====
    query_file <- build_input_path(gene_group, PROCESSING_LEVELS[1], 
                                   COUNT_TYPES[1], GENE_TYPES[1],
                                   matrices_dir, config$master_reference)
    
    cat("  Looking for query file:", query_file, "\n")
    
    if (!file.exists(query_file)) {
      cat("  Skipped: Query gene file not found\n")
      cat("  Trying alternative path...\n")
      # Alternative: try with Gene_ID type
      query_file <- build_input_path(gene_group, PROCESSING_LEVELS[1], 
                                     COUNT_TYPES[1], "Gene_ID",
                                     matrices_dir, config$master_reference)
      cat("  Alternative path:", query_file, "\n")
      if (!file.exists(query_file)) {
        cat("  Skipped: Neither path found\n")
        next
      }
    }
    
    query_matrix <- read_count_matrix(query_file)
    query_genes <- rownames(query_matrix)
    cat("  Query genes (to highlight):", length(query_genes), "\n")
    cat("    Genes:", paste(head(query_genes, 5), collapse = ", "), 
        if (length(query_genes) > 5) "..." else "", "\n")
    
    # ===== STEP 2: Read the FULL TRANSCRIPTOME matrix =====
    full_file <- build_input_path(config$master_reference, PROCESSING_LEVELS[1], 
                                  COUNT_TYPES[1], GENE_TYPES[1],
                                  matrices_dir, config$master_reference)
    
    validation <- validate_and_read_matrix(full_file, MIN_GENES_WGCNA)
    if (!validation$success) {
      cat("  Skipped: Full transcriptome -", validation$reason, "\n")
      next
    }
    
    cat("  Full transcriptome genes:", nrow(validation$data), "\n")
    
    # ===== STEP 3: Ensure query genes are in the full transcriptome =====
    query_genes_matched <- query_genes[query_genes %in% rownames(validation$data)]
    if (length(query_genes_matched) == 0) {
      cat("  Warning: No query genes found in transcriptome, trying Gene_ID format...\n")
      # Try with Gene_ID format instead
      query_file_geneid <- build_input_path(gene_group, PROCESSING_LEVELS[1], 
                                            COUNT_TYPES[1], "Gene_ID",
                                            matrices_dir, config$master_reference)
      if (file.exists(query_file_geneid)) {
        query_matrix_geneid <- read_count_matrix(query_file_geneid)
        query_genes <- rownames(query_matrix_geneid)
        query_genes_matched <- query_genes[query_genes %in% rownames(validation$data)]
      }
    }
    
    cat("  Query genes matched in transcriptome:", length(query_genes_matched), "\n")
    
    # ===== STEP 4: Variance stabilization and filtering =====
    data_log <- log2(validation$data + 1)
    
    # Filter genes with low variance
    gene_vars <- apply(data_log, 1, var, na.rm = TRUE)
    var_threshold <- quantile(gene_vars, 0.25, na.rm = TRUE)
    
    # IMPORTANT: Always keep query genes even if low variance
    keep_genes <- gene_vars > var_threshold | rownames(data_log) %in% query_genes_matched
    data_filtered <- data_log[keep_genes, , drop = FALSE]
    
    cat("  Genes after variance filtering:", nrow(data_filtered), "\n")
    cat("  Query genes retained:", sum(query_genes_matched %in% rownames(data_filtered)), "\n")
    
    # ===== STEP 5: Subsample for speed if too many genes =====
    if (nrow(data_filtered) > TOP_VAR_GENES) {
      # Keep top variable genes BUT always include query genes
      gene_vars_filtered <- apply(data_filtered, 1, var, na.rm = TRUE)
      top_var_genes <- names(sort(gene_vars_filtered, decreasing = TRUE)[1:TOP_VAR_GENES])
      
      # Ensure query genes are included
      genes_to_keep <- unique(c(top_var_genes, query_genes_matched))
      genes_to_keep <- genes_to_keep[genes_to_keep %in% rownames(data_filtered)]
      
      data_filtered <- data_filtered[genes_to_keep, , drop = FALSE]
      cat("  Subsampled to top variable genes:", nrow(data_filtered), "\n")
    }
    
    data_matrix <- t(data_filtered)  # WGCNA expects samples as rows
    
    # ===== STEP 6: Run WGCNA pipeline =====
    cat("  Running soft threshold selection...\n")
    soft_power <- pick_soft_threshold(data_matrix, output_dir, gene_group)
    cat("  Soft power:", soft_power, "\n")
    
    cat("  Building network and detecting modules...\n")
    network <- build_network_and_detect_modules(data_matrix, soft_power, output_dir, 
                                                 gene_group, query_genes_matched)
    cat("  Modules detected:", length(unique(network$module_colors)), "\n")
    
    MEs <- calculate_module_eigengenes(data_matrix, network$module_colors, output_dir, gene_group)
    
    # ===== STEP 7: Export gene-module assignments (with query gene flag) =====
    kME <- signedKME(data_matrix, MEs, outputColumnName = "kME")
    gene_info <- data.frame(
      Gene = colnames(data_matrix),
      Module = network$module_colors,
      Is_Query_Gene = colnames(data_matrix) %in% query_genes_matched,
      stringsAsFactors = FALSE
    )
    gene_info <- cbind(gene_info, kME)
    
    # Sort to show query genes first
    gene_info <- gene_info[order(-gene_info$Is_Query_Gene, gene_info$Module), ]
    
    write.table(gene_info, file.path(output_dir, paste0(gene_group, "_module_assignments.tsv")),
                sep = "\t", row.names = FALSE, quote = FALSE)
    
    # ===== STEP 8: Query gene module summary =====
    query_gene_info <- gene_info[gene_info$Is_Query_Gene, ]
    cat("  Query genes module distribution:\n")
    print(table(query_gene_info$Module))
    
    write.table(query_gene_info, 
                file.path(output_dir, paste0(gene_group, "_query_genes_modules.tsv")),
                sep = "\t", row.names = FALSE, quote = FALSE)
    
    # ===== STEP 9: Hub genes =====
    hubs <- identify_hub_genes(gene_info)
    # Mark query genes in hub list
    hubs$Is_Query_Gene <- hubs$Gene %in% query_genes_matched
    write.table(hubs, file.path(output_dir, paste0(gene_group, "_hub_genes.tsv")),
                sep = "\t", row.names = FALSE, quote = FALSE)
    
    # ===== STEP 10: Correlation network with query genes bold =====
    cat("  Creating correlation network...\n")
    create_correlation_network(t(data_matrix), query_genes_matched, output_dir, gene_group)
    
    # ===== STEP 11: Find genes co-expressed with query genes =====
    cat("  Finding genes co-expressed with query genes...\n")
    cor_matrix <- gpu_cor(data_filtered)
    
    coexpr_results <- data.frame()
    for (qg in query_genes_matched) {
      if (qg %in% rownames(cor_matrix)) {
        cors <- cor_matrix[qg, ]
        cors <- cors[names(cors) != qg]  # Remove self-correlation
        top_coexpr <- sort(cors, decreasing = TRUE)[1:min(N_COEXPRESSED_GENES, length(cors))]
        
        qg_coexpr <- data.frame(
          Query_Gene = qg,
          Coexpressed_Gene = names(top_coexpr),
          Correlation = as.numeric(top_coexpr),
          Is_Also_Query = names(top_coexpr) %in% query_genes_matched,
          stringsAsFactors = FALSE
        )
        coexpr_results <- rbind(coexpr_results, qg_coexpr)
      }
    }
    
    if (nrow(coexpr_results) > 0) {
      write.table(coexpr_results, 
                  file.path(output_dir, paste0(gene_group, "_coexpressed_genes.tsv")),
                  sep = "\t", row.names = FALSE, quote = FALSE)
    }
    
    # ===== STEP 12: Export RAW RESULTS for downstream analysis =====
    if (OUTPUT_RAW) {
    cat("  Exporting raw results...\n")
    raw_results_dir <- file.path(output_dir, "raw_results")
    if (!dir.exists(raw_results_dir)) dir.create(raw_results_dir, recursive = TRUE)
    
    # 12a: Save full correlation matrix as TSV
    cat("    Saving correlation matrix...\n")
    cor_df <- as.data.frame(cor_matrix)
    cor_df$Gene <- rownames(cor_matrix)
    cor_df <- cor_df[, c("Gene", setdiff(names(cor_df), "Gene"))]  # Move Gene to first column
    write.table(cor_df, file.path(raw_results_dir, paste0(gene_group, "_correlation_matrix.tsv")),
                sep = "\t", row.names = FALSE, quote = FALSE)
    
    # 12b: Save module eigengenes
    cat("    Saving module eigengenes...\n")
    me_df <- as.data.frame(MEs)
    me_df$Sample <- rownames(MEs)
    me_df <- me_df[, c("Sample", setdiff(names(me_df), "Sample"))]
    write.table(me_df, file.path(raw_results_dir, paste0(gene_group, "_module_eigengenes.tsv")),
                sep = "\t", row.names = FALSE, quote = FALSE)
    
    # 12c: Save network as RDS for complete reproducibility
    cat("    Saving network object (RDS)...\n")
    saveRDS(network, file.path(raw_results_dir, paste0(gene_group, "_network.rds")))
    
    # 12d: Save gene info with kME as RDS
    saveRDS(gene_info, file.path(raw_results_dir, paste0(gene_group, "_gene_info.rds")))
    
    # 12e: Save soft threshold analysis
    sft_file <- file.path(output_dir, paste0(gene_group, "_soft_threshold.tsv"))
    if (file.exists(sft_file)) {
      file.copy(sft_file, file.path(raw_results_dir, paste0(gene_group, "_soft_threshold.tsv")))
    }
    
    # 12f: Create a summary file with all parameters used
    params_summary <- data.frame(
      Parameter = c("TOP_VAR_GENES", "MIN_MODULE_SIZE_DEFAULT", "MIN_MODULE_SIZE_SMALL",
                    "COR_THRESHOLD", "N_HUB_GENES", "N_COEXPRESSED_GENES", "N_NETWORK_GENES",
                    "MERGE_CUT_HEIGHT", "NETWORK_TYPE", "Soft_Power_Used",
                    "Total_Genes_Analyzed", "Query_Genes_Matched"),
      Value = c(TOP_VAR_GENES, MIN_MODULE_SIZE_DEFAULT, MIN_MODULE_SIZE_SMALL,
                COR_THRESHOLD, N_HUB_GENES, N_COEXPRESSED_GENES, N_NETWORK_GENES,
                MERGE_CUT_HEIGHT, NETWORK_TYPE, soft_power,
                ncol(data_matrix), length(query_genes_matched)),
      stringsAsFactors = FALSE
    )
    write.table(params_summary, file.path(raw_results_dir, paste0(gene_group, "_parameters.tsv")),
                sep = "\t", row.names = FALSE, quote = FALSE)
    
    # 12g: Save data matrix used (expression values)
    cat("    Saving expression matrix...\n")
    expr_df <- as.data.frame(t(data_matrix))  # Genes as rows
    expr_df$Gene <- rownames(expr_df)
    expr_df <- expr_df[, c("Gene", setdiff(names(expr_df), "Gene"))]
    write.table(expr_df, file.path(raw_results_dir, paste0(gene_group, "_expression_matrix.tsv")),
                sep = "\t", row.names = FALSE, quote = FALSE)
    
    cat("  Raw results saved to:", raw_results_dir, "\n")
    }  # End OUTPUT_RAW block
    
    successful <- successful + 1
    cat("  Complete!\n")
  }
  
  print_summary(successful, total)
}

if (!interactive() && identical(environment(), globalenv())) {
  run_wgcna()
}
