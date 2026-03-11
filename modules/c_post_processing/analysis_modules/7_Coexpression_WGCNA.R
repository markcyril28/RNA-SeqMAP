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
})

allowWGCNAThreads()

SCRIPT_DIR <- Sys.getenv("ANALYSIS_MODULES_DIR", ".")
source(file.path(SCRIPT_DIR, "0_shared_config.R"))
source(file.path(SCRIPT_DIR, "1_utility_functions.R"))

# ===============================================
# CONFIGURATION
# ===============================================

WGCNA_OUT_DIR <- file.path(CONSOLIDATED_BASE_DIR, OUTPUT_SUBDIRS$WGCNA)

# WGCNA parameters (STRINGENT settings for high-confidence results)
SOFT_POWER_RANGE <- 1:30           # Extended range for optimal fit
# MIN_MODULE_SIZE: Set dynamically based on gene count (default 30, min 10)
# For small gene groups (<100), use smaller module size
MIN_MODULE_SIZE_DEFAULT <- 50      # STRINGENT: Larger modules only (was 30)
MIN_MODULE_SIZE_SMALL <- 10        # STRINGENT: Min 10 genes per module (was 5)
MERGE_CUT_HEIGHT <- 0.15           # STRINGENT: More conservative module merging (was 0.10)
NETWORK_TYPE <- "signed"           # "signed" preserves biological interpretation
DEEP_SPLIT <- 2                    # STRINGENT: Moderate sensitivity (was 3, range 0-4)
PAM_STAGE <- TRUE                  # PAM for more accurate module assignment
# Note: MIN_GENES_WGCNA is defined in 0_shared_config.R (default: 20)

# Coexpression parameters (STRINGENT - high-confidence connections only)
N_HUB_GENES <- 5                   # STRINGENT: Top 5 hub genes per module (was 10)
N_COEXPRESSED_GENES <- 50          # STRINGENT: Top 50 co-expressed genes (was 200)
N_NETWORK_GENES <- 30              # STRINGENT: Smaller focused network (was 100)
COR_THRESHOLD <- 0.90              # STRINGENT: High correlation only (was 0.80)
TOP_VAR_GENES <- 5000              # Use top N variable genes from transcriptome for speed

# Low-expression gene filtering (removes unreliable near-zero expression genes)
# Genes must have expression >= MIN_EXPR_THRESHOLD in at least MIN_EXPR_SAMPLES samples
# Set MIN_EXPR_THRESHOLD to 0 to disable this filter
MIN_EXPR_THRESHOLD <- 1            # Minimum expression value (TPM/count) to consider "expressed"
MIN_EXPR_SAMPLES <- 2              # Minimum number of samples that must meet the threshold

# Figure toggles (all enabled for comprehensive output)
GENERATE_WGCNA_FIGURES <- list(
  soft_threshold_plot = TRUE,      # Scale-free topology fit
  module_dendrogram = TRUE,        # Gene clustering dendrogram with query genes marked
  eigengene_adjacency = TRUE,      # Module relationships
  correlation_network = TRUE,      # Network visualization with query genes bold
  query_gene_focus = TRUE          # Focused network around query genes only
)

# Output format toggles
OUTPUT_FORMATS <- list(
  png = TRUE,                      # Static PNG images
  html = TRUE                      # Interactive HTML (zoomable, hoverable)
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
    if (OUTPUT_FORMATS$png) {
      png(file.path(output_dir, paste0(gene_group, "_module_dendrogram.png")), 
          width = 1600, height = 800, res = 100)
      
      plotDendroAndColors(net$dendrograms[[1]], color_matrix,
                          dendroLabels = FALSE, hang = 0.03,
                          addGuide = TRUE, guideHang = 0.05,
                          main = paste0(gene_group, " - Module Dendrogram (Query genes in RED)"))
      dev.off()
    }
    
    # ===== HTML OUTPUT: INTERACTIVE MODULE SUMMARY =====
    if (OUTPUT_FORMATS$html) {
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
      
      for (i in seq_len(nrow(module_counts))) {
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
      for (i in seq_len(nrow(query_subset))) {
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
    if (OUTPUT_FORMATS$png) {
      png(file.path(output_dir, paste0(gene_group, "_eigengene_adjacency.png")), 
          width = 800, height = 800, res = 100)
      plotEigengeneNetworks(MEs, "", marDendro = c(0, 4, 1, 2), marHeatmap = c(3, 4, 1, 2),
                            plotDendrograms = TRUE, xLabelsAngle = 90)
      dev.off()
    }
    
    # ===== HTML OUTPUT: INTERACTIVE EIGENGENE HEATMAP =====
    if (OUTPUT_FORMATS$html) {
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
          for (i in seq_len(nrow(ME_cor))) {
            html_content <- paste0(html_content, "<tr><th>", rownames(ME_cor)[i], "</th>")
            for (j in seq_len(ncol(ME_cor))) {
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
  
  # Get top correlated genes for each query gene
  all_genes <- matched
  for (qg in matched) {
    cors <- cor_matrix[qg, ]
    cors <- cors[!names(cors) %in% matched]
    top_cors <- names(sort(abs(cors), decreasing = TRUE)[1:min(N_NETWORK_GENES, length(cors))])
    all_genes <- unique(c(all_genes, top_cors))
  }
  
  # Create adjacency for network
  sub_cor <- cor_matrix[all_genes, all_genes]
  adj <- abs(sub_cor)
  adj[adj < COR_THRESHOLD] <- 0
  diag(adj) <- 0
  
  # Create igraph network
  g <- graph_from_adjacency_matrix(adj, mode = "undirected", weighted = TRUE)
  V(g)$is_query <- V(g)$name %in% matched
  
  # Query genes: RED, BOLD, LARGER | Others: blue, normal, smaller
  colors <- ifelse(V(g)$is_query, "#B2182B", "#4393C3")
  sizes <- ifelse(V(g)$is_query, 18, 6)
  label_sizes <- ifelse(V(g)$is_query, 1.2, 0.5)
  label_colors <- ifelse(V(g)$is_query, "black", "gray40")
  frame_widths <- ifelse(V(g)$is_query, 3, 1)
  
  # ===== PNG OUTPUT: FULL NETWORK =====
  if (OUTPUT_FORMATS$png) {
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
         vertex.label.font = ifelse(V(g)$is_query, 2, 1),  # 2 = BOLD
         vertex.frame.width = frame_widths,
         vertex.frame.color = ifelse(V(g)$is_query, "black", NA),
         edge.width = E(g)$weight * 3,
         edge.color = adjustcolor("gray50", alpha.f = 0.5),
         main = paste0(gene_group, " - Correlation Network\n(Query genes in RED/BOLD)"))
    
    legend("bottomleft", 
           legend = c(paste0("Query genes (n=", sum(V(g)$is_query), ")"), 
                      paste0("Correlated genes (n=", sum(!V(g)$is_query), ")")),
           col = c("#B2182B", "#4393C3"), 
           pch = 19, 
           pt.cex = c(2, 1.2),
           bty = "n",
           cex = 0.9)
    
    dev.off()
  }
  
  # ===== HTML OUTPUT: INTERACTIVE FULL NETWORK =====
  if (OUTPUT_FORMATS$html) {
    # Create nodes dataframe for visNetwork
    nodes <- data.frame(
      id = V(g)$name,
      label = V(g)$name,
      group = ifelse(V(g)$is_query, "Query Gene", "Correlated Gene"),
      color = ifelse(V(g)$is_query, "#B2182B", "#4393C3"),
      size = ifelse(V(g)$is_query, 30, 15),
      font.size = ifelse(V(g)$is_query, 18, 12),
      font.bold = V(g)$is_query,
      borderWidth = ifelse(V(g)$is_query, 3, 1),
      title = paste0("<b>", V(g)$name, "</b><br>",
                     ifelse(V(g)$is_query, "Query Gene", "Correlated Gene")),
      stringsAsFactors = FALSE
    )
    
    # Create edges dataframe
    edge_list <- as_data_frame(g, what = "edges")
    if (nrow(edge_list) > 0) {
      edges <- data.frame(
        from = edge_list$from,
        to = edge_list$to,
        width = edge_list$weight * 5,
        color = "rgba(100, 100, 100, 0.5)",
        title = paste0("Correlation: ", round(edge_list$weight, 3)),
        stringsAsFactors = FALSE
      )
    } else {
      edges <- data.frame(from = character(), to = character(), 
                          width = numeric(), color = character(), 
                          title = character(), stringsAsFactors = FALSE)
    }
    
    # Create interactive network
    vis_net <- visNetwork(nodes, edges, 
                          main = paste0(gene_group, " - Correlation Network (Interactive)"),
                          submain = "Query genes in RED (larger). Hover for details. Zoom/pan enabled.") %>%
      visGroups(groupname = "Query Gene", color = "#B2182B", 
                font = list(size = 18, bold = TRUE)) %>%
      visGroups(groupname = "Correlated Gene", color = "#4393C3",
                font = list(size = 12)) %>%
      visLegend(addNodes = list(
        list(label = "Query Gene", shape = "dot", color = "#B2182B", size = 20),
        list(label = "Correlated Gene", shape = "dot", color = "#4393C3", size = 12)
      ), useGroups = FALSE) %>%
      visOptions(highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE),
                 nodesIdSelection = TRUE,
                 selectedBy = "group") %>%
      visInteraction(navigationButtons = TRUE,
                     keyboard = TRUE,
                     dragNodes = TRUE,
                     dragView = TRUE,
                     zoomView = TRUE) %>%
      visPhysics(solver = "forceAtlas2Based",
                 forceAtlas2Based = list(gravitationalConstant = -50))
    
    # Save HTML
    saveWidget(vis_net, 
               file.path(output_dir, paste0(gene_group, "_correlation_network.html")),
               selfcontained = TRUE)
    cat("  Saved interactive HTML:", paste0(gene_group, "_correlation_network.html"), "\n")
  }
  
  # ===== QUERY GENE FOCUSED NETWORK =====
  if (GENERATE_WGCNA_FIGURES$query_gene_focus && length(matched) >= 2) {
    # Create network showing only query genes and direct connections between them
    query_sub_cor <- cor_matrix[matched, matched]
    query_adj <- abs(query_sub_cor)
    query_adj[query_adj < 0.7] <- 0  # STRINGENT: Higher threshold for query-only network
    diag(query_adj) <- 0
    
    g_query <- graph_from_adjacency_matrix(query_adj, mode = "undirected", weighted = TRUE)
    
    # ===== PNG OUTPUT: QUERY GENES NETWORK =====
    if (OUTPUT_FORMATS$png) {
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
           vertex.label.font = 2,  # BOLD
           vertex.frame.width = 3,
           vertex.frame.color = "black",
           edge.width = E(g_query)$weight * 5,
           edge.color = adjustcolor("darkred", alpha.f = 0.6),
           edge.label = round(E(g_query)$weight, 2),
           edge.label.cex = 0.7,
           main = paste0(gene_group, " - Query Genes Co-expression\n(Edge weights = correlation)"))
      
      dev.off()
    }
    
    # ===== HTML OUTPUT: INTERACTIVE QUERY GENES NETWORK =====
    if (OUTPUT_FORMATS$html) {
      # Create nodes dataframe
      nodes_q <- data.frame(
        id = V(g_query)$name,
        label = V(g_query)$name,
        color = "#B2182B",
        size = 35,
        font.size = 16,
        font.bold = TRUE,
        borderWidth = 3,
        title = paste0("<b>", V(g_query)$name, "</b><br>Query Gene"),
        stringsAsFactors = FALSE
      )
      
      # Create edges dataframe
      edge_list_q <- as_data_frame(g_query, what = "edges")
      if (nrow(edge_list_q) > 0) {
        edges_q <- data.frame(
          from = edge_list_q$from,
          to = edge_list_q$to,
          width = edge_list_q$weight * 8,
          label = as.character(round(edge_list_q$weight, 2)),
          color = "rgba(139, 0, 0, 0.7)",
          title = paste0("Correlation: ", round(edge_list_q$weight, 3)),
          font.size = 14,
          stringsAsFactors = FALSE
        )
      } else {
        edges_q <- data.frame(from = character(), to = character(),
                              width = numeric(), label = character(),
                              color = character(), title = character(),
                              font.size = numeric(), stringsAsFactors = FALSE)
      }
      
      # Create interactive network
      vis_query <- visNetwork(nodes_q, edges_q,
                              main = paste0(gene_group, " - Query Genes Co-expression (Interactive)"),
                              submain = "Edge labels show correlation values. Drag nodes to rearrange.") %>%
        visOptions(highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE)) %>%
        visInteraction(navigationButtons = TRUE,
                       keyboard = TRUE,
                       dragNodes = TRUE,
                       dragView = TRUE,
                       zoomView = TRUE) %>%
        visPhysics(solver = "forceAtlas2Based",
                   forceAtlas2Based = list(gravitationalConstant = -100))
      
      # Save HTML
      saveWidget(vis_query,
                 file.path(output_dir, paste0(gene_group, "_query_genes_network.html")),
                 selfcontained = TRUE)
      cat("  Saved interactive HTML:", paste0(gene_group, "_query_genes_network.html"), "\n")
    }
  }
  
  return(g)
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
    
    output_folder_name <- get_output_folder_name(gene_group, CURRENT_DATASET)
    output_dir <- file.path(WGCNA_OUT_DIR, output_folder_name)
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
    
    # ===== STEP 4: Low-expression filtering (pre-log transformation) =====
    # Filter out genes with near-zero expression that produce unreliable correlations
    if (MIN_EXPR_THRESHOLD > 0) {
      raw_data <- validation$data
      # Count how many samples have expression >= threshold for each gene
      samples_expressed <- apply(raw_data, 1, function(x) sum(x >= MIN_EXPR_THRESHOLD, na.rm = TRUE))
      
      # Keep genes expressed in at least MIN_EXPR_SAMPLES samples
      # BUT always keep query genes regardless of expression level
      expr_pass <- samples_expressed >= MIN_EXPR_SAMPLES | rownames(raw_data) %in% query_genes_matched
      
      # Track which query genes have low expression (for warning)
      low_expr_query <- query_genes_matched[!query_genes_matched %in% rownames(raw_data)[samples_expressed >= MIN_EXPR_SAMPLES]]
      if (length(low_expr_query) > 0) {
        cat("  WARNING: Low-expression query genes (kept but may have unreliable correlations):\n")
        cat("    ", paste(low_expr_query, collapse = ", "), "\n")
      }
      
      validation$data <- raw_data[expr_pass, , drop = FALSE]
      cat("  Genes after expression filter (>=", MIN_EXPR_THRESHOLD, " in >=", MIN_EXPR_SAMPLES, " samples):", 
          nrow(validation$data), "\n")
    }
    
    # ===== STEP 5: Variance stabilization and filtering =====
    data_log <- log2(validation$data + 1)
    
    # Filter genes with low variance
    gene_vars <- apply(data_log, 1, var, na.rm = TRUE)
    var_threshold <- quantile(gene_vars, 0.25, na.rm = TRUE)
    
    # IMPORTANT: Always keep query genes even if low variance
    keep_genes <- gene_vars > var_threshold | rownames(data_log) %in% query_genes_matched
    data_filtered <- data_log[keep_genes, , drop = FALSE]
    
    cat("  Genes after variance filtering:", nrow(data_filtered), "\n")
    cat("  Query genes retained:", sum(query_genes_matched %in% rownames(data_filtered)), "\n")
    
    # ===== STEP 6: Subsample for speed if too many genes =====
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
    
    # ===== STEP 7: Run WGCNA pipeline =====
    cat("  Running soft threshold selection...\n")
    soft_power <- pick_soft_threshold(data_matrix, output_dir, gene_group)
    cat("  Soft power:", soft_power, "\n")
    
    cat("  Building network and detecting modules...\n")
    network <- build_network_and_detect_modules(data_matrix, soft_power, output_dir, 
                                                 gene_group, query_genes_matched)
    cat("  Modules detected:", length(unique(network$module_colors)), "\n")
    
    MEs <- calculate_module_eigengenes(data_matrix, network$module_colors, output_dir, gene_group)
    
    # ===== STEP 8: Export gene-module assignments (with query gene flag) =====
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
    
    # ===== STEP 9: Query gene module summary =====
    query_gene_info <- gene_info[gene_info$Is_Query_Gene, ]
    cat("  Query genes module distribution:\n")
    print(table(query_gene_info$Module))
    
    write.table(query_gene_info, 
                file.path(output_dir, paste0(gene_group, "_query_genes_modules.tsv")),
                sep = "\t", row.names = FALSE, quote = FALSE)
    
    # ===== STEP 10: Hub genes =====
    hubs <- identify_hub_genes(gene_info)
    # Mark query genes in hub list
    hubs$Is_Query_Gene <- hubs$Gene %in% query_genes_matched
    write.table(hubs, file.path(output_dir, paste0(gene_group, "_hub_genes.tsv")),
                sep = "\t", row.names = FALSE, quote = FALSE)
    
    # ===== STEP 11: Correlation network with query genes bold =====
    cat("  Creating correlation network...\n")
    create_correlation_network(t(data_matrix), query_genes_matched, output_dir, gene_group)
    
    # ===== STEP 12: Find genes co-expressed with query genes =====
    cat("  Finding genes co-expressed with query genes...\n")
    cor_matrix <- gpu_cor(t(data_filtered))
    
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
    
    # ===== STEP 13: Export RAW RESULTS for downstream analysis =====
    cat("  Exporting raw results...\n")
    raw_results_dir <- file.path(output_dir, "raw_results")
    if (!dir.exists(raw_results_dir)) dir.create(raw_results_dir, recursive = TRUE)
    
    # 13a: Save full correlation matrix as TSV
    cat("    Saving correlation matrix...\n")
    cor_df <- as.data.frame(cor_matrix)
    cor_df$Gene <- rownames(cor_matrix)
    cor_df <- cor_df[, c("Gene", setdiff(names(cor_df), "Gene"))]  # Move Gene to first column
    write.table(cor_df, file.path(raw_results_dir, paste0(gene_group, "_correlation_matrix.tsv")),
                sep = "\t", row.names = FALSE, quote = FALSE)
    
    # 13b: Save module eigengenes
    cat("    Saving module eigengenes...\n")
    me_df <- as.data.frame(MEs)
    me_df$Sample <- rownames(MEs)
    me_df <- me_df[, c("Sample", setdiff(names(me_df), "Sample"))]
    write.table(me_df, file.path(raw_results_dir, paste0(gene_group, "_module_eigengenes.tsv")),
                sep = "\t", row.names = FALSE, quote = FALSE)
    
    # 13c: Save network as RDS for complete reproducibility
    cat("    Saving network object (RDS)...\n")
    saveRDS(network, file.path(raw_results_dir, paste0(gene_group, "_network.rds")))
    
    # 13d: Save gene info with kME as RDS
    saveRDS(gene_info, file.path(raw_results_dir, paste0(gene_group, "_gene_info.rds")))
    
    # 13e: Save soft threshold plot (PNG saved by pick_soft_threshold())
    sft_file <- file.path(output_dir, paste0(gene_group, "_soft_threshold.png"))
    if (file.exists(sft_file)) {
      file.copy(sft_file, file.path(raw_results_dir, paste0(gene_group, "_soft_threshold.png")))
    }
    
    # 13f: Create a summary file with all parameters used
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
    
    # 13g: Save data matrix used (expression values)
    cat("    Saving expression matrix...\n")
    expr_df <- as.data.frame(t(data_matrix))  # Genes as rows
    expr_df$Gene <- rownames(expr_df)
    expr_df <- expr_df[, c("Gene", setdiff(names(expr_df), "Gene"))]
    write.table(expr_df, file.path(raw_results_dir, paste0(gene_group, "_expression_matrix.tsv")),
                sep = "\t", row.names = FALSE, quote = FALSE)
    
    cat("  Raw results saved to:", raw_results_dir, "\n")
    
    successful <- successful + 1
    cat("  Complete!\n")
  }
  
  print_summary(successful, total)
}

if (!interactive() && identical(environment(), globalenv())) {
  run_wgcna()
}
