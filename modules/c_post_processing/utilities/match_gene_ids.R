#!/usr/bin/env Rscript

# ===============================================
# GENE ID MATCHING UTILITY
# ===============================================
# Shared logic for matching a gene list to matrix row names, handling
# version suffixes common in eggplant IDs (e.g., SMEL4.1_06g023900.1.01).
#
# Used by: tximport_salmon_to_matrices.R, tximport_star_to_matrices.R,
#          tximport_rsem_to_matrices.R, 3_ranking_stability.R,
#          3_Matrix_Creation_utils.R

match_gene_ids <- function(gene_list, data_rownames) {
  if (is.null(gene_list) || length(gene_list) == 0) return(character(0))
  if (is.null(data_rownames) || length(data_rownames) == 0) return(character(0))

  # Single-pass regex handles both .X and .X.XX suffixes: O(n) vs 2 × O(n)
  base_ids <- sub("(\\.[0-9]+){1,2}$", "", data_rownames)

  # Build lookup: base_id -> data_rownames indices (vectorized, O(n) via split)
  # Use environment as hash map for O(1) lookups
  base_to_rows <- new.env(hash = TRUE, parent = emptyenv(), size = length(data_rownames))
  # split() + seq_along avoids O(n^2) c() concatenation that occurs with incremental appends
  # list2env is a single C-level call vs O(n) R-level loop
  idx_groups <- split(seq_along(data_rownames), base_ids)
  list2env(idx_groups, envir = base_to_rows)
  # Resolve indices to row names (deferred to avoid repeated string concatenation)
  .resolve_rows <- function(key) {
    idx <- base_to_rows[[key]]
    if (is.null(idx)) return(NULL)
    data_rownames[idx]
  }
  # Use %in% on data_rownames directly instead of building a second env
  # (match_gene_ids.R:59 was the only user of rowname_set)

  # Vectorized: exact matches first
  exact_mask <- gene_list %in% data_rownames
  matched_list <- list(gene_list[exact_mask])

  # Non-exact: try forward match (gene_list ID is base -> find suffixed data rows)
  # Pre-allocate list to avoid O(n^2) c() concatenation
  non_exact <- gene_list[!exact_mask]
  if (length(non_exact) > 0) {
    ne_results <- vector("list", length(non_exact))
    for (i in seq_along(non_exact)) {
      gene <- non_exact[i]
      hits <- .resolve_rows(gene)
      if (!is.null(hits)) {
        ne_results[[i]] <- hits
      } else {
        # Reverse: strip suffix from gene_list ID to match base-level row IDs
        # Use same double-suffix regex as line 18 for consistency
        gene_base <- sub("(\\.[0-9]+){1,2}$", "", gene)
        if (gene_base != gene) {
          hits2 <- .resolve_rows(gene_base)
          if (!is.null(hits2)) {
            ne_results[[i]] <- hits2
          } else if (gene_base %in% data_rownames) {
            ne_results[[i]] <- gene_base
          }
        }
      }
    }
    matched_list <- c(matched_list, ne_results)
  }
  unique(unlist(matched_list, use.names = FALSE))
}
