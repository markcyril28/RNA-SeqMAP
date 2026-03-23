#!/usr/bin/env Rscript

# ===============================================
# SHARED METHOD LOADERS FOR CONCORDANCE ANALYSIS
# ===============================================
# Per-method TPM loading functions used by cross-genome, cross-gene-group,
# and cross-method concordance loaders.
#
# Prerequisites: source 0_concordance_config.R first (defines ALIGNMENT_BASE,
# POST_PROC_BASE, get_method_ref_dir, etc.)
#
# Usage:
#   source("0_concordance_config.R")
#   source("0_method_loaders.R")
#   mat <- load_method_for_ref(method, ref_dir)

if (exists(".METHOD_LOADERS_SOURCED") && .METHOD_LOADERS_SOURCED) {
  # Already sourced — skip
} else {
  .METHOD_LOADERS_SOURCED <- TRUE

  # -----------------------------------------------
  # Dependencies: data.table, parallel
  # -----------------------------------------------

  if (!exists(".use_dt")) .use_dt <- requireNamespace("data.table", quietly = TRUE)
  if (!exists(".use_parallel")) .use_parallel <- requireNamespace("parallel", quietly = TRUE)
  if (!exists(".n_cores")) {
    .n_cores <- if (.use_parallel) {
      .threads_env <- as.integer(Sys.getenv("THREADS", unset = "0"))
      if (!is.na(.threads_env) && .threads_env > 1) min(.threads_env, 8L) else {
        max(1L, parallel::detectCores(logical = FALSE) %/% 2L)
      }
    } else 1L
  }
  if (!exists(".par_lapply")) {
    # Threshold: mclapply fork+merge overhead (~5-10ms per worker) exceeds benefit
    # for fewer than 5 items. For 1-4 samples, sequential lapply is faster.
    .par_lapply <- function(X, FUN, ...) {
      if (.use_parallel && .n_cores > 1L && length(X) > 4L) {
        tryCatch(
          parallel::mclapply(X, FUN, ..., mc.cores = .n_cores),
          error = function(e) { message("[CONCORDANCE] mclapply failed: ", e$message); lapply(X, FUN, ...) }
        )
      } else {
        lapply(X, FUN, ...)
      }
    }
  }

  # -----------------------------------------------
  # Shared helpers
  # -----------------------------------------------

  if (!exists(".fast_read_tsv")) {
    .fast_read_tsv <- function(path, ...) {
      if (.use_dt) data.table::fread(path, data.table = FALSE, ...)
      else read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                      check.names = FALSE, comment.char = "", quote = "")
    }
  }

  if (!exists(".assemble_tpm_matrix")) {
    .assemble_tpm_matrix <- function(tpm_list, filter_genes = TRUE) {
      if (length(tpm_list) == 0) return(NULL)
      if (.use_dt) {
        dt_list <- lapply(names(tpm_list), function(srr) {
          data.table::data.table(gene = names(tpm_list[[srr]]),
                                 tpm  = as.numeric(tpm_list[[srr]]),
                                 srr  = srr)
        })
        long <- data.table::rbindlist(dt_list)
        if (filter_genes) long <- long[nzchar(gene) & !is.na(gene)]
        wide <- data.table::dcast(long, gene ~ srr, value.var = "tpm", fill = 0)
        mat <- as.matrix(wide[, -1, with = FALSE])
        rownames(mat) <- wide$gene
        if (nrow(mat) == 0) { warning("Matrix assembly produced 0 rows"); return(NULL) }
      } else {
        # use.names=FALSE avoids allocating name attributes on the O(G×S) intermediate vector.
        # Memory peak: O(G×S) for the unlisted vector + O(G) for unique result.
        # Acceptable for typical RNA-seq (20K genes × 50 samples); data.table path preferred for larger datasets.
        all_genes <- unique(unlist(lapply(tpm_list, names), use.names = FALSE))
        if (filter_genes) all_genes <- all_genes[nzchar(all_genes) & !is.na(all_genes)]
        mat <- matrix(0, nrow = length(all_genes), ncol = length(tpm_list),
                      dimnames = list(all_genes, names(tpm_list)))
        # Pre-compute gene->row index map: O(G) once, then O(1) named lookup per sample.
        # Replaces O(G) intersect() per sample → total O(G + S×g) vs prior O(S×G).
        gene_idx <- setNames(seq_along(all_genes), all_genes)
        for (srr in names(tpm_list)) {
          sv <- tpm_list[[srr]]
          idx <- gene_idx[names(sv)]
          valid <- !is.na(idx)
          mat[idx[valid], srr] <- sv[valid]
        }
        if (nrow(mat) == 0) { warning("Matrix assembly produced 0 rows"); return(NULL) }
      }
      return(mat)
    }
  }

  # -----------------------------------------------
  # Reference directory resolver
  # -----------------------------------------------
  # Handles naming inconsistencies (e.g., "GPE001970_genome" vs "GPE001970",
  # "Eggplant_V4.1_genome" vs "Eggplant_V4.1") by trying suffix variants.

  .resolve_ref_dir <- function(parent_dir, ref_dir) {
    if (!dir.exists(parent_dir)) return(NULL)
    # Single list.dirs() up front replaces up to 4× dir.exists() + 4× list.dirs()
    # + 1× list.dirs(parent_dir) = 9 syscalls → 1 listing + O(1) match() lookups.
    # Called 1-2× per method × 5 methods per concordance run.
    .all_children <- list.dirs(parent_dir, recursive = FALSE, full.names = TRUE)
    if (length(.all_children) == 0L) return(NULL)
    .child_names <- basename(.all_children)
    # Helper: TRUE if dir has at least one subdirectory (content check)
    .has_content <- function(d) {
      length(list.dirs(d, recursive = FALSE, full.names = FALSE)) > 0
    }

    # Try exact match via O(1) hash lookup instead of filesystem stat
    .idx <- match(ref_dir, .child_names)
    if (!is.na(.idx) && .has_content(.all_children[.idx])) return(.all_children[.idx])
    # Try without _genome / _transcripts suffix
    stripped <- sub("_(genome|transcripts)(\\..+)?$", "", ref_dir)
    if (stripped != ref_dir) {
      .idx <- match(stripped, .child_names)
      if (!is.na(.idx) && .has_content(.all_children[.idx])) return(.all_children[.idx])
    }
    # Try adding _genome suffix
    .idx <- match(paste0(ref_dir, "_genome"), .child_names)
    if (!is.na(.idx) && .has_content(.all_children[.idx])) return(.all_children[.idx])
    # Glob: match any dir starting with the base name (reuses cached listing)
    base_pattern <- paste0("^", gsub("([.()])", "\\\\\\1", stripped))
    matches <- .all_children[grepl(base_pattern, .child_names)]
    # Prefer directories with content; fall back to any match
    with_content <- matches[vapply(matches, .has_content, logical(1))]
    if (length(with_content) >= 1L) return(with_content[1L])
    if (length(matches) == 1L) return(matches[1L])
    return(NULL)
  }

  # -----------------------------------------------
  # Per-method TPM loaders
  # -----------------------------------------------
  # Each takes a reference directory name and returns a genes x samples matrix.

  .load_stringtie_tpm <- function(method, ref_dir) {
    base_dir <- .resolve_ref_dir(file.path(ALIGNMENT_BASE, method, "stringtie_WD"), ref_dir)
    if (is.null(base_dir)) { cat(" [", ref_dir, "] Not found in:", file.path(ALIGNMENT_BASE, method, "stringtie_WD"), "\n"); return(NULL) }
    sample_dirs <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)
    sample_dirs <- sample_dirs[startsWith(basename(sample_dirs), "SRR")]
    if (length(sample_dirs) == 0) { cat(" [", ref_dir, "] No samples\n"); return(NULL) }

    tpm_list <- setNames(.par_lapply(sample_dirs, function(sdir) {
      abundance_files <- list.files(sdir, pattern = "gene_abundances.*\\.tsv$", full.names = TRUE)
      if (length(abundance_files) == 0) return(NULL)
      df <- tryCatch(.fast_read_tsv(abundance_files[1]), error = function(e) NULL)
      if (is.null(df) || !"TPM" %in% colnames(df)) return(NULL)
      setNames(df$TPM, df[[1]])
    }), basename(sample_dirs))
    .assemble_tpm_matrix(tpm_list[lengths(tpm_list) > 0L], filter_genes = TRUE)
  }

  .load_stringtie_denovo_tpm <- function(method, ref_dir) {
    base_dir <- .resolve_ref_dir(file.path(ALIGNMENT_BASE, method, "stringtie_WD"), ref_dir)
    if (is.null(base_dir)) { cat(" [", ref_dir, "] Not found in:", file.path(ALIGNMENT_BASE, method, "stringtie_WD"), "\n"); return(NULL) }
    sample_dirs <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)
    sample_dirs <- sample_dirs[startsWith(basename(sample_dirs), "SRR")]
    if (length(sample_dirs) == 0) { cat(" [", ref_dir, "] No samples\n"); return(NULL) }

    tpm_list <- setNames(.par_lapply(sample_dirs, function(sdir) {
      abundance_files <- list.files(sdir, pattern = "gene_abundances.*\\.tsv$", full.names = TRUE)
      if (length(abundance_files) == 0) return(NULL)
      df <- .fast_read_tsv(abundance_files[1])
      if (!all(c("TPM", "Reference") %in% colnames(df))) return(NULL)
      gene_ids <- sub("(\\.[0-9]+){1,2}$", "", df$Reference)
      valid <- nzchar(gene_ids) & !gene_ids %in% c(".", "-")
      # Use data.table grouped-max when available (~2-5x faster than tapply for large gene sets)
      if (.use_dt) {
        dt <- data.table::data.table(gene = gene_ids[valid], tpm = df$TPM[valid])
        res <- dt[, .(tpm = max(tpm, na.rm = TRUE)), by = gene]
        setNames(res$tpm, res$gene)
      } else {
        tapply(df$TPM[valid], gene_ids[valid], max, na.rm = TRUE)
      }
    }), basename(sample_dirs))
    .assemble_tpm_matrix(tpm_list[lengths(tpm_list) > 0L], filter_genes = TRUE)
  }

  .load_star_salmon_tpm <- function(method, ref_dir) {
    # Resolve ref_dir against alignment directory (e.g., GPE001970 vs GPE001970_genome)
    resolved_align <- .resolve_ref_dir(file.path(ALIGNMENT_BASE, method), ref_dir)
    actual_ref_name <- if (!is.null(resolved_align)) basename(resolved_align) else ref_dir
    quant_base <- file.path(ALIGNMENT_BASE, method, actual_ref_name, "6_salmon", "quant")

    # Resolve ref_dir against post-proc directory too
    resolved_pp <- .resolve_ref_dir(file.path(POST_PROC_BASE, method, "count_matrices_from_STAR"), ref_dir)
    pp_ref_name <- if (!is.null(resolved_pp)) basename(resolved_pp) else ref_dir
    tx2gene_file <- file.path(POST_PROC_BASE, method, "count_matrices_from_STAR", pp_ref_name,
                              paste0("tx2gene_", pp_ref_name, ".tsv"))
    tx2gene <- NULL
    if (file.exists(tx2gene_file)) {
      tx2gene <- .fast_read_tsv(tx2gene_file)
      if (ncol(tx2gene) > 2) tx2gene <- tx2gene[, 1:2, drop = FALSE]
      colnames(tx2gene) <- c("transcript_id", "gene_id")
    }

    sample_dirs <- if (dir.exists(quant_base)) {
      list.dirs(quant_base, recursive = FALSE, full.names = TRUE)
    } else character(0)
    sample_dirs <- sample_dirs[startsWith(basename(sample_dirs), "SRR")]

    # Tissue-specific fallback
    if (length(sample_dirs) == 0 && dir.exists(quant_base)) {
      tissue_dirs <- list.dirs(quant_base, recursive = FALSE, full.names = TRUE)
      # unlist(lapply()) avoids O(T²) c() accumulation overhead
      sample_dirs <- unlist(lapply(tissue_dirs, function(td) {
        sub_dirs <- list.dirs(td, recursive = FALSE, full.names = TRUE)
        sub_dirs[startsWith(basename(sub_dirs), "SRR")]
      }), use.names = FALSE)
    }

    if (length(sample_dirs) == 0) {
      tpm_search_base <- if (!is.null(resolved_pp)) resolved_pp else {
        file.path(POST_PROC_BASE, method, "count_matrices_from_STAR", ref_dir)
      }
      tpm_files <- list.files(tpm_search_base, pattern = "_tpm_Gene_ID_.*\\.csv$",
                              recursive = TRUE, full.names = TRUE)
      if (length(tpm_files) > 0) {
        cat(" [", ref_dir, "] Using pre-built TPM:", tpm_files[1], "\n")
        df <- if (.use_dt) data.table::fread(tpm_files[1], data.table = FALSE) else {
          read.csv(tpm_files[1], header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
        }
        rownames(df) <- df[[1]]; df <- df[, -1, drop = FALSE]
        return(as.matrix(df))
      }
      cat(" [", ref_dir, "] No quant data for M3\n")
      return(NULL)
    }

    # Pre-build named vector for O(1) hash lookup inside parallel workers.
    # Replaces per-worker match() which is O(G_transcripts) linear scan per sample.
    .tx2gene_map <- if (!is.null(tx2gene)) setNames(tx2gene$gene_id, tx2gene$transcript_id) else NULL

    tpm_list <- setNames(.par_lapply(sample_dirs, function(sdir) {
      qsf <- file.path(sdir, "quant.sf")
      if (!file.exists(qsf)) return(NULL)
      df <- .fast_read_tsv(qsf)
      if (!is.null(.tx2gene_map)) {
        df$gene_id <- .tx2gene_map[df$Name]
        unmapped <- is.na(df$gene_id)
        if (any(unmapped)) df$gene_id[unmapped] <- sub("(\\.[0-9]+){1,2}$", "", df$Name[unmapped])
      } else {
        df$gene_id <- sub("(\\.[0-9]+){1,2}$", "", df$Name)
      }
      rs <- rowsum(df$TPM, df$gene_id, reorder = FALSE, na.rm = TRUE)
      setNames(rs[, 1], rownames(rs))
    }), basename(sample_dirs))
    .assemble_tpm_matrix(tpm_list[lengths(tpm_list) > 0L], filter_genes = TRUE)
  }

  .load_salmon_saf_tpm <- function(method, ref_dir) {
    base_dir <- .resolve_ref_dir(file.path(ALIGNMENT_BASE, method, "Salmon_Quant"), ref_dir)
    if (is.null(base_dir)) {
      resolved_pp <- .resolve_ref_dir(file.path(POST_PROC_BASE, method, "count_matrices_from_Salmon_Quant"), ref_dir)
      tpm_search_base <- if (!is.null(resolved_pp)) resolved_pp else {
        file.path(POST_PROC_BASE, method, "count_matrices_from_Salmon_Quant", ref_dir)
      }
      tpm_files <- list.files(tpm_search_base, pattern = "_tpm_Gene_ID_.*\\.csv$",
                              recursive = TRUE, full.names = TRUE)
      if (length(tpm_files) > 0) {
        cat(" [", ref_dir, "] Using pre-built TPM:", tpm_files[1], "\n")
        df <- if (.use_dt) data.table::fread(tpm_files[1], data.table = FALSE) else {
          read.csv(tpm_files[1], header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
        }
        rownames(df) <- df[[1]]; df <- df[, -1, drop = FALSE]
        return(as.matrix(df))
      }
      cat(" [", ref_dir, "] Salmon quant not found\n")
      return(NULL)
    }

    sample_dirs <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)
    sample_dirs <- sample_dirs[startsWith(basename(sample_dirs), "SRR")]
    if (length(sample_dirs) == 0) { cat(" [", ref_dir, "] No samples\n"); return(NULL) }

    tpm_list <- setNames(.par_lapply(sample_dirs, function(sdir) {
      qsf <- file.path(sdir, "quant.sf")
      if (!file.exists(qsf)) return(NULL)
      df <- .fast_read_tsv(qsf)
      df$gene_id <- sub("(\\.[0-9]+){1,2}$", "", df$Name)
      rs <- rowsum(df$TPM, df$gene_id, reorder = FALSE, na.rm = TRUE)
      setNames(rs[, 1], rownames(rs))
    }), basename(sample_dirs))
    .assemble_tpm_matrix(tpm_list[lengths(tpm_list) > 0L], filter_genes = TRUE)
  }

  .load_rsem_tpm <- function(method, ref_dir) {
    resolved <- .resolve_ref_dir(file.path(ALIGNMENT_BASE, method, "RSEM_Quant_WD"), ref_dir)
    rsem_quant_base <- if (!is.null(resolved)) resolved else {
      file.path(ALIGNMENT_BASE, method, "RSEM_Quant_WD", ref_dir)
    }

    if (dir.exists(rsem_quant_base)) {
      sample_dirs <- list.dirs(rsem_quant_base, recursive = FALSE, full.names = TRUE)
      sample_dirs <- sample_dirs[startsWith(basename(sample_dirs), "SRR")]
      if (length(sample_dirs) > 0) {
        tpm_list <- setNames(.par_lapply(sample_dirs, function(sdir) {
          srr <- basename(sdir)
          results_file <- file.path(sdir, paste0(srr, ".genes.results"))
          if (!file.exists(results_file)) return(NULL)
          df <- .fast_read_tsv(results_file)
          gene_ids <- sub("(\\.[0-9]+){1,2}$", "", df$gene_id)
          rs <- rowsum(df$TPM, gene_ids, reorder = FALSE, na.rm = TRUE)
          setNames(rs[, 1], rownames(rs))
        }), basename(sample_dirs))
        tpm_list <- tpm_list[lengths(tpm_list) > 0L]
        if (length(tpm_list) > 0) {
          mat <- .assemble_tpm_matrix(tpm_list, filter_genes = TRUE)
          if (!is.null(mat)) return(mat)
        }
      }
    }

    resolved_pp <- .resolve_ref_dir(file.path(POST_PROC_BASE, method, "count_matrices_from_RSEM_Quant"), ref_dir)
    tpm_search_base <- if (!is.null(resolved_pp)) resolved_pp else {
      file.path(POST_PROC_BASE, method, "count_matrices_from_RSEM_Quant", ref_dir)
    }
    tpm_files <- list.files(tpm_search_base, pattern = "_tpm_Gene_ID_.*\\.csv$",
                            recursive = TRUE, full.names = TRUE)
    if (length(tpm_files) > 0) {
      cat(" [", ref_dir, "] Using pre-built TPM:", tpm_files[1], "\n")
      df <- if (.use_dt) data.table::fread(tpm_files[1], data.table = FALSE) else {
        read.csv(tpm_files[1], header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
      }
      rownames(df) <- df[[1]]; df <- df[, -1, drop = FALSE]
      return(as.matrix(df))
    }

    cat(" [", ref_dir, "] RSEM TPM not found\n")
    return(NULL)
  }

  # -----------------------------------------------
  # Unified loader: load one method for one reference
  # -----------------------------------------------

  load_method_for_ref <- function(method, ref_dir) {
    if (method == "M1_HISAT2_RefGuided") return(.load_stringtie_tpm(method, ref_dir))
    if (method == "M2_HISAT2_DeNovo")    return(.load_stringtie_denovo_tpm(method, ref_dir))
    if (method == "M3_STAR_Align")       return(.load_star_salmon_tpm(method, ref_dir))
    if (method == "M4_Salmon_Saf")       return(.load_salmon_saf_tpm(method, ref_dir))
    if (method == "M5_RSEM_Bowtie2")     return(.load_rsem_tpm(method, ref_dir))
    cat("  [WARN] Unknown method:", method, "\n")
    return(NULL)
  }
}
