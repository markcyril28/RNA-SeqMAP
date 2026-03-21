#!/usr/bin/env Rscript

# ===============================================
# CROSS-METHOD CONCORDANCE - LOAD & HARMONIZE MATRICES
# ===============================================
# Loads gene-level TPM matrices from all 5 methods for a given reference,
# harmonizes gene IDs and sample sets, and saves the result as an RDS file.
#
# Data sources per method:
#   M1 (HISAT2+StringTie ref-guided): per-sample gene_abundances.tsv -> TPM column
#   M2 (HISAT2+StringTie de novo):    per-sample gene_abundances.tsv -> TPM column (STRG.N -> Reference mapping)
#   M3 (STAR+Salmon):                 per-sample quant.sf -> TPM column (transcript -> gene aggregation)
#   M4 (Salmon pseudo-align):         pre-built gene_count_matrix.csv or genes.counts.matrix + tximport for TPM
#   M5 (Bowtie2+RSEM):                per-sample .genes.results -> TPM column (fallback: tximport _tpm_Gene_ID_*.csv)
#
# Output: HARMONIZED_RDS containing a list with:
#   $tpm_matrices  - named list of gene x sample TPM matrices (one per method)
#   $common_genes  - character vector of genes present in all methods
#   $common_samples - character vector of samples present in all methods
#   $method_stats  - data.frame with per-method stats (n_genes, n_samples, etc.)

source(file.path(Sys.getenv("CONCORDANCE_SCRIPT_DIR", "."), "0_concordance_config.R"))

# Use data.table for fast file I/O when available
.use_dt <- requireNamespace("data.table", quietly = TRUE)

# Use parallel::mclapply for concurrent per-sample file reads when available.
# On HPC (Linux), mclapply uses fork-based parallelism with zero serialization overhead.
# On Windows/WSL with fork issues, falls back gracefully to sequential lapply.
# Big O: reduces wall-clock from O(S × read_time) to O(S × read_time / ncores).
.use_parallel <- requireNamespace("parallel", quietly = TRUE)
.n_cores <- if (.use_parallel) {
  # Use THREADS env var if set (consistent with pipeline config), else detectCores/2
  .threads_env <- as.integer(Sys.getenv("THREADS", unset = "0"))
  if (!is.na(.threads_env) && .threads_env > 1) min(.threads_env, 8L) else {
    max(1L, parallel::detectCores(logical = FALSE) %/% 2L)
  }
} else 1L
# Wrapper: mclapply with safe fallback
.par_lapply <- function(X, FUN, ...) {
  if (.use_parallel && .n_cores > 1L && length(X) > 1L) {
    tryCatch(
      parallel::mclapply(X, FUN, ..., mc.cores = .n_cores),
      error = function(e) { message("[CONCORDANCE] mclapply failed, using sequential: ", e$message); lapply(X, FUN, ...) }
    )
  } else {
    lapply(X, FUN, ...)
  }
}

# -----------------------------------------------
# Helper: assemble gene x sample matrix from a named list of named vectors
# Replaces the repeated unique(unlist(lapply())) + sparse-assignment loop
# that appeared 5 times (once per method) with a single vectorized merge.
# -----------------------------------------------
# Big O: data.table path: O(total_entries) via rbindlist + dcast (hash-based reshape).
# Fallback path: O(genes × samples) sparse assignment with pre-allocated matrix.
.assemble_tpm_matrix <- function(tpm_list, filter_genes = TRUE) {
  if (length(tpm_list) == 0) return(NULL)

  if (.use_dt) {
    # Build long-form data.table and reshape — avoids O(genes * samples) sparse assignment
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
    if (nrow(mat) == 0) { warning("Matrix assembly produced 0 rows — check input data"); return(NULL) }
  } else {
    # Fallback: original approach
    all_genes <- unique(unlist(lapply(tpm_list, names)))
    if (filter_genes) all_genes <- all_genes[nzchar(all_genes) & !is.na(all_genes)]
    mat <- matrix(0, nrow = length(all_genes), ncol = length(tpm_list),
                  dimnames = list(all_genes, names(tpm_list)))
    for (srr in names(tpm_list)) {
      genes <- intersect(names(tpm_list[[srr]]), all_genes)
      mat[genes, srr] <- tpm_list[[srr]][genes]
    }
    if (nrow(mat) == 0) { warning("Matrix assembly produced 0 rows — check input data"); return(NULL) }
  }
  return(mat)
}

# Helper: fast file read (data.table::fread when available, else read.table)
.fast_read_tsv <- function(path, ...) {
  if (.use_dt) {
    data.table::fread(path, data.table = FALSE, ...)
  } else {
    read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
               check.names = FALSE, comment.char = "", quote = "")
  }
}

cat("\n=== STEP 1: Loading Expression Matrices ===\n\n")

# -----------------------------------------------
# M1: HISAT2 + StringTie (ref-guided)
# -----------------------------------------------
# Per-sample gene_abundances.tsv files with columns:
#   Gene ID | Gene Name | Reference | Strand | Start | End | Coverage | FPKM | TPM

load_m1_tpm <- function() {
  method <- "M1_HISAT2_RefGuided"
  ref_dir <- get_method_ref_dir(method)
  stringtie_base <- file.path(ALIGNMENT_BASE, method, "stringtie_WD", ref_dir)

  if (!dir.exists(stringtie_base)) {
    cat("[M1] StringTie directory not found:", stringtie_base, "\n")
    return(NULL)
  }

  sample_dirs <- list.dirs(stringtie_base, recursive = FALSE, full.names = TRUE)
  # Filter to SRR directories (exclude deseq2_input, etc.)
  sample_dirs <- sample_dirs[grepl("^SRR", basename(sample_dirs))]

  if (length(sample_dirs) == 0) {
    cat("[M1] No sample directories found\n")
    return(NULL)
  }

  # Parallel sample loading via mclapply (I/O-bound: benefits from concurrent reads)
  tpm_list <- setNames(.par_lapply(sample_dirs, function(sdir) {
    abundance_files <- list.files(sdir, pattern = "gene_abundances.*\\.tsv$", full.names = TRUE)
    if (length(abundance_files) == 0) return(NULL)
    df <- tryCatch(.fast_read_tsv(abundance_files[1]), error = function(e) {
      cat("  Warning: failed to read", abundance_files[1], ":", e$message, "\n")
      return(NULL)
    })
    if (is.null(df) || !"TPM" %in% colnames(df)) return(NULL)
    setNames(df$TPM, df[[1]])
  }), basename(sample_dirs))
  tpm_list <- Filter(Negate(is.null), tpm_list)

  tpm_matrix <- .assemble_tpm_matrix(tpm_list, filter_genes = TRUE)
  if (!is.null(tpm_matrix))
    cat("[M1] Loaded:", nrow(tpm_matrix), "genes x", ncol(tpm_matrix), "samples\n")
  return(tpm_matrix)
}

# -----------------------------------------------
# M2: HISAT2 + StringTie (de novo)
# -----------------------------------------------
# De novo assembly uses STRG.N gene IDs. The "Reference" column in the
# abundance file maps to the transcript ID from the reference. We extract
# the gene-level ID by stripping the transcript suffix.

load_m2_tpm <- function() {
  method <- "M2_HISAT2_DeNovo"
  ref_dir <- get_method_ref_dir(method)
  stringtie_base <- file.path(ALIGNMENT_BASE, method, "stringtie_WD", ref_dir)

  if (!dir.exists(stringtie_base)) {
    cat("[M2] StringTie directory not found:", stringtie_base, "\n")
    return(NULL)
  }

  sample_dirs <- list.dirs(stringtie_base, recursive = FALSE, full.names = TRUE)
  sample_dirs <- sample_dirs[grepl("^SRR", basename(sample_dirs))]

  if (length(sample_dirs) == 0) {
    cat("[M2] No sample directories found\n")
    return(NULL)
  }

  # Parallel sample loading via mclapply (matches M1/M4/M5 pattern)
  tpm_list <- setNames(.par_lapply(sample_dirs, function(sdir) {
    srr <- basename(sdir)
    abundance_files <- list.files(sdir, pattern = "gene_abundances.*\\.tsv$", full.names = TRUE)
    if (length(abundance_files) == 0) return(NULL)

    df <- .fast_read_tsv(abundance_files[1])

    if (!all(c("TPM", "Reference") %in% colnames(df))) {
      cat("[M2] Warning: Missing TPM/Reference column in", basename(abundance_files[1]), "for", srr, "\n")
      return(NULL)
    }
    # Map STRG.N -> reference gene ID via the Reference column
    # Single-pass suffix stripping for double-suffixed transcript IDs
    # (e.g., Sme2.5_01g005840.1.01 -> gene-level). O(n) regex vs 2 × O(n).
    gene_ids <- sub("(\\.[0-9]+){1,2}$", "", df$Reference)
    valid <- nzchar(gene_ids) & !gene_ids %in% c(".", "-")
    n_unmapped <- sum(!valid)
    if (n_unmapped > 0) {
      pct_unmapped <- round(100 * n_unmapped / length(valid), 1)
      cat("[M2]", srr, ":", n_unmapped, "of", length(valid), "transcripts unmapped (",
          pct_unmapped, "%)\n")
      if (pct_unmapped > 50)
        cat("[M2] WARNING: >50% unmapped transcripts in", srr, "— check assembly quality\n")
    }
    # Use MAX (not SUM) to aggregate multiple STRG.N entries mapping to the same
    # reference gene. MAX matches matrix_builder.py's visualization aggregation.
    tapply(df$TPM[valid], gene_ids[valid], max, na.rm = TRUE)
  }), basename(sample_dirs))
  tpm_list <- Filter(Negate(is.null), tpm_list)

  tpm_matrix <- .assemble_tpm_matrix(tpm_list, filter_genes = TRUE)
  if (!is.null(tpm_matrix))
    cat("[M2] Loaded:", nrow(tpm_matrix), "genes x", ncol(tpm_matrix), "samples\n")
  return(tpm_matrix)
}

# -----------------------------------------------
# M3: STAR + Salmon
# -----------------------------------------------
# Per-sample quant.sf files from Salmon quantification after STAR alignment.
# quant.sf columns: Name | Length | EffectiveLength | TPM | NumReads
# Transcript-level -> gene-level aggregation via tx2gene mapping.

load_m3_tpm <- function() {
  method <- "M3_STAR_Align"
  ref_dir <- get_method_ref_dir(method)
  quant_base <- file.path(ALIGNMENT_BASE, method, ref_dir, "6_salmon", "quant")

  if (!dir.exists(quant_base)) {
    cat("[M3] Salmon quant directory not found:", quant_base, "\n")
    return(NULL)
  }

  # Load tx2gene mapping
  tx2gene_file <- file.path(POST_PROC_BASE, method, "count_matrices_from_STAR", ref_dir,
                            paste0("tx2gene_", ref_dir, ".tsv"))
  tx2gene <- NULL
  if (file.exists(tx2gene_file)) {
    # Detect header: check if first line contains known column-name keywords.
    # Use word-boundary anchors to avoid matching real transcript IDs that happen
    # to start with "tx" (e.g., "tx_12345" or "txSMEL_001").
    first_line <- readLines(tx2gene_file, n = 1)
    has_header <- grepl("(^|\\t)(transcript_id|tx_id|TXNAME)(\\t|$)", first_line, ignore.case = TRUE)
    tx2gene <- read.table(tx2gene_file, header = has_header, sep = "\t",
                          stringsAsFactors = FALSE)
    # Keep only first 2 columns (transcript_id, gene_id)
    if (ncol(tx2gene) > 2) tx2gene <- tx2gene[, 1:2, drop = FALSE]
    colnames(tx2gene) <- c("transcript_id", "gene_id")
  }

  sample_dirs <- list.dirs(quant_base, recursive = FALSE, full.names = TRUE)
  sample_dirs <- sample_dirs[grepl("^SRR", basename(sample_dirs))]

  # Tissue-specific fallback: quant/{tissue}/{SRR}/quant.sf layout
  if (length(sample_dirs) == 0) {
    cat("[M3] No SRR directories at top level; trying tissue-specific layout...\n")
    tissue_dirs <- list.dirs(quant_base, recursive = FALSE, full.names = TRUE)
    for (td in tissue_dirs) {
      sub_dirs <- list.dirs(td, recursive = FALSE, full.names = TRUE)
      sub_dirs <- sub_dirs[grepl("^SRR", basename(sub_dirs))]
      sample_dirs <- c(sample_dirs, sub_dirs)
    }
    if (length(sample_dirs) > 0) {
      cat("[M3] Found", length(sample_dirs), "samples via tissue-specific layout\n")
    }
  }

  if (length(sample_dirs) == 0) {
    cat("[M3] No sample directories found under", quant_base, "\n")
    # Fallback: try pre-built TPM matrix from tximport (matches M4 fallback pattern)
    tpm_search_base <- file.path(POST_PROC_BASE, method,
                                 "count_matrices_from_STAR", ref_dir)
    tpm_files <- list.files(tpm_search_base, pattern = "_tpm_Gene_ID_.*\\.csv$",
                            recursive = TRUE, full.names = TRUE)
    if (length(tpm_files) > 0) {
      tpm_file <- tpm_files[1]
      cat("[M3] Using pre-built TPM matrix:", tpm_file, "\n")
      df <- if (.use_dt) data.table::fread(tpm_file, data.table = FALSE) else {
        read.csv(tpm_file, header = TRUE, stringsAsFactors = FALSE,
                 check.names = FALSE)
      }
      rownames(df) <- df[[1]]; df <- df[, -1, drop = FALSE]
      return(as.matrix(df))
    }
    return(NULL)
  }

  # Parallel sample loading via mclapply (matches M1/M4/M5 pattern)
  tpm_list <- setNames(.par_lapply(sample_dirs, function(sdir) {
    qsf <- file.path(sdir, "quant.sf")
    if (!file.exists(qsf)) return(NULL)

    df <- .fast_read_tsv(qsf)

    if (!is.null(tx2gene)) {
      df$gene_id <- tx2gene$gene_id[match(df$Name, tx2gene$transcript_id)]
      unmapped <- is.na(df$gene_id)
      if (any(unmapped)) {
        # Single-pass suffix stripping for unmapped transcripts. O(n) vs 2 × O(n).
        df$gene_id[unmapped] <- sub("(\\.[0-9]+){1,2}$", "", df$Name[unmapped])
      }
    } else {
      # Single-pass suffix stripping. O(n) vs 2 × O(n).
      df$gene_id <- sub("(\\.[0-9]+){1,2}$", "", df$Name)
    }

    # rowsum() is a base-R C routine — ~30% faster than tapply for grouped sums. O(n).
    rs <- rowsum(df$TPM, df$gene_id, reorder = FALSE, na.rm = TRUE)
    setNames(rs[, 1], rownames(rs))
  }), basename(sample_dirs))
  tpm_list <- Filter(Negate(is.null), tpm_list)

  tpm_matrix <- .assemble_tpm_matrix(tpm_list, filter_genes = TRUE)
  if (!is.null(tpm_matrix))
    cat("[M3] Loaded:", nrow(tpm_matrix), "genes x", ncol(tpm_matrix), "samples\n")
  return(tpm_matrix)
}

# -----------------------------------------------
# M4: Salmon (pseudo-alignment)
# -----------------------------------------------
# Pre-built genes.counts.matrix exists but that's counts, not TPM.
# Load per-sample quant.sf and aggregate to gene-level TPM.

load_m4_tpm <- function() {
  method <- "M4_Salmon_Saf"
  ref_dir <- get_method_ref_dir(method)
  quant_base <- file.path(ALIGNMENT_BASE, method, "Salmon_Quant", ref_dir)

  if (!dir.exists(quant_base)) {
    cat("[M4] Salmon quant directory not found:", quant_base, "\n")
    return(NULL)
  }

  sample_dirs <- list.dirs(quant_base, recursive = FALSE, full.names = TRUE)
  sample_dirs <- sample_dirs[grepl("^SRR", basename(sample_dirs))]

  if (length(sample_dirs) == 0) {
    cat("[M4] No sample directories found\n")
    return(NULL)
  }

  # Parallel sample loading via mclapply (replaces sequential for-loop)
  tpm_list <- setNames(.par_lapply(sample_dirs, function(sdir) {
    qsf <- file.path(sdir, "quant.sf")
    if (!file.exists(qsf)) return(NULL)
    df <- .fast_read_tsv(qsf)
    # Single-pass suffix stripping for double-suffixed transcript IDs. O(n) vs 2 × O(n).
    df$gene_id <- sub("(\\.[0-9]+){1,2}$", "", df$Name)
    # rowsum() is a base-R C routine — ~30% faster than tapply for grouped sums. O(n).
    rs <- rowsum(df$TPM, df$gene_id, reorder = FALSE, na.rm = TRUE)
    setNames(rs[, 1], rownames(rs))
  }), basename(sample_dirs))
  tpm_list <- Filter(Negate(is.null), tpm_list)

  if (length(tpm_list) == 0) {
    # Fallback: try pre-built TPM matrix from tximport
    # tximport_salmon_to_matrices.R saves as {prefix}_tpm_Gene_ID_from_{ref}_gene_level.csv
    tpm_search_base <- file.path(POST_PROC_BASE, method,
                                 "count_matrices_from_Salmon_Quant", ref_dir)
    tpm_files <- list.files(tpm_search_base, pattern = "_tpm_Gene_ID_.*\\.csv$",
                            recursive = TRUE, full.names = TRUE)
    if (length(tpm_files) > 0) {
      tpm_file <- tpm_files[1]
      cat("[M4] Using pre-built TPM matrix:", tpm_file, "\n")
      df <- if (.use_dt) data.table::fread(tpm_file, data.table = FALSE) else {
        read.csv(tpm_file, header = TRUE, stringsAsFactors = FALSE,
                 check.names = FALSE)
      }
      rownames(df) <- df[[1]]; df <- df[, -1, drop = FALSE]
      return(as.matrix(df))
    }
    return(NULL)
  }

  tpm_matrix <- .assemble_tpm_matrix(tpm_list, filter_genes = TRUE)
  if (!is.null(tpm_matrix))
    cat("[M4] Loaded:", nrow(tpm_matrix), "genes x", ncol(tpm_matrix), "samples\n")
  return(tpm_matrix)
}

# -----------------------------------------------
# M5: Bowtie2 + RSEM
# -----------------------------------------------
# Primary: per-sample .genes.results from RSEM alignment output
# Fallback: pre-built TPM matrix from tximport (_tpm_Gene_ID_*.csv)

load_m5_tpm <- function() {
  method <- "M5_RSEM_Bowtie2"
  ref_dir <- get_method_ref_dir(method)

  # Primary: load per-sample .genes.results from RSEM alignment output
  rsem_quant_base <- file.path(ALIGNMENT_BASE, method, "RSEM_Quant_WD", ref_dir)

  if (dir.exists(rsem_quant_base)) {
    sample_dirs <- list.dirs(rsem_quant_base, recursive = FALSE, full.names = TRUE)
    sample_dirs <- sample_dirs[grepl("^SRR", basename(sample_dirs))]

    if (length(sample_dirs) > 0) {
      # Parallel sample loading via mclapply (replaces sequential for-loop)
      tpm_list <- setNames(.par_lapply(sample_dirs, function(sdir) {
        srr <- basename(sdir)
        results_file <- file.path(sdir, paste0(srr, ".genes.results"))
        if (!file.exists(results_file)) return(NULL)
        df <- .fast_read_tsv(results_file)
        # Single-pass suffix stripping. O(n) vs 2 × O(n).
        gene_ids <- sub("(\\.[0-9]+){1,2}$", "", df$gene_id)
        # rowsum() is a base-R C routine — ~30% faster than tapply for grouped sums. O(n).
        rs <- rowsum(df$TPM, gene_ids, reorder = FALSE, na.rm = TRUE)
        setNames(rs[, 1], rownames(rs))
      }), basename(sample_dirs))
      tpm_list <- Filter(Negate(is.null), tpm_list)

      if (length(tpm_list) > 0) {
        tpm_matrix <- .assemble_tpm_matrix(tpm_list, filter_genes = TRUE)
        if (!is.null(tpm_matrix)) {
          cat("[M5] Loaded:", nrow(tpm_matrix), "genes x", ncol(tpm_matrix),
              "samples (per-sample .genes.results)\n")
          return(tpm_matrix)
        }
      }
    }
  }

  # Fallback: try pre-built TPM matrix from tximport
  # tximport_rsem_to_matrices.R saves as {prefix}_tpm_Gene_ID_from_{ref}_gene_level.csv
  # under count_matrices_from_RSEM_Quant/{ref}/gene_level/
  tpm_search_base <- file.path(POST_PROC_BASE, method,
                               "count_matrices_from_RSEM_Quant", ref_dir)
  tpm_files <- list.files(tpm_search_base, pattern = "_tpm_Gene_ID_.*\\.csv$",
                          recursive = TRUE, full.names = TRUE)
  if (length(tpm_files) > 0) {
    tpm_file <- tpm_files[1]
    cat("[M5] Using pre-built TPM matrix:", tpm_file, "\n")
    df <- if (.use_dt) data.table::fread(tpm_file, data.table = FALSE) else {
      read.csv(tpm_file, header = TRUE, stringsAsFactors = FALSE,
               check.names = FALSE)
    }
    rownames(df) <- df[[1]]; df <- df[, -1, drop = FALSE]
    return(as.matrix(df))
  }

  cat("[M5] TPM matrix not found\n")
  return(NULL)
}

# -----------------------------------------------
# Load all methods
# -----------------------------------------------

loader_map <- list(
  M1_HISAT2_RefGuided = load_m1_tpm,
  M2_HISAT2_DeNovo    = load_m2_tpm,
  M3_STAR_Align       = load_m3_tpm,
  M4_Salmon_Saf       = load_m4_tpm,
  M5_RSEM_Bowtie2     = load_m5_tpm
)

tpm_matrices <- list()
# Pre-allocate stats list to avoid O(n^2) list growth
stats_list <- vector("list", length(CONCORDANCE_METHODS))
stats_idx <- 0L

for (method in CONCORDANCE_METHODS) {
  cat("\nLoading", method, "...\n")
  loader_fn <- loader_map[[method]]
  if (is.null(loader_fn)) {
    cat("  [WARN] No loader for method:", method, "\n")
    next
  }

  mat <- tryCatch(loader_fn(), error = function(e) {
    cat("  [ERROR]", e$message, "\n")
    NULL
  })

  if (!is.null(mat) && nrow(mat) > 0 && ncol(mat) > 0) {
    tpm_matrices[[method]] <- mat
    stats_idx <- stats_idx + 1L
    stats_list[[stats_idx]] <- data.frame(
      method = method,
      short_name = get_short_name(method),
      n_genes_raw = nrow(mat),
      n_samples_raw = ncol(mat),
      stringsAsFactors = FALSE
    )
  } else {
    cat("  [WARN] No data loaded for", method, "\n")
  }
}
stats_list <- stats_list[seq_len(stats_idx)]
method_stats <- if (.use_dt) {
  data.table::setDF(data.table::rbindlist(stats_list, use.names = TRUE, fill = TRUE))
} else {
  do.call(rbind, stats_list)
}

if (length(tpm_matrices) < 2) {
  stop("Need at least 2 methods with data for concordance analysis. Found: ",
       length(tpm_matrices))
}

cat("\n--- Loaded", length(tpm_matrices), "methods ---\n")

# -----------------------------------------------
# Harmonize gene IDs across methods
# -----------------------------------------------
# All methods for GPE001970 use SMEL5_XXgXXXXXX gene IDs (after suffix stripping).
# M1 uses the same IDs natively. M2 maps STRG.N -> Reference transcript -> gene.
# We strip transcript suffixes (.N) to get gene-level IDs for all methods.

cat("\n--- Harmonizing gene IDs ---\n")

# Ensure all gene IDs are at gene level (strip .N suffix if present)
for (method in names(tpm_matrices)) {
  mat <- tpm_matrices[[method]]
  rn <- rownames(mat)
  # Only strip if IDs look like they have transcript suffixes
  # Handles both GPE001970 (SMEL5_XXgXXXXXX.N) and Eggplant_V4.1 (Sme2.5_XXgXXXXXX.N)
  if (any(grepl("^(SMEL|Sme)[0-9].*\\.[0-9]+$", rn))) {
    # Single-pass suffix stripping for double-suffixed IDs
    # (e.g., SMEL4.1_06g023900.1.01 -> SMEL4.1_06g023900). O(n) vs 2 × O(n).
    new_rn <- sub("(\\.[0-9]+){1,2}$", "", rn)
    if (any(duplicated(new_rn))) {
      # rowsum() is a base-R C routine optimized for grouped column sums on matrices
      # — faster than data.table for this use case (avoids matrix→DT→matrix round-trip)
      mat_agg <- rowsum(mat, group = new_rn, reorder = FALSE)
      tpm_matrices[[method]] <- mat_agg
    } else {
      rownames(mat) <- new_rn
      tpm_matrices[[method]] <- mat
    }
  }
}

# Find common genes across all methods — O(sum(|G_i|)) single-pass table count
# replaces Reduce(intersect,...) which allocates M-1 temporary vectors
gene_sets <- lapply(tpm_matrices, rownames)
n_methods <- length(gene_sets)
if (n_methods == 1L) {
  common_genes <- gene_sets[[1L]]
} else {
  .gene_counts <- table(unlist(gene_sets, use.names = FALSE))
  common_genes <- names(.gene_counts[.gene_counts == n_methods])
  rm(.gene_counts)
}
cat("  Common genes across all methods:", length(common_genes), "\n")

if (length(common_genes) == 0) {
  cat("\n  [ERROR] Zero common genes across methods. Per-method gene counts:\n")
  for (m in names(gene_sets)) {
    cat("    ", get_short_name(m), ":", length(gene_sets[[m]]), "genes",
        "(sample IDs:", paste(head(gene_sets[[m]], 3), collapse = ", "), "...)\n")
  }
  stop("No common genes found across methods. Check gene ID formats (suffix stripping may be too aggressive).")
}

# Pairwise gene overlaps for reporting
# Pre-compute short names — avoids O(M²) calls to get_short_name()
n_mat <- length(tpm_matrices)
mat_names <- names(tpm_matrices)
short_names_map <- setNames(vapply(mat_names, get_short_name, character(1)), mat_names)
for (i in seq_len(n_mat - 1)) {
  for (j in seq(i + 1, n_mat)) {
    overlap <- length(intersect(gene_sets[[mat_names[i]]], gene_sets[[mat_names[j]]]))
    cat("  ", short_names_map[mat_names[i]], " & ", short_names_map[mat_names[j]], ":", overlap, "shared genes\n")
  }
}

# -----------------------------------------------
# Harmonize sample IDs across methods
# -----------------------------------------------

cat("\n--- Harmonizing sample IDs ---\n")

# Find common samples across all methods — O(sum(|S_i|)) single-pass table count
# (matches gene intersection approach above; avoids M-1 intermediate vectors from Reduce)
sample_sets <- lapply(tpm_matrices, colnames)
if (length(sample_sets) == 1L) {
  common_samples <- sample_sets[[1L]]
} else {
  .sample_counts <- table(unlist(sample_sets, use.names = FALSE))
  common_samples <- names(.sample_counts[.sample_counts == length(sample_sets)])
  rm(.sample_counts)
}
cat("  Common samples across all methods:", length(common_samples), "\n")

if (length(common_samples) < CONCORDANCE_MIN_SAMPLES) {
  stop("Too few common samples (", length(common_samples), "). Need at least ",
       CONCORDANCE_MIN_SAMPLES)
}

# -----------------------------------------------
# Subset to common genes and samples
# -----------------------------------------------

cat("\n--- Subsetting to common features ---\n")

for (method in names(tpm_matrices)) {
  tpm_matrices[[method]] <- tpm_matrices[[method]][common_genes, common_samples, drop = FALSE]
  cat("  ", get_short_name(method), ":", nrow(tpm_matrices[[method]]), "x",
      ncol(tpm_matrices[[method]]), "\n")
}

# Update stats with harmonized counts — vectorized via precomputed maps.
# O(M) vapply replaces O(M) sapply with per-element function dispatch.
.nrow_map <- vapply(tpm_matrices, nrow, integer(1))
.ncol_map <- vapply(tpm_matrices, ncol, integer(1))
method_stats$n_genes_harmonized <- .nrow_map[method_stats$method]
method_stats$n_samples_harmonized <- .ncol_map[method_stats$method]

# -----------------------------------------------
# Filter lowly-expressed genes
# -----------------------------------------------

cat("\n--- Filtering lowly-expressed genes ---\n")
cat("  Criteria: TPM >=", CONCORDANCE_MIN_EXPR, "in at least", CONCORDANCE_MIN_SAMPLES, "samples\n")
cat("  (Applied per method; keeping genes that pass in ANY method)\n")

expressed_genes_per_method <- lapply(tpm_matrices, function(mat) {
  rownames(mat)[rowSums(mat >= CONCORDANCE_MIN_EXPR) >= CONCORDANCE_MIN_SAMPLES]
})

# Keep genes expressed in at least one method (union)
expressed_union <- unique(unlist(expressed_genes_per_method))
# Also track genes expressed in ALL methods (intersection) for stricter analysis
expressed_intersect <- Reduce(intersect, expressed_genes_per_method)

cat("  Genes expressed in ANY method:", length(expressed_union), "\n")
cat("  Genes expressed in ALL methods:", length(expressed_intersect), "\n")

# Use union for general concordance (more inclusive)
filtered_genes <- intersect(common_genes, expressed_union)
for (method in names(tpm_matrices)) {
  tpm_matrices[[method]] <- tpm_matrices[[method]][filtered_genes, , drop = FALSE]
}

cat("  Final gene count:", length(filtered_genes), "\n")

# Free intermediate objects before saving (prevents memory bloat during concordance analysis)
# NOTE: gene_sets and sample_sets are kept — they are referenced in the result list below.
rm(expressed_genes_per_method, expressed_union)
gc(verbose = FALSE)

# -----------------------------------------------
# Save harmonized data
# -----------------------------------------------

result <- list(
  tpm_matrices = tpm_matrices,
  common_genes = filtered_genes,
  common_samples = common_samples,
  expressed_in_all = expressed_intersect,
  method_stats = method_stats,
  gene_sets_raw = gene_sets,
  sample_sets_raw = sample_sets
)

saveRDS(result, HARMONIZED_RDS)
cat("\n[DONE] Harmonized data saved to:", HARMONIZED_RDS, "\n")

# Also save method stats table
write.csv(method_stats, file.path(TABLES_DIR, "method_stats.csv"), row.names = FALSE)
