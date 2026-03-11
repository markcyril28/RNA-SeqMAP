#!/usr/bin/env Rscript
# ==============================================================================
# TXIMPORT STAR HELPER - STAR+Salmon quantification import for DESeq2
# ==============================================================================
# Usage: Rscript tximport_star_helper.R <quant_dir> <metadata_file> <tx2gene_file> <output_dir> [master_ref]
#
# Outputs (all in <output_dir>):
#   gene_level/
#     {master_ref}_NumReads_Gene_ID_from_{master_ref}_gene_level.tsv
#     {master_ref}_tpm_Gene_ID_from_{master_ref}_gene_level.tsv
#   tximport_star_salmon.rds   (txi object - for downstream R scripts)
#   deseq2_dataset_star.rds    (DESeqDataSet - ready for DESeq())
# ==============================================================================

suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript tximport_star_helper.R <quant_dir> <metadata_file> <tx2gene_file> <output_dir> [master_ref]")
}

quant_dir     <- args[1]
metadata_file <- args[2]
tx2gene_file  <- args[3]
output_dir    <- args[4]
# master_ref defaults to basename(output_dir) when not supplied explicitly.
# This matches the convention STAR_MATRIX_ROOT = .../count_matrices_from_STAR/{fasta_tag}/
master_ref    <- if (length(args) >= 5 && nzchar(args[5])) args[5] else basename(output_dir)

cat("STAR + Salmon tximport for DESeq2\n")
cat("Quant dir:   ", quant_dir, "\n")
cat("Metadata:    ", metadata_file, "\n")
cat("TX2Gene:     ", tx2gene_file, "\n")
cat("Output dir:  ", output_dir, "\n")
cat("Master ref:  ", master_ref, "\n\n")

# ---------------------------------------------------------------------------
# Load metadata
# ---------------------------------------------------------------------------
coldata <- read.delim(metadata_file, header = TRUE, stringsAsFactors = FALSE)

# Accept "sample" or "SampleID" as the sample column name
sample_col <- intersect(c("sample", "SampleID", "sampleID"), colnames(coldata))[1]
if (is.na(sample_col)) {
  stop("metadata file must have a 'sample' or 'SampleID' column. Found: ",
       paste(colnames(coldata), collapse = ", "))
}
samples <- coldata[[sample_col]]
rownames(coldata) <- coldata[[sample_col]]

# Accept "condition" or "Condition" as the design column
cond_col <- intersect(c("condition", "Condition"), colnames(coldata))[1]
if (is.na(cond_col)) {
  stop("metadata file must have a 'condition' or 'Condition' column. Found: ",
       paste(colnames(coldata), collapse = ", "))
}
# Standardise column name to 'condition' for DESeq2 formula
coldata$condition <- coldata[[cond_col]]

# ---------------------------------------------------------------------------
# Locate quant.sf files
# ---------------------------------------------------------------------------
files <- file.path(quant_dir, samples, "quant.sf")
names(files) <- samples

missing <- files[!file.exists(files)]
if (length(missing) > 0) {
  stop("Missing quant.sf files:\n", paste(missing, collapse = "\n"))
}
cat("Found", length(files), "quant.sf files\n")

# ---------------------------------------------------------------------------
# Load tx2gene mapping (transcript_id -> gene_id, created from GTF)
# ---------------------------------------------------------------------------
tx2gene <- read.delim(tx2gene_file, header = FALSE, stringsAsFactors = FALSE,
                      col.names = c("TXNAME", "GENEID"))
tx2gene$TXNAME <- trimws(tx2gene$TXNAME)
tx2gene$GENEID <- trimws(tx2gene$GENEID)
cat("Loaded tx2gene:", nrow(tx2gene), "entries\n\n")

# ---------------------------------------------------------------------------
# Import with tximport (gene-level via tx2gene)
# ---------------------------------------------------------------------------
cat("Running tximport...\n")
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
cat("Imported:", ncol(txi$counts), "samples,", nrow(txi$counts), "genes\n\n")

# ---------------------------------------------------------------------------
# Create DESeq2 object
# ---------------------------------------------------------------------------
dds <- DESeqDataSetFromTximport(txi, colData = coldata, design = ~ condition)

# ---------------------------------------------------------------------------
# Save RDS outputs (for direct DESeq2 use by downstream scripts)
# ---------------------------------------------------------------------------
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(txi, file.path(output_dir, "tximport_star_salmon.rds"))
saveRDS(dds, file.path(output_dir, "deseq2_dataset_star.rds"))
cat("Saved RDS objects to:", output_dir, "\n")

# ---------------------------------------------------------------------------
# Save TSV matrices following the pipeline naming convention:
#   {master_ref}_{count_type}_{gene_type}_from_{master_ref}_{level}.tsv
# GeneID is stored as a proper first column (not row names with col.names=NA)
# ---------------------------------------------------------------------------
save_matrix_tsv <- function(mat, count_type, gene_type, out_dir, mr, level) {
  fname <- paste0(mr, "_", count_type, "_", gene_type, "_from_", mr, "_", level, ".tsv")
  fpath <- file.path(out_dir, fname)
  df <- as.data.frame(mat, check.names = FALSE)
  df <- cbind(GeneID = rownames(mat), df)
  rownames(df) <- NULL
  write.table(df, fpath, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("Saved:", fname, "\n")
}

gene_level_dir <- file.path(output_dir, "gene_level")
dir.create(gene_level_dir, showWarnings = FALSE, recursive = TRUE)

save_matrix_tsv(txi$counts,    "NumReads", "Gene_ID", gene_level_dir, master_ref, "gene_level")
save_matrix_tsv(txi$abundance, "tpm",      "Gene_ID", gene_level_dir, master_ref, "gene_level")

cat("\nTximport completed successfully!\n")
cat("TSV matrices saved to:", gene_level_dir, "\n")
