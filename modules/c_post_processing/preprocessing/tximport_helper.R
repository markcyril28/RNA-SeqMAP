#!/usr/bin/env Rscript
# ==============================================================================
# TXIMPORT HELPER - Salmon/RSEM quantification import for DESeq2
# ==============================================================================
# Usage: Rscript tximport_helper.R <method> <quant_dir> <metadata_file> [output_dir]
# Methods: salmon, rsem
# ==============================================================================

suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript tximport_helper.R <method> <quant_dir> <metadata_file> [output_dir]")
}

method <- args[1]
quant_dir <- args[2]
metadata_file <- args[3]
output_dir <- if (length(args) > 3) args[4] else dirname(metadata_file)

cat(method, "tximport for DESeq2\n")
cat("Quant dir:", quant_dir, "\n")
cat("Metadata:", metadata_file, "\n")

# Load metadata
coldata <- read.delim(metadata_file, header = TRUE, stringsAsFactors = FALSE, sep = ",")
if (ncol(coldata) == 1) {
  coldata <- read.delim(metadata_file, header = TRUE, stringsAsFactors = FALSE)
}
# Accept "sample", "SampleID", or "sampleID" as the sample column name
sample_col <- intersect(c("sample", "SampleID", "sampleID"), colnames(coldata))[1]
if (is.na(sample_col)) {
  stop("metadata file must have a 'sample', 'SampleID', or 'sampleID' column. Found: ",
       paste(colnames(coldata), collapse = ", "))
}
samples <- coldata[[sample_col]]

# Find quantification files
if (tolower(method) == "salmon") {
  files <- file.path(quant_dir, samples, "quant.sf")
} else if (tolower(method) == "rsem") {
  files <- file.path(quant_dir, samples, paste0(samples, ".genes.results"))
} else {
  stop("Unknown method: ", method, ". Use 'salmon' or 'rsem'.")
}
names(files) <- samples

if (!all(file.exists(files))) {
  missing <- files[!file.exists(files)]
  stop("Missing quantification files:\n", paste(missing, collapse = "\n"))
}

# Import with tximport
if (tolower(method) == "salmon") {
  # Check for tx2gene mapping file (enables gene-level aggregation)
  # Look for GENE_TRANS_MAP_FILE env var, or auto-discover near quant_dir
  tx2gene_file <- Sys.getenv("GENE_TRANS_MAP_FILE", unset = "")
  if (!nzchar(tx2gene_file)) {
    parent_dir <- dirname(quant_dir)
    candidates <- list.files(parent_dir, pattern = "\\.gene_trans_map$",
                             recursive = TRUE, full.names = TRUE)
    if (length(candidates) > 0) tx2gene_file <- candidates[1]
  }
  if (nzchar(tx2gene_file) && file.exists(tx2gene_file)) {
    cat("Using tx2gene mapping for gene-level aggregation:", tx2gene_file, "\n")
    tx2gene <- read.table(tx2gene_file, header = FALSE, sep = "\t",
                          stringsAsFactors = FALSE, colClasses = c("character", "character"),
                          strip.white = TRUE)
    colnames(tx2gene) <- c("GENEID", "TXNAME")
    tx2gene$GENEID <- trimws(tx2gene$GENEID)
    tx2gene$TXNAME <- trimws(tx2gene$TXNAME)
    tx2gene <- tx2gene[, c("TXNAME", "GENEID")]
    txi <- tximport(files, type = "salmon", txIn = TRUE, txOut = FALSE,
                    tx2gene = tx2gene, ignoreTxVersion = FALSE, ignoreAfterBar = FALSE)
  } else {
    cat("No tx2gene mapping found — producing transcript-level output\n")
    cat("Set GENE_TRANS_MAP_FILE env var for gene-level aggregation\n")
    txi <- tximport(files, type = "salmon", txOut = TRUE)
  }
} else if (tolower(method) == "rsem") {
  txi <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)
}

# Accept "condition" or "Condition" as the design column
cond_col <- intersect(c("condition", "Condition"), colnames(coldata))[1]
if (is.na(cond_col)) {
  stop("metadata file must have a 'condition' or 'Condition' column. Found: ",
       paste(colnames(coldata), collapse = ", "))
}
# Standardise column name to 'condition' for DESeq2 formula
coldata$condition <- coldata[[cond_col]]

# Create DESeq2 object
dds <- DESeqDataSetFromTximport(txi, colData = coldata, design = ~condition)

# Save outputs
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(txi, file.path(output_dir, paste0("tximport_", tolower(method), ".rds")))
saveRDS(dds, file.path(output_dir, paste0("deseq2_dataset_", tolower(method), ".rds")))

# Export count matrices with explicit GeneID column (consistent with other preprocessing scripts)
# Salmon with tx2gene produces gene-level data; Salmon without produces transcript-level
# RSEM always produces gene-level data
is_gene_level <- tolower(method) == "rsem" || exists("tx2gene", inherits = FALSE)
id_label <- if (is_gene_level) "GeneID" else "TranscriptID"
level <- if (is_gene_level) "gene" else "transcript"

save_tximport_tsv <- function(mat, filename) {
  df <- as.data.frame(mat, check.names = FALSE)
  df <- cbind(setNames(data.frame(rownames(mat), stringsAsFactors = FALSE), id_label), df)
  write.table(df, file.path(output_dir, filename),
              sep = "\t", quote = FALSE, row.names = FALSE)
}
save_tximport_tsv(txi$counts,    paste0(level, "_counts_tximport.tsv"))
save_tximport_tsv(txi$abundance, paste0(level, "_tpm_tximport.tsv"))

cat("\nTximport completed successfully!\n")
cat("Outputs saved to:", output_dir, "\n")
