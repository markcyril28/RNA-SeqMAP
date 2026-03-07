#!/usr/bin/env Rscript
# ==============================================================================
# TXIMPORT STAR HELPER - STAR+Salmon quantification import for DESeq2
# ==============================================================================
# Usage: Rscript tximport_star_helper.R <quant_dir> <metadata_file> <tx2gene_file> [output_dir]
# ==============================================================================

suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript tximport_star_helper.R <quant_dir> <metadata_file> <tx2gene_file> [output_dir]")
}

quant_dir <- args[1]
metadata_file <- args[2]
tx2gene_file <- args[3]
output_dir <- if (length(args) > 3) args[4] else dirname(metadata_file)

cat("STAR + Salmon tximport for DESeq2\n")
cat("Quant dir:", quant_dir, "\n")
cat("Metadata:", metadata_file, "\n")
cat("TX2Gene:", tx2gene_file, "\n")

# Load metadata
coldata <- read.delim(metadata_file, header = TRUE, stringsAsFactors = FALSE)
samples <- coldata$sample

# Find quantification files
files <- file.path(quant_dir, samples, "quant.sf")
names(files) <- samples

if (!all(file.exists(files))) {
  missing <- files[!file.exists(files)]
  stop("Missing quant.sf files:\n", paste(missing, collapse = "\n"))
}

# Load tx2gene mapping
tx2gene <- read.delim(tx2gene_file, header = FALSE, col.names = c("TXNAME", "GENEID"))

# Import with tximport
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)

# Create DESeq2 object
dds <- DESeqDataSetFromTximport(txi, colData = coldata, design = ~condition)

# Save outputs
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(txi, file.path(output_dir, "tximport_star_salmon.rds"))
saveRDS(dds, file.path(output_dir, "deseq2_dataset_star.rds"))

# Export count matrices
write.table(txi$counts, file.path(output_dir, "gene_counts_tximport.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)
write.table(txi$abundance, file.path(output_dir, "gene_tpm_tximport.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)

cat("\nTximport completed successfully!\n")
cat("Outputs saved to:", output_dir, "\n")
