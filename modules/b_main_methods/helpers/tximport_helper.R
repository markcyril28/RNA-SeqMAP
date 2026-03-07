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
  library(readr)
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
samples <- coldata$sample

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
  txi <- tximport(files, type = "salmon", txOut = TRUE)
} else if (tolower(method) == "rsem") {
  txi <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)
}

# Create DESeq2 object
dds <- DESeqDataSetFromTximport(txi, colData = coldata, design = ~condition)

# Save outputs
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(txi, file.path(output_dir, paste0("tximport_", tolower(method), ".rds")))
saveRDS(dds, file.path(output_dir, paste0("deseq2_dataset_", tolower(method), ".rds")))

# Export count matrices
write.table(txi$counts, file.path(output_dir, "gene_counts_tximport.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)
write.table(txi$abundance, file.path(output_dir, "gene_tpm_tximport.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)

cat("\nTximport completed successfully!\n")
cat("Outputs saved to:", output_dir, "\n")
