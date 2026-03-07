#!/usr/bin/env Rscript
# ==============================================================================
# tximport Script for Proper Import of Transcript Quantifications into DESeq2
# ==============================================================================
# This script properly handles transcript-level uncertainty when importing
# RSEM or Salmon quantifications into DESeq2, which is the recommended approach
# over using expected_count values directly.
#
# WHY USE TXIMPORT?
# - RSEM produces 'expected counts' that are probabilistically estimated
# - Direct use in DESeq2 can violate negative binomial assumptions
# - tximport properly propagates quantification uncertainty to DESeq2
# - This is the gold standard for publication-quality RNA-seq analysis
#
# REFERENCE:
# Soneson et al. (2015) "Differential analyses for RNA-seq: transcript-level
# estimates improve gene-level inferences" F1000Research
# ==============================================================================

# Load required libraries
suppressPackageStartupMessages({
    if (!requireNamespace("tximport", quietly = TRUE)) {
        stop("\n❌ Package 'tximport' is required.\n",
             "   Install with: BiocManager::install('tximport')\n",
             "   Or: install.packages('BiocManager'); BiocManager::install('tximport')")
    }
    if (!requireNamespace("DESeq2", quietly = TRUE)) {
        stop("\n❌ Package 'DESeq2' is required.\n",
             "   Install with: BiocManager::install('DESeq2')")
    }

    library(tximport)
    library(DESeq2)
})

# Configuration
METHOD <- "METHOD_PLACEHOLDER"  # Will be replaced: "rsem" or "salmon"
QUANT_DIR <- "QUANT_DIR_PLACEHOLDER"
METADATA_FILE <- "METADATA_FILE_PLACEHOLDER"
OUTPUT_DIR <- dirname(METADATA_FILE)

cat("\n")
cat("==============================================================================\n")
cat("  tximport Pipeline for", toupper(METHOD), "Quantifications\n")
cat("  Publication-Ready RNA-seq Quantification Import\n")
cat("==============================================================================\n\n")
cat("📖 This script imports RSEM quantifications using tximport, which:\n")
cat("   • Properly handles transcript-level uncertainty\n")
cat("   • Preserves statistical validity for DESeq2 analysis\n")
cat("   • Is the recommended approach for publication\n\n")

# Load sample metadata
cat("Loading sample metadata from:", METADATA_FILE, "\n")
metadata <- read.csv(METADATA_FILE, stringsAsFactors = FALSE)
cat("  Samples:", nrow(metadata), "\n")
cat("  Conditions:", paste(unique(metadata$condition), collapse=", "), "\n\n")

# Find quantification files
cat("Locating quantification files in:", QUANT_DIR, "\n")
if (METHOD == "rsem") {
    files <- list.files(QUANT_DIR, pattern = "\\.genes\\.results$", 
                       full.names = TRUE, recursive = TRUE)
    names(files) <- gsub("\\.genes\\.results$", "", basename(files))
} else if (METHOD == "salmon") {
    # Salmon: look for quant.sf files
    sample_dirs <- list.dirs(QUANT_DIR, recursive = FALSE)
    files <- file.path(sample_dirs, "quant.sf")
    names(files) <- basename(sample_dirs)
    files <- files[file.exists(files)]
} else {
    stop("Unknown method: ", METHOD)
}

cat("  Found", length(files), "quantification files\n\n")

# Verify all samples in metadata have quantification files
missing <- setdiff(metadata$sample, names(files))
if (length(missing) > 0) {
    warning("Samples in metadata missing quantification files: ", 
           paste(missing, collapse=", "))
}

# Import with tximport
cat("Importing quantifications with tximport...\n")
if (METHOD == "rsem") {
    txi <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)
} else if (METHOD == "salmon") {
    # For Salmon, may need tx2gene mapping for gene-level summarization
    # If using transcript-level index, provide tx2gene
    # Otherwise, use txOut=FALSE for gene-level quantifications
    txi <- tximport(files, type = "salmon", txIn = FALSE, txOut = FALSE)
}

cat("  Imported", nrow(txi$counts), "genes\n")
cat("  Across", ncol(txi$counts), "samples\n\n")

# Create DESeq2 dataset
cat("Creating DESeq2 dataset...\n")
# Ensure metadata sample order matches count matrix columns
metadata <- metadata[match(colnames(txi$counts), metadata$sample), ]

# Check if we have proper experimental design
if (length(unique(metadata$condition)) < 2) {
    stop("ERROR: Need at least 2 different conditions for differential expression analysis!\n",
         "       Current conditions: ", paste(unique(metadata$condition), collapse=", "), "\n",
         "       Please edit the metadata file to assign proper experimental conditions.")
}

# Create DESeq2 object with proper design
dds <- DESeqDataSetFromTximport(txi, 
                                colData = metadata, 
                                design = ~ condition)

cat("  DESeq2 dataset created successfully\n")
cat("  Design formula: ~ condition\n\n")

# Save the DESeq2 object for downstream analysis
output_rds <- file.path(OUTPUT_DIR, paste0(METHOD, "_deseq2_tximport.rds"))
saveRDS(dds, output_rds)
cat("Saved DESeq2 object to:", output_rds, "\n\n")

# Save normalized counts for inspection
cat("Saving normalized count tables...\n")
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized = TRUE)
write.csv(normalized_counts, 
          file.path(OUTPUT_DIR, paste0(METHOD, "_normalized_counts_tximport.csv")),
          row.names = TRUE)

# Basic quality control plots
cat("Generating QC plots...\n")
pdf(file.path(OUTPUT_DIR, paste0(METHOD, "_tximport_qc.pdf")), width = 10, height = 8)

# Sample correlation heatmap
library(pheatmap)
sample_cor <- cor(normalized_counts)
pheatmap(sample_cor, 
         annotation_col = data.frame(Condition = metadata$condition, row.names = metadata$sample),
         main = "Sample Correlation (tximport normalized counts)")

# PCA plot
vsd <- vst(dds, blind = FALSE)
pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
library(ggplot2)
ggplot(pcaData, aes(PC1, PC2, color = condition, label = name)) +
    geom_point(size = 3) +
    geom_text(vjust = -0.5, hjust = 0.5, size = 3) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    ggtitle("PCA Plot (tximport data)") +
    theme_minimal()

dev.off()

cat("\n")
cat("==============================================================================\n")
cat("✅ tximport Pipeline Completed Successfully!\n")
cat("==============================================================================\n\n")
cat("📁 OUTPUT FILES:\n")
cat("  1. DESeq2 object (RDS):", basename(output_rds), "\n")
cat("  2. Normalized counts (CSV):", basename(file.path(OUTPUT_DIR, paste0(METHOD, "_normalized_counts_tximport.csv"))), "\n")
cat("  3. QC plots (PDF):", basename(file.path(OUTPUT_DIR, paste0(METHOD, "_tximport_qc.pdf"))), "\n\n")
cat("📂 All files saved to:", OUTPUT_DIR, "\n\n")
cat("==============================================================================\n")
cat("🔬 NEXT STEPS: Differential Expression Analysis\n")
cat("==============================================================================\n\n")
cat("1️⃣  Load the DESeq2 object:\n")
cat("   dds <- readRDS('", output_rds, "')\n\n", sep="")
cat("2️⃣  Run differential expression:\n")
cat("   dds <- DESeq(dds)\n\n")
cat("3️⃣  Extract results (adjust conditions as needed):\n")
cat("   # List available conditions:\n")
cat("   levels(dds$condition)\n\n")
cat("   # Example: compare 'treatment' vs 'control'\n")
cat("   res <- results(dds, contrast=c('condition', 'treatment', 'control'))\n\n")
cat("4️⃣  Visualize results:\n")
cat("   # MA plot\n")
cat("   plotMA(res)\n\n")
cat("   # Significant genes\n")
cat("   summary(res)\n")
cat("   sig_genes <- subset(res, padj < 0.05)\n\n")
cat("==============================================================================\n")
cat("📚 CITATION:\n")
cat("   Soneson et al. (2015) 'Differential analyses for RNA-seq:\n")
cat("   transcript-level estimates improve gene-level inferences'\n")
cat("   F1000Research. doi: 10.12688/f1000research.7563.1\n")
cat("==============================================================================\n\n")