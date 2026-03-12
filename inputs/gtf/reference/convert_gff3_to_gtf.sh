#!/bin/bash
# ==============================================================================
# GFF3 TO GTF CONVERTER FOR HISAT2 REFERENCE-GUIDED ALIGNMENT
# ==============================================================================
# Converts GPE001970.gff (GFF3 format) to GTF format compatible with:
#   - hisat2_extract_splice_sites.py
#   - hisat2_extract_exons.py
#   - StringTie
#
# Key transformations:
#   1. GFF3 attributes (ID=, Parent=) → GTF attributes (transcript_id, gene_id)
#   2. Seqname (col1) set to transcript/mRNA ID (matching FASTA headers)
#   3. Feature types: gene→skipped, mRNA→transcript, exon/CDS→kept
#   4. Genomic coordinates preserved from GFF3
#
# Usage: ./convert_gff3_to_gtf.sh [input.gff] [output.gtf]
# ==============================================================================

set -euo pipefail

# ==============================================================================
# CONFIGURATION
# ==============================================================================

INPUT_GFF="${1:-GPE001970.gff}"
OUTPUT_GTF="${2:-GPE001970_transcripts.gtf}"

# ==============================================================================
# VALIDATION
# ==============================================================================

if [[ ! -f "$INPUT_GFF" ]]; then
    echo "ERROR: Input file not found: $INPUT_GFF"
    exit 1
fi

echo "Converting GFF3 to GTF..."
echo "  Input:  $INPUT_GFF"
echo "  Output: $OUTPUT_GTF"
echo ""

# ==============================================================================
# CONVERSION (two-pass awk)
# ==============================================================================
# Pass 1: Build mRNA→gene mapping from mRNA lines
# Pass 2: Output transcript/exon/CDS lines in GTF format
#
# Output col1 = mRNA ID (matches transcript FASTA headers)
# Output col9 = transcript_id "mRNA_ID"; gene_id "gene_ID";
# ==============================================================================

awk -F'\t' -v OFS='\t' '
# Skip comment lines
/^#/ { next }

# --- Parse attribute helpers ---
function get_attr(attrs, key,    parts, i, kv) {
    split(attrs, parts, ";")
    for (i in parts) {
        # Trim leading whitespace
        gsub(/^[ \t]+/, "", parts[i])
        if (parts[i] ~ "^" key "=") {
            sub("^" key "=", "", parts[i])
            return parts[i]
        }
    }
    return ""
}

# --- Pass 1 & 2 combined via NR==FNR ---
# First pass: collect mRNA→gene_id mapping
NR == FNR {
    if ($3 == "mRNA") {
        mrna_id = get_attr($9, "ID")
        parent  = get_attr($9, "Parent")
        if (mrna_id != "" && parent != "") {
            gene_of[mrna_id] = parent
        }
    }
    next
}

# --- Second pass: output GTF lines ---
/^#/ { next }

$3 == "mRNA" {
    mrna_id = get_attr($9, "ID")
    gene_id = gene_of[mrna_id]
    if (mrna_id == "" || gene_id == "") next

    # transcript line
    print mrna_id, $2, "transcript", $4, $5, $6, $7, ".", \
        "transcript_id \"" mrna_id "\"; gene_id \"" gene_id "\";"
    next
}

$3 == "exon" || $3 == "CDS" {
    parent = get_attr($9, "Parent")
    if (parent == "") next
    gene_id = gene_of[parent]
    if (gene_id == "") next

    # Use frame from GFF for CDS, "." for exon
    frame = ($3 == "CDS") ? $8 : "."

    print parent, $2, $3, $4, $5, $6, $7, frame, \
        "transcript_id \"" parent "\"; gene_id \"" gene_id "\";"
    next
}
' "$INPUT_GFF" "$INPUT_GFF" > "$OUTPUT_GTF"

# ==============================================================================
# VALIDATION
# ==============================================================================

INPUT_GENES=$(awk -F'\t' '$3=="gene"' "$INPUT_GFF" | wc -l)
INPUT_MRNAS=$(awk -F'\t' '$3=="mRNA"' "$INPUT_GFF" | wc -l)
OUTPUT_LINES=$(wc -l < "$OUTPUT_GTF")
OUTPUT_TRANSCRIPTS=$(awk -F'\t' '$3=="transcript"' "$OUTPUT_GTF" | wc -l)
OUTPUT_EXONS=$(awk -F'\t' '$3=="exon"' "$OUTPUT_GTF" | wc -l)
OUTPUT_CDS=$(awk -F'\t' '$3=="CDS"' "$OUTPUT_GTF" | wc -l)

echo "Input  (GFF3):  $INPUT_GENES genes, $INPUT_MRNAS mRNAs"
echo "Output (GTF):   $OUTPUT_TRANSCRIPTS transcripts, $OUTPUT_EXONS exons, $OUTPUT_CDS CDS ($OUTPUT_LINES total lines)"
echo ""

echo "Sample output (first 5 lines):"
head -5 "$OUTPUT_GTF"
echo ""

echo "Validating GTF format..."
if head -1 "$OUTPUT_GTF" | grep -qE 'transcript_id "[^"]+"; gene_id "[^"]+";'; then
    echo "  OK: Attribute format is correct"
else
    echo "  WARNING: Attribute format may need review"
fi

# Verify all transcripts have matching FASTA entries (if FASTA exists nearby)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FASTA_FILE="$SCRIPT_DIR/../../fasta/reference_genomes/GPE001970_transcripts.fa"
if [[ -f "$FASTA_FILE" ]]; then
    FASTA_IDS=$(grep -c '^>' "$FASTA_FILE")
    echo "  FASTA transcripts: $FASTA_IDS"
    echo "  GTF transcripts:   $OUTPUT_TRANSCRIPTS"
    if [[ "$FASTA_IDS" -eq "$OUTPUT_TRANSCRIPTS" ]]; then
        echo "  OK: Transcript counts match"
    else
        echo "  WARNING: Transcript count mismatch"
    fi
fi

echo ""
echo "Done! Output saved to: $OUTPUT_GTF"
