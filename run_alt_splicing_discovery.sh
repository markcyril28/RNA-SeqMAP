#!/usr/bin/env bash
set -euo pipefail

# Reproducible orchestrator:
# 1) Discover alternatively spliced genes from a single GTF.
# 2) Extract GTF entries for those genes.
# 3) Extract transcript nucleotide FASTA from genome FASTA.
# 4) Create a small test gene-group CSV for downstream post-processing.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DISCOVER_PY="$SCRIPT_DIR/modules/c_post_processing/utilities/discover_alt_splice_from_gtf.py"

GTF_PATH=""
GENOME_FASTA=""
OUTPUT_DIR=""
GENE_INFO_CSV="$SCRIPT_DIR/inputs/mapping/Eggplant_V4.1_transcripts.function.gene_info.csv"
GENE_GROUP_CSV=""
N_TEST_GENES=8
PYTHON_BIN="python3"

# User-editable defaults (used when a CLI argument is not provided)
DEFAULT_GTF_PATH="inputs/gtf/reference/Eggplant_V4.1_function_IPR_final_stringtie.gtf"
DEFAULT_GENOME_FASTA="inputs/fasta/reference_genome/Eggplant_V4.1.fa"
DEFAULT_OUTPUT_DIR="inputs/alt_splicing_Eggplant_V4.1"
DEFAULT_GENE_INFO_CSV="inputs/mapping/Eggplant_V4.1_transcripts.function.gene_info.csv"
DEFAULT_GENE_GROUP_CSV="inputs/gene_groups_csv/experimental/Eggplant_V4.1/AltSplice_Test_Genes.csv"
DEFAULT_N_TEST_GENES=8

usage() {
    cat <<EOF
Usage:
  bash run_alt_splicing_discovery.sh \\
    --gtf <path/to/annotation.gtf> \\
    --genome-fasta <path/to/genome.fa> \\
    --output-dir <path/to/output_dir> \\
    [--gene-info-csv <path/to/gene_info.csv>] \\
    [--gene-group-csv <path/to/test_gene_group.csv>] \\
    [--n-test-genes 8] \\
    [--python python3]

Notes:
- Alternative splicing is inferred from exon-structure differences between transcripts within each gene.
- gene_info.csv is used to assign readable names in the test gene-group CSV.
- If --gene-group-csv is not supplied, it is written under output-dir.
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --gtf)
            GTF_PATH="$2"; shift 2 ;;
        --genome-fasta)
            GENOME_FASTA="$2"; shift 2 ;;
        --output-dir)
            OUTPUT_DIR="$2"; shift 2 ;;
        --gene-info-csv)
            GENE_INFO_CSV="$2"; shift 2 ;;
        --gene-group-csv)
            GENE_GROUP_CSV="$2"; shift 2 ;;
        --n-test-genes)
            N_TEST_GENES="$2"; shift 2 ;;
        --python)
            PYTHON_BIN="$2"; shift 2 ;;
        -h|--help)
            usage; exit 0 ;;
        *)
            echo "ERROR: Unknown argument: $1" >&2
            usage
            exit 1 ;;
    esac
done

# Fill from top-of-file defaults when missing from CLI.
if [[ -z "$GTF_PATH" ]]; then GTF_PATH="$DEFAULT_GTF_PATH"; fi
if [[ -z "$GENOME_FASTA" ]]; then GENOME_FASTA="$DEFAULT_GENOME_FASTA"; fi
if [[ -z "$OUTPUT_DIR" ]]; then OUTPUT_DIR="$DEFAULT_OUTPUT_DIR"; fi
if [[ -z "$GENE_INFO_CSV" ]]; then GENE_INFO_CSV="$DEFAULT_GENE_INFO_CSV"; fi
if [[ -z "$GENE_GROUP_CSV" ]]; then GENE_GROUP_CSV="$DEFAULT_GENE_GROUP_CSV"; fi
if [[ -z "${N_TEST_GENES:-}" ]]; then N_TEST_GENES="$DEFAULT_N_TEST_GENES"; fi

if [[ "$GTF_PATH" != /* ]]; then GTF_PATH="$SCRIPT_DIR/$GTF_PATH"; fi
if [[ "$GENOME_FASTA" != /* ]]; then GENOME_FASTA="$SCRIPT_DIR/$GENOME_FASTA"; fi
if [[ "$OUTPUT_DIR" != /* ]]; then OUTPUT_DIR="$SCRIPT_DIR/$OUTPUT_DIR"; fi
if [[ "$GENE_INFO_CSV" != /* ]]; then GENE_INFO_CSV="$SCRIPT_DIR/$GENE_INFO_CSV"; fi

mkdir -p "$OUTPUT_DIR"

if [[ "$GENE_GROUP_CSV" != /* ]]; then
    GENE_GROUP_CSV="$SCRIPT_DIR/$GENE_GROUP_CSV"
fi

if [[ ! -f "$DISCOVER_PY" ]]; then
    echo "ERROR: Missing discovery script: $DISCOVER_PY" >&2
    exit 1
fi

if [[ ! -f "$GTF_PATH" ]]; then
    echo "ERROR: GTF not found: $GTF_PATH" >&2
    exit 1
fi

if [[ ! -f "$GENOME_FASTA" ]]; then
    echo "ERROR: Genome FASTA not found: $GENOME_FASTA" >&2
    exit 1
fi

echo "[1/3] Discovering alternatively spliced genes from GTF..."
"$PYTHON_BIN" "$DISCOVER_PY" \
    --gtf "$GTF_PATH" \
    --genome-fasta "$GENOME_FASTA" \
    --output-dir "$OUTPUT_DIR"

ALT_LIST="$OUTPUT_DIR/alternatively_spliced_genes.txt"
if [[ ! -s "$ALT_LIST" ]]; then
    echo "ERROR: Alternative spliced gene list missing or empty: $ALT_LIST" >&2
    exit 1
fi

echo "[2/3] Building a small test gene-group CSV (${N_TEST_GENES} genes)..."
"$PYTHON_BIN" - "$ALT_LIST" "$GENE_INFO_CSV" "$GENE_GROUP_CSV" "$N_TEST_GENES" <<'PY'
import csv
import re
import sys
from pathlib import Path

alt_list = Path(sys.argv[1])
gene_info = Path(sys.argv[2])
out_csv = Path(sys.argv[3])
n = int(sys.argv[4])

genes = []
with alt_list.open('r', encoding='utf-8') as fh:
    for line in fh:
        g = line.strip()
        if g:
            genes.append(g)

if not genes:
    raise SystemExit(f"No genes found in {alt_list}")

selected = genes[:n]

name_map = {}
if gene_info.exists():
    with gene_info.open('r', encoding='utf-8', newline='') as fh:
        reader = csv.DictReader(fh)
        # Expected columns: Gene_ID,Name where Gene_ID may be transcript-level (e.g., .01)
        for row in reader:
            tid = (row.get('Gene_ID') or '').strip()
            nm = (row.get('Name') or '').strip()
            if not tid:
                continue
            gid = re.sub(r'\.\d+$', '', tid)
            if gid and gid not in name_map:
                name_map[gid] = nm

out_csv.parent.mkdir(parents=True, exist_ok=True)
with out_csv.open('w', encoding='utf-8', newline='') as fh:
    w = csv.writer(fh)
    w.writerow(['Gene_ID', 'Shortened_Name'])
    for i, gid in enumerate(selected, start=1):
        name = name_map.get(gid, f'AS_Gene_{i:02d}')
        short = re.sub(r'[^A-Za-z0-9_\-]+', '_', name).strip('_')
        if not short:
            short = f'AS_Gene_{i:02d}'
        short = short[:40]
        w.writerow([gid, short])

print(f"Selected genes: {len(selected)}")
print(f"Gene-group CSV: {out_csv}")
PY

echo "[3/3] Done. Output summary:"
echo "- Gene summary TSV: $OUTPUT_DIR/gene_splicing_summary.tsv"
echo "- Alt gene list:    $OUTPUT_DIR/alternatively_spliced_genes.txt"
echo "- Alt-only GTF:     $OUTPUT_DIR/alternatively_spliced_genes.gtf"
echo "- Transcript FASTA: $OUTPUT_DIR/alternatively_spliced_transcripts.fa"
echo "- Test gene-group:  $GENE_GROUP_CSV"
