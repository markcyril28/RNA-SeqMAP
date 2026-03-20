#!/usr/bin/env python3
"""
Discover alternatively spliced genes from a single GTF and optionally extract
nucleotide transcript FASTA from a reference genome FASTA.

Definition used here:
- A gene is flagged as alternatively spliced when it has >=2 transcripts and
  >=2 unique exon structures (exon coordinate chains).
"""

from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple


@dataclass
class ExonRec:
    seqname: str
    start: int
    end: int
    strand: str


def parse_gtf_attributes(attr: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for part in attr.strip().split(";"):
        part = part.strip()
        if not part:
            continue
        if " " not in part:
            continue
        key, val = part.split(" ", 1)
        out[key.strip()] = val.strip().strip('"')
    return out


def reverse_complement(seq: str) -> str:
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]


def read_fasta(path: Path) -> Dict[str, str]:
    seqs: Dict[str, List[str]] = {}
    curr = None
    with path.open("r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                curr = line[1:].split()[0]
                if curr not in seqs:
                    seqs[curr] = []
            else:
                if curr is None:
                    raise ValueError(f"Invalid FASTA format in {path}")
                seqs[curr].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def wrap_fasta(seq: str, width: int = 80) -> str:
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width))


def discover(gtf_path: Path):
    # gene -> transcript -> exon records
    gene_tx_exons: Dict[str, Dict[str, List[ExonRec]]] = defaultdict(lambda: defaultdict(list))

    with gtf_path.open("r", encoding="utf-8") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            seqname, _source, feature, start, end, _score, strand, _frame, attrs = fields
            if feature != "exon":
                continue
            ad = parse_gtf_attributes(attrs)
            gene_id = ad.get("gene_id")
            tx_id = ad.get("transcript_id")
            if not gene_id or not tx_id:
                continue
            try:
                s, e = int(start), int(end)
            except ValueError:
                continue
            if s > e:
                continue
            gene_tx_exons[gene_id][tx_id].append(
                ExonRec(seqname=seqname, start=s, end=e, strand=strand)
            )

    alt_genes = []
    all_gene_stats = []

    # O(G × T × E) where G=genes, T=transcripts/gene, E=exons/transcript
    for gene_id, tx_map in gene_tx_exons.items():
        tx_structures = {}
        for tx_id, exons in tx_map.items():
            exons_sorted = sorted(exons, key=lambda x: (x.start, x.end))
            # exon chain includes coordinates only; seqname and strand should match within transcript
            chain = tuple((e.start, e.end) for e in exons_sorted)
            tx_structures[tx_id] = chain

        tx_count = len(tx_structures)
        unique_structures = len(set(tx_structures.values()))
        is_alt = tx_count >= 2 and unique_structures >= 2

        all_gene_stats.append((gene_id, tx_count, unique_structures, is_alt))
        if is_alt:
            alt_genes.append(gene_id)

    return gene_tx_exons, all_gene_stats, set(alt_genes)


def write_gene_stats(path: Path, all_gene_stats) -> None:
    with path.open("w", encoding="utf-8") as out:
        out.write("gene_id\ttranscript_count\tunique_exon_structures\tis_alternative_spliced\n")
        for gene_id, tx_count, unique_structures, is_alt in all_gene_stats:
            out.write(f"{gene_id}\t{tx_count}\t{unique_structures}\t{str(is_alt).upper()}\n")


def write_alt_gene_list(path: Path, alt_genes) -> None:
    with path.open("w", encoding="utf-8") as out:
        for g in sorted(alt_genes):
            out.write(f"{g}\n")


def filter_gtf_by_genes(gtf_in: Path, gtf_out: Path, alt_genes: set) -> int:
    kept = 0
    with gtf_in.open("r", encoding="utf-8") as fin, gtf_out.open("w", encoding="utf-8") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            attrs = parse_gtf_attributes(fields[8])
            gene_id = attrs.get("gene_id")
            if gene_id in alt_genes:
                fout.write(line)
                kept += 1
    return kept


def extract_transcript_fasta(
    out_fa: Path,
    gene_tx_exons: Dict[str, Dict[str, List[ExonRec]]],
    alt_genes: set,
    genome_fasta: Path,
) -> Tuple[int, int]:
    genome = read_fasta(genome_fasta)
    written = 0
    skipped = 0
    skipped_seqnames: set = set()

    # O(A × T × E) where A=alt-spliced genes, T=transcripts/gene, E=exons/transcript
    with out_fa.open("w", encoding="utf-8") as out:
        for gene_id in sorted(alt_genes):
            for tx_id, exons in sorted(gene_tx_exons[gene_id].items()):
                if not exons:
                    skipped += 1
                    continue

                seqname = exons[0].seqname
                strand = exons[0].strand
                if seqname not in genome:
                    skipped += 1
                    skipped_seqnames.add(seqname)
                    continue

                chrom_seq = genome[seqname]
                exons_sorted = sorted(exons, key=lambda e: (e.start, e.end))
                pieces = []
                for e in exons_sorted:
                    # GTF is 1-based inclusive; Python slice is 0-based, end-exclusive.
                    s = max(1, e.start)
                    t = min(len(chrom_seq), e.end)
                    if t < s:
                        continue
                    pieces.append(chrom_seq[s - 1 : t])
                tx_seq = "".join(pieces)
                if strand == "-":
                    tx_seq = reverse_complement(tx_seq)

                if not tx_seq:
                    skipped += 1
                    continue

                header = f">{tx_id} gene={gene_id} seq={seqname} strand={strand} len={len(tx_seq)}"
                out.write(header + "\n")
                out.write(wrap_fasta(tx_seq) + "\n")
                written += 1

    if skipped_seqnames:
        print(
            f"WARNING: {len(skipped_seqnames)} GTF seqname(s) not found in genome FASTA: "
            f"{', '.join(sorted(skipped_seqnames)[:5])}"
            f"{'...' if len(skipped_seqnames) > 5 else ''}",
            file=sys.stderr,
        )
    return written, skipped


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Discover alternatively spliced genes from one GTF and optionally extract transcript FASTA."
    )
    parser.add_argument("--gtf", required=True, help="Input GTF path")
    parser.add_argument("--output-dir", required=True, help="Output directory")
    parser.add_argument(
        "--genome-fasta",
        default=None,
        help="Optional genome FASTA path to extract transcript nucleotide sequences",
    )
    args = parser.parse_args()

    gtf_path = Path(args.gtf)
    if not gtf_path.is_file():
        print(f"ERROR: GTF file not found: {gtf_path}", file=sys.stderr)
        sys.exit(1)
    if args.genome_fasta and not Path(args.genome_fasta).is_file():
        print(f"ERROR: Genome FASTA not found: {args.genome_fasta}", file=sys.stderr)
        sys.exit(1)

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    gene_tx_exons, all_gene_stats, alt_genes = discover(gtf_path)

    if not all_gene_stats:
        print("WARNING: No valid exon records with gene_id and transcript_id found in GTF.", file=sys.stderr)
        print("Check that the GTF contains exon features with gene_id and transcript_id attributes.", file=sys.stderr)

    stats_tsv = out_dir / "gene_splicing_summary.tsv"
    alt_genes_txt = out_dir / "alternatively_spliced_genes.txt"
    alt_gtf = out_dir / "alternatively_spliced_genes.gtf"

    write_gene_stats(stats_tsv, all_gene_stats)
    write_alt_gene_list(alt_genes_txt, alt_genes)
    kept_lines = filter_gtf_by_genes(gtf_path, alt_gtf, alt_genes)

    print(f"Parsed genes: {len(all_gene_stats)}")
    print(f"Alternatively spliced genes: {len(alt_genes)}")
    print(f"Filtered GTF lines written: {kept_lines}")
    print(f"Summary TSV: {stats_tsv}")
    print(f"Gene list: {alt_genes_txt}")
    print(f"Filtered GTF: {alt_gtf}")

    if args.genome_fasta:
        out_fa = out_dir / "alternatively_spliced_transcripts.fa"
        written, skipped = extract_transcript_fasta(
            out_fa=out_fa,
            gene_tx_exons=gene_tx_exons,
            alt_genes=alt_genes,
            genome_fasta=Path(args.genome_fasta),
        )
        print(f"Transcript FASTA: {out_fa}")
        print(f"Transcript FASTA written: {written}")
        print(f"Transcript FASTA skipped: {skipped}")


if __name__ == "__main__":
    main()
