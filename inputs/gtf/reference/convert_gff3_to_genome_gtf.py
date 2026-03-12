#!/usr/bin/env python3
"""
GFF3 to GTF Converter for HISAT2 Ref-Guided and STAR Alignment
===============================================================
Converts a genome-level GFF3 file to GTF format compatible with:
  - HISAT2 reference-guided alignment (hisat2_extract_splice_sites.py,
    hisat2_extract_exons.py)
  - STAR genome index generation (--sjdbGTFfile)
  - StringTie quantification

Key transformations:
  1. GFF3 attributes (ID=, Parent=) → GTF attributes (gene_id "", transcript_id "")
  2. Column 1 (seqname) kept as genomic chromosome/scaffold ID
  3. Feature types: gene→gene, mRNA→transcript, exon→exon, CDS→CDS
  4. Genomic coordinates preserved (1-based, as in source GFF3)

Usage:
  python3 convert_gff3_to_genome_gtf.py <input.gff> <output.gtf>

Example:
  python3 convert_gff3_to_genome_gtf.py Eggplant_V4.1_function_IPR_final.gff \\
                                        Eggplant_V4.1_function_IPR_final.gtf
"""

import sys
import re


def get_attr(attr_string, key):
    """Extract a value from a GFF3 attribute string (key=value pairs)."""
    match = re.search(r'(?:^|;)\s*' + re.escape(key) + r'=([^;]+)', attr_string)
    return match.group(1).strip() if match else ""


def gff3_to_gtf(input_gff, output_gtf):
    # --- Pass 1: Build gene ID lookup (mRNA_id → gene_id) ---
    mrna_to_gene = {}
    with open(input_gff) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 9:
                continue
            feat = cols[2]
            if feat == 'mRNA':
                mrna_id = get_attr(cols[8], 'ID')
                parent  = get_attr(cols[8], 'Parent')
                if mrna_id and parent:
                    mrna_to_gene[mrna_id] = parent

    # --- Pass 2: Write GTF ---
    gene_count = trans_count = exon_count = cds_count = 0

    with open(input_gff) as fh, open(output_gtf, 'w') as out:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 9:
                continue

            seqname = cols[0]
            source  = cols[1]
            feat    = cols[2]
            start   = cols[3]
            end     = cols[4]
            score   = cols[5]
            strand  = cols[6]
            frame   = cols[7]
            attrs   = cols[8]

            if feat == 'gene':
                gene_id = get_attr(attrs, 'ID')
                if not gene_id:
                    continue
                gtf_attr = f'gene_id "{gene_id}"; transcript_id "";'
                out.write('\t'.join([seqname, source, 'gene',
                                     start, end, score, strand, '.',
                                     gtf_attr]) + '\n')
                gene_count += 1

            elif feat == 'mRNA':
                mrna_id = get_attr(attrs, 'ID')
                gene_id = get_attr(attrs, 'Parent')
                if not mrna_id or not gene_id:
                    continue
                gtf_attr = f'gene_id "{gene_id}"; transcript_id "{mrna_id}";'
                out.write('\t'.join([seqname, source, 'transcript',
                                     start, end, score, strand, '.',
                                     gtf_attr]) + '\n')
                trans_count += 1

            elif feat == 'exon':
                parent = get_attr(attrs, 'Parent')
                gene_id = mrna_to_gene.get(parent, '')
                if not parent or not gene_id:
                    continue
                gtf_attr = f'gene_id "{gene_id}"; transcript_id "{parent}";'
                out.write('\t'.join([seqname, source, 'exon',
                                     start, end, score, strand, '.',
                                     gtf_attr]) + '\n')
                exon_count += 1

            elif feat == 'CDS':
                parent = get_attr(attrs, 'Parent')
                gene_id = mrna_to_gene.get(parent, '')
                if not parent or not gene_id:
                    continue
                # Use frame from GFF3 for CDS; warn if unknown and defaulting to 0
                if frame not in ('0', '1', '2'):
                    print(f"  Warning: CDS with unknown frame '{frame}' for {parent}, defaulting to 0", file=sys.stderr)
                cds_frame = frame if frame in ('0', '1', '2') else '0'
                gtf_attr = f'gene_id "{gene_id}"; transcript_id "{parent}";'
                out.write('\t'.join([seqname, source, 'CDS',
                                     start, end, score, strand, cds_frame,
                                     gtf_attr]) + '\n')
                cds_count += 1

    print(f"Conversion complete: {input_gff} → {output_gtf}")
    print(f"  Genes:       {gene_count}")
    print(f"  Transcripts: {trans_count}")
    print(f"  Exons:       {exon_count}")
    print(f"  CDS:         {cds_count}")
    print(f"  Total lines: {gene_count + trans_count + exon_count + cds_count}")


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print(f"Usage: python3 {sys.argv[0]} <input.gff> <output.gtf>")
        sys.exit(1)
    gff3_to_gtf(sys.argv[1], sys.argv[2])
