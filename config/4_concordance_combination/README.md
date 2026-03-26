# Concordance Combination Config Crosses

This folder defines reusable config crosses for `run_concordance_analysis.sh`.

## How To Use

Run with a single cross:

```bash
bash run_concordance_analysis.sh config/4_concordance_combination/defaults.toml
```

Run with cross basename:

```bash
bash run_concordance_analysis.sh defaults
```

Chain multiple crosses using env var (later cross overrides earlier values):

```bash
CONCORDANCE_CONFIG_CROSSES="defaults,cross_methods_vs_methods" bash run_concordance_analysis.sh
```

## Cross Types

- `defaults.toml`: baseline defaults (single-run behavior)
- `cross_genomes_vs_genomes.toml`: compare TPM across reference genomes listed in `concordance_genomes`
- `cross_equivalent_gene_between_genomes.toml`: per-gene concordance of positionally equivalent genes across genomes
- `cross_methods_vs_methods.toml`: iterate across `METHOD_COMBINATIONS`
- `cross_within_gene_group_vs_within_gene_group.toml`: iterate across `GENE_GROUP_COMBINATIONS`
- `cross_full_factorial_example.toml`: nested combinations (genomes x methods x genes)

## Combination Flags

- `RUN_ALL_MASTER_REFERENCES=TRUE`
- `RUN_ALL_METHOD_COMBINATIONS=TRUE`
- `RUN_ALL_GENE_GROUP_COMBINATIONS=TRUE`

## Combination Arrays

- `MASTER_REFERENCES=(...)`
- `METHOD_COMBINATIONS=("M1_HISAT2_RefGuided,M3_STAR_Align" ... )`
- `GENE_GROUP_COMBINATIONS=("GroupA,GroupB" ... )`

## Output Layout

All outputs default under `4_CONCORDANCE_ANALYSIS/`.
Nested combination mode writes to per-combination subfolders to avoid overwrite.

## TPM Comparability Note

StringTie TPM (M1/M2) and tximport TPM (M3/M4/M5) use different normalization algorithms:

- **StringTie TPM**: Per-base coverage normalized to sum = 1 million. Tolerates genes without exon annotations.
- **tximport TPM**: Length-normalized, bias-corrected read counts normalized to sum = 1 million. Requires exon annotations; genes with zero effective length are excluded.

Empirical differences: up to ~30% for high-expression genes; order-of-magnitude or presence/absence for annotation-edge-case genes (e.g., genes lacking exon annotations). Rank-order comparisons (Spearman) are appropriate; absolute TPM value comparisons across methods are not.

Per-method gene counts (GPE001970): M1 = 32,892; M3/M4 = 32,862 (-30); M5 = 32,444 (-448). The concordance pipeline intersects gene sets before comparison.

## M1 Per-Sample Strandedness Caveat

The GPE001970 M1 outputs (2026-03-10/11 run) were generated with `--fr` (forward-stranded) applied to all 32 samples. Eight PRJNA328564 samples are not forward-stranded, causing StringTie to discard 80-99% of reads. These samples show M1 Spearman correlations of 0.03-0.40 with M2-M5 (which agree at >0.94 among themselves):

- SRR3884679 (Pistils), SRR3884686 (Buds_0.7cm), SRR3884684 (Senescent_leaves), SRR3884620 (Fruits_Stage_1), SRR3884675 (Roots), SRR3884677 (Cotyledons), SRR3884631 (Fruits_6cm), SRR3884608 (Fruits_1cm)

M1 TPM, coverage, and FPKM values for these 8 samples are unreliable. The config has been updated to `hisat2_strandness = ""` (unstranded). A re-run with the current config will fix this issue.

## Count-Model Notes

- M1 prepDE integer counts: valid for DESeq2 via `DESeqDataSetFromMatrix()`.
- M3/M4 "NumReads" and M5 "expected_count": fractional EM-estimated values from Salmon/RSEM. Valid for DESeq2 via `DESeqDataSetFromTximport()`. Do not round to integers.
