# Concordance Combination Config Crosses

This folder defines reusable config crosses for `run_concordance.sh`.

## How To Use

Run with a single cross:

```bash
bash run_concordance.sh config/4_concordance_combination/defaults.toml
```

Run with cross basename:

```bash
bash run_concordance.sh defaults
```

Chain multiple crosses using env var (later cross overrides earlier values):

```bash
CONCORDANCE_CONFIG_CROSSES="defaults,cross_methods_vs_methods" bash run_concordance.sh
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
