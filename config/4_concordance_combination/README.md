# Concordance Combination Config Classes

This folder defines reusable config classes for `run_concordance.sh`.

## How To Use

Run with a single class:

```bash
bash run_concordance.sh config/4_concordance_combination/defaults.toml
```

Run with class basename:

```bash
bash run_concordance.sh defaults
```

Chain multiple classes using env var (later class overrides earlier values):

```bash
CONCORDANCE_CONFIG_CLASSES="defaults,class_methods_vs_methods" bash run_concordance.sh
```

## Class Types

- `defaults.toml`: baseline defaults (single-run behavior)
- `class_genomes_vs_genomes.toml`: iterate across `MASTER_REFERENCES`
- `class_methods_vs_methods.toml`: iterate across `METHOD_COMBINATIONS`
- `class_genes_vs_genes.toml`: iterate across `GENE_GROUP_COMBINATIONS`
- `class_full_factorial_example.toml`: nested combinations (genomes x methods x genes)

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
