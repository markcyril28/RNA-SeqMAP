# 0_INPUT_FASTAs Directory Organization

This directory contains all input reference files, gene sequences, and annotations used in the HeatSeq pipeline.

## Directory Structure

```
0_INPUT_FASTAs/
├── fasta/
│   ├── reference_genomes/     # Main reference genome FASTA files
│   │   ├── Eggplant_V4.1_transcripts.function.fa
│   │   ├── All_Smel_Genes.fasta
│   │   └── All_Smel_Genes_Full_Name.fasta
│   ├── gene_sets/             # Specific gene sets (DMP, GIF, GRF)
│   │   ├── SmelDMPs.fa
│   │   ├── SmelGIF.fasta
│   │   └── SmelGRF.fasta
│   ├── control_genes/         # Control and reference genes
│   │   ├── Best_Control_Genes.fa
│   │   ├── Best_Cell_Cycle_Associated_Control_Genes.fa
│   │   ├── Putative_Control_Genes.fasta
│   │   └── Putative_Cell_Cycle_Associated_Control_Genes.fasta
│   └── experimental/          # Experimental gene combinations
│       ├── SmelDMP_CDS_with_Best_Control_Genes.fasta
│       ├── SmelDMP_CDS_with_Putative_Control_Genes.fasta
│       ├── SmelGIF_with_Best_Control_Cyclo.fasta
│       ├── SmelGIF_with_Cell_Cycle_Control_genes.fasta
│       ├── SmelGRF_with_Best_Control_Cyclo.fasta
│       ├── SmelGRF_with_Cell_Cycle_Control_genes.fasta
│       ├── SmelGRF-GIF_with_Best_Control_Cyclo.fasta
│       ├── SmelGRF-GIF_with_Best_Control_Cyclo_Experiment.fasta
│       └── TEST.fasta
├── gtf/
│   └── reference/             # Gene annotation GTF files
│       ├── All_Smel_Genes_Full_Name.gtf
│       ├── All_Smel_Genes_Full_Name_FIXED.gtf
│       ├── All_Smel_Genes_Full_Name_reformatted.gtf
│       ├── All_Smel_Genes_transcript_based.gtf
│       └── Eggplant_V4.1_function_IPR_final_formatted_v3.gtf
├── mappings/                  # Gene mapping and info files
│   ├── All_Smel_Genes.fasta.gene_trans_map
│   ├── Eggplant_V4.1_transcripts.function.fa.gene_trans_map
│   ├── chr_to_sequence_mapping.tsv
│   └── Eggplant_V4.1_transcripts.function.gene_info.csv
├── scripts/                   # Utility scripts for processing
│   └── extract_gene_info.py
├── docs/                      # Documentation and notes
│   ├── All_Genes_Info.md
│   └── GTF_FIX_README.md
└── temp/                      # Temporary files
    ├── d2utmpCYrzNM
    └── d2utmpPW39OU
```

## File Descriptions

### Reference Genomes
- **Eggplant_V4.1_transcripts.function.fa**: Primary eggplant reference transcriptome
- **All_Smel_Genes.fasta**: Complete gene set for Solanum melongena
- **All_Smel_Genes_Full_Name.fasta**: Same as above with full gene names

### Gene Sets
- **SmelDMPs.fa**: DNA Methylation Polymorphism genes
- **SmelGIF.fa**: Growth-Regulating Factor-Interacting Factor genes
- **SmelGRF.fa**: Growth-Regulating Factor genes

### Control Genes
- Reference genes used for normalization and validation
- Includes both validated and putative control genes
- Cell cycle associated genes for temporal analysis

### GTF Files
- **All_Smel_Genes_Full_Name_reformatted.gtf**: PRIMARY GTF file used in most analyses
- **All_Smel_Genes_Full_Name_FIXED.gtf**: Fixed version resolving StringTie issues
- Other GTF variants for specific use cases

## Usage

Reference these files in config files using the new paths:
```bash
# Example paths
"0_INPUT_FASTAs/fasta/reference_genomes/Eggplant_V4.1_transcripts.function.fa"
"0_INPUT_FASTAs/gtf/reference/All_Smel_Genes_Full_Name_reformatted.gtf"
"0_INPUT_FASTAs/mappings/Eggplant_V4.1_transcripts.function.gene_info.csv"
```

## Notes

- All FASTA files have been processed with dos2unix during init_setup.sh
- GTF files are formatted for StringTie compatibility
- See docs/GTF_FIX_README.md for details on GTF processing issues and solutions
