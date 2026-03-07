# 0_INPUT_FASTAs Directory Reorganization - Migration Guide

**Date**: January 18, 2026  
**Status**: Completed

## Summary

The `0_INPUT_FASTAs` directory has been reorganized into logical subdirectories to improve clarity and maintainability. All scripts referencing these files have been updated.

## Directory Structure Changes

### Previous Structure
```
0_INPUT_FASTAs/
├── All_Smel_Genes.fasta
├── Eggplant_V4.1_transcripts.function.fa
├── *.gtf files
├── *.gene_trans_map files
├── Control genes
├── Experimental files
└── Documentation
```

### New Structure
```
0_INPUT_FASTAs/
├── fasta/
│   ├── reference_genomes/     # Main reference FASTA files
│   ├── gene_sets/             # Specific gene sets (DMP, GIF, GRF)
│   ├── control_genes/         # Control and reference genes
│   └── experimental/          # Experimental gene combinations
├── gtf/
│   └── reference/             # Gene annotation GTF files
├── mappings/                  # Gene mapping files (.gene_trans_map, .csv)
├── scripts/                   # Utility scripts
├── docs/                      # Documentation
└── temp/                      # Temporary files
```

## Files Updated

### Configuration Files
1. **config/HPC_full_run_config.sh**
   - Updated `Eggplant_V4_1_transcripts_function_FASTA_FILE` path
   - Updated `ALL_Smel_Genes_Full_Name_reformatted_GTF_FILE` path
   - Updated `decoy` file path
   - Updated all paths in `ALL_FASTA_FILES` array

2. **config/local_test_config.sh**
   - Same updates as HPC_full_run_config.sh

### R Scripts
3. **4_POST_PROC_v1_postproc_test/M3_STAR_Align/B_M3_modules/6_coexpression_wgcna.R**
   - Updated `gene_info_paths` to include `/mappings/` subdirectory

4. **4_POST_PROC_v1_postproc_test/main_modules/preprocessing/Salmon/tximport_salmon_to_matrices.R**
   - Updated `tx2gene_file` path to look in `/mappings/` subdirectory

5. **4_POST_PROC_v1_postproc_test/main_modules/analysis_modules/1_utility_functions.R**
   - Updated gene_info.csv lookup to check `/mappings/` subdirectory

6. **4_POST_PROC_v1_postproc_test/M3_STAR_Align/B_M3_modules/8_gene_set_enrichment.R**
   - Updated comment example for custom gene sets to use `/docs/` subdirectory

### Shell Scripts
7. **init_setup.sh**
   - Updated find command to search `./0_INPUT_FASTAs/fasta` for FASTA files

## Path Translation Table

| Old Path | New Path | File Type |
|----------|----------|-----------|
| `0_INPUT_FASTAs/Eggplant_V4.1_transcripts.function.fa` | `0_INPUT_FASTAs/fasta/reference_genomes/Eggplant_V4.1_transcripts.function.fa` | Reference FASTA |
| `0_INPUT_FASTAs/All_Smel_Genes.fasta` | `0_INPUT_FASTAs/fasta/reference_genomes/All_Smel_Genes.fasta` | Reference FASTA |
| `0_INPUT_FASTAs/All_Smel_Genes_Full_Name_reformatted.gtf` | `0_INPUT_FASTAs/gtf/reference/All_Smel_Genes_Full_Name_reformatted.gtf` | GTF Annotation |
| `0_INPUT_FASTAs/TEST.fasta` | `0_INPUT_FASTAs/fasta/experimental/TEST.fasta` | Experimental FASTA |
| `0_INPUT_FASTAs/SmelGIF.fasta` | `0_INPUT_FASTAs/fasta/gene_sets/SmelGIF.fasta` | Gene Set |
| `0_INPUT_FASTAs/Best_Control_Genes.fa` | `0_INPUT_FASTAs/fasta/control_genes/Best_Control_Genes.fa` | Control Genes |
| `0_INPUT_FASTAs/Eggplant_V4.1_transcripts.function.gene_info.csv` | `0_INPUT_FASTAs/mappings/Eggplant_V4.1_transcripts.function.gene_info.csv` | Mapping CSV |
| `0_INPUT_FASTAs/*.gene_trans_map` | `0_INPUT_FASTAs/mappings/*.gene_trans_map` | Mapping Files |
| `0_INPUT_FASTAs/extract_gene_info.py` | `0_INPUT_FASTAs/scripts/extract_gene_info.py` | Python Script |
| `0_INPUT_FASTAs/All_Genes_Info.md` | `0_INPUT_FASTAs/docs/All_Genes_Info.md` | Documentation |

## Files Moved

### FASTA Files - Reference Genomes
- Eggplant_V4.1_transcripts.function.fa → fasta/reference_genomes/
- All_Smel_Genes.fasta → fasta/reference_genomes/
- All_Smel_Genes_Full_Name.fasta → fasta/reference_genomes/

### FASTA Files - Gene Sets
- SmelDMPs.fa → fasta/gene_sets/
- SmelGIF.fasta → fasta/gene_sets/
- SmelGRF.fasta → fasta/gene_sets/

### FASTA Files - Control Genes
- Best_Control_Genes.fa → fasta/control_genes/
- Best_Cell_Cycle_Associated_Control_Genes.fa → fasta/control_genes/
- Putative_Control_Genes.fasta → fasta/control_genes/
- Putative_Cell_Cycle_Associated_Control_Genes.fasta → fasta/control_genes/

### FASTA Files - Experimental
- SmelDMP_CDS_with_Best_Control_Genes.fasta → fasta/experimental/
- SmelDMP_CDS_with_Putative_Control_Genes.fasta → fasta/experimental/
- SmelGIF_with_Best_Control_Cyclo.fasta → fasta/experimental/
- SmelGIF_with_Cell_Cycle_Control_genes.fasta → fasta/experimental/
- SmelGRF_with_Best_Control_Cyclo.fasta → fasta/experimental/
- SmelGRF_with_Cell_Cycle_Control_genes.fasta → fasta/experimental/
- SmelGRF-GIF_with_Best_Control_Cyclo.fasta → fasta/experimental/
- SmelGRF-GIF_with_Best_Control_Cyclo_Experiment.fasta → fasta/experimental/
- TEST.fasta → fasta/experimental/

### GTF Files
- All_Smel_Genes_Full_Name.gtf → gtf/reference/
- All_Smel_Genes_Full_Name_FIXED.gtf → gtf/reference/
- All_Smel_Genes_Full_Name_reformatted.gtf → gtf/reference/
- All_Smel_Genes_transcript_based.gtf → gtf/reference/
- Eggplant_V4.1_function_IPR_final_formatted_v3.gtf → gtf/reference/

### Mapping Files
- All_Smel_Genes.fasta.gene_trans_map → mappings/
- Eggplant_V4.1_transcripts.function.fa.gene_trans_map → mappings/
- chr_to_sequence_mapping.tsv → mappings/
- Eggplant_V4.1_transcripts.function.gene_info.csv → mappings/

### Scripts
- extract_gene_info.py → scripts/

### Documentation
- All_Genes_Info.md → docs/
- GTF_FIX_README.md → docs/

### Temporary Files
- d2utmpCYrzNM → temp/
- d2utmpPW39OU → temp/

## Testing Recommendations

After migration, verify:

1. **Configuration Loading**: Ensure config files load correctly
   ```bash
   source config/HPC_full_run_config.sh
   echo $Eggplant_V4_1_transcripts_function_FASTA_FILE
   ```

2. **R Script Paths**: Verify R scripts can find gene_info.csv
   ```R
   # Test from within R scripts
   gene_info_file <- "0_INPUT_FASTAs/mappings/Eggplant_V4.1_transcripts.function.gene_info.csv"
   file.exists(gene_info_file)
   ```

3. **Pipeline Runs**: Test a small pipeline run with updated paths

## Backward Compatibility

**Note**: Old paths will no longer work. Any custom scripts or notebooks referencing the old paths must be updated manually using the Path Translation Table above.

## Documentation

New README created at:
- `0_INPUT_FASTAs/README.md` - Complete directory structure documentation

## Contact

For issues or questions about this migration, check:
- `0_INPUT_FASTAs/README.md` - Directory organization
- `0_INPUT_FASTAs/docs/GTF_FIX_README.md` - GTF file processing details
- `0_INPUT_FASTAs/docs/All_Genes_Info.md` - Gene information
