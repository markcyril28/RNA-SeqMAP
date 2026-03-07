# STAR Alignment Preprocessing

This folder will contain STAR-specific matrix building scripts.

## Status: TODO

STAR alignment uses:
- **Aligner**: STAR (Spliced Transcripts Alignment to a Reference)
- **Quantification Options**: 
  - Built-in `--quantMode GeneCounts`
  - featureCounts (subread)
  - RSEM (via STAR-RSEM pipeline)

## Expected Scripts

- `star_matrix_builder.sh` - Build count matrices from STAR output
- Or use featureCounts/RSEM output with appropriate tximport

## Count Types

- Gene counts (from `--quantMode GeneCounts`)
- TPM (if using RSEM post-alignment)
