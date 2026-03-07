# Quick Start Guide - Using the Reorganized 0_INPUT_FASTAs

## What Changed?

The `0_INPUT_FASTAs` directory is now organized into subdirectories by file type. **All scripts have been updated to use the new paths.**

## Common File Locations

### For Most Users

If you're using the standard pipeline, **no action needed** - all config files are already updated!

### If You're Writing Custom Scripts

Use these paths:

```bash
# Primary reference genome
REFERENCE_FASTA="0_INPUT_FASTAs/fasta/reference_genomes/Eggplant_V4.1_transcripts.function.fa"

# Primary GTF annotation
REFERENCE_GTF="0_INPUT_FASTAs/gtf/reference/All_Smel_Genes_Full_Name_reformatted.gtf"

# Gene information CSV
GENE_INFO="0_INPUT_FASTAs/mappings/Eggplant_V4.1_transcripts.function.gene_info.csv"

# Gene-transcript mapping
GENE_TRANS_MAP="0_INPUT_FASTAs/mappings/Eggplant_V4.1_transcripts.function.fa.gene_trans_map"
```

### Directory Quick Reference

| Need | Look In |
|------|---------|
| Reference genomes (FASTA) | `fasta/reference_genomes/` |
| GTF annotation files | `gtf/reference/` |
| Gene mapping files (.csv, .tsv, .gene_trans_map) | `mappings/` |
| Control genes | `fasta/control_genes/` |
| Experimental gene sets | `fasta/experimental/` |
| Utility scripts | `scripts/` |
| Documentation | `docs/` |

## Updating Old Scripts

Replace old paths with new ones:

```bash
# OLD:
INPUT="0_INPUT_FASTAs/Eggplant_V4.1_transcripts.function.fa"

# NEW:
INPUT="0_INPUT_FASTAs/fasta/reference_genomes/Eggplant_V4.1_transcripts.function.fa"
```

See [MIGRATION_GUIDE.md](MIGRATION_GUIDE.md) for a complete path translation table.

## Need Help?

- 📖 **Full documentation**: See [README.md](README.md)
- 🔄 **Migration details**: See [MIGRATION_GUIDE.md](MIGRATION_GUIDE.md)
- 📊 **Organization summary**: See [ORGANIZATION_SUMMARY.md](ORGANIZATION_SUMMARY.md)

## Already Updated Scripts

These scripts are already updated and will work without changes:
- ✅ All config files (`config/*.sh`)
- ✅ All R analysis scripts
- ✅ `init_setup.sh`

You can continue using the pipeline as normal!
