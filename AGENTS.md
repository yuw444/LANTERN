# PROJECT KNOWLEDGE BASE

**Generated:** 2026-03-17
**Commit:** f06e101
**Branch:** main

## OVERVIEW

R bioinformatics pipeline for ancestry-specific rare variant association analysis. Splits VCF by ancestry (African/European), performs association testing with ancestry-specific effects, and finds optimal weights.

## STRUCTURE

```
./
├── src/                    # R pipeline scripts
│   ├── step1_vcf_split_by_ancestry.R
│   ├── step2_association_detection.R
│   └── step3_weight_finding.R
├── test/                   # Test data + shell runners
│   ├── *.sh                # Pipeline step runners
│   ├── data/               # Input test data
│   ├── output/             # Pipeline output + cache
│   └── src/data_sim.R      # Data simulation
└── README.md               # Full documentation
```

## WHERE TO LOOK

| Task | Location | Notes |
|------|----------|-------|
| Run full pipeline | test/*.sh | Shell scripts for each step |
| Modify Step 1 | src/step1_vcf_split_by_ancestry.R | VCF splitting logic (501 lines) |
| Modify Step 2 | src/step2_association_detection.R | Association testing |
| Modify Step 3 | src/step3_weight_finding.R | Weight finding |
| Input format docs | README.md | Detailed pipeline documentation |
| Test data | test/data/ | GDS, VCF, gene group files |

## CODE MAP

| Symbol | Type | Location | Refs | Role |
|--------|------|----------|------|------|
| opt_parser | OptionParser | step1 | 1 | CLI arg parsing |
| parse_args | function | step1 | 1 | Parse CLI args |
| data.table | library | step1,2 | - | Data manipulation |
| vcfR | library | step1 | - | VCF parsing |
| SeqArray | library | step2 | - | GDS operations |

## CONVENTIONS

- **CLI**: Uses optparse for Rscript argument parsing
- **Paths**: Absolute paths hardcoded as defaults (see README for /scratch/g/pauer/Yu/...)
- **Output**: Writes to out_path directory, creates if needed
- **Parallel**: Uses doParallel for parallel processing

## ANTI-PATTERNS (THIS PROJECT)

- No unit tests (test/ only contains integration/data files)
- No CI/CD configuration
- No package management (no DESCRIPTION, no renv)
- Hardcoded absolute paths in default arguments
- No error handling for missing input files

## UNIQUE STYLES

- Bioinformatics workflow with ancestry-specific genetics
- Uses PLINK .bed/.bim/.fam files for ancestry matrix
- VCF → GDS conversion for efficient genotype storage
- Cauchy-weighted combination for multi-ancestry p-values

## COMMANDS

```bash
# Step 1: Split VCF by ancestry
Rscript src/step1_vcf_split_by_ancestry.R --bed <plink.bed> --bim <plink.bim> --fam <plink.fam> --vcf_path <vcf> --out_path <dir> --chr_id <chr>

# Step 2: Run association
Rscript src/step2_association_detection.R --african_gds <gds> --european_gds <gds> --data_file <tsv> --gene_group_file <tsv> --response_type <type> --kinship_rds <rds> --out_file <rds>

# Step 3: Find weights
Rscript src/step3_weight_finding.R --pt <pt_matrix.tsv> --gene_group <gene_group.tsv> --out_file <tsv>
```

## NOTES

- Requires R > 4.3.3 with: data.table, dplyr, tidyr, doParallel, vcfR, stringr, SeqArray, SeqVarTools, snpStats, GMMAT
- Requires: bcftools > 1.20, bgzip, tabix
- Input VCF should be chromosome-split for memory efficiency
- Kinship matrix requires column named "id"
