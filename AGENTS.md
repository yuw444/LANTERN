# LANTERN

R/C package for **ancestry-specific rare variant association analysis** (AFR vs EUR local-ancestry split). The C backend (`lantern/src/ancestry.c`) is the source of truth for the p1/p2 mixed-ancestry splitting algorithm; the standalone R scripts in `src/` are the original (pre-package) pipeline.

## STRUCTURE

```
./
├── lantern/                       # R package (current focus of development)
│   ├── DESCRIPTION                # Package: lantern, v0.1.0
│   ├── NAMESPACE                  # roxygen-generated, useDynLib(lantern)
│   ├── R/
│   │   ├── lantern.R              # C wrappers + exported count_/split_/write_/subset_
│   │   └── workflow.R             # run_ancestry_pipeline() (overlap handling, VCF I/O)
│   ├── src/
│   │   ├── ancestry.c             # C backend (split_by_ancestry, count_ancestry_codes)
│   │   ├── ancestry.h             # C function declarations
│   │   ├── init.c                 # R_CallMethodDef table
│   │   └── Makevars               # PKG_CFLAGS = -O3 -march=native -Wall
│   └── tests/testthat/
│       ├── test-functions.R       # Tests for count_ancestry_codes, split_by_ancestry
│       └── test-workflow.R        # Tests for run_ancestry_pipeline (overlap, monomorphic)
├── src/                           # Original standalone R pipeline (legacy)
│   ├── step1_vcf_split_by_ancestry.R    # VCF → AA/EE GDS + DS annotation
│   ├── step2_association_detection.R    # GMMAT burden tests via SMMAT
│   └── step3_weight_finding.R           # Gene-level AFR/EUR ancestry weights
├── test/
│   ├── src/data_sim.R             # Simulates phenotype + kinship for testing
│   ├── data/genes_oi_group.tsv    # gene_group example (no header: gene chr pos ref alt weight)
│   ├── data/1000g/                # 1000G chr22 test data (small subset tracked; large VCFs gitignored by *.gz)
│   ├── output/                    # Pipeline outputs (gitignored)
│   ├── step1.sh / step2.sh / step3.sh   # SLURM runners (paths hardcoded to Tractor-RVA)
│   └── *.out                      # SBATCH log files
├── download_1000g.sh              # Fetches 1000G chr22 VCF + panel, builds AFR/EUR subset
├── prepare_1000g.R                # Generates simulated PT matrix for 1000G subset
├── test_lantern_full.R            # Integration test: end-to-end pipeline on 1000G data
├── pixi.toml                      # Pixi manifest: R 4.5.3 + Bioconductor + CRAN deps + bcftools
├── R.Makevars.local               # Forces CC=x86_64-conda-linux-gnu-cc (bypasses ~/.R/Makevars gcc 9.3)
├── raw/                           # Local raw data (gitignored)
└── .vscode/                       # Points to pixi-bundled R + radian
```

## WHERE TO LOOK

| Task | Location | Notes |
|------|----------|-------|
| Modify splitting algorithm | `lantern/src/ancestry.c` | Source of truth for p1/p2; C, not R |
| Add exported R function | `lantern/R/lantern.R` or `workflow.R` | Then re-run roxygen |
| Register new C entry point | `lantern/src/init.c` (CallEntries) AND `ancestry.h` | Both must agree |
| Add testthat case | `lantern/tests/testthat/test-*.R` | Existing pattern: one file per module |
| End-to-end smoke test | `test_lantern_full.R` | Reads 1000G subset, exercises full API |
| Get 1000G test data | `download_1000g.sh` | Downloads to `test/data/1000g/` |
| Run legacy 3-step pipeline | `test/step{1,2,3}.sh` | SLURM only; hardcoded paths |
| Algorithm reference | `README.md` §3 | Math (p1, p2 formulas) |

## CORE ALGORITHM (p1/p2 split for mixed het)

Codes: `pt=1`=EUR/EUR, `pt=2`=AFR/EUR (mixed), `pt=3`=AFR/AFR. Genotypes: 0/1/2 alt-dosage.

For each variant, count per-row combinations N1..N8 (see `README.md` table and `ancestry.c:14-30`). Then:
```
p1 = (2*N1 + N2 + N4) / (sum(gt) - N5)
p2 = (N4 + 2*N7 + N8) / (sum(gt) - N5)
```
**Singleton edge case**: if `sum(gt) == N5` (only mixed hets carry alt alleles), `p1 = p2 = 0.5`. This is hard-coded in both C (`ancestry.c:93-96`) and R (`step1_vcf_split_by_ancestry.R:215-219`, `241-245`).

**`p * 0.9999` quirk**: the R script in `src/step1_*` multiplies p1/p2 by 0.9999 to avoid boundary issues. The C backend does NOT do this. If you change one, change both to stay consistent.

## R PACKAGE BUILD & TEST

All R commands should run through pixi (system R 4.3 + system gcc 9.3 are too old):

```bash
# Install the lantern package (R 4.5 + pixi gcc 15, C23-compatible)
pixi run install-lantern

# Run the testthat suite
pixi run test

# Or directly:
pixi run Rscript -e 'testthat::test_dir("lantern/tests/testthat")'
```

C build flags: `-O3 -march=native -Wall` (see `lantern/src/Makevars`). Compiles standalone; no external libs.

## END-TO-END 1000G TEST WORKFLOW

```bash
./download_1000g.sh                          # ~250 MB; chr22 only
Rscript prepare_1000g.R                      # builds simulated PT matrix
Rscript test_lantern_full.R                  # exercises count_/split_/run_ancestry_pipeline
```

AFR pops: `YRI,LWK,MSL,GWD,ESN,ACB,ASW`. EUR pops: `CEU,TSI,FIN,GBR,IBS`. PT matrix in `test/data/1000g/subset/pt_matrix.tsv` is **synthetic** (random ancestry flip, seed=42, 15% admixed) — not real local-ancestry inference.

## LEGACY PIPELINE (src/step*.R)

Three-step pipeline, run via SLURM:

| Step | Script | Input | Output |
|------|--------|-------|--------|
| 1 | `src/step1_vcf_split_by_ancestry.R` | PLINK .bed/.bim/.fam + VCF | `aa_ds.gds`, `ee_ds.gds`, `cache/pt_matrix_chr*.tsv` |
| 2 | `src/step2_association_detection.R` | AA/EE GDS + phenotype + kinship + gene_group | `results.rds` (SMMAT burden tests) |
| 3 | `src/step3_weight_finding.R` | `pt_matrix_chr*.tsv` + gene_group | `weights_*.tsv` (median AFR/EUR per gene) |

**Path bug in test/step*.sh**: all three reference `/scratch/g/pauer/Yu/Tractor-RVA/...` but the repo lives at `/scratch/g/pauer/Yu/LANTERN/`. Update before running, or symlink.

External tools required by step 1: `bcftools`, `bgzip`, `tabix`. Step 2 requires: `GMMAT`, `MASS`, `SeqArray`, `future`, `future.apply`. PLINK files are read with `snpStats::read.plink`.

## VS CODE / R ENVIRONMENT

**This project uses pixi.** R 4.5.3 + all required packages + `bcftools`/`htslib`/`samtools` are declared in `pixi.toml` and resolved into `.pixi/envs/default/`. `.vscode/settings.json` points `r.rpath.linux` to `${workspaceFolder}/.pixi/envs/default/bin/R` (the pixi-bundled R) and `r.rterm.linux` to `${workspaceFolder}/.pixi/envs/default/bin/radian` (the pixi-bundled radian).

**Required VS Code settings (R extension)**: `r.sessionWatcher: true` is mandatory for `r.plot.useHttpgd: true` to work. The pixi env pre-installs `r-httpgd` so plot viewing works out of the box.

**Tasks** (run with `pixi run <task>`):
- `pixi run install-lantern` — `R CMD INSTALL lantern/` (sets `R_MAKEVARS_USER=R.Makevars.local` to use the pixi-bundled gcc 15; this is needed because the user's `~/.R/Makevars` points to system gcc 9.3.0, which is too old for R 4.5's C23 enum syntax)
- `pixi run test` — `devtools::test("lantern")` (also uses `R.Makevars.local`)
- `pixi run install-cran-gmmat` — installs `GMMAT` from CRAN (GMMAT is not in conda-forge or bioconda; only on CRAN)
- `pixi run r` — launch pixi-bundled `radian`

**Known env quirk**: the system `~/.R/Makevars` hardcodes `CC=/hpc/apps/gcc/9.3.0/bin/gcc`. Without the project-local `R.Makevars.local` override, R 4.5 fails to compile any C code that includes `R_ext/Boolean.h` (C23 `enum :int` syntax). Always run R-related tasks through pixi (which sets `R_MAKEVARS_USER` automatically).

## CONVENTIONS

- **Ancestry codes**: integer 1/2/3 in R; same in C (`int` comparisons). Do not use strings or factor levels.
- **PT matrix orientation**: rows = samples, cols = variants/regions. GT matrix is the opposite (rows = variants). `run_ancestry_pipeline` handles this internally.
- **Output Rds**: step 2 writes `list(aa=..., ee=..., observed=...)` — three GMMAT/SMMAT results.
- **C entry points** (must agree across `init.c`, `ancestry.h`, and `ancestry.c`): `count_ancestry_codes_C`, `split_by_ancestry_C`, `read_bed_file_C`, `write_vcf_with_ancestry_C`, `subset_vcf_by_range_C`. **Header mismatch**: `ancestry.h` declares `read_ancestry_plink`/`read_pt_matrix`/`write_ancestry_vcf` but `init.c` registers different names (`read_bed_file`/`write_vcf_with_ancestry`/`subset_vcf_by_range`). The C source implements the init.c names; the header is out of sync. Only `count_ancestry_codes` and `split_by_ancestry` are wrapped/exported from R.
- **Hardcoded absolute paths**: `src/step1_vcf_split_by_ancestry.R` defaults all `--bed/--bim/--fam/--vcf_path/--out_path` to `/scratch/g/pauer/Yu/smmat/...` (wrong machine). Always pass explicit args via CLI.

## ANTI-PATTERNS

- **No CI/CD.** No `.github/workflows/`, no pre-commit, no R CMD check config.
- **No DESCRIPTION for the repo root.** Only `lantern/DESCRIPTION`. The repo is a package + scripts, not a meta-package.
- **Step 1 defaults are wrong for any machine other than the original.** Always pass `--bed/--bim/--fam/--vcf_path/--out_path/--chr_id` explicitly.
- **R and C implementations of split diverge slightly** (the `* 0.9999` factor). If you change one, change the other; tests live in C only.
- **Hardcoded email in SLURM scripts**: `ywang@mcw.edu`. Update before sharing runs.
- **monomorphic variants**: step 1 in R drops them; the C package's `run_ancestry_pipeline` also filters (returns `n_monomorphic_filtered` in overlap). Don't try to "fix" a 0-row output by suppressing the filter.

## DEPENDENCIES

All deps are managed by `pixi.toml` and installed into `.pixi/envs/default/` via `pixi install`. Do not use the system R 4.3 / system gcc 9.3 — the lantern C backend requires R 4.5's C23 support and the pixi-bundled gcc 15.

**R packages (lantern package)**: `data.table`, `dplyr`, `SeqArray`, `SeqVarTools`, `snpStats`. `testthat` (>= 3.0.0) for tests. `devtools`/`roxygen2` for development.

**R packages (legacy pipeline, src/)**: `optparse`, `data.table`, `dplyr`, `tidyr`, `doParallel`, `vcfR`, `stringr`, `SeqArray`, `SeqVarTools`, `snpStats`, `MASS`, `future`, `future.apply`. `GMMAT` is **not in conda** — install it once via `pixi run install-cran-gmmat` (downloads from CRAN).

**System tools** (all on bioconda/conda-forge): `bcftools` (>= 1.23), `htslib` (for `bgzip`/`tabix`), `samtools`.

**VS Code IDE**: `r-languageserver`, `r-httpgd`, `radian` — all pre-installed in the pixi env.

## NOTES

- Branch: `r-package-c-backend` is the active development branch. `main` is the older standalone-script version.
- The `raw/` directory holds local-only data and is gitignored.
- `test/output/` and `test/*.out` are gitignored — they're regenerated by runs.
- `seqOpen`/`seqClose` for GDS files: in step 2, only the AA GDS is opened to fetch sample IDs; both GDS files are then read by SMMAT directly.
- The "Cauchy-weighted combination" mentioned in legacy AGENTS.md is computed downstream of step 2, not inside the package.
