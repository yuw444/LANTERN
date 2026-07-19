# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Branch Warning

**Work on `main`.** The `package` branch (R package work) was merged into `main`; `r-package-c-backend` and `multi-ancestry` are also fully merged (0 commits ahead of `main`) and are stale. CI (`R-CMD-check.yaml`, `pkgdown.yaml`) now triggers off `main`.

## Build and Test

All R work requires pixi (system R 4.3 and gcc 9.3 are too old):

```bash
pixi run install-lantern          # R CMD INSTALL lantern/ (applies R.Makevars.local)
pixi run test                     # run testthat suite
pixi run install-cran-gmmat       # one-time: GMMAT is not in conda
pixi run r                        # launch radian REPL
```

End-to-end 1000G smoke test:
```bash
./download_1000g.sh               # ~250 MB; chr22 only; do NOT commit the output
pixi run Rscript prepare_1000g.R
pixi run Rscript test_lantern_full.R
```

Run a single testthat file:
```bash
pixi run Rscript -e 'devtools::test("lantern", filter="phased")'
```

## Architecture

LANTERN splits genotype dosages by inferred local ancestry (African vs European) to enable ancestry-stratified rare-variant burden tests.

```
lantern/                  R package (the target product)
  R/ancestry_split.R      Step 1: ancestry_split(), split_diploid()/split_haplotype() C wrappers, read_bed_file()
  R/write_ancestry_gds.R  Step 2: write_ancestry_gds(), write_dosage_gds()
  R/ancestry_smmat.R      Step 3: ancestry_smmat() (SMMAT + per-gene ancestry weights + Cauchy combination)
  R/utils.R               cauchy_combine()
  src/ancestry.c          ALL core algorithms: p1/p2 split, phased split, bed reader, vcf writer
  src/init.c              .Call registration — must stay in sync with ancestry.h
src/              legacy standalone SLURM-CLI pipeline (thin wrappers around the package)
  step1_*.R       VCF + RFMix MSP file → ancestry-split GDS, via lantern::ancestry_split()/write_ancestry_gds()
  step2_*.R       GMMAT/SMMAT association tests, via lantern::ancestry_smmat() (also computes Cauchy-combination weights)
test/             SLURM runners + test data (paths are hardcoded to /scratch, don't run directly)
simulation/       earlier power/type-I simulations
raw/              real JHS cohort data (chr19, ~3000 samples)
```

`src/step1_*.R` takes the same RFMix-MSP-based input as
`lantern::ancestry_split()` (Step 1) — it's just an `optparse` CLI wrapper
around it plus `write_ancestry_gds()`, for standalone SLURM-per-chromosome
use. See [src/AGENTS.md](src/AGENTS.md) for usage and flags.
`lantern::read_bed_file()` (PLINK-BED-encoded ancestry, a different input
format) is still a package function but no longer used by `src/step1_*.R`.

## Core Algorithm

Ancestry codes: `1`=EUR/EUR, `2`=AFR/EUR (mixed het), `3`=AFR/AFR. Genotype: 0/1/2 alt-dosage.

```
p1 (AFR proportion) = (2*N1 + N2 + N4) / (sum(gt) - N5)
p2 (EUR proportion) = (N4 + 2*N7 + N8) / (sum(gt) - N5)
```

Singleton edge case (only mixed hets carry alt allele): `p1 = p2 = 0.5`. Implemented at `ancestry.c:93-96`.

**Do not "fix" the `* 0.9999` divergence** between legacy R and C without updating both; it's an intentional quirk. (This applied to the *old* pure-R `src/step1_*.R` math, which no longer exists — that script is now a thin CLI wrapper around `lantern::ancestry_split(mode = "dosage")`, with no `0.9999` factor.)

## Key Conventions

- Ancestry codes are always integers 1/2/3 (never strings or factors).
- PT matrix shape: rows=samples, cols=variants. GT matrix shape: rows=variants, cols=samples. `ancestry_split_dosage()` handles the transpose.
- New C entry points must be registered in both `init.c` AND declared in `ancestry.h`.
- After editing `lantern/R/lantern.R` exported functions, re-run roxygen: `pixi run Rscript -e 'devtools::document("lantern")'`

## Simulation Infrastructure

All phased GDS outputs are pre-built at `simulation/gds/` (44271 variants × 3313 samples). Gene group: `simulation/data/genes_chr19_group.tsv` (1470 genes, 50336 variants).

**Type I error** (null phenotype, all 1470 genes):
```bash
# Submit 10-job SLURM array (20,000 iters total, 40 cores/job)
sbatch simulation/slurm/slurm_alpha.sh
# Output: simulation/output/alpha/alpha_task{01..10}.rds
# Re-run one failed task: sbatch --array=3 simulation/slurm/slurm_alpha.sh
```

Key design in `simulation/01_simulate_alpha.R`:
- Cholesky of `sigma_g * K + sigma_e * I` is pre-computed once in the parent process before `foreach` (avoids re-decomposing per iteration)
- Each iteration: draw null phenotype → `glmmkin` → 3× `SMMAT` (aa/ee/obs) → collect result
- Each task writes one RDS (`alpha_task<NN>.rds`) with 2,000 iteration results; resume-safe (skips if output exists)

Timing from `simulation/archive/small_batch_timing.R`: ~255s/iteration with all 1470 genes. At 200 SLURM CPUs (40 R workers × 10 jobs): estimated **6–8 hours total**.

## Vignettes

`lantern/vignettes/real-data-chr19.Rmd` needs the gitignored `raw/` JHS
cohort data, so it can never be built by CI (`R-CMD-check.yaml`,
`pkgdown.yaml`) or on a machine without that data. It follows the standard
`.Rmd.orig` pattern: edit `real-data-chr19.Rmd.orig` (unguarded, real code),
then re-knit it locally into the committed `real-data-chr19.Rmd` (static,
no live code, safe for CI):

```bash
pixi run Rscript -e '
  setwd("lantern/vignettes")
  knitr::knit("real-data-chr19.Rmd.orig", output = "real-data-chr19.Rmd")
'
```

Commit both files together. The other two vignettes (`split-intuition.Rmd`,
`toy-vcf-msp-example.Rmd`) use only synthetic data and build normally.

## Known Issues

- 5 pre-existing test failures in `lantern/tests/testthat/` (wrong expected values + monomorphic-filter logic); not environment bugs.
- `lantern/src/ancestry.o`, `init.o`, `lantern.so` are incorrectly committed; `.gitignore` has unresolved conflict markers (`<<<<<<<`/`=======`/`>>>>>>>`) still present on `main` — needs manual resolution.
- `~/.R/Makevars` may pin gcc 9.3 system-wide; `R.Makevars.local` (passed via `R_MAKEVARS_USER`) overrides this for pixi tasks.
- `src/step1_*.R` defaults to hardcoded `/scratch/g/pauer/Yu/smmat/...` paths — always pass CLI flags explicitly.
- PLINK `.bed` files are SNP-major; byte offset for variant `v` is `3 + v * ceil(N_samples / 4)` — fully random-access, no need to load the whole file.
