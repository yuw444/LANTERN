# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Branch Warning

**Always work on `r-package-c-backend`** — `origin/main` is a reverted branch pointing at a different codebase (scRNA-seq / data portal). Run `git checkout r-package-c-backend` if you're ever on the wrong branch.

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
lantern/          R package (the target product)
  R/lantern.R     exported R wrappers; calls C via .Call()
  R/workflow.R    high-level run_ancestry_pipeline()
  src/ancestry.c  ALL core algorithms: p1/p2 split, phased split, bed reader, vcf writer
  src/init.c      .Call registration — must stay in sync with ancestry.h
src/              legacy standalone pipeline (pre-package, kept for reference)
  step1_*.R       PLINK+VCF → ancestry-split GDS (uses snpStats + bcftools shellouts)
  step2_*.R       GMMAT/SMMAT association tests
  step3_*.R       Cauchy-combination weight finding
test/             SLURM runners + test data (paths are hardcoded to /scratch, don't run directly)
simulation/       earlier power/type-I simulations
raw/              real JHS cohort data (chr19, ~3000 samples)
```

**Main objective**: wire the `read_bed_file_c` function in `lantern/src/ancestry.c:166-227` through `R/lantern.R` and replace the `snpStats::read.plink()` + bcftools shellouts in `src/step1_*.R` with package calls. See [AGENTS.md](AGENTS.md) for full detail.

## Core Algorithm

Ancestry codes: `1`=EUR/EUR, `2`=AFR/EUR (mixed het), `3`=AFR/AFR. Genotype: 0/1/2 alt-dosage.

```
p1 (AFR proportion) = (2*N1 + N2 + N4) / (sum(gt) - N5)
p2 (EUR proportion) = (N4 + 2*N7 + N8) / (sum(gt) - N5)
```

Singleton edge case (only mixed hets carry alt allele): `p1 = p2 = 0.5`. Implemented at `ancestry.c:93-96`.

**Do not "fix" the `* 0.9999` divergence** between legacy R and C without updating both; it's an intentional quirk.

## Key Conventions

- Ancestry codes are always integers 1/2/3 (never strings or factors).
- PT matrix shape: rows=samples, cols=variants. GT matrix shape: rows=variants, cols=samples. `run_ancestry_pipeline()` handles the transpose.
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

## Known Issues

- 5 pre-existing test failures in `lantern/tests/testthat/` (wrong expected values + monomorphic-filter logic); not environment bugs.
- `lantern/src/ancestry.o`, `init.o`, `lantern.so` are incorrectly committed; `.gitignore` has unresolved conflict markers on HEAD — the `r-package-c-backend` branch has a clean version.
- `~/.R/Makevars` may pin gcc 9.3 system-wide; `R.Makevars.local` (passed via `R_MAKEVARS_USER`) overrides this for pixi tasks.
- `src/step1_*.R` defaults to hardcoded `/scratch/g/pauer/Yu/smmat/...` paths — always pass CLI flags explicitly.
- PLINK `.bed` files are SNP-major; byte offset for variant `v` is `3 + v * ceil(N_samples / 4)` — fully random-access, no need to load the whole file.
