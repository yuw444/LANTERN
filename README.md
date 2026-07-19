# LANTERN <img src="lantern/man/figures/logo.png" align="right" height="139" alt="LANTERN logo" />

**L**everaging Local **AN**cestry **T**racts to **E**nhance **R**are-Varia**N**t Aggregate Association Testing

## 1. Background

* LANTERN is a method for conducting aggregate rare-variant association tests using inferred local ancestry as additional information. 
* For simplicity, in the following we demonstrate LANTERN on two ancestries, **african(AFR)** and **european(EUR)**.
* LANTERN can be implemented on up to 5 ancestries. 



## 2. Data

### Ancestry Matrix

| Genome track        | Subject_S1(AN) | Subject_S2(AN) | Subject_S3(AN) |
|---------------------|:------------:|:------------:|:------------:|
| chr1:1,000-2,000    | 03         | 02         | 01         |
| chr1:2,000-3,000    | 02         | 02         | 01         |
| chr2:5,000-6,000    | 03         | 03         | 02         |

Explanation:
- Track: Half-open [start, end) (end excluded)
- Rows: contiguous genome tracks or windows used in analysis
- Columns: individual subjects in the cohort.
- Entries: Inferred ancestry (03 = 2 AFR chromosomes, 02 = 1 AFR and 1 EUR chromosome , 01 = 2 EUR chromosomes).

### Variant Matrix

| Variant (CHROM:POS REF>ALT) | Subject_S1 (GT) | Subject_S2 (GT) | Subject_S3 (GT) |
|-----------------------------|:-----------------:|:-----------------:|:-----------------:|
| chr1:1,234 A>G              | 0/1             | 1/1             | 0/0             |
| chr1:2,345 T>C              | 1/0             | 0/1             | 0/0             |
| chr2:5,678 G>A              | 0/0             | 0/1             | 1/1             |

Notes:
- Rows: variant loci (chromosome:position, REF>ALT).
- Columns: per-sample genotype (GT) field from VCF.
- GT meanings: 1/1 = homozygous alternate, 0/1 or 1/0 = heterozygous (phase unknown), 0/0 = homozygous reference.

### Variant Matrix with Ancestry Annotation

Combine GT with the ancestry matrix to map alternate alleles to ancestry-of-origin where possible.

| Variant (CHROM:POS REF>ALT) | Subject_S1 (GT:AN) | Subject_S2 (GT:AN) | Subject_S3 (GT:AN) |
|-----------------------------|:--------------------:|:--------------------:|:--------------------:|
| chr1:1,234 A>G              | 0/1:01             | 1/1:02             | 0/0:02             |
| chr1:2,345 T>C              | 1/0:03             | 0/1:03             | 0/0:03             |
| chr2:5,678 G>A              | 0/0:02             | 0/1:01             | 1/1:01             |

## 3. Method

### Ancestry Specific Variant Matrix Split

#### 2-Ancestry Case (AFR + EUR)

All counts are per variant.

| Ancestry         | Genotype         | Freq  | Condition           | $\mathbf{x}_{AFR}$ | $\mathbf{x}_{EUR}$ |
|------------------|------------------|-------|---------------------|:------------------:|:------------------:|
| **AFR/AFR (03)** | 1/1 (2)          | $N_1$ | pt = 3 & gt = 2     | 2                  | 0                  |
|                  | 0/1 (1)          | $N_2$ | pt = 3 & gt = 1     | 1                  | 0                  |
|                  | 0/0 (0)          | $N_3$ | pt = 3 & gt = 0     | 0                  | 0                  |
| **AFR/EUR (02)** | 1/1 (2)          | $N_4$ | pt = 2 & gt = 2     | 1                  | 1                  |
|                  | 0/1 (1)          | $N_5$ | pt = 2 & gt = 1     | $p_1$              | $p_2$              |
|                  | 0/0 (0)          | $N_6$ | pt = 2 & gt = 0     | 0                  | 0                  |
| **EUR/EUR (01)** | 1/1 (2)          | $N_7$ | pt = 1 & gt = 2     | 0                  | 2                  |
|                  | 0/1 (1)          | $N_8$ | pt = 1 & gt = 1     | 0                  | 1                  |
|                  | 0/0 (0)          | $N_9$ | pt = 1 & gt = 0     | 0                  | 0                  |

$$p_1 = \frac{2N_1 + N_2 + N_4}{D}, \qquad p_2 = \frac{N_4 + 2N_7 + N_8}{D}, \qquad D = 2N_1 + N_2 + 2N_4 + 2N_7 + N_8$$

Special case: if the only alt carriers are AFR/EUR heterozygotes ($N_5 \ge 1$, all other $N = 0$), then $p_1 = p_2 = 0.5$.

---

#### 3-Ancestry Case (AFR + EUR + NAT)

With a third ancestry (NAT = Native/other), there are 3 pure types and 3 mixed pairs.  
All counts are per variant; 0/0 rows (gt = 0) contribute nothing and are omitted for brevity.

| Ancestry         | Genotype | Freq       | $\mathbf{x}_{AFR}$               | $\mathbf{x}_{EUR}$               | $\mathbf{x}_{NAT}$               |
|------------------|----------|------------|:--------------------------------:|:--------------------------------:|:--------------------------------:|
| **AFR/AFR**      | 1/1 (2)  | $N_1$      | 2                                | 0                                | 0                                |
|                  | 0/1 (1)  | $N_2$      | 1                                | 0                                | 0                                |
| **EUR/EUR**      | 1/1 (2)  | $N_3$      | 0                                | 2                                | 0                                |
|                  | 0/1 (1)  | $N_4$      | 0                                | 1                                | 0                                |
| **NAT/NAT**      | 1/1 (2)  | $N_5$      | 0                                | 0                                | 2                                |
|                  | 0/1 (1)  | $N_6$      | 0                                | 0                                | 1                                |
| **AFR/EUR**      | 1/1 (2)  | $N_7$      | 1                                | 1                                | 0                                |
|                  | 0/1 (1)  | $N_8$      | $\dfrac{p_1}{p_1+p_2}$           | $\dfrac{p_2}{p_1+p_2}$           | 0                                |
| **AFR/NAT**      | 1/1 (2)  | $N_9$      | 1                                | 0                                | 1                                |
|                  | 0/1 (1)  | $N_{10}$   | $\dfrac{p_1}{p_1+p_3}$           | 0                                | $\dfrac{p_3}{p_1+p_3}$           |
| **EUR/NAT**      | 1/1 (2)  | $N_{11}$   | 0                                | 1                                | 1                                |
|                  | 0/1 (1)  | $N_{12}$   | 0                                | $\dfrac{p_2}{p_2+p_3}$           | $\dfrac{p_3}{p_2+p_3}$           |

$$D = 2N_1 + N_2 + 2N_3 + N_4 + 2N_5 + N_6 + 2N_7 + 2N_9 + 2N_{11}$$

$$p_1 = \frac{2N_1 + N_2 + N_7 + N_9}{D}, \quad p_2 = \frac{2N_3 + N_4 + N_7 + N_{11}}{D}, \quad p_3 = \frac{2N_5 + N_6 + N_9 + N_{11}}{D}$$

By construction $p_1 + p_2 + p_3 = 1$.  Each mixed hom-alt individual (e.g. $N_7$, AFR/EUR 1/1) contributes one alt allele to each parent population's numerator and two to the shared denominator $D$.  Ambiguous hets in each mixed pair are then split by the **pairwise ratio** of the two parent populations' proportions.

---

#### General K-Ancestry Case

With $K$ ancestries labelled $1, \ldots, K$ there are $K$ pure types and $\binom{K}{2}$ mixed pairs.

**Counts per variant:**

| Symbol | Meaning |
|--------|---------|
| $N_{kk}^{(2)}$ | count of pure pop-$k$ samples with gt = 2 |
| $N_{kk}^{(1)}$ | count of pure pop-$k$ samples with gt = 1 |
| $N_{ij}^{(2)}$ | count of pop-$i$/pop-$j$ mixed samples with gt = 2, $i < j$ |
| $N_{ij}^{(1)}$ | count of pop-$i$/pop-$j$ mixed samples with gt = 1 (ambiguous het), $i < j$ |

**Splitting rules:**

| Ancestry type | gt | $\mathbf{x}_k$ |
|---|---|---|
| pure pop-$k$ | 2 | 2 for $k$, 0 for all others |
| pure pop-$k$ | 1 | 1 for $k$, 0 for all others |
| mixed pop-$i$/pop-$j$ | 2 | 1 for $i$, 1 for $j$, 0 for all others |
| mixed pop-$i$/pop-$j$ | 1 | $\dfrac{p_i}{p_i + p_j}$ for $i$, $\dfrac{p_j}{p_i + p_j}$ for $j$, 0 for others |
| any | 0 | 0 for all |

**Denominator** (total unambiguous alt alleles):

$$D = \sum_{k=1}^{K} \left(2N_{kk}^{(2)} + N_{kk}^{(1)}\right) + 2\sum_{i < j} N_{ij}^{(2)}$$

**Population proportions** ($k = 1, \ldots, K$):

$$p_k = \frac{2N_{kk}^{(2)} + N_{kk}^{(1)} + \displaystyle\sum_{j \neq k} N_{\langle kj \rangle}^{(2)}}{D}$$

where $N_{\langle kj \rangle}^{(2)} = N_{\min(k,j),\max(k,j)}^{(2)}$.

It follows that $\sum_{k=1}^{K} p_k = 1$.

**Interpretation:** $p_k$ estimates the fraction of population-level alt alleles attributable to ancestry $k$, using only observations with unambiguous ancestry.  The mixed-pair hom-alts ($N_{ij}^{(2)}$) contribute one allele to each parent population.  Ambiguous hets in pair $(i,j)$ are then split by the conditional pairwise ratio $p_i / (p_i + p_j)$, i.e. the AFR fraction among AFR and EUR only for an AFR/EUR het, regardless of how many other populations exist.

**Singleton special case:** if the only alt carriers in a mixed pair $(i,j)$ are heterozygotes and all pure-ancestry and hom-alt counts are zero, set $p_i = p_j = 0.5$ for that pair.

**Number of ancestry types by K:**

| K | Pure types | Mixed pairs | Total types |
|---|-----------|-------------|-------------|
| 2 | 2 | 1 | 3 |
| 3 | 3 | 3 | 6 |
| 4 | 4 | 6 | 10 |
| 5 | 5 | 10 | 15 |

---

### Association Detection

* Per-ancestry association: $\mathbf{x}_k \sim \text{Response} \to p_k^{\text{assoc}}$ for each $k = 1,\ldots,K$
* Observed association: $\mathbf{x}_{obs} \sim \text{Response} \to p_{obs}$
* Cauchy combination of per-ancestry p-values $\to p_c$

For the 2-ancestry case: pure African $\to p_a$, pure European $\to p_e$, Cauchy combination $\to p_c$.


### 4. R Package Installation

The `lantern` R package is on GitHub. It has two categories of dependencies:

- **Bioconductor R packages** — [SeqArray](https://bioconductor.org/packages/SeqArray/) and [SeqVarTools](https://bioconductor.org/packages/SeqVarTools/) are not on CRAN; `devtools::install_github()` will not find them automatically.
- **System binary** — the phased-haplotype pipeline shells out to [`bcftools`](https://samtools.github.io/bcftools/) (≥ 1.10) for VCF reading and sample subsetting. Install it via your system package manager (`brew install bcftools` on macOS, `apt install bcftools` on Debian/Ubuntu) before using phased-mode functions.

**Recommended installation** (handles Bioconductor dependencies):

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("yuw444/LANTERN/lantern")
```

**Alternative**: pre-install Bioconductor dependencies first, then use devtools:

```r
BiocManager::install(c("SeqArray", "SeqVarTools"))
devtools::install_github("yuw444/LANTERN", subdir = "lantern")
```

### 5. Pipeline (legacy `./src`)

Two SLURM-friendly CLI scripts (`step1_vcf_split_by_ancestry.R`,
`step2_association_detection.R`) that wrap the `lantern` R package —
`optparse` argument parsing and file I/O around
`lantern::ancestry_split()` / `write_ancestry_gds()` / `ancestry_smmat()`.
Use these if you want each step as its own SLURM job; if you're scripting
in R directly, just call the package functions yourself (see
`vignette("split-intuition")`, `vignette("toy-vcf-msp-example")`,
`vignette("real-data-chr19")`). There is no separate weight-finding
script — per-gene ancestry weights for the Cauchy combination are
computed automatically inside Step 2.

#### Prerequisites

* The `lantern` R package installed (Section 4 above)
* `GMMAT` (`pixi run install-cran-gmmat` — not on conda/CRAN mirrors used by pixi)
* `bcftools` ≥ 1.20 on `PATH`
* An RFMix2-produced `.msp.tsv` file — see `vignette("generating-local-ancestry")`
  if you need to produce one from a phased VCF

No installation of `./src` itself is needed — just run the scripts with `Rscript`.

#### Step 1 — split a VCF by local ancestry

```bash
Rscript src/step1_vcf_split_by_ancestry.R \
  --vcf_path  cohort.phased.vcf.gz \
  --msp_path  cohort.msp.tsv \
  --out_path  out/ \
  --chr_id    22 \
  --mode      dosage
```

| Flag | Required | Meaning |
|------|----------|---------|
| `--vcf_path` / `-i` | yes | Phased VCF/BCF (bgzipped). May contain other chromosomes besides `--chr_id` — see "Multi-chromosome input" below. |
| `--msp_path` | yes | RFMix `.msp.tsv` (plain text or gzipped). May also contain other chromosomes' tracts; only `--chr_id`'s are used. |
| `--out_path` | yes | Output directory, created if missing. |
| `--chr_id` | yes | Chromosome to process, e.g. `22` or `chr22` — must identify the same chromosome in both `--vcf_path` and `--msp_path` (a `chr` prefix mismatch between the two is handled automatically). |
| `--mode` | no (default `dosage`) | `dosage` = proportional p1/p2 split (unphased-friendly) or `haplotype` = deterministic per-haplotype split (needs a truly phased VCF). See `vignette("split-intuition")` for the difference. |

**Multi-chromosome input**: `--vcf_path`/`--msp_path` don't need to be
pre-split per chromosome — `--chr_id` selects one chromosome out of a
genome-wide file. When possible (the VCF's contigs can be resolved),
`bcftools query` itself is restricted to `--chr_id`, so the other
chromosomes' genotypes are never read into R — this keeps memory bounded
to one chromosome's variant count even if you hand it a whole-genome VCF.
Run one Step 1 job per chromosome (e.g. a SLURM array over `--chr_id`)
rather than combining chromosomes into a single call.

**Output** (in `--out_path`):

| File | Contents |
|------|----------|
| `<POP>.gds` (one per population named in the MSP header, e.g. `AFR.gds`, `EUR.gds`) | Ancestry-specific dosage GDS, ready for `GMMAT::SMMAT()` (`is.dosage = TRUE`) |
| `split_meta_chr<chr_id>.rds` | `list(gds_paths, variant_info, ancestry_counts, sample_ids, mode, chr_id)` — bundles the GDS paths above plus per-variant ancestry counts. This whole file is Step 2's `--split_meta` input. |

#### Step 2 — ancestry-stratified association testing

```bash
Rscript src/step2_association_detection.R \
  --split_meta      out/split_meta_chr22.rds \
  --data_file       pheno.tsv \
  --gene_group_file genes.tsv \
  --kinship_rds     kinship.rds \
  --response_type   continuous \
  --out_file        results.rds \
  --ncores          8
```

| Flag | Required | Meaning |
|------|----------|---------|
| `--split_meta` | yes | `split_meta_chr<chr_id>.rds` from Step 1 — supplies the GDS paths for every population plus per-variant ancestry counts (used to weight the Cauchy combination; equal weights if missing). |
| `--data_file` | yes | Phenotype file, csv/tsv/rds — see format below. |
| `--gene_group_file` | yes | Gene-group file passed straight to `GMMAT::SMMAT()` — see format below. |
| `--kinship_rds` | yes | Kinship matrix RDS — see format below. |
| `--response_type` | no (default `continuous`) | `continuous` (`gaussian`), `binary` (`binomial`), or `count` (`poisson`). |
| `--out_file` | yes | Output RDS path. |
| `--ncores` | no | Cores for `GMMAT::SMMAT()`. Defaults to `$SLURM_CPUS_PER_TASK` when set (i.e. automatically picks up `--cpus-per-task` in a SLURM job), else `1`. Pass explicitly to override. |

**Input file formats**:

* **Phenotype** (`--data_file`, header row required): column 1 must be
  `id`, column 2 is the response, remaining columns are covariates (all
  used — the formula is built as `<col2> ~ <col3> + <col4> + ...`).
  ```
  id	response	age	sex	PC1	PC2
  sample_001	1	63	M	-0.012	0.034
  sample_002	0	57	F	0.104	-0.021
  ```
* **Gene group** (`--gene_group_file`, no header): `gene, chr, pos, ref,
  alt, weight`, one row per variant per gene.
  ```
  GOLGA6L22	15	22460882	G	T	1
  GOLGA6L22	15	22462401	G	C	1
  HERC2P2	15	22554572	G	A	1
  ```
* **Kinship** (`--kinship_rds`): an RDS of a square numeric matrix with
  row and column names equal to sample IDs.

**Output** (`--out_file`): an RDS of `ancestry_smmat()`'s return value —
`list(results, smmat_results)`:
* `results` — data.frame, one row per gene: `gene`, one `p_<POP>` and one
  `w_<POP>` column per population (e.g. `p_AFR`, `w_AFR`, `p_EUR`,
  `w_EUR`), and `p_cauchy` (the combined p-value).
* `smmat_results` — named list (one entry per population) of the raw
  `GMMAT::SMMAT()` output data.frames, for anyone who needs the full detail.
