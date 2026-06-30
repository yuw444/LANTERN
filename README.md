# LANTERN (Leveraging local ANcestry Tracts to Enhance Rare variaNt aggregate associations)

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

* **Prerequiste**:
    * R > 4.3.3
      * data.table
      * dplyr
      * tidyr
      * doParallel
      * vcfR
      * stingr
      * SeqArray
      * SeqVarTools
      * snpStats
      * GMMAT
    * bcftools > 1.20
    * bgzip
    * tabix
* **Installation**:
  * No installation is needed, just copy the `./src` to your local machine

* **Input**: 
  * plink file that contains inferred local ancestry matrix
  * vcf file that contains variant matrix, <u>we highly suggest doing  this by chrosomes</u>.
  * kinship matrix rds `kinship_rds` file, with column name as `id`.
  * gene group tsv file without header `gene_group.tsv`, columns are gene, chr, pos, ref, alt, weight. For example,
```
GOLGA6L22	15	22460882	G	T	1
GOLGA6L22	15	22462401	G	C	1
GOLGA6L22	15	22464252	G	A	1
GOLGA6L22	15	22465059	G	A	1
GOLGA6L22	15	22465078	C	T	1
GOLGA6L22	15	22466304	A	G	1
HERC2P2	15	22554572	G	A	1
HERC2P2	15	22554572	G	A	1
```
  * covariate tsv(csv, rds) file `data_file` with header, columns are id, response, var1, var2, ...
```
id	response	age	sex	PC1	PC2
sample_001	1	63	M	-0.012	0.034
sample_002	0	57	F	0.104	-0.021
sample_003	0	45	M	-0.045	0.110
sample_004	1	52	F	0.008	-0.076
sample_005	0	39	M	0.212	0.003
sample_006	1	71	F	-0.131	-0.044
```
  
* **Step 1**: Split VCF by Ancestry
```
Rscript /path/to/step1_vcf_split_by_ancestry.R \
  --bed /path/to/plink/bed \
  --bim /path/to/plink/bim \
  --fam /path/to/plink/fam \
  --vcf_path /path/to/vcf \
  --out_path /path/to/output/dir \
  --chr_id 15
```

* **Step 2**: Model the Association
  * use the `african_gds`, `european_gds` generated from **Step1**
  * `response_type` could be one of *continous*, *binary*, or *count*
```
Rscript /path/to/step2_association_detection.R \
  --african_gds /path/to/gds \
  --european_gds /path/to/gds \
  --data_file /path/to/data/file \
  --gene_group_file /path/to/gene_group_file
  --response_type type \
  --kinship_rds /path/to/kinship/rds \
  --out_file /path/to/rds/file
```

* **Step 3**: Get the Weights for Ancestry
  * use the `pt_matrix_chr*.tsv` from **Step1**
```
Rscript /path/to/step3_weight_finding.R \
  --pt /path/to/cache/pt_matrix_chr*.tsv \
  --gene_group /path/to/gene_group.tsv \
  --out_file /path/to/tsv/file \
  --chr_id 15
```
