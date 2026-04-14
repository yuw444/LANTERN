# LANTERN

## 1. Background

* Homo sapiens(diploid), the genome regions on homogolous pair of chromsomes could originate from different ancestries. 
* We currently focus on two ancestries, **african(AFR)** and **european(EUR)**.
* Similarly, two alleles of each variant could come from **pure african(03)**, **mixed ancestry(02)**, and **pure european(01)** in each subject.
* The B allele count could be **1/1**, **0/1**, **1/0**, and **0/0**, as the phase information is unkown.
* We assume that the effect of allele could be ancestry specific(**pure african**, **pure european**) or not(**mixed ancestry**).


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
- Entries: ancestry-of-origin code for the track in that subject (03 = pure African, 02 = mixed, 01 = pure European).
- Interpretation example: for chr1:1,000-2,000, Subject_S1 has African-origin ancestry (03), Subject_S3 has European-origin (01), and Subject_S2 has mixed ancestry (02). This matrix is used to assign ancestry-specific allele effects or to filter variants by ancestry background.


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

| Ancestry            | Genotype               | Frequency | Condition              | $\mathbf x_{AFR,AFR}$ | $\mathbf x_{EUR,EUR}$ |
|---------------------|------------------------|-----------|------------------------|:-----------------:|:-----------------:|
| **AFR,AFR (03)**    | 1/1 (2)                | $N_1$     | (pt = 3 & gt = 2)      | 2                | 0                |
|                     | 1/0 (1), 0/1 (1)       | $N_2$     | (pt = 3 & gt = 1)      | 1                | 0                |
|                     | 0/0 (0)                | $N_3$     | (pt = 3 & gt = 0)      | 0                | 0                |
| **AFR,EUR (02)**    | 1/1 (2)                | $N_4$     | (pt = 2 & gt = 2)      | 1                | 1                |
|                     | 1/0 (1), 0/1 (1)       | $N_5$     | (pt = 2 & gt = 1)      | $p_1$            | $p_2$            |
|                     | 0/0 (0)                | $N_6$     | (pt = 2 & gt = 0)      | 0                | 0                |
| **EUR,EUR (01)**    | 1/1 (2)                | $N_7$     | (pt = 1 & gt = 2)      | 0                | 2                |
|                     | 1/0 (1), 0/1 (1)       | $N_8$     | (pt = 1 & gt = 1)      | 0                | 1                |
|                     | 0/0 (0)                | $N_9$     | (pt = 1 & gt = 0)      | 0                | 0                |

To obtain the $p_1$ and $p_2$, we use the following formulas. 

$p_1 = \frac{2*N_1 + N_2 + N_4}{2*N_1 + N_2 + 2*N_4 + 2*N_7 + N_8}$

$p_2 = \frac{N_4 + 2*N_7 + N_8} {2*N_1 + N_2 + 2*N_4 + 2*N_7 + N_8}$.

Special case: When singleton SNP has an ancestry of AFR,EUR, i.e. , $N_5 = 1$ and $N_1..N_9 = 0$, then $p_1 = p_2 = 0.5$.


### Association Detection

* Pure African association: $\mathbf x_{AFR,AFR} \sim \text{Response} \to p_{a}$
* Pure European association: $\mathbf x_{EUR,EUR} \sim \text{Response} \to p_{e}$
* Observation association: $\mathbf x_{obs} \sim \text{Response} \to p_{obs}$
* Cauchy-Weighted combination between $p_a$ and $p_e$, result $p_c$.


### 4. Pipeline

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
  * plink file that contains ancestry matrix
  * vcf file that contains variant matrix, <u>highly suggest do it by chrosomes for the sake of memory</u>.
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
  * clinical covariate tsv(csv, rds) file `data_file` with header, columns are id, response, var1, var2, ...
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
