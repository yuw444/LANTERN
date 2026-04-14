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


### Association Detection

* For testing association with rare alleles on AFR haplotypes we obtain p-value: $\mathbf p_{a}$
* For testing association with rare alleles on EUR haplotypes we obtain p-value: $\mathbf p_{e}$
* For testing association with rare alleles ignoring ancestry we obtain p-value: Observation association: $\mathbf p_{obs}$
* The Cauchy-weighted combination between $p_a$ and $p_e$ results in $p_c$


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
