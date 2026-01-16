#!/bin/bash
#SBATCH --job-name=step1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=6gb
#SBATCH --time=7:00:00
#SBATCH --account=pauer
#SBATCH --partition=normal
#SBATCH --output=%x.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ywang@mcw.edu

ml R/4.3.3

Rscript /scratch/g/pauer/Yu/Tractor-RVA/src/step1_vcf_split_by_ancestry.R \
  --bed /scratch/g/pauer/Yu/smmat/rawdata/rfmix_merged.bed \
  --bim /scratch/g/pauer/Yu/smmat/rawdata/rfmix_merged.bim \
  --fam /scratch/g/pauer/Yu/smmat/rawdata/rfmix_merged.fam \
  --vcf_path /scratch/g/pauer/Yu/smmat/src/python_split/output/chr15.maf0.01.intersected.vcf.gz \
  --out_path /scratch/g/pauer/Yu/Tractor-RVA/test/output \
  --chr_id 15