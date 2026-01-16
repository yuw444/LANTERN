#!/bin/bash
#SBATCH --job-name=step2
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

Rscript /scratch/g/pauer/Yu/Tractor-RVA/src/step2_association_detection.R \
  --african_gds /scratch/g/pauer/Yu/Tractor-RVA/test/output/aa_ds.gds \
  --european_gds /scratch/g/pauer/Yu/Tractor-RVA/test/output/ee_ds.gds \
  --data_file /scratch/g/pauer/Yu/Tractor-RVA/test/output/data_file.tsv \
  --gene_group_file /scratch/g/pauer/Yu/Tractor-RVA/test/data/genes_oi_group.tsv \
  --response_type count \
  --kinship_rds /scratch/g/pauer/Yu/Tractor-RVA/test/output/kinship.rds \
  --out_file /scratch/g/pauer/Yu/Tractor-RVA/test/output/results.rds