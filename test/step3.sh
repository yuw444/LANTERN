#!/bin/bash
#SBATCH --job-name=step3
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

Rscript /scratch/g/pauer/Yu/Tractor-RVA/src/step3_weight_finding.R \
  --pt /scratch/g/pauer/Yu/Tractor-RVA/test/output/cache/pt_matrix_chr15.tsv \
  --gene_group /scratch/g/pauer/Yu/Tractor-RVA/test/data/genes_oi_group.tsv \
  --out_file /scratch/g/pauer/Yu/Tractor-RVA/test/output/weights_genes_oi.tsv
