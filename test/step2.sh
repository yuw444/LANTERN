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

REPO=/scratch/g/pauer/Yu/LANTERN
cd "$REPO"
export PATH="$HOME/.pixi/bin:$PATH"

# --data_file/--gene_group_file/--kinship_rds below still point at the old
# Tractor-RVA test data location; update before running (see src/AGENTS.md).
pixi run Rscript src/step2_association_detection.R \
  --split_meta test/output/split_meta_chr15.rds \
  --data_file /scratch/g/pauer/Yu/Tractor-RVA/test/output/data_file.tsv \
  --gene_group_file /scratch/g/pauer/Yu/Tractor-RVA/test/data/genes_oi_group.tsv \
  --response_type count \
  --kinship_rds /scratch/g/pauer/Yu/Tractor-RVA/test/output/kinship.rds \
  --out_file test/output/results.rds \
  --ncores "${SLURM_CPUS_PER_TASK}"
