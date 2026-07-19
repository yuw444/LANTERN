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

REPO=/scratch/g/pauer/Yu/LANTERN
cd "$REPO"
export PATH="$HOME/.pixi/bin:$PATH"

# --msp_path below is a placeholder following the old --bed/--bim/--fam
# trio's naming convention -- update to the actual RFMix .msp.tsv path
# before running (see src/AGENTS.md).
pixi run Rscript src/step1_vcf_split_by_ancestry.R \
  --vcf_path /scratch/g/pauer/Yu/smmat/src/python_split/output/chr15.maf0.01.intersected.vcf.gz \
  --msp_path /scratch/g/pauer/Yu/smmat/rawdata/rfmix_merged.msp.tsv \
  --out_path test/output \
  --chr_id 15
