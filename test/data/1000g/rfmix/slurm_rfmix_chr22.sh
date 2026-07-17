#!/bin/bash
#SBATCH --job-name=rfmix_chr22
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=2gb
#SBATCH --time=02:00:00
#SBATCH --account=pauer
#SBATCH --partition=normal
#SBATCH --output=test/data/1000g/rfmix/logs/%x_%j.out
#SBATCH --error=test/data/1000g/rfmix/logs/%x_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ywang@mcw.edu

# Local ancestry inference on 1000G chr22: 157 admixed query samples
# (ASW+ACB) against a 1,007-sample AFR/EUR reference panel (~1M variants,
# unpruned). See vignette("generating-local-ancestry") for how these
# inputs were built.
#
# Submit: sbatch test/data/1000g/rfmix/slurm_rfmix_chr22.sh
# Output: test/data/1000g/rfmix/output/chr22.msp.tsv (+ .fb.tsv, .rfmix.Q, .sis.tsv)

set -euo pipefail

REPO=/scratch/g/pauer/Yu/LANTERN
cd "$REPO"

mkdir -p test/data/1000g/rfmix/output test/data/1000g/rfmix/logs

export PATH="$HOME/.pixi/bin:$PATH"

echo "Start: $(date)"
t0=$(date +%s)

pixi run rfmix \
  -f test/data/1000g/rfmix/query.vcf.gz \
  -r test/data/1000g/rfmix/reference.vcf.gz \
  -m test/data/1000g/rfmix/sample_map.tsv \
  -g test/data/1000g/genetic_map/chr22_rfmix.gmap \
  -o test/data/1000g/rfmix/output/chr22 \
  --chromosome=22 \
  --n-threads="${SLURM_CPUS_PER_TASK}"

t1=$(date +%s)
echo "End: $(date)"
echo "RFMIX_ELAPSED_SECONDS: $((t1 - t0))"
