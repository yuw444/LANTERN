#!/bin/bash
# Download and prepare 1000 Genomes data for lantern testing
# Usage: ./download_1000g.sh [output_dir]

set -e

OUTPUT_DIR="${1:-test/data/1000g}"
SUBSET_DIR="${OUTPUT_DIR}/subset"
BASE_URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502"

echo "============================================"
echo "  Downloading 1000 Genomes Phase 3 Data"
echo "============================================"
echo ""

mkdir -p "${SUBSET_DIR}"

# Download chr22 VCF (smallest chromosome for testing)
echo "[1/4] Downloading chr22 VCF..."
if [ ! -f "${OUTPUT_DIR}/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz" ]; then
    curl -L -C - -o "${OUTPUT_DIR}/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz" \
        "${BASE_URL}/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
else
    echo "   Already exists, skipping..."
fi

# Download VCF index
echo "[2/4] Downloading chr22 VCF index..."
if [ ! -f "${OUTPUT_DIR}/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi" ]; then
    curl -L -C - -o "${OUTPUT_DIR}/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi" \
        "${BASE_URL}/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi"
else
    echo "   Already exists, skipping..."
fi

# Download sample panel
echo "[3/4] Downloading sample panel..."
if [ ! -f "${OUTPUT_DIR}/integrated_call_samples_v3.20130502.ALL.panel" ]; then
    curl -L -o "${OUTPUT_DIR}/integrated_call_samples_v3.20130502.ALL.panel" \
        "${BASE_URL}/integrated_call_samples_v3.20130502.ALL.panel"
else
    echo "   Already exists, skipping..."
fi

# Check for bcftools
echo "[4/4] Preparing subset..."
if command -v bcftools &> /dev/null; then
    echo "   bcftools found, creating AFR/EUR subset..."
    
    # Get AFR and EUR sample IDs
    AFR_POPS="YRI,LWK,MSL,GWD,ESN,ACB,ASW"
    EUR_POPS="CEU,TSI,FIN,GBR,IBS"
    
    # Create keep file for AFR+EUR
    grep -E ",(${AFR_POPS}|${EUR_POPS})," "${OUTPUT_DIR}/integrated_call_samples_v3.20130502.ALL.panel" | \
        cut -f1 > "${SUBSET_DIR}/keep_afr_eur.txt"
    
    # Subset VCF
    bcftools view -S "${SUBSET_DIR}/keep_afr_eur.txt" \
        -Oz -o "${SUBSET_DIR}/chr22_subset.vcf.gz" \
        "${OUTPUT_DIR}/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
    
    # Index
    bcftools index -t "${SUBSET_DIR}/chr22_subset.vcf.gz"
    
    echo "   Created: ${SUBSET_DIR}/chr22_subset.vcf.gz"
else
    echo "   bcftools not found, skipping VCF subset."
    echo "   Install bcftools to create subset: conda install -c bioconda bcftools"
fi

echo ""
echo "============================================"
echo "  Download Complete!"
echo "============================================"
echo ""
echo "Files downloaded:"
ls -lh "${OUTPUT_DIR}"/*.gz "${OUTPUT_DIR}"/*.tbi "${OUTPUT_DIR}"/*.panel 2>/dev/null || true
echo ""
echo "Next steps:"
echo "  1. Run R script to create PT matrix:"
echo "     Rscript prepare_1000g.R"
echo ""
echo "  2. Test lantern package:"
echo "     Rscript test_lantern_full.R"
