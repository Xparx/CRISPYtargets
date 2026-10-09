#!/bin/bash

# ==============================================================================
# CRISPYtargets design pipeline example
# ==============================================================================

set -euo pipefail

# 1. SETUP ENVIRONMENT
echo "Setting up environment..."
# Install the package and its Python dependencies
pip install -e .

# Download the required FlashFry JAR if not present, and point the pipeline at it
JAR_FILE="FlashFry-assembly-1.15.jar"
if [ ! -f "tools/$JAR_FILE" ]; then
    echo "Downloading FlashFry JAR..."
    wget -P tools/ "https://github.com/mckennalab/FlashFry/releases/download/1.15/$JAR_FILE"
fi
export CRISPY_FLASHFRY_JAR="tools/$JAR_FILE"

# 2. DEFINE FILE PATHS
# Update these paths to point to your actual data files
GENES_CSV="data/target_genes.csv"
GFF_FILE="data/genome_annotation.gff3"
REF_FASTA="data/genome_reference.fa"
ACC_BED="data/accessomeKY.bed"
SNP_BED="data/cint.snps.bed"      # Optional
FILTER_BED="data/cint.bed"        # Optional (e.g., conserved regions)

# 3. RUN THE PIPELINE
echo "Starting the design pipeline..."

crispytargets design \
    --genes "$GENES_CSV" \
    --gff "$GFF_FILE" \
    --ref "$REF_FASTA" \
    --accessome "$ACC_BED" \
    --snps "$SNP_BED" \
    --filter-bed "$FILTER_BED" \
    --cutoff 0.5

echo "Pipeline finished. Check output/final_sgrnas.tsv for results."
