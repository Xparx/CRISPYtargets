#!/bin/bash

# ==============================================================================
# Crispy Pipeline Execution Script
# ==============================================================================

# 1. SETUP ENVIRONMENT
echo "Setting up environment..."
# Install python dependencies
pip install -r crispy/requirements.txt

# Download the required FlashFry JAR if not present
JAR_FILE="FlashFry-assembly-1.15.jar"
if [ ! -f "$JAR_FILE" ]; then
    echo "Downloading FlashFry JAR..."
    wget https://github.com/mckennalab/FlashFry/releases/download/1.15/$JAR_FILE
fi

# 2. DEFINE FILE PATHS
# Update these paths to point to your actual data files
GENES_CSV="data/target_genes.csv"
GFF_FILE="data/genome_annotation.gff3"
REF_FASTA="data/genome_reference.fa"
ACC_BED="data/accessomeKY.bed"
SNP_BED="data/cint.snps.bed"      # Optional
FILTER_BED="data/cint.bed"        # Optional (e.g., conserved regions)

# 3. RUN THE PIPELINE
echo "Starting Crispy Pipeline..."

python crispy/pipeline.py \
    --genes "$GENES_CSV" \
    --gff "$GFF_FILE" \
    --ref "$REF_FASTA" \
    --accessome "$ACC_BED" \
    --snps "$SNP_BED" \
    --filter-bed "$FILTER_BED" \
    --cutoff 0.5

echo "Pipeline finished. Check final_sgrnas.csv for results."
