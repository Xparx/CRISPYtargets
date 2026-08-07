# Crispy Pipeline Documentation

## Overview
**Crispy** is a specialized bioinformatics pipeline designed to identify, score, and design CRISPR guide RNAs (sgRNAs). It automates the process of moving from a genomic reference (GFF3 and FASTA) to a ranked list of target sequences with ready-to-order oligonucleotide designs.

The pipeline integrates **FlashFry** (a high-performance CRISPR target discovery tool) with custom Python logic to handle gene hierarchy resolution, promoter region definition, accessibility scoring, and SNP-based entropy calculation.

It is inspired by and aims to, in part, replicate the output of the [CRISPRtargets](https://github.com/ChristiaenLab/CRISPRtargets) package but in python.

---

## How the Code Works

### 1. Target Preparation (`prepare_targets`)
The pipeline identifies potential target regions based on the GFF3 annotation:
*   **Exons/CDS**: Directly extracts all coding sequences.
*   **Promoters**: Automatically defines promoter regions by using transcript-level features as anchors for the **Transcription Start Site (TSS)**. The pipeline then defines a window extending primarily upstream of this anchor to capture the regulatory promoter region.
    *   *Default Window*: Positive strand (TSS - 1107bp to TSS + 107bp), Negative strand (TSS - 107bp to TSS + 1107bp).
    *   *Customization*: Users can specify the anchor features via `--tss-features` and the window size via `--promoter-upstream` and `--promoter-downstream`.
*   **Recursive Gene Filtering**: If a gene list is provided, the pipeline doesn't just look for the gene ID; it recursively finds all children (e.g., Gene $\rightarrow$ mRNA $\rightarrow$ Exon) to ensure no target regions are missed.

### 2. Sequence Extraction (`extract_sequences`)
Using `pyfaidx`, the pipeline extracts the actual DNA sequences from the reference FASTA. It formats the FASTA headers in a specific way (`GeneID|Feature|Coords|Strand`) so that the origin of every guide can be traced back after scoring.

### 3. FlashFry Integration (`run_flashfry`)
The pipeline wraps the FlashFry Java JAR to perform three main tasks:
1.  **Indexing**: Creates a genomic database for fast off-target searching.
2.  **Discovery**: Finds all possible sgRNA targets within the extracted sequences.
3.  **Scoring**: Applies metrics (default: `doench2014ontarget`) to predict guide efficiency.

### 4. Post-Processing & Ranking (`post_process`)
This is where the raw FlashFry data is refined:
*   **Accessibility**: If an "accessome" BED file is provided, it calculates the overlap between the target and open chromatin regions.
*   **SNP Entropy**: If a SNP BED file is provided (from the `lift_snps.py` liftover pipeline), each overlapping SNP is converted to its **IUPAC ambiguity code** (e.g., A/G → R, C/T → Y). These codes replace bases in the guide sequence, and Shannon entropy is calculated per guide. Guides at polymorphic sites receive higher entropy scores and rank lower. Individual SNP positions are annotated in the `snp_overlaps` column.
*   **SNP Re-scoring** (optional): If `--rescore-snps` is enabled, FlashFry is re-run on SNP-modified sequences for guides that directly overlap SNPs. Results appear in `{metric}_snp` columns alongside reference scores, allowing direct comparison of reference vs. variant guide performance.
*   **Oligo Design**: Generates forward and reverse primer sequences for the guides, including a 'G' prefix for T7 promoter efficiency and specific overhangs (`agat` / `aaac`).
*   **Final Ranking**: Sorts guides based on a hierarchy: `Cutoff` $\rightarrow$ `Doench Score` $\rightarrow$ `Accessibility` $\rightarrow$ `Entropy`.

---


## Local Setup and Execution

### 1. Prerequisites
1.  **Java Runtime Environment (JRE)**: Required to run FlashFry.
2.  **FlashFry JAR**: You must have `FlashFry-assembly-1.15.jar` in the same directory as the script.
3.  **Python 3.x**

### 2. Environment Setup
It is highly recommended to use a virtual environment to ensure dependency isolation:

```bash
# Create and activate a virtual environment
python -m venv venv
source venv/bin/activate  # Windows: venv\Scripts\activate

# Install required Python packages
pip install -r crispy/requirements.txt
```

### 3. Acquiring FlashFry
The pipeline requires the FlashFry assembly JAR. Download it from the [FlashFry releases page](https://github.com/mckennalab/FlashFry/releases):

```bash
# Option A: Download to crispy/ directory (default location)
wget -P crispy/ https://github.com/mckennalab/FlashFry/releases/download/1.15/FlashFry-assembly-1.15.jar

# Option B: Download to a custom location
wget -P /path/to/tools/ https://github.com/mckennalab/FlashFry/releases/download/1.15/FlashFry-assembly-1.15.jar
```

By default, the pipeline looks for `FlashFry-assembly-1.15.jar` in the `crispy/` directory. Override with `--flashfry-jar`:
```bash
python3 crispy/pipeline.py --flashfry-jar /path/to/FlashFry-assembly-1.15.jar ...
```

### Configuration and Execution

The pipeline can be run using command-line arguments or a YAML configuration file for better reproducibility.

#### Option A: Using a Configuration File (Recommended)
Create a `config.yaml` file:
```yaml
gff: "genome_annotation.gff3"
ref: "genome_reference.fasta"
genes: "target_genes.csv"
accessome: "chromatin_accessibility.bed"
snps: "known_snps.bed"
tss_features: "transcript,mRNA,noncoding_transcript,lncRNA,sRNA,tRNA,rRNA"
promoter_upstream: 1107
promoter_downstream: 107
cutoff: 0.6
mem: "16g"
metrics: "doench2014ontarget"
```
Run the pipeline with:
```bash
python crispy/pipeline.py --config config.yaml
```

#### Option B: Using Command Line Arguments
```bash
python crispy/pipeline.py \
    --gff genome_annotation.gff3 \
    --ref genome_reference.fasta \
    --genes target_genes.csv \
    --accessome chromatin_accessibility.bed \
    --snps known_snps.bed \
    --tss-features "transcript,mRNA,lncRNA" \
    --promoter-upstream 1200 \
    --promoter-downstream 100 \
    --cutoff 0.6 \
    --mem 16g \
    --metrics doench2014ontarget,mitchell2018
```

#### Argument Breakdown:
| Argument | Required | Description |
| :--- | :--- | :--- |
| `--config` | No | Path to a YAML configuration file. |
| `--gff` | **Yes** | Path to the GFF3 annotation file. |
| `--ref` | **Yes** | Path to the reference genome FASTA. |
| `--genes` | No | CSV containing a `Gene.ID` column to limit targets to specific genes. |
| `--filter-bed` | No | BED file used to strictly limit the search space to specific genomic regions. |
| `--accessome`| No | BED file of accessible regions to prioritize "open" DNA in the final ranking. |
| `--snps` | No | BED file of SNPs to calculate sequence entropy. |
| `--tss-features` | No | Comma-separated list of GFF3 features to use as TSS anchors. |
| `--promoter-upstream` | No | Distance upstream of TSS for promoter window (Default: 1107). |
| `--promoter-downstream` | No | Distance downstream of TSS for promoter window (Default: 107). |
| `--cutoff` | No | Minimum Doench score to be considered "over cutoff" (Default: 0.65). |
| `--mem` | No | Memory allocated to the Java VM (Default: 8g). |
| `--metrics` | No | Comma-separated FlashFry scoring metrics. |
| `--rescore-snps` | No | Re-run FlashFry scoring on SNP-modified sequences |
| `--browser-bed` | No | Output BED file for genome browser visualization |
| `--verbose` | No | Show FlashFry Java output (suppressed by default) |


### Difference between `--filter-bed` and `--accessome`

While both arguments take a BED file as input, they serve completely different purposes in the pipeline: one is a **hard filter** and the other is a **soft ranking metric**.

#### 1. `--filter-bed` (The "Hard Filter")
*   **When it happens:** This occurs at the very beginning of the pipeline during the `prepare_targets` phase.
*   **What it does:** It performs a geometric intersection. If a target region (exon or promoter) does not overlap with the regions defined in the `filter-bed` file, it is **deleted** from the list of targets.
*   **Result:** These regions are never sent to FlashFry. They are completely excluded from the results.
*   **Use Case:** Use this when you only care about a specific chromosome, a specific genomic window, or a set of known "safe harbor" regions.

#### 2. `--accessome` (The "Soft Ranker")
*   **When it happens:** This occurs at the end of the pipeline during the `post_process` phase.
*   **What it does:** It calculates the total number of base pairs (bp) that a target overlaps with the accessible regions defined in the BED file. This value is stored as an `accessibility` score.
*   **Result:** No targets are deleted. Instead, the `accessibility` score is used as a **sorting priority**. Targets with higher accessibility are pushed to the top of the `final_sgrnas.tsv` list.
*   **Use Case:** Use this when you want to see all possible guides, but you want the ones most likely to be biologically active (due to open chromatin) to be your top recommendations.


### Output Files

| File | Description |
|------|-------------|
| `targets.bed` | Genomic regions selected for targeting |
| `targets.fasta` | Sequences sent to FlashFry |
| `targets.output.scored` | Raw scored output from FlashFry |
| **`final_sgrnas.tsv`** | Final ranked list with scores, annotations, and oligos |
| `final_sgrnas.bed` | Genome browser BED file (with `--browser-bed`) |
| `targets_snp.fasta` | SNP-modified sequences (with `--rescore-snps`) |
| `targets_snp.output.scored` | Scored output with SNP variants (with `--rescore-snps`) |

### Key Columns in `final_sgrnas.tsv`

| Column | Description |
|--------|-------------|
| `gene_id` | Resolved root gene ID from GFF3 (first column) |
| `feature` | Feature type (exon, CDS, promoter) |
| `Chromosome`, `Start`, `End`, `Strand` | Genomic coordinates |
| `target` | 23bp guide + PAM sequence |
| `Doench2014OnTarget` | On-target efficiency score (reference) |
| `entropy` | Shannon entropy based on overlapping SNPs (higher = more polymorphic) |
| `snp_overlaps` | List of overlapping SNPs (e.g., `5:A>G\|12:C>T`) |
| `num_snps` | Number of overlapping SNPs |
| `{metric}_snp` | Re-scored metric with SNP variants (only if `--rescore-snps`) |
| `fwd`, `rev` | Designed oligos for cloning |
| `over_cutoff` | Whether guide passes the Doench cutoff |

### Genome Browser Visualization

With `--browser-bed`, a BED file is generated for IGV/UCSC:

| Color | Meaning |
|-------|---------|
| Green | Passes cutoff, no SNPs |
| Orange | Passes cutoff, but overlaps SNPs |
| Red | Fails cutoff |

Load in IGV: **File → Load from File → final_sgrnas.bed**

### On-Target Scoring

The pipeline computes two on-target efficiency scores: **Doench 2016 (Rule Set 2)** and **Rule Set 3**.

#### Rule Set 2 -- Doench 2016 (`doench2016.py`)

Pure numpy GradientBoostingRegressor (100 trees, lr=0.1, max_depth=3) using 627 features extracted from each 30mer context sequence. No sklearn dependency -- prediction is pure numpy tree traversal.

| Column | Scale | Description |
|--------|-------|-------------|
| `Doench2016OnTarget` | 0 to 1 | On-target efficiency score (higher = better) |

The model (`doench2016_model.json`) was extracted from the azimuth package bundled with CRISPOR. Scores may differ from the CRISPOR website by ~10 points on average, as the production server appears to use a different model than what is distributed in the repository.

#### Rule Set 3 (`rs3score.py`)

LightGBM model (2,279 trees) from the [GPP/rs3](https://github.com/gpp-rnd/rs3) package, using `sglearn` for feature extraction. Requires `sglearn` and `lightgbm` (see requirements.txt). Uses tracrRNA identity features (Hsu2013 default).

| Column | Scale | Description |
|--------|-------|-------------|
| `RuleSet3` | unbounded | On-target efficiency score (higher = better) |

Enable with `--rs3` flag. Spearman correlation with Doench 2016: ~0.72.

#### Setup for RS3

One-time model conversion (requires `rs3test` environment):

```bash
conda create -n rs3test python=3.10 scikit-learn==1.0.2 lightgbm joblib -c conda-forge -y
conda activate rs3test
cd crispy
python convert_rs3_model.py --rs3-dir ~/projects/src/rs3
```

#### Output columns

| Column | Description |
|--------|-------------|
| `Doench2014OnTarget` | On-target score from FlashFry (Rule Set 1, deprecated) |
| `Doench2016OnTarget` | On-target score from Doench 2016 (Rule Set 2) |
| `RuleSet3` | On-target score from Rule Set 3 (requires `--rs3`) |

Ranking uses: cutoff -> Doench 2016 -> Rule Set 3 -> accessibility -> entropy.

The `Doench2016OnTarget` column is used for cutoff filtering (`over_cutoff`).

### Running with SNP Input

To use SNP data, first run the liftover pipeline (`lift_snps.py`) to convert SNP positions to your reference coordinates, then pass the resulting BED file to Crispy:

```bash
# Step 1: Lift SNPs to new assembly coordinates
python3 idmapping/lift_snps.py \
    --old-genome data/old_assembly.fa \
    --new-genome data/new_assembly.fa \
    --vcf data/snps.vcf.gz \
    --output output/liftover \
    --bed \
    --threads 8 \
    --preset asm10

# Step 2: Run Crispy with lifted SNPs and RS3 scoring
python3 crispy/pipeline.py \
    --gff data/genome_annotation.gff3 \
    --ref data/new_assembly.fa \
    --snps output/liftover/snps_lifted.bed \
    --genes data/target_genes.csv \
    --outdir output/sgRNAs \
    --rescore-snps \
    --browser-bed \
    --rs3
```
