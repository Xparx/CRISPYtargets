# Design: sgRNA discovery and scoring

## Overview
The **design** capability identifies, scores, and designs CRISPR guide RNAs (sgRNAs). It automates the process of moving from a genomic reference (GFF3 and FASTA) to a ranked list of target sequences with ready-to-order oligonucleotide designs.

Run it as `crispytargets design` (see the [top-level README](../README.md) for installation).

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
3.  **Scoring**: Applies the `--metrics` list (default: `hsu2013`, the MIT/Hsu off-target specificity score; `otCount` is always emitted alongside). On-target efficiency metrics from FlashFry are opt-in, since our own RS2/RS3 provide those.

### 4. Post-Processing & Ranking (`post_process`)
This is where the raw FlashFry data is refined:
*   **Accessibility**: If an "accessome" BED file is provided, it calculates the overlap between the guide's own 23bp window and open chromatin regions.
*   **SNP Entropy**: If a SNP BED file is provided (from the `crispytargets lift` tool), a SNP is assigned to a guide only when it falls inside that guide's 23bp window (computed in the guide's own orientation, so a guide inside a large exon is not credited with SNPs elsewhere in the exon). Each overlapping SNP is converted to its **IUPAC ambiguity code** (e.g., A/G → R, C/T → Y). These codes replace bases in the guide sequence, and Shannon entropy is calculated per guide. Guides at polymorphic sites receive higher entropy scores and rank lower. Individual SNP positions are annotated in the `snp_overlaps` column, together with the minimum distance from a SNP to the PAM (`snp_pam_dist`, positive = in guide, <=0 = inside PAM) and to the cut site (`snp_cut_dist`). These annotations need only `--snps`, no scoring.
*   **SNP Re-scoring** (optional): If `--rescore-snps` is enabled, guides that overlap known SNPs are re-scored with Doench 2016 (RS2) and, with `--rs3`, Rule Set 3, using the concrete variant (ALT) allele at each SNP position and the guide's real flanks and PAM. For each guide this yields three variant scores, `avg1snp` (mean over single-SNP variants, the population-average effect), `worst1snp` (most damaging single variant, the carrier worst case), and `allsnps` (every overlapping SNP substituted at once), each reported as a raw score and as a delta versus the reference score. Unlike the reference scores, these come from our own scorers, so they are directly comparable to `Doench2016OnTarget` and `RuleSet3`.
*   **Oligo Design**: Generates forward and reverse primer sequences for the guides, including a 'G' prefix for T7 promoter efficiency and specific overhangs (`agat` / `aaac`).
*   **Final Ranking**: Sorts guides based on a hierarchy: `Cutoff` $\rightarrow$ `Doench Score` $\rightarrow$ `Accessibility` $\rightarrow$ `Entropy`.

---


## Local Setup and Execution

### 1. Prerequisites
1.  **Java Runtime Environment (JRE)**: Required to run FlashFry.
2.  **FlashFry JAR**: `FlashFry-assembly-1.15.jar`, resolved via `--flashfry-jar` or `$CRISPY_FLASHFRY_JAR`.
3.  **Python 3.x**

### 2. Environment Setup
It is highly recommended to use a virtual environment to ensure dependency isolation:

```bash
# Create and activate a virtual environment
python -m venv venv
source venv/bin/activate  # Windows: venv\Scripts\activate

# Install the package and its Python dependencies
pip install -e .
```

### 3. Acquiring FlashFry
The pipeline requires the FlashFry assembly JAR. Download it from the [FlashFry releases page](https://github.com/mckennalab/FlashFry/releases):

```bash
# Option A: keep it with the other external tools in this repo
wget -P tools/ https://github.com/mckennalab/FlashFry/releases/download/1.15/FlashFry-assembly-1.15.jar

# Option B: download to any custom location
wget -P /path/to/tools/ https://github.com/mckennalab/FlashFry/releases/download/1.15/FlashFry-assembly-1.15.jar
```

The JAR is an untracked external binary and is not shipped inside the
package, so the pipeline resolves it from `--flashfry-jar` or the
`CRISPY_FLASHFRY_JAR` environment variable:

```bash
export CRISPY_FLASHFRY_JAR=tools/FlashFry-assembly-1.15.jar
crispytargets design --gff genome.gff3 --ref genome.fa ...
# or per run:
crispytargets design --flashfry-jar tools/FlashFry-assembly-1.15.jar ...
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
metrics: "hsu2013"
```
Run the pipeline with:
```bash
crispytargets design --config config.yaml
```

#### Option B: Using Command Line Arguments
```bash
crispytargets design \
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
    --metrics hsu2013,doench2014ontarget,moreno2015
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
| `--snps` | No | BED file of SNPs; gives `num_snps`, `snp_overlaps`, `entropy`, `snp_pam_dist` and `snp_cut_dist`. |
| `--tss-features` | No | Comma-separated list of GFF3 features to use as TSS anchors. |
| `--promoter-upstream` | No | Distance upstream of TSS for promoter window (Default: 1107). |
| `--promoter-downstream` | No | Distance downstream of TSS for promoter window (Default: 107). |
| `--cutoff` | No | Minimum Doench score to be considered "over cutoff" (Default: 0.65). |
| `--mem` | No | Memory allocated to the Java VM (Default: 8g). |
| `--metrics` | No | Comma-separated FlashFry scoring metrics (Default: `hsu2013`). On-target metrics (`doench2014ontarget`, `doench2016cfd`, `moreno2015`, ...) are opt-in and are written as trailing columns of `final_sgrnas.tsv`. |
| `--rescore-snps` | No | Add RS2/RS3 variant scores for SNP-overlapping guides (avg/worst single-SNP and all-SNP scores plus deltas). This is the expensive step; SNP annotation needs only `--snps`. |
| `--strict` | No | Write only guides over `--cutoff` to `final_sgrnas.tsv` (default keeps all guides, ranked, with `over_cutoff` flag) |
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
| `design.log` | Per-run log: executed commands with timestamps, plus a final line with stage timings (prepare targets, FlashFry, post-process, total) |

### Key Columns in `final_sgrnas.tsv`

Columns are grouped so related scores sit next to each other, in this order:

| # | Column(s) | Description |
|---|-----------|-------------|
| 1-3 | `gene_id`, `feature`, `over_cutoff` | Identity: resolved root gene ID, feature type (exon, CDS, promoter), passes the Doench cutoff |
| 4-10 | `Chromosome`, `Start`, `End`, `Strand`, `start`, `stop`, `contig` | Coordinates: gene location plus the guide's position within it |
| 11-15 | `Doench2016OnTarget`, `RuleSet3`, `Hsu2013`, `otCount`, `accessibility` | Reference scores: our on-target scores (RS2, RS3), the default off-target metrics from FlashFry, and open-chromatin overlap |
| 16-21 | `target`, `context`, `orientation`, `overflow`, `fwd`, `rev` | Guide sequence: 23bp guide + PAM, flanking context, strand, and the designed oligos |
| 22-26 | `num_snps`, `snp_overlaps`, `entropy`, `snp_pam_dist`, `snp_cut_dist` | SNP annotation, grouped with the variant scores below because both are absent when no SNP data is given |
| 27-38 | `Doench2016OnTarget_avg1snp` / `_worst1snp` / `_allsnps` + `Doench2016_delta_avg` / `_worst` / `_all`, then the same six for `RuleSet3` | Variant scores (only with `--rescore-snps`, RS3 columns need `--rs3`): raw variant scores plus deltas versus the reference score (variant minus reference; positive = improvement) |
| trailing | any further `--metrics` column, e.g. `Doench2014OnTarget` | FlashFry on-target scores, written last as pass-through whenever requested via `--metrics` |

Notes on individual columns:

| Column | Description |
|--------|-------------|
| `entropy` | Shannon entropy based on overlapping SNPs (higher = more polymorphic) |
| `Strand`, `orientation` | `Strand` is the target region's strand from the GFF3; `orientation` is the guide's own strand (FWD/RVS). Window overlaps, flanks, variant alleles and the browser BED all follow `orientation` |
| `snp_overlaps` | List of overlapping SNPs (e.g., `5:A>G\|12:C>T`) |
| `snp_pam_dist` | Min distance from a SNP to the PAM boundary (positive = in guide, <=0 = inside PAM) |
| `snp_cut_dist` | Min distance from a SNP to the cut site (between guide positions 17/18) |
| `Hsu2013` | MIT/Hsu off-target specificity score from FlashFry (0-100, higher = more specific; default `--metrics`) |
| `otCount` | Number of off-target hits found by FlashFry (always written) |

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

The 30mer uses the Azimuth layout: **4bp 5' flank + 20bp guide + 3bp PAM + 3bp 3' flank**, so the guide sits at offset 4..23 and the PAM at 24..26. Feature slices (GC from `s[4:24]`, NGGX from `s[24]+s[27]`, Tm segments) follow Azimuth's `featurization.py` exactly.

| Column | Scale | Description |
|--------|-------|-------------|
| `Doench2016OnTarget` | 0 to 1 | On-target efficiency score (higher = better) |

The model (`doench2016_model.json`) is a **retrain**: a GradientBoostingRegressor fit from scratch on Azimuth's `FC_plus_RES_withPredictions.csv` training data (see `tools/train_doench2016.py`), converted to pure-numpy JSON. It is *not* the production Azimuth V3 model, so scores differ from the CRISPOR website by ~9-10 points on average. Exact CRISPOR `fusi` reproduction requires the original `V3_model_nopos.pickle`, which needs scikit-learn 1.0.x to unpickle (see `reproduce_crispor.py` for a reference attempt).

#### Rule Set 3 (`rs3score.py`)

LightGBM model (2,279 trees) from the [GPP/rs3](https://github.com/gpp-rnd/rs3) package, using `sglearn` for feature extraction. Requires `sglearn` and `lightgbm` (install with `pip install -e ".[rs3]"`). Uses tracrRNA identity features (Hsu2013 default).

| Column | Scale | Description |
|--------|-------|-------------|
| `RuleSet3` | unbounded | On-target efficiency score (higher = better) |

Enable with `--rs3` flag. Spearman correlation with Doench 2016: ~0.72.

The 30mers fed to RS3 use the same Azimuth layout (4+20+3+3, guide at offset 4), matching what CRISPOR feeds its `rs3` column, so these scores reproduce the CRISPOR "Doench-RuleSet3" column. The model itself is the genuine published Rule Set 3, converted faithfully from `gpp-rnd/rs3`'s `RuleSet3.pkl`.

#### Setup for RS3

One-time model conversion (requires `rs3test` environment):

```bash
conda create -n rs3test python=3.10 scikit-learn==1.0.2 lightgbm joblib -c conda-forge -y
conda activate rs3test
python tools/convert_rs3_model.py --rs3-dir ~/projects/src/rs3
```

#### Output columns

| Column | Description |
|--------|-------------|
| `Doench2016OnTarget` | On-target score from Doench 2016 (Rule Set 2) |
| `RuleSet3` | On-target score from Rule Set 3 (requires `--rs3`) |
| `Doench2014OnTarget` | On-target score from FlashFry (Rule Set 1), only when requested via `--metrics`; written as a trailing column |

Ranking uses: cutoff -> Doench 2016 -> Rule Set 3 -> accessibility -> entropy.

The `Doench2016OnTarget` column is used for cutoff filtering (`over_cutoff`).

### Running with SNP Input

To use SNP data, first run `crispytargets lift` to convert SNP positions to your reference coordinates, then pass the resulting BED file to the design pipeline:

```bash
# Step 1: Lift SNPs to new assembly coordinates
crispytargets lift \
    --old-genome data/old_assembly.fa \
    --new-genome data/new_assembly.fa \
    --vcf data/snps.vcf.gz \
    --output output/liftover \
    --bed \
    --threads 8 \
    --preset asm10

# Step 2: Run Crispy with lifted SNPs and RS3 scoring
crispytargets design \
    --gff data/genome_annotation.gff3 \
    --ref data/new_assembly.fa \
    --snps output/liftover/snps_lifted.bed \
    --genes data/target_genes.csv \
    --outdir output/sgRNAs \
    --rescore-snps \
    --browser-bed \
    --rs3
```
