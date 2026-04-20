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
*   **Promoters**: Automatically defines promoter regions based on the Transcription Start Site (TSS). 
    *   *Positive strand*: TSS - 1107bp to TSS + 107bp.
    *   *Negative strand*: TSS - 107bp to TSS + 1107bp.
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
*   **SNP Entropy**: If a SNP BED file is provided, it replaces bases with IUPAC ambiguity codes and calculates an entropy score to penalize guides that hit highly polymorphic sites.
*   **Oligo Design**: Generates forward and reverse primer sequences for the guides, including a 'G' prefix for T7 promoter efficiency and specific overhangs (`agat` / `aaac`).
*   **Final Ranking**: Sorts guides based on a hierarchy: `Cutoff` $\rightarrow$ `Doench Score` $\rightarrow$ `Accessibility` $\rightarrow$ `Entropy`.

---

## Local Setup and Execution

### Prerequisites
1.  **Java Runtime Environment (JRE)**: Required to run FlashFry.
2.  **FlashFry JAR**: You must have `FlashFry-assembly-1.15.jar` in the same directory as the script.
3.  **Python 3.x**

### Installation
```bash
# 1. Create and activate a virtual environment
python -m venv venv
source venv/bin/activate  # Windows: venv\Scripts\activate

# 2. Install required Python packages
pip install -r crispy/requirements.txt
```

### Running Example
To run the pipeline, execute `pipeline.py` from the command line. Below is a comprehensive example including optional filters.

#### Example Command:
```bash
python crispy/pipeline.py \
    --gff genome_annotation.gff3 \
    --ref genome_reference.fasta \
    --genes target_genes.csv \
    --accessome chromatin_accessibility.bed \
    --snps known_snps.bed \
    --cutoff 0.6 \
    --mem 16g \
    --metrics doench2014ontarget,mitchell2018
```

#### Argument Breakdown:
| Argument | Required | Description |
| :--- | :--- | :--- |
| `--gff` | **Yes** | Path to the GFF3 annotation file. |
| `--ref` | **Yes** | Path to the reference genome FASTA. |
| `--genes` | No | CSV containing a `Gene.ID` column to limit targets to specific genes. |
| `--filter-bed` | No | BED file used to strictly limit the search space to specific genomic regions. |
| `--accessome`| No | BED file of accessible regions to prioritize "open" DNA in the final ranking. |
| `--snps` | No | BED file of SNPs to calculate sequence entropy. |
| `--cutoff` | No | Minimum Doench score to be considered "over cutoff" (Default: 0.65). |
| `--mem` | No | Memory allocated to the Java VM (Default: 8g). |
| `--metrics` | No | Comma-separated FlashFry scoring metrics. |


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
After a successful run, the following files will be generated in your working directory:
*   `targets.bed`: The genomic regions selected for targeting.
*   `targets.fasta`: The sequences sent to FlashFry.
*   `targets.output.scored`: The raw scored output from FlashFry.
*   **`final_sgrnas.tsv`**: The final ranked list containing the best guides, their scores, and the designed `fwd` and `rev` oligos.
