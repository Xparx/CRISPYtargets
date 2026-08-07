# IDMapping: Incremental Orthology Mapping Engine

## Overview
The **IDMapping** module is a professional bioinformatics utility designed to map genomic feature IDs from a source organism to their corresponding homologs in any number of target species.

Unlike simple BLAST wrappers, this engine implements an **Incremental Backbone** strategy. It allows researchers to build a master orthology table one species at a time, resolving complex GFF3 hierarchies to ensure that mappings are attributed to the root **Gene** rather than individual isoforms or exons.

---

## How the Code Works

### 1. GFF3 Hierarchy Resolution (`GFF3Resolver`)
The `GFF3Resolver` ensures biological accuracy by navigating the genomic hierarchy:
*   **Normalization**: Automatically removes prefixes (e.g., `KY21:`) to match BLAST query IDs.
*   **BFS Resolution**: Uses a Breadth-First Search (BFS) to climb from any feature (exon, CDS, mRNA) to its closest ancestor labeled as a `'gene'`. This correctly handles multi-parent features and alternative splicing.

### 2. Optimal Hit Selection Logic
Instead of simply taking the top BLAST hit, the engine identifies the **Optimal Hit** for each gene:
*   **Robustness Check**: If the top $N$ hits all resolve to the same target gene, that mapping is marked as high-confidence.
*   **Significance Gap**: The engine evaluates the bitscore gap between the top hit and subsequent hits to distinguish between a clear ortholog and a family of paralogs.
*   **Primary vs. Secondary**: The best biological candidate is assigned as the `Primary` hit. All other hits passing the E-value threshold are concatenated into a `Secondary` column, separated by a pipe (`|`).

### 3. Incremental Backbone Workflow
The module maintains a master table (the "Backbone") of all genes in the reference GFF3:
*   **Initialization**: If no output file exists, the backbone is initialized with all genes from the reference GFF3.
*   **Additive Mapping**: Each new species is joined to the backbone via a left-merge.
*   **Collision Handling**: If a species is re-mapped, the existing columns for that species are overwritten with the new results.

---

## 🛠️ Setup and Data Preparation

The `idmapping` module does not run BLAST internally; it processes the results of a BLAST search. To generate the necessary input files, follow this end-to-end workflow.

### 1. Install NCBI BLAST+
You must have the **NCBI BLAST+** suite installed on your system to generate the homology results.
*   **Ubuntu/Debian:** `sudo apt-get install ncbi-blast+`
*   **MacOS (Homebrew):** `brew install blast`
*   **Conda:** `conda install -c bioconda blast`

### 2. Prepare Target Proteomes (Databases)
To find homologs, you need a high-quality reference proteome for the target species (e.g., Human and Mouse).
1.  Go to [UniProt.org](https://www.uniprot.org/).
2.  Search for the target species (e.g., **Human** [Taxon: 9606] or **Mouse** [Taxon: 10090]).
3.  **Crucial Filters:** 
    *   Filter by **"Reviewed (Swiss-Prot)"** to ensure curated, high-quality sequences.
    *   Download in **FASTA (canonical)** format.
4.  Save these as `human_proteome.fasta` and `mouse_proteome.fasta`.

### 3. Create BLAST Databases
BLAST requires FASTA files to be indexed into a binary database format before they can be searched.
```bash
# Create Human database
makeblastdb -in human_proteome.fasta -dbtype prot -out human_db

# Create Mouse database
makeblastdb -in mouse_proteome.fasta -dbtype prot -out mouse_db
```

### 4. Generate Homology Results (`blastp`)
Run a protein-protein BLAST search. 

**⚠️ Important:** The `query_proteins.fasta` file must contain the protein sequences of your source organism, and the headers must match the IDs used in your **GFF3 file**.

```bash
# Search against Human
blastp -query query_proteins.fasta \
       -db human_db \
       -out human_results.tsv \
       -outfmt 6 \
       -evalue 1e-5 \
       -max_target_seqs 5 \
       -num_threads 4

# Search against Mouse
blastp -query query_proteins.fasta \
       -db mouse_db \
       -out mouse_results.tsv \
       -outfmt 6 \
       -evalue 1e-5 \
       -max_target_seqs 5 \
       -num_threads 4
```
*Note: `-outfmt 6` is mandatory as the Python mapper expects a tab-separated table.*

---

## 🚀 Execution

### Prerequisites
*   **Python 3.x**
*   **Pandas & PyYAML**: `pip install pandas PyYAML`

### Running the Incremental Mapper
The mapping is performed one species at a time. Run the script for each species you wish to add to your master table.

#### Step 1: Map the first species (e.g., Human)
```bash
python idmapping/run_mapping.py \
    --gff source_organism.gff3 \
    --blast-file human_results.tsv \
    --species Human \
    --output orthologs.tsv
```
*This initializes the table with all genes from your GFF and adds `Human_primary` and `Human_secondary` columns.*

#### Step 2: Map subsequent species (e.g., Mouse)
```bash
python idmapping/run_mapping.py \
    --gff source_organism.gff3 \
    --blast-file mouse_results.tsv \
    --species Mouse \
    --output orthologs.tsv
```
*This loads the existing `orthologs.tsv` and appends `Mouse_primary` and `Mouse_secondary` columns.*

#### Step 3: Using a Configuration File
Create a `mapping_config.yaml`:
```yaml
gff: "source_organism.gff3"
output: "orthologs.tsv"
evalue: 1e-50
include_top: true
```
Run:
```bash
python idmapping/run_mapping.py --config mapping_config.yaml --blast-file zfish.tsv --species Zebrafish
```

#### Argument Breakdown:
| Argument | Required | Description |
| :--- | :--- | :--- |
| `--config` | No | Path to a YAML configuration file. |
| `--gff` | **Yes** | Path to the source organism's GFF3 file. |
| `--blast-file` | **Yes** | Path to the BLAST results TSV for the species. |
| `--species` | **Yes** | Label for the species (e.g., 'Human', 'Mouse'). |
| `--output` | **Yes** | Path to the output TSV (will be updated incrementally). |
| `--evalue` | No | E-value threshold for significance (Default: `1e-50`). |
| `--include-top`| No | If present, keeps the top hit even if it's above the E-value threshold. |

### Output Format
The resulting `orthologs.tsv` is a wide table where each row is a unique gene from your reference genome.

| gene_id | {Species}_primary | {Species}_secondary | ... |
| :--- | :--- | :--- | :--- |
| GENE_001 | H_GENE_A | H_GENE_B$\vert$H_GENE_C | ... |
| GENE_002 | $\text{NaN}$ | $\text{NaN}$ | ... |

*   **`gene_id`**: The resolved root gene ID from your reference GFF3.
*   **`{Species}_primary`**: The optimal hit based on bitscore and robustness.
*   **`{Species}_secondary`**: All other hits passing the E-value threshold, separated by a pipe (`|`).


# SNP Liftover Pipeline

Map SNP positions from an old genome assembly to a new genome assembly using whole-genome alignment.

## Background

This pipeline converts SNP positions from one genome assembly to another by aligning the two genomes with minimap2 and lifting coordinates through the alignment blocks.

- **Input**: VCF with old-genome coordinates
- **Output**: VCF with new-genome coordinates
- **Method**: PAF-based whole-genome alignment (no BAM/CIGAR intermediate)

## Requirements

Install all dependencies via conda:

```bash
conda create -n liftover python=3.10 pysam minimap2 bcftools tabix ucsc-twobittofa -c bioconda -c conda-forge
conda activate liftover
```

Or install into an existing environment:

```bash
conda install -c bioconda -c conda-forge pysam minimap2 bcftools tabix ucsc-twobittofa
```

### Tool summary

| Tool | Purpose |
|------|---------|
| `pysam` | Read/write VCF files in Python |
| `minimap2` | Whole-genome alignment between assemblies (PAF output) |
| `bcftools` | Filter and manipulate VCF files |
| `tabix` | Index VCF files |
| `ucsc-twobittofa` | Convert `.2bit` genome to FASTA (optional) |

## Quick Start

### Basic usage

```bash
python3 lift_snps.py \
    --old-genome old_assembly.fa \
    --new-genome new_assembly.fa \
    --vcf input_snps.vcf.gz \
    --output output_snps.vcf.gz \
    --threads 8 \
    -v
```

### With all options

```bash
python3 lift_snps.py \
    --old-genome old_assembly.fa \
    --new-genome new_assembly.fa \
    --vcf input_snps.vcf.gz \
    --output output_snps.vcf.gz \
    --bed \
    --threads 8 \
    --preset asm10 \
    --keep-intermediates \
    -v
```

### From a .2bit file (auto-converts)

```bash
python3 lift_snps.py \
    --old-genome old_assembly.2bit \
    --new-genome new_assembly.fa \
    --vcf input_snps.vcf.gz \
    --output output_snps.vcf.gz \
    --threads 8 \
    -v
```

## CLI Arguments

| Argument | Required | Default | Description |
|----------|----------|---------|-------------|
| `--old-genome` | Yes | — | Old assembly (FASTA or .2bit) |
| `--new-genome` | Yes | — | New assembly FASTA |
| `--vcf` | Yes | — | Input VCF(.gz) with old-genome coordinates |
| `--output` | Yes | — | Output VCF (use `.vcf.gz` for compressed output) |
| `--bed` | No | off | Also output a BED file with lifted positions |
| `--threads` | No | 8 | Number of threads for minimap2 |
| `--preset` | No | `asm5` | Minimap2 preset: `asm5`, `asm10`, `asm20`, `map-ont`, `map-hifi` |
| `--keep-intermediates` | No | off | Keep the intermediate PAF alignment file |
| `-v, --verbose` | No | off | Enable debug logging |

## What the Pipeline Does

1. **Whole-genome alignment** — Aligns old genome (query) to new genome (reference) using minimap2 with PAF output
2. **Block extraction** — Parses PAF records to extract aligned blocks with coordinate correspondence between genomes
3. **Position mapping** — For each SNP position, binary-searches the aligned blocks to find the corresponding position in the new genome
4. **VCF output** — Writes a new VCF with lifted coordinates, preserving REF, ALT, QUAL, FILTER, and INFO fields
5. **Indexing** — Creates a `.tbi` index for the output VCF

### Why PAF instead of BAM?

The pipeline uses minimap2's PAF output directly rather than converting to BAM. This avoids:
- Complex CIGAR string parsing
- Supplementary alignment handling (PAF records directly provide query start/end coordinates)
- SAM/BAM sorting and indexing overhead

Each PAF record provides clean coordinate ranges:
- `qstart`–`qend`: position in the old genome
- `tstart`–`tend`: position in the new genome

## Preset Selection

| Preset | Use case |
|--------|----------|
| `asm5` | Same species, <5% divergence (default) |
| `asm10` | Same genus, 5–10% divergence |
| `asm20` | Different species, 10–20% divergence |
| `map-ont` | Nanopore reads |
| `map-hifi` | PacBio HiFi reads |

## Output

- **VCF** (`output.vcf.gz`) — Lifted VCF with new-genome coordinates, tabix-indexed
- **BED** (optional, `--bed`) — Tab-separated file with `chrom`, `start`, `end`, `name`, `score`, `strand`
- **PAF** (intermediate) — Whole-genome alignment (deleted unless `--keep-intermediates`)

The log reports how many variants were successfully lifted and how many failed, with per-chromosome failure counts.

## Troubleshooting

### Many variants fail to map

- Run with `-v` to see per-chromosome failure counts
- Check that the old genome FASTA contains the expected chromosomes
- Try a more relaxed preset (`asm10` or `asm20`) for divergent assemblies
- Inspect the PAF file to check alignment coverage: `wc -l alignment.paf`

### Alignment is slow

- Increase `--threads` for parallelism
- Only align the main chromosomes (filter the old genome FASTA before running)
- Reduce memory by using `--preset asm5` for closely related assemblies

### Chromosome names don't match expected mapping

- The mapping is determined by the alignment, not hardcoded
- Check PAF output to see which old-genome chromosomes map to which new-genome chromosomes: `cut -f1,6 alignment.paf | sort -u`

## SNP Source and Preprocessing

The input SNPs were obtained from [CRISPOR's ci3 genome page](https://crispor.gi.ucsc.edu/genomes/ci3), which provides variant annotations for the *Ciona robusta* (ci3) assembly.

### Obtaining the SNP VCF

The SNP data comes from the KH2013 genome build (ci3 assembly) as provided by CRISPOR:

1. Navigate to [https://crispor.gi.ucsc.edu/genomes/ci3](https://crispor.gi.ucsc.edu/genomes/ci3)
2. Download the variant file `Ghost_SNPs.vcf.gz` — this contains SNP positions on the KH2013 (ci3) assembly
3. Filter to main chromosomes (KhC1–KhC14) if the VCF contains scaffolds:

```bash
bcftools view \
    -r KhC1,KhC2,KhC3,KhC4,KhC5,KhC6,KhC7,KhC8,KhC9,KhC10,KhC11,KhC12,KhC13,KhC14 \
    Ghost_SNPs.vcf.gz | bgzip > Ghost_SNPs_main.vcf.gz
tabix -p vcf Ghost_SNPs_main.vcf.gz
```

### Lifting SNPs to the new assembly

The SNPs are in KH2013 (ci3/old) coordinates and must be lifted to the new KY2021 assembly (`HT.Ref.fasta`) before they can be used by the Crispy pipeline:

```bash
python3 idmapping/lift_snps.py \
    --old-genome data/KH2013_main.fa \
    --new-genome data/HT.Ref.fasta \
    --vcf data/Ghost_SNPs_main.vcf.gz \
    --output data/Ghost_SNPs_lifted.vcf.gz \
    --bed \
    --threads 8 \
    --preset asm10 \
    --outdir data/
```

### Output files

| File | Description |
|------|-------------|
| `Ghost_SNPs_lifted_kh2013.vcf.gz` | Lifted VCF with new-genome coordinates |
| `Ghost_SNPs_lifted_kh2013.bed` | BED file for Crispy pipeline input |

The BED file format (6 columns):

| Column | Example | Description |
|--------|---------|-------------|
| chrom | Chr1 | New-genome chromosome |
| start | 108229 | 0-based start position |
| end | 108230 | 1-based end position (1bp SNP) |
| name | KhC1:211:C/A | Original SNP ID from ci3 |
| score | 1000 | QUAL from VCF (or 1000 if missing) |
| strand | . | Always . for SNPs |

### Feeding SNPs into Crispy

Use the lifted BED file with the `--snps` argument:

```bash
python3 crispy/pipeline.py \
    --gff data/HT.KY21Gene.2.gff3 \
    --ref data/HT.Ref.fasta \
    --snps data/Ghost_SNPs_lifted_kh2013.bed \
    --genes data/genelist.tsv \
    --rescore-snps \
    --browser-bed
```

## Notes

- Chromosome name correspondence is derived automatically from the alignment
- Small scaffolds may not map reliably — filter to main chromosomes for best results
- Output VCF is sorted by chromosome and position
- The pipeline handles both gzipped and ungzipped VCF input/output
