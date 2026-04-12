# IDMapping Module Documentation

## Overview
The **IDMapping** module is a bioinformatics utility designed to map genomic feature IDs from a source organism to their corresponding homologs in other species (specifically Human and Mouse). 

It solves a common problem in comparative genomics: BLAST results often return hits for specific transcripts or exons (children), but researchers typically need to know which top-level **Gene** (parent) those features belong to. This module parses GFF3 annotation files to resolve these hierarchical relationships and filters BLAST results based on statistical significance.

---

## How the Code Works

### 1. GFF3 Hierarchy Resolution (`GFF3Resolver`)
The `GFF3Resolver` class acts as a map of the genome's organization.
*   **Normalization**: It automatically removes prefixes (like `KY21:`) from IDs to ensure they match the query IDs (`qseqid`) found in BLAST output files.
*   **The "Climb"**: The `resolve_to_gene` method implements a recursive lookup. If a feature ID is provided (e.g., an exon), the resolver looks up its parent, then that parent's parent, and so on, until it reaches a feature labeled as a `'gene'`. This ensures that all sub-features are correctly attributed to their root gene.

### 2. BLAST Result Processing (`create_mapping`)
The module processes BLAST output files for multiple species using the following logic:
*   **Significance Filtering**: By default, only hits with an E-value below a specific threshold (e.g., $1e-50$) are kept.
*   **Top-Hit Fallback**: If the `--include-top` flag is enabled, the pipeline will keep the best possible hit for a query even if it fails to meet the E-value threshold.
*   **Data Enrichment**: For every hit, the module attaches the resolved **Gene ID**, the **Feature Type**, and the **Species Gene Symbol**.

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
*   **Pandas**: `pip install pandas`

### Running the Resolver
Once you have your `.tsv` results and your GFF3 annotation, run the mapping script to resolve the feature IDs to top-level genes.

#### Example Command:
```bash
python idmapping/run_mapping.py \
    --gff source_organism.gff3 \
    --human human_results.tsv \
    --mouse mouse_results.tsv \
    --output final_ortholog_map.tsv \
    --evalue 1e-50 \
    --include-top
```

#### Argument Breakdown:
| Argument | Required | Description |
| :--- | :--- | :--- |
| `--gff` | **Yes** | Path to the source organism's GFF3 file. |
| `--human` | **Yes** | Path to the BLAST results against the Human genome. |
| `--mouse` | **Yes** | Path to the BLAST results against the Mouse genome. |
| `--output` | **Yes** | Path where the final mapping TSV will be saved. |
| `--evalue` | No | E-value threshold for significance (Default: `1e-50`). |
| `--include-top`| No | If present, keeps the top hit even if it's above the E-value threshold. |

### Output Format
The resulting `final_ortholog_map.tsv` is a tab-separated file with the following columns:
*   `qseqid`: The original query ID from the BLAST search.
*   `feature_type`: The type of feature (e.g., exon, CDS).
*   `gene_id`: The resolved top-level Gene ID.
*   `species_gene_symbol`: The gene symbol of the homolog in the target species.
*   `species`: Either "Human" or "Mouse".
*   `evalue`: The BLAST E-value.
*   `bitscore`: The BLAST bitscore.
