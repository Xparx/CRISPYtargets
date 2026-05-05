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
