# CRISPYtargets

CRISPR guide RNA design pipeline for non-model organisms, with genomic ID mapping tools.

## What's in this repo

| Directory | Purpose |
|-----------|---------|
| `crispy/` | sgRNA discovery, scoring, and oligo design |
| `idmapping/` | SNP liftover between genome assemblies + ortholog mapping |

## Requirements

```bash
conda create -n crispy python=3.10 pandas pyranges biopython pyfaidx pyyaml pysam minimap2 bcftools tabix numpy -c conda-forge -c bioconda
conda activate crispy
pip install -r crispy/requirements.txt
```

FlashFry (already in `crispy/`):
```bash
ls crispy/FlashFry-assembly-1.15.jar   # should exist
```

## Quick Start

### Workflow 1: Design sgRNAs (Crispy)

Find and rank CRISPR guides for your target genes.

```bash
python3 crispy/pipeline.py \
    --gff /path/to/genome.gff3 \
    --ref /path/to/genome.fa \
    --outdir output/sgRNAs
```

**With a gene list:**
```bash
python3 crispy/pipeline.py \
    --gff /path/to/genome.gff3 \
    --ref /path/to/genome.fa \
    --genes /path/to/genes.csv \
    --outdir output/sgRNAs
```

**With SNPs and browser output:**
```bash
python3 crispy/pipeline.py \
    --gff /path/to/genome.gff3 \
    --ref /path/to/genome.fa \
    --snps /path/to/lifted_snps.bed \
    --outdir output/sgRNAs \
    --rescore-snps --browser-bed
```

**With chromatin accessibility:**
```bash
python3 crispy/pipeline.py \
    --gff /path/to/genome.gff3 \
    --ref /path/to/genome.fa \
    --accessome /path/to/open_chromatin.bed \
    --outdir output/sgRNAs
```

### Workflow 2: Lift SNPs between assemblies

Map SNP positions from one genome assembly to another.

```bash
python3 idmapping/lift_snps.py \
    --old-genome /path/to/old_assembly.fa \
    --new-genome /path/to/new_assembly.fa \
    --vcf /path/to/snps.vcf.gz \
    --output output/liftover \
    --bed --preset asm10
```

### Workflow 3: Map orthologs between species

Find human/mouse orthologs via BLAST.

```bash
# Create BLAST databases (one-time)
makeblastdb -in /path/to/human_proteome.fasta -dbtype prot -out output/human_db
makeblastdb -in /path/to/mouse_proteome.fasta -dbtype prot -out output/mouse_db

# Run BLAST (outfmt 6 is required)
blastp -query /path/to/your_proteins.fasta -db output/human_db \
    -out output/human_results.tsv -outfmt 6 -evalue 1e-5
blastp -query /path/to/your_proteins.fasta -db output/mouse_db \
    -out output/mouse_results.tsv -outfmt 6 -evalue 1e-5

# Map orthologs (add species one at a time)
python3 idmapping/run_mapping.py \
    --gff /path/to/genome.gff3 \
    --blast-file output/human_results.tsv \
    --species Human \
    --output output/orthologs

python3 idmapping/run_mapping.py \
    --gff /path/to/genome.gff3 \
    --blast-file output/mouse_results.tsv \
    --species Mouse \
    --output output/orthologs
```

## Full pipeline example

```bash
# 1. Lift SNPs
python3 idmapping/lift_snps.py \
    --old-genome /path/to/old_assembly.fa \
    --new-genome /path/to/new_assembly.fa \
    --vcf /path/to/snps.vcf.gz \
    --output output/liftover \
    --bed --preset asm10

# 2. Design sgRNAs with SNP awareness
python3 crispy/pipeline.py \
    --gff /path/to/genome.gff3 \
    --ref /path/to/new_assembly.fa \
    --snps output/liftover/*_lifted.bed \
    --outdir output/sgRNAs \
    --rescore-snps --browser-bed
```

## All CLI arguments

<details>
<summary>crispy/pipeline.py</summary>

| Argument | Default | Description |
|----------|---------|-------------|
| `--gff` | — | Reference GFF3 |
| `--ref` | — | Reference FASTA |
| `--outdir` | `./output` | Output directory |
| `--genes` | — | Gene list CSV (`Gene.ID` column) |
| `--snps` | — | SNP BED file |
| `--accessome` | — | Accessible regions BED |
| `--filter-bed` | — | Hard region filter BED |
| `--cutoff` | 0.65 | Doench 2016 score cutoff |
| `--metrics` | doench2014ontarget | FlashFry scoring metrics |
| `--rescore-snps` | off | Re-score with SNP variants |
| `--browser-bed` | off | Output genome browser BED |
| `--verbose` | off | Show FlashFry output |

</details>

<details>
<summary>idmapping/lift_snps.py</summary>

| Argument | Default | Description |
|----------|---------|-------------|
| `--old-genome` | — | Old assembly (FASTA or .2bit) |
| `--new-genome` | — | New assembly FASTA |
| `--vcf` | — | Input VCF |
| `--output` | — | Output directory |
| `--bed` | off | Also output BED |
| `--preset` | asm5 | minimap2 preset |
| `--threads` | 8 | Threads |
| `--keep-intermediates` | off | Keep PAF file |
| `-v, --verbose` | off | Debug logging |

</details>

<details>
<summary>idmapping/run_mapping.py</summary>

| Argument | Default | Description |
|----------|---------|-------------|
| `--gff` | — | Source organism GFF3 |
| `--blast-file` | — | BLAST results TSV (outfmt 6) |
| `--species` | — | Species label |
| `--output` | — | Output directory |
| `--evalue` | 1e-50 | E-value threshold |
| `--include-top` | off | Keep top hit even if above threshold |

</details>

## Citations

This pipeline uses software and models from the following sources:

- **FlashFry** — McDonel & King. *Nucleic Acids Research*, 2015. [doi:10.1093/nar/gkw379](https://doi.org/10.1093/nar/gkw379)
- **Doench 2016 (Rule Set 2 / Azimuth)** — Fusi, Listgarten et al. *bioRxiv*, 2015. [doi:10.1101/021568](https://doi.org/10.1101/021568). Model extracted from the [crisporWebsite](https://github.com/maximilianh/crisporWebsite) repository.
- **Rule Set 3** — Doench et al. [gpp-rnd/rs3](https://github.com/gpp-rnd/rs3). Model converted from `RuleSet3.pkl`.
- **sglearn** — Feature extraction for Rule Set 3. [PyPI](https://pypi.org/project/sglearn/)
- **LightGBM** — Gradient boosting framework used for Rule Set 3 inference.
- **CRISPOR** — Haeussler et al. *Nucleic Acids Research*, 2018. [doi:10.1093/nar/gky406](https://doi.org/10.1093/nar/gky406). Model files sourced from the [crisporWebsite](https://github.com/maximilianh/crisporWebsite) repository.
- **BioPython** — Cock et al. *Bioinformatics*, 2009. [doi:10.1093/bioinformatics/btp163](https://doi.org/10.1093/bioinformatics/btp163). Used for melting temperature calculations (Tm_NN).

## Detailed documentation

- [crispy/README.md](crispy/README.md) — full Crispy documentation, scoring details, Doench 2016 model extraction
- [idmapping/README.md](idmapping/README.md) — liftover pipeline, SNP source, ortholog mapping details
