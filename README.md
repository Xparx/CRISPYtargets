# CRISPYtargets

CRISPR guide RNA design pipeline for non-model organisms, with genomic ID mapping tools.

Two largely independent capabilities, packaged together because both exist to
support work in organisms without established tooling:

- **Design**: discover, score, and rank sgRNAs from a GFF3 + FASTA pair
- **Mapping**: lift SNP positions between assemblies, and map gene IDs to
  orthologs in a well-annotated reference

## Installation

```bash
conda create -n crispy python=3.10 minimap2 bcftools tabix -c conda-forge -c bioconda
conda activate crispy
pip install -e .            # core: design + mapping
pip install -e ".[rs3]"     # optional Rule Set 3 scoring (--rs3)
pip install -e ".[test]"    # the pytest suite
```

Python dependencies are declared in `pyproject.toml`. Three external tools are
required at runtime and are not installed by pip:

| Tool | Used by | Install |
|------|---------|---------|
| Java (JRE) | FlashFry, used by `design` | system package |
| FlashFry JAR | `design` | [releases](https://github.com/mckennalab/FlashFry/releases) |
| minimap2, tabix | `lift` | conda (above) |
| BLAST+ (`makeblastdb`, `blastp`) | `map` | `conda install -c bioconda blast` |

The FlashFry JAR is an untracked binary and is **not** shipped inside the
package. Tell the pipeline where it is with either:

```bash
export CRISPY_FLASHFRY_JAR=/path/to/FlashFry-assembly-1.15.jar
# or pass --flashfry-jar /path/to/FlashFry-assembly-1.15.jar per run
```

## Repository layout

```
src/crispytargets/
    design/     sgRNA discovery and scoring (FlashFry, Doench 2016 RS2/RS3)
    mapping/    SNP liftover and BLAST ortholog mapping
    cli/        the `crispytargets` subcommand dispatcher
tests/          pytest suite
tools/          developer-only scripts (model conversion, training)
examples/       run_crispy.sh workflow template
docs/           detailed per-capability documentation
```

## Command line

One command with a subcommand per capability:

```
crispytargets design         find and score sgRNAs
crispytargets map            build an ortholog table from BLAST results
crispytargets lift           lift SNP positions between assemblies
crispytargets clean-mapping  collapse isoform IDs to locus IDs
```

Every subcommand can also be run without installing, as
`python -m crispytargets <command> ...`.

### Workflow 1: Design sgRNAs

```bash
crispytargets design \
    --gff /path/to/genome.gff3 \
    --ref /path/to/genome.fa \
    --outdir output/sgRNAs
```

With a gene list, SNPs, and browser output:

```bash
crispytargets design \
    --gff /path/to/genome.gff3 \
    --ref /path/to/genome.fa \
    --genes /path/to/genes.csv \
    --snps /path/to/lifted_snps.bed \
    --outdir output/sgRNAs \
    --rescore-snps --browser-bed
```

### Workflow 2: Lift SNPs between assemblies

```bash
crispytargets lift \
    --old-genome /path/to/old_assembly.fa \
    --new-genome /path/to/new_assembly.fa \
    --vcf /path/to/snps.vcf.gz \
    --output output/liftover \
    --bed --preset asm10
```

### Workflow 3: Map orthologs

```bash
# One-time: build BLAST databases
makeblastdb -in human_proteome.fasta -dbtype prot -out output/human_db

# Search (outfmt 6 is required)
blastp -query your_proteins.fasta -db output/human_db \
    -out output/human_results.tsv -outfmt 6 -evalue 1e-5

# Map one species at a time into a master table
crispytargets map \
    --gff /path/to/genome.gff3 \
    --blast-file output/human_results.tsv \
    --species Human \
    --output output/orthologs.tsv
```

### Full pipeline example

```bash
# 1. Lift SNPs onto the assembly you will design against
crispytargets lift \
    --old-genome /path/to/old_assembly.fa \
    --new-genome /path/to/new_assembly.fa \
    --vcf /path/to/snps.vcf.gz \
    --output output/liftover \
    --bed --preset asm10

# 2. Design guides with SNP awareness
crispytargets design \
    --gff /path/to/genome.gff3 \
    --ref /path/to/new_assembly.fa \
    --snps output/liftover/*_lifted.bed \
    --outdir output/sgRNAs \
    --rescore-snps --browser-bed
```

## Detailed documentation

- [docs/tutorial_design.md](docs/tutorial_design.md),
  [docs/tutorial_snps.md](docs/tutorial_snps.md),
  [docs/tutorial_mapping.md](docs/tutorial_mapping.md) — start-to-finish
  tutorials, each from downloads in your own working directory to one
  deliverable (guide table, lifted SNP track, ortholog table)
- [docs/provenance.md](docs/provenance.md) — checksums and the verification
  log behind the Ciona worked example
- [docs/design.md](docs/design.md) — target preparation, promoter windows,
  accessibility vs. hard filtering, output files and columns, Doench 2016
  Rule Set 2 and Rule Set 3 scoring details
- [docs/mapping.md](docs/mapping.md) — BLAST ortholog mapping and its
  incremental backbone, SNP liftover internals, preset selection, SNP source
  and preprocessing for *Ciona robusta*
- [docs/cluster.md](docs/cluster.md) — running the package on the Unity
  cluster: layout, the scratch rule, the workspace, and the job pattern
- [docs/kh2013_dictionary_update.md](docs/kh2013_dictionary_update.md) —
  the KH2013 columns added to the KY21 dictionary: method, checksums,
  quality figures, and how to read them

## Testing

```bash
pip install -e ".[test]"
python -m pytest tests/
```

The smoke and liftover tests skip themselves when java/FlashFry or minimap2
are absent. Golden Doench 2016 scores live in
`tests/fixtures/doench2016_scores.tsv` and should only be regenerated when a
scoring change is intentional
(`python tests/fixtures/regenerate_doench2016_scores.py`).

## Citations

This pipeline uses software and models from the following sources:

- **FlashFry** — McDonel & King. *Nucleic Acids Research*, 2015. [doi:10.1093/nar/gkw379](https://doi.org/10.1093/nar/gkw379)
- **Doench 2016 (Rule Set 2 / Azimuth)** — Fusi, Listgarten et al. *bioRxiv*, 2015. [doi:10.1101/021568](https://doi.org/10.1101/021568). Model extracted from the [crisporWebsite](https://github.com/maximilianh/crisporWebsite) repository.
- **Rule Set 3** — Doench et al. [gpp-rnd/rs3](https://github.com/gpp-rnd/rs3). Model converted from `RuleSet3.pkl`.
- **sglearn** — Feature extraction for Rule Set 3. [PyPI](https://pypi.org/project/sglearn/)
- **LightGBM** — Gradient boosting framework used for Rule Set 3 inference.
- **CRISPOR** — Haeussler et al. *Nucleic Acids Research*, 2018. [doi:10.1093/nar/gky406](https://doi.org/10.1093/nar/gky406). Model files sourced from the [crisporWebsite](https://github.com/maximilianh/crisporWebsite) repository.
- **BioPython** — Cock et al. *Bioinformatics*, 2009. [doi:10.1093/bioinformatics/btp163](https://doi.org/10.1093/bioinformatics/btp163). Used for melting temperature calculations (Tm_NN).
- **CRISPRtargets** — [ChristiaenLab/CRISPRtargets](https://github.com/ChristiaenLab/CRISPRtargets). This project was originally a Python port of it.
