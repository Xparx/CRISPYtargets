# Tutorial: map your gene IDs to human and mouse orthologs

**Deliverable:** `orthologs.tsv`, a table with one row per gene in your
organism and columns naming its human and mouse orthologs (with gene symbols
and a combined `<symbol>-<gene_id>` label). For Ciona this table *is* the
`KY21_dictionary.tsv` that [cionalib](https://gitlab.com/Xparx/cionalib)
ships behind `cionalib.ky21_to_human()` and `cionalib.ky21_to_mouse()`.

## What you need

- `crispytargets` installed:
  `pip install git+https://github.com/Xparx/CRISPYtargets.git`
- BLAST+ (`conda install -c bioconda blast`)
- A source genome annotation and its protein FASTA (your organism)
- Target proteomes, e.g. from https://www.uniprot.org: search
  `reviewed:true AND organism_name:"Homo sapiens"`, Download → FASTA; same
  for mouse. Keep the downloaded file names, they are referenced below.

Set up your workspace:

```bash
working_dir=~/my_mapping                      # change to taste
mkdir -p "$working_dir"/inputs "$working_dir"/outputs "$working_dir"/blast
```

Example inputs for Ciona (put them in `"$working_dir"/inputs/`):

| Input | Where it comes from |
|-------|---------------------|
| `HT.KY21Gene.2.gff3` (source annotation) | https://ghost.zool.kyoto-u.ac.jp/datas/HT.KY21Gene.2.gff3.zip |
| `HT.KY21Gene.protein.2.fasta` (source proteins, the BLAST query) | https://ghost.zool.kyoto-u.ac.jp/datas/HT.KY21Gene.protein.2.fasta.zip |
| `uniprotkb_Human_...fasta`, `uniprotkb_Mouse_...fasta` | UniProt, as above |
| optional: old-build proteins, to map across builds of your own genome | for Ciona: https://ghost.zool.kyoto-u.ac.jp/datas/HT.KYGene.protein.fasta.zip |

## Step 1: build BLAST databases and search

One database per target species, then one search per target species. The
**query is always your source proteins**, so the hits carry your organism's
gene IDs:

```bash
cd "$working_dir"
makeblastdb -in inputs/uniprotkb_Human_*.fasta -dbtype prot -out blast/human_db
makeblastdb -in inputs/uniprotkb_Mouse_*.fasta -dbtype prot -out blast/mouse_db

blastp -query inputs/HT.KY21Gene.protein.2.fasta -db blast/human_db \
       -evalue 1e-5 -max_target_seqs 20 -outfmt 6 -out blast/human_vs_source.tsv
blastp -query inputs/HT.KY21Gene.protein.2.fasta -db blast/mouse_db \
       -evalue 1e-5 -max_target_seqs 20 -outfmt 6 -out blast/mouse_vs_source.tsv
```

`-outfmt 6` is required (tabular, twelve columns). The `map` step below
applies the stricter significance filter, so the search settings here are
forgiving. Runtime is dominated by these searches; the mapping steps after
them take seconds.

## Step 2: map, one run per species

The table is incremental: each run adds columns for one species to the same
file. **Run order sets the column order of the result.**

```bash
crispytargets map --gff inputs/HT.KY21Gene.2.gff3 \
    --blast-file blast/human_vs_source.tsv --species Human \
    --proteome inputs/uniprotkb_Human_*.fasta \
    --output outputs/orthologs.tsv --evalue 1e-50

crispytargets map --gff inputs/HT.KY21Gene.2.gff3 \
    --blast-file blast/mouse_vs_source.tsv --species Mouse \
    --proteome inputs/uniprotkb_Mouse_*.fasta \
    --output outputs/orthologs.tsv --evalue 1e-50
```

`--proteome` adds the `{Species}_symbol` and `{Species}_uniq` columns from the
UniProt `GN=` fields; drop it if you only want UniProt entry names. Mapping to
an older build of your own genome (the Ciona KY19 column) is identical: make a
database from the old proteins, blast against it, and run a third `map` with
`--species KY19`.

If one of your input columns carries isoform-level IDs, collapse them to gene
level:

```bash
crispytargets clean-mapping -i outputs/orthologs.tsv -o outputs/orthologs.tsv \
    -c "KY19_primary,KY19_secondary"
```

## What you get

`outputs/orthologs.tsv`, one row per source gene:

| Column | Meaning |
|--------|---------|
| `gene_id` | your organism's gene ID, from the GFF3 |
| `{Species}_primary` | the best ortholog hit |
| `{Species}_secondary` | every other hit that passes the threshold, `\|`-separated (often the recent duplicates) |
| `{Species}_symbol` | gene symbol, with `--proteome` (e.g. `Fgf13`) |
| `{Species}_uniq` | `<symbol>-<gene_id>`, a label that is unique and readable at once (e.g. `Fgf13-KY21.Chr1.1`) |

Two things worth knowing when you read the results: several of your genes
pointing at one human gene usually means your genome has recent duplicates,
not a mapping error; and assigning orthologs across species is a hypothesis,
so treat `primary` as the candidate and `secondary` as the alternatives.