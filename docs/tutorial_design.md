# Tutorial: design guide RNAs

**Deliverable:** `final_sgrnas.tsv`, a ranked table of candidate sgRNAs with
efficiency scores, off-target counts, oligos for cloning, and (optionally)
SNP-aware columns.

## What you need

- `crispytargets` installed:
  `pip install "git+https://github.com/Xparx/CRISPYtargets.git"`
  (add Rule Set 3 scoring with
  `pip install "crispytargets[rs3] @ git+https://github.com/Xparx/CRISPYtargets.git"`)
- Java (any recent JRE), for FlashFry:
  `conda install -c conda-forge openjdk` or your system package manager
- The FlashFry jar, downloaded once from
  https://github.com/mckennalab/FlashFry/releases/download/1.15/FlashFry-assembly-1.15.jar
  (it is not shipped with the package; keep it anywhere and pass
  `--flashfry-jar`, or export `CRISPY_FLASHFRY_JAR=/path/to/it`)

Set up your workspace:

```bash
working_dir=~/ciona_design                     # change to taste
mkdir -p "$working_dir"/inputs "$working_dir"/outputs
```

Download into `"$working_dir"/inputs/`:

| Input | Where it comes from |
|-------|---------------------|
| Genome FASTA + gene annotation (GFF3) | for Ciona: https://ghost.zool.kyoto-u.ac.jp/download_ht.html (the direct files are `HT.Ref.fasta.zip` and `HT.KY21Gene.2.gff3.zip`); for any other organism, its own assembly and annotation |
| Optional: a gene list, a one-column file with the header `Gene.ID` and one gene ID per line, to design for specific genes only | you write this; it can live anywhere |
| Optional: a SNP BED, from [tutorial_snps.md](tutorial_snps.md) | SNP data is optional and slow, see below |

## Step 1: the run

Without SNP data this takes about 3 minutes for the Ciona example (24,800
candidate guides found and scored). Discovery and scoring scale with the size
of the region you target: a whole-genome run, with no `--genes` filter, works
through vastly more sequence and takes proportionally longer.

```bash
crispytargets design \
    --gff "$working_dir"/inputs/HT.KY21Gene.2.gff3 \
    --ref "$working_dir"/inputs/HT.Ref.fasta \
    --outdir "$working_dir"/outputs/guides \
    --flashfry-jar /path/to/FlashFry-assembly-1.15.jar
```

Useful options:

- `--genes "$working_dir"/inputs/my_genes.csv` to restrict to specific genes
- `--strict` to keep only guides above the score cutoff (default: keep
  everything, ranked, with an `over_cutoff` flag)
- `--browser-bed` to additionally write `final_sgrnas.bed` for IGV/UCSC

## Step 2: SNP awareness (optional, slow)

Add the SNP BED from [tutorial_snps.md](tutorial_snps.md) to get variant-aware
columns:

```bash
crispytargets design \
    ... \
    --snps "$working_dir"/outputs/snps_lifted/Ghost_SNPs_main.bed_lifted.bed \
    --rescore-snps \
    --rs3
```

What the two flags cost:

- `--snps` alone is cheap (seconds): it adds `num_snps`, `snp_overlaps`,
  `entropy`, `snp_pam_dist`, `snp_cut_dist`.
- `--rescore-snps` is the expensive part: every SNP-overlapping guide is
  re-scored with each variant allele present, adding the twelve variant
  columns (average / worst single-SNP / all-SNP scores and their deltas). It
  added about 18 minutes to a run over **5 genes** (24,800 candidate guides).
  Cost scales with how many guides you scan, so budget accordingly: a
  **whole-genome run (10,000 genes) is a different order of effort** and the
  rescore can go from minutes to days. For genome-wide work either leave
  `--rescore-snps` off and keep the cheap `--snps` annotations, or narrow the
  run with `--genes` to the loci you actually intend to edit.
- `--rs3` adds the Rule Set 3 score (needs the `[rs3]` extra installed).

If your species has no variant data, leave all three out; nothing else in the
output changes.

## What you get

In `"$working_dir"/outputs/guides/`:

| File | What it is |
|------|-----------|
| **`final_sgrnas.tsv`** | the deliverable: one row per candidate guide |
| `final_sgrnas.bed` | the same guides coloured by status, for a genome browser (`--browser-bed`) |
| `design.log` | the run log, ending with a timing line (prepare / FlashFry / post-process / total) |

The columns are grouped: identity, coordinates, reference scores
(`Doench2016OnTarget`, `RuleSet3`, `Hsu2013`, `otCount`, `accessibility`),
guide sequence and oligos, SNP annotation, variant scores, and finally any
extra FlashFry metrics you requested with `--metrics`. Full column reference:
[design.md](design.md).

Scoring, in one paragraph: efficiency comes from our own Doench 2016
Rule Set 2 (and optionally Rule Set 3) implementations; off-target reporting
comes from FlashFry (`Hsu2013` specificity score and `otCount`, always on).
FlashFry's own on-target metrics are opt-in via `--metrics` and appended as
trailing columns.

## Optional: attach ortholog IDs to the table

If you are working on Ciona and want a mouse ortholog ID as column 2 (for
example `Depdc1b-KY21.Chr2.1411`), join the mapping dictionary that cionalib
ships:

```bash
python - <<'EOF'
import pandas as pd
p = "outputs/guides/final_sgrnas.tsv"
df = pd.read_csv(p, sep="\t")
m = pd.read_csv("/path/to/KY21_dictionary.tsv", sep="\t").set_index("gene_id")["Mouse_uniq"]
df.insert(1, "Mouse_uniq", df["gene_id"].str.split(":", n=1).str[-1].map(m).fillna(""))
df.to_csv(p, sep="\t", index=False)
EOF
```

The dictionary itself comes from [tutorial_mapping.md](tutorial_mapping.md)
or pre-built from [cionalib](https://gitlab.com/Xparx/cionalib).