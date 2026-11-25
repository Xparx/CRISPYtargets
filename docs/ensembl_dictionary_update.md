# Ensembl gene mappings for the KY21 dictionary

Description of a derived data product: Ensembl gene identifiers joined onto
the KY21 gene dictionary. Self-contained; nothing here assumes a particular
consumer.

- Produced 2026-11-24, Unity job 2233878, 15 minutes.
- Code: `crispytargets` commit `fa4194b`, conda env `mycrispytargets`.
- Job script: `tools/map_species.sbatch` with `SPECIES=ENSEMBL`; QA gate:
  `tools/qa_species_mapping.py`.

## The data

Two columns appended after the KH2013 pair, taking the dictionary from 13 to
15 columns:

| Column | Meaning |
|---|---|
| `ENSEMBL_primary` | Best Ensembl gene for the KY21 gene, or empty |
| `ENSEMBL_secondary` | Other Ensembl genes above the threshold, pipe-separated, or empty |

- Values are gene-level Ensembl identifiers (`ENSCING00000004107`), not
  peptide identifiers. The source proteome is protein-level (`ENSCINP`); a
  subject map collapses peptides to their gene after BLAST.
- `ENSEMBL_secondary` never repeats the primary. It is empty for 6,934 of
  the 12,378 assigned rows; the remaining 5,444 list at least one
  alternative gene.
- There are no `_symbol` or `_uniq` companions, matching the KH2013 columns.
- The 13 pre-existing columns are unchanged row by row, asserted by the QA
  gate.

## Files

| File | Rows | sha256 |
|---|---|---|
| `KY21_dictionary.extended.tsv` (15 columns) | 18,788 | `e477bb5e4c3f9d8f08dc98032afa12a8a92bbbed3d56365d5536809c6bcf622a` |
| `KY21_to_ENSEMBL.tsv` (`gene_id` plus the two columns) | 18,788 | `7fca09eca467ab70530d4a377a7123aa87b1fbe9f502c5b5e78ce1c13c168815` |
| `qa_report.tsv` (metric/value TSV, 19 metrics) | 19 metrics | `cf60b6fb2e9075d811d6dae8d5c612335c49045f0cc5b54316ccbb5c08ce0f89` |

Verify with `sha256sum <file>`; the values above are the contract.

Copies: `output/ensembl_mapping/` in the CRISPYtargets working tree (that
directory is not tracked by git), and on Unity under
`/work/pi_ybernadskaya_umass_edu/dtjaernberg_umass_edu/projects/crispytargets/results/ensembl_mapping/`.
Staged inputs and the regenerable BLAST artefacts live in the project
scratch workspace, not in those directories.

## How it was produced

Same direction as the other species columns: KY21 proteins are the BLAST
query, the Ensembl proteome is the database.

```
blastp -query HT.KY21Gene.protein.2.fasta        # 55,505 sequences
       -db    Ciona_intestinalis.KH.pep.all      # 17,302 peptides, 16,671 genes
       -outfmt 6 -evalue 1e-5 -max_target_seqs 20
crispytargets map --species ENSEMBL --evalue 1e-50 \
    --subject-map ensembl_protein_to_gene.tsv --output <dictionary copy>
crispytargets clean-mapping -c ENSEMBL_primary,ENSEMBL_secondary
```

The proteome is Ensembl release 116 for *Ciona intestinalis*, assembly KH
(GCA_000224145.1), `pep.all`: 17,302 peptides over 16,671 protein-coding
genes, one to a few isoforms per gene. Its alphabet is the 20 standard
residues plus X, so unlike the KH2013 file no stop markers had to be
stripped. The subject map (`protein_id`, `gene_id`) comes from the FASTA
headers and is what turns peptide hits into gene identifiers; without it the
`ENSCINP` accessions would pass through untouched.

Both thresholds are the same as the KH2013 run: BLAST reports hits down to
`1e-5`, and the mapper keeps those at `1e-50` or better.

## Quality

- 12,378 of 18,788 KY21 genes (65.9%) receive an Ensembl primary, over
  11,145 distinct Ensembl genes; 13,973 distinct genes appear in the primary
  or secondary column.
- This is 1,743 rows more than the existing KY19-chained bridge (10,635,
  56.6%), the only Ensembl mapping available before this run.
- Against that bridge: 90.8% of its 10,635 pairs appear as our primary or
  secondary; 80.1% of the 10,282 genes mapped by both agree exactly; and
  90.4% of our primaries fall inside the bridge's candidate set, which lists
  several genes for 4,110 KY21 genes.
- Coverage is below the KH2013 column's 75.2% because the Ensembl proteome
  carries about one protein per gene (17,302 for 16,671 genes), while the
  KH2013 proteome carries 55,017 isoforms. BLAST has less sequence to find,
  not fewer genes to compare.

## Why this could not be a straight join

Ensembl serves the same assembly as KH: 1,128 of the 1,272 sequences are
md5-identical to `KH2013.fa`, including chromosome 1 at 10,041,005 bp. It is
still not a rename of the KH gene models.

- The GFF3, GTF and peptide headers carry no KH or KY identifiers, only
  `ENSCING` and `ENSCINP`; the ENA cross-reference links genes to KH
  scaffolds and to old `FAA` protein accessions, not to KH gene IDs.
- The genebuild has 16,671 genes against 15,288 in KH2013, with a different
  merge and split structure: 11,012 KH genes overlap an Ensembl gene, 4,479
  of them overlap more than one, and only 378 of 17,302 Ensembl peptides are
  identical to a KH2013 protein.
- A coordinate-only chain (KY21 to KH2013 to Ensembl, reciprocal and
  unambiguous) reached 5,923 KY21 genes (31.5%) and agreed with the bridge
  only 54% where it did map, so it is not a substitute for the BLAST.

## Reading the result correctly

- The table is keyed on KY21 genes. The reverse direction is one-to-many.
- The KH2013 and Ensembl columns are independent protein-similarity maps,
  so where both carry a gene they need not name the same locus. Of the
  10,140 such pairs whose Ensembl gene sits on a scaffold shared with the
  KH2013 file, 94.1% are on the same chromosome or scaffold and 72.2% lie
  within 10 kb; the two builds resolve tandem families differently and KY21
  carries more copies in several of them. Neither column is ground truth,
  and the secondary column is informational.
- `ENSEMBL_secondary` lists genes, not peptides, and is not a confidence
  measure.

## Verification

Checked 2026-11-24, after the run:

- The QA gate reports 18,788 rows before and after and zero mismatches
  across the 13 pre-existing columns, so the two Ensembl columns are the
  only change.
- Columns 1-13 of the extended dictionary hash to `aa05c786…`, the same
  value as the dictionary the run started from.
- The local copies in `output/ensembl_mapping/` and the cluster copies in
  `/work/.../results/ensembl_mapping/` hash to the values in the table
  above.
- `sha256sum -c CHECKSUMS.sha256` passes in
  `output/ensembl_mapping/handoff/`.
- The 63 tests pass locally and on the cluster at commit `fa4194b`, the
  commit the run used; only documentation changed afterwards.

## Handover bundle

The three files in the table above, plus this document. A TSV join on
`gene_id` is enough to use them, and the QA report carries the exact
coverage and agreement figures.
