# KH2013 columns for the KY21 dictionary

Run date 2026-11-15. Unity job 1792936, 57 minutes, `cpu` partition.
Code: `crispytargets` at commit `2ed538a`, env `mycrispytargets`.
Job script: `tools/kh2013_map.sbatch`. QA gate: `tools/qa_kh2013_mapping.py`.

## What was added

Two columns appended to the existing 11, in this order:

| Column | Meaning |
|---|---|
| `KH2013_primary` | Best KH2013 gene ID for the KY21 gene, or empty |
| `KH2013_secondary` | Other KH2013 genes passing the threshold, pipe-separated, or empty |

Values are gene-level KH IDs (`KH.C1.1`), not isoform IDs. The 11 existing
columns are unchanged, which the QA gate asserts row by row.

## Files

| File | sha256 |
|---|---|
| `KY21_dictionary.extended.tsv` (13 columns, 18,788 rows) | `aa05c786e9fe7d163f0885e374a153f009d98f3417825f7afac093283296b3a2` |
| `KY21_to_KH2013.tsv` (`gene_id`, the two columns) | `b5684115bea2ed05b6c827d97a69abb4162b4061e955b8bc1c79999eb48cb820` |
| `qa_report.tsv` | metrics below |

Cluster: `/work/pi_ybernadskaya_umass_edu/dtjaernberg_umass_edu/projects/crispytargets/results/kh2013_mapping/`
Bulk (BLAST database and table) and staged inputs:
`/scratch3/workspace/dtjaernberg_umass_edu-crispytargets/kh2013/`

## Method

Same shape as the KY19, Human and Mouse columns: KY21 proteins are the
BLAST query, the species proteome is the database.

```
blastp -query HT.KY21Gene.protein.2.fasta          # 55,505 sequences
       -db    KH.KHGene.2013.Longest.protein      # 55,017 sequences
       -outfmt 6 -evalue 1e-5 -max_target_seqs 20
crispytargets map --species KH2013 --evalue 1e-50 --output <dictionary copy>
crispytargets clean-mapping -c KH2013_primary,KH2013_secondary
```

The KH2013 proteome is GHOST's `KH.KHGene.2013.Longest.protein.zip`, the
same longest-ORF convention as the KY21 file. Its sequences carry `+` as a
stop marker (1,357 lines), which `makeblastdb` rejects, so the staged copy
has `+` stripped; nothing else is altered.

## Results

- 14,119 of 18,788 KY21 genes (75.2%) receive a KH2013 primary.
- 12,805 distinct KH2013 genes appear as a primary; 1,314 KH2013 genes are
  primary for more than one KY21 gene, which is expected where the older
  build merged recent duplicates.
- 4,948 assigned rows carry a genuine alternative gene in
  `KH2013_secondary`; 9,171 have no other gene passing the threshold.
- Against GHOST's own `KHtoKY21.txt.gz` (13,906 pairs, update 2026.11):
  90.9% of its pairs appear as our primary or secondary, and 89.9% of our
  primaries agree with GHOST.
- 4,669 KY21 genes (24.8%) have no KH2013 assignment at this threshold.

Secondary hits are gene-level and never repeat the primary. Two causes
were fixed before these files were finalised: the mapper could list the
same accession twice under different subject formats (`sp|P1|GENE` and
`tr|P1|GENE`), and collapsing isoform IDs to locus IDs could turn a list of
alternative isoforms into the primary locus itself. The mapper now
deduplicates stripped accessions, and `clean-mapping` drops the primary
from the secondary column.

## Why there is no KH2012 column

GHOST also serves the 2012 release, and it was evaluated and rejected: its
longest-ORF protein set is a strict subset of 2013's.

- 2012 has 54,733 sequences; every ID also exists in 2013 and the shared
  sequences are byte-identical. 2013 adds 284 sequences.
- At gene level 15,285 IDs are shared, 4 genes exist only in 2013 and none
  only in 2012. Both releases use the same `KH.C1.1` namespace.
- Measured against this run's BLAST table, the 284 missing proteins appear
  in 4,275 hit rows and were the top hit for 223 KY21 isoforms, which
  belong to 13 of the 18,788 genes. Fourteen of those queries keep the same
  KH gene through the next-best hit, 184 lose their hit entirely, and 25
  switch gene.

A KH2012 column would therefore differ from `KH2013_primary` on roughly 13
rows, under the same identifiers, so it is not built. The consequence for
consumers: a KH2012 gene ID can be answered from the KH2013 column, since
the namespaces coincide.

The numeric identifiers in older tabulations (`KH.100170058` and similar,
as in `ciona/data/genomic_resources/KH2012_KY2019.tsv`) are a different
scheme, not the `KH.C` gene names. Mapping those to KY21 is a join through
KY19 (`KY19_primary`), not a new BLAST.

## What this means for cionalib

The current `ky21_to_kh2013()` walks KY21 -> KY19 -> KH2013 through
`feature.KH.KY.anno.tsv` and resolves 14,781 genes (78.7%), slightly more
than the 75.2% here. The two are different evidence, not a correction of
one by the other: the annotation chain follows gene-model curation, this
follows protein similarity. Options:

1. Replace `cionalib/data/KY21_dictionary.tsv` with the extended file and
   leave `ky21_to_kh2013()` on the two-hop path. Nothing else changes; the
   new columns become available for direct joins.
2. Switch `ky21_to_kh2013()` to read `KH2013_primary`. This changes
   coverage from 78.7% to 75.2% and makes `annotate_ghost_ids()`
   inconsistent with it, because that function still maps through
   `feature.KH.KY.anno.tsv`; the test `annotate_ghost_ids() disagrees with
   ky21_to_kh2013()` will fail unless both paths move together.
3. Keep both: rename the existing function to `ky21_to_kh2013_curated()`
   and add `ky21_to_kh2013_blast()`, so callers choose the evidence.

Option 3 is the honest one if downstream analyses depend on the current
coverage; option 1 is the smallest change.
