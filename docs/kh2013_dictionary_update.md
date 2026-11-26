# KH2013 gene mappings for the KY21 dictionary

Description of a derived data product: KH2013 gene identifiers joined onto
the KY21 gene dictionary. Self-contained; nothing here assumes a particular
consumer.

- Produced 2026-11-15, Unity job 1792936, 57 minutes.
- Code: `crispytargets` commit `2ed538a`, conda env `mycrispytargets`.
- Job script and QA gate: `tools/map_species.sbatch` and
  `tools/qa_species_mapping.py` (named `kh2013_map.sbatch` and
  `qa_kh2013_mapping.py` when this run was made; the KH2013 defaults
  reproduce it).

## The data

Two columns appended to the existing 11 of the KY21 dictionary, in this
order:

| Column | Meaning |
|---|---|
| `KH2013_primary` | Best KH2013 gene for the KY21 gene, or empty |
| `KH2013_secondary` | Other KH2013 genes above the threshold, pipe-separated, or empty |

- Values are gene-level KH identifiers in the `KH.C1.1`, `KH.L28.8` form,
  not isoform identifiers. The source proteome is isoform-level; the
  collapse happens after mapping.
- `KH2013_secondary` never repeats the primary. It is empty for 9,171 of
  the 14,119 assigned rows; the remaining 4,948 list at least one
  alternative gene.
- Unlike the Human and Mouse columns there are no `_symbol` or `_uniq`
  companions: KH identifiers are already the stable gene names, and GHOST
  supplies no separate symbol per gene.
- The 11 pre-existing columns are unchanged row by row, asserted by the QA
  gate.

## Files

| File | Rows | sha256 |
|---|---|---|
| `KY21_dictionary.extended.tsv` (13 columns) | 18,788 | `aa05c786e9fe7d163f0885e374a153f009d98f3417825f7afac093283296b3a2` |
| `KY21_to_KH2013.tsv` (`gene_id` plus the two columns) | 18,788 | `b5684115bea2ed05b6c827d97a69abb4162b4061e955b8bc1c79999eb48cb820` |
| `qa_report.tsv` (metric/value TSV) | 19 metrics | `cdb61e7e6147cd718f8d5e1f672d4d6dc03d1c9f56757d1fd29328f1bfbaa33b` |

Verify with `sha256sum <file>`; the values above are the contract.

Copies: `output/kh2013_mapping/` in the CRISPYtargets working tree (that
directory is not tracked by git), and on Unity under
`/work/pi_ybernadskaya_umass_edu/dtjaernberg_umass_edu/projects/crispytargets/results/kh2013_mapping/`.
Staged inputs and the regenerable BLAST artefacts lived in the
project scratch workspace during the run, not in those directories;
the workspace was released on 2026-11-26.

## How it was produced

Same shape as the existing KY19, Human and Mouse columns: KY21 proteins are
the BLAST query, the species proteome is the database.

```
blastp -query HT.KY21Gene.protein.2.fasta        # 55,505 sequences
       -db    KH.KHGene.2013.Longest.protein     # 55,017 sequences
       -outfmt 6 -evalue 1e-5 -max_target_seqs 20
crispytargets map --species KH2013 --evalue 1e-50 --output <dictionary copy>
crispytargets clean-mapping -c KH2013_primary,KH2013_secondary
```

The proteome is GHOST's `KH.KHGene.2013.Longest.protein.zip`, the same
longest-ORF convention as the KY21 file. Its sequences use `+` as a stop
marker in 1,357 lines, which `makeblastdb` rejects, so the staged copy has
`+` stripped; nothing else is altered. Input checksums are recorded in the
job's `PROVENANCE.sha256` and verified before the run starts.

Two thresholds matter when reading the result. BLAST reports hits down to
`1e-5`, and the mapper keeps only those at `1e-50` or better for this
column, matching how the earlier species columns were built.

## Quality

- 14,119 of 18,788 KY21 genes (75.2%) receive a KH2013 primary.
- 12,805 distinct KH2013 genes appear as a primary. 685 of them are the
  primary for more than one KY21 gene, one for 46, which is the expected
  consequence of recent duplication resolved in KY21 and merged in the
  older build.
- 4,669 KY21 genes (24.8%) have no KH2013 assignment at this threshold.
- Independent check against GHOST's own `KHtoKY21.txt.gz` mapping (13,906
  pairs, update 2026.11): 90.9% of its pairs appear as our primary or
  secondary, and 89.9% of our primaries agree with GHOST.

## Reading the result correctly

- The table is keyed on KY21 genes, one row each. It is not a
  KH2013-to-KY21 lookup: because 685 KH2013 genes are shared, the reverse
  direction is one-to-many and has to be built as long-form pairs.
- An existing curated cross-build chain (KY21 to KY19 to KH2013 through the
  legacy annotation table) resolves 78.7% of KY21 genes. It is different
  evidence, gene-model curation rather than protein similarity, and it is
  not a superset or a correction of this column. Mixing the two silently
  would make coverage and assignments depend on which code path ran, so a
  consumer should pick one source per use and state it.
- `KH2013_secondary` lists genes, not isoforms, and is informational. It is
  not a measure of confidence.

## Why there is no KH2012 column

GHOST also serves the 2012 release. It was evaluated and rejected because
its longest-ORF protein set is a strict subset of 2013's.

- 2012 has 54,733 sequences; every ID also exists in 2013 and the shared
  sequences are identical. 2013 adds 284 sequences.
- At gene level 15,285 identifiers are shared, 4 genes exist only in 2013
  and none only in 2012. Both releases use the same `KH.C1.1` namespace.
- Measured against this run's BLAST table, the 284 missing proteins appear
  in 4,275 hit rows and were the top hit for 223 KY21 isoforms, across 13
  of the 18,788 genes. Fourteen of those queries keep the same KH gene via
  the next-best hit, 184 lose their hit, 25 switch gene.

A KH2012 column would therefore differ from `KH2013_primary` on about 13
rows under the same identifiers. A KH2012 `KH.C` identifier can be answered
from the KH2013 column.

The numeric identifiers in older tabulations (`KH.100170058` and similar,
as in `ciona/data/genomic_resources/KH2012_KY2019.tsv`) are a different
scheme, not `KH.C` gene names. Resolving those to KY21 is a join through
KY19 (`KY19_primary`), not a new BLAST.

## Handover bundle

The three files in the table above, plus this document. They are
independent of the pipeline: a TSV join on `gene_id` is enough to use them.
The QA report carries the exact coverage and agreement figures so a
consumer can re-verify without rerunning anything.
