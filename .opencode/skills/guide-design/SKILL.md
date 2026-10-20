---
name: Guide Design
description: sgRNA design and scoring runs — prepare targets, run the crispy pipeline, Doench 2016 Rule Set 3 scoring, FlashFry off-target analysis, SNP/flank/accessibility overlays, ranked candidate output. Use when designing, scoring, or ranking CRISPR guides or changing the design pipeline.
---

# Guide Design

## Roles

Full role definitions in `roles.org`. Read them if not already
in context. The active roles for this workflow are:

- Expert Bioinformatics Researcher
- Expert Software Developer

Single unified voice, expertise internal, no role attribution
in output.

## When Active

- Preparing targets from GFF3/FASTA (CDS, promoters, TSS
  features) with the =crispytargets design= subcommand
- On-target scoring: Doench 2016 Rule Set 3 feature extraction,
  GC normalization, model artifacts
- Off-target analysis through FlashFry (Java JAR, .dict/.fasta
  companion files)
- Overlaying guide-level annotations: SNPs, strand-aware 5'
  flanks, accessibility, PAM geometry
- Ranking and exporting candidate guides (TSV outputs under
  =output/=)
- Changing scoring behavior, feature order, or pipeline
  contracts in =crispy/=

## Key Rules

- Biology first: a guide set where every gene scores near
  maximum means a feature-extraction bug, not a good day
- Model-transfer honesty: scores come from human-trained
  models; state that and suggest orthogonal checks (GC,
  specificity, manual inspection)
- RULES.org "Research Integrity Rules": every reported number
  traces to a run recorded in =docs/provenance.md= (inputs,
  parameters, code version)
- Regression tests travel with behavior changes: feature order,
  strand-aware flanks, and PAM handling have fixed real bugs
  before (=tests/=, run with =python -m pytest tests/= in the
  =mainscworkspace= conda env)
- New tabular writes use TSV per RULES.org "Data Format Rules"
- Propose before modifying: state what needs to change and why,
  confirm before writing

