---
name: Variant Mapping
description: Assembly liftover and ID mapping — minimap2 chain lifting, indel-aware SNP position transfer, BLAST ortholog inference, KY21/KH2013 cross-references via the map and lift subcommands. Use when lifting variants between assemblies or mapping gene IDs across sources.
---

# Variant Mapping

## Roles

Full role definitions in `roles.org`. Read them if not already
in context. The active roles for this workflow are:

- Expert Bioinformatics Researcher
- Expert Software Developer

Single unified voice, expertise internal, no role attribution
in output.

## When Active

- Reference-guided liftover with minimap2 (=asm5=, =asm10=,
  =asm20= presets, =-c= for CIGAR) via =crispytargets lift=
- Transferring SNP positions across assemblies with
  =idmapping/lift_snps.py= (indel-aware block splitting)
- Ortholog inference from BLAST results (evalue thresholds,
  top-hit retention) and KY21/KH2013 ID cross-references via
  =crispytargets map= and =clean-mapping=
- Validating lifted coordinates against variants near indels
  and structural differences
- Changing chain/Paf parsing, block logic, or mapping output
  contracts

## Key Rules

- Positions after an indel belong to a new collinear block:
  split records on CIGAR, never treat a whole Paf record as
  one block (this bug shipped once)
- minimap2 must run with =-c=; without CIGAR output the lift
  cannot be validated
- Model/assembly honesty: state which assembly and annotation
  version every lifted position came from and what was lost in
  the lift
- RULES.org "Research Integrity Rules": lifted outputs land in
  =docs/provenance.md= with inputs, parameters, and code
  version before any number from them is reported
- Round-trip tests (=tests/test_lift_snps= style) accompany
  changes to lifting or mapping logic
- Propose before modifying: state what needs to change and why,
  confirm before writing

