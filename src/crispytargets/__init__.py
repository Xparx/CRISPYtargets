"""CRISPYtargets: CRISPR guide RNA design and genomic ID mapping.

Two largely independent capabilities, packaged together because both exist
to support work in organisms without established tooling:

- :mod:`crispytargets.design`  find and score sgRNAs from a GFF3 + FASTA pair
- :mod:`crispytargets.mapping` lift SNPs between assemblies and map gene IDs
  to orthologs in a well-annotated reference

The command-line entry point is :mod:`crispytargets.cli`.
"""
