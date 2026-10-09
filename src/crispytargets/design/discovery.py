"""FlashFry integration: index, discover, and score sgRNA targets.

Wraps the three FlashFry stages the design pipeline runs: building the
genomic index, discovering targets in extracted sequences, and scoring
them. Sequence extraction from the reference lives here too, since it
produces the FASTA that discovery consumes.
"""

import logging
import os

from pyfaidx import Fasta

from ..common.exec import run_command


def extract_sequences(targets_gr, outdir, ref, gene_map):
    """Write target sequences to FASTA for FlashFry discovery.

    Headers are ``gene_id|feature|chrom:start-end|strand`` so every guide
    can be traced back to its origin after scoring. gene_map resolves each
    feature ID to its root gene.
    """
    logging.info("Extracting sequences for FlashFry...")
    fasta_out = os.path.join(outdir, "targets.fasta")
    genome = Fasta(ref)

    with open(fasta_out, 'w') as f_out:
        for row in targets_gr.df.itertuples():
            # 1. Identify the feature's own ID or its immediate parent
            fid = getattr(row, 'ID', None)
            if fid is None or str(fid) == 'nan':
                fid = getattr(row, 'Parent', None)

            # 2. Use the gene_map to resolve this ID to the top-level Gene ID
            gene_id = gene_map.get(str(fid), "unknown")

            seq = genome[row.Chromosome][row.Start:row.End].seq
            # Include row.Feature in the header to track origin (CDS, exon, promoter)
            f_out.write(f">{gene_id}|{row.Feature}|{row.Chromosome}:{row.Start}-{row.End}|{row.Strand}\n{seq}\n")

    return fasta_out


def index_database(ff_jar, ref, outdir):
    """Build the FlashFry index, or reuse it if already present."""
    db = os.path.join(outdir, "htky_db")
    if os.path.exists(db):
        logging.info("FlashFry database already exists: %s", db)
        return db
    # Create a temporary directory for FlashFry indexing
    os.makedirs(os.path.join(outdir, "build"), exist_ok=True)
    run_command(f"java -Xmx4g -jar {ff_jar} index --database {db} --reference {ref} "
                f"--enzyme spcas9ngg --tmpLocation {os.path.join(outdir, 'build')}",
                verbose=True)
    return db


def discover_and_score(ff_jar, db, ref, outdir, targets_gr, args, gene_map):
    """Run FlashFry discovery and scoring, returning the scored file path."""
    logging.info("Running FlashFry...")
    index_database(ff_jar, ref, outdir)

    # Extract sequences from PyRanges to FASTA
    targets_fasta = extract_sequences(targets_gr, outdir, ref, gene_map)

    # Discover
    run_command(f"java -Xmx{args.mem} -jar {ff_jar} discover --database {db} "
                f"--fasta {targets_fasta} --positionOutput "
                f"--output {os.path.join(outdir, 'targets.output')}",
                verbose=getattr(args, 'verbose', False))

    # Score
    run_command(f"java -Xmx{args.mem} -jar {ff_jar} score "
                f"--input {os.path.join(outdir, 'targets.output')} "
                f"--output {os.path.join(outdir, 'targets.output.scored')} "
                f"--database {db} --scoringMetrics={args.metrics}",
                verbose=getattr(args, 'verbose', False))

    return os.path.join(outdir, "targets.output.scored")


def rescore_with_snps(ff_jar, db, outdir, args, snp_fasta):
    """Re-run discovery and scoring on SNP-modified sequences.

    Returns the path to the scored output, which the caller merges back
    for guides that directly overlap SNPs.
    """
    logging.info("Running FlashFry on SNP-modified sequences...")
    snp_output = os.path.join(outdir, 'targets_snp.output')
    snp_scored = os.path.join(outdir, 'targets_snp.output.scored')

    run_command(f"java -Xmx{args.mem} -jar {ff_jar} discover --database {db} "
                f"--fasta {snp_fasta} --positionOutput --output {snp_output}",
                verbose=getattr(args, 'verbose', False))
    run_command(f"java -Xmx{args.mem} -jar {ff_jar} score --input {snp_output} "
                f"--output {snp_scored} --database {db} --scoringMetrics={args.metrics}",
                verbose=getattr(args, 'verbose', False))

    return snp_scored