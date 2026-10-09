"""CRISPYtargets command-line interface.

One command with a subcommand per capability:

    crispytargets design         find and score sgRNAs
    crispytargets map            build an ortholog table from BLAST results
    crispytargets lift           lift SNP positions between assemblies
    crispytargets clean-mapping  collapse isoform IDs to locus IDs

Argument parsing stays with each tool; this module only routes, so every
flag keeps the behaviour it had when the tools were standalone scripts.
"""

import argparse
import sys

from ..design.pipeline import main as design_main
from ..mapping.clean_mapping import main as clean_main
from ..mapping.lift_snps import main as lift_main
from ..mapping.run_mapping import main as map_main

PROG = "crispytargets"

SUBCOMMANDS = (
    ("design", "Find and score sgRNAs from a GFF3 + FASTA pair", design_main),
    ("map", "Build an ortholog table from BLAST outfmt 6 results", map_main),
    ("lift", "Lift SNP positions between genome assemblies", lift_main),
    ("clean-mapping", "Collapse isoform IDs to locus IDs in a mapping table", clean_main),
)


def build_root_parser():
    """Root parser, used only to render the subcommand list in --help."""
    parser = argparse.ArgumentParser(
        prog=PROG,
        description="CRISPR guide RNA design and genomic ID mapping "
                    "for non-model organisms.",
    )
    subparsers = parser.add_subparsers(dest="command", metavar="<command>")
    for name, help_text, _ in SUBCOMMANDS:
        subparsers.add_parser(name, help=help_text, description=help_text, add_help=False)
    return parser


def main(argv=None):
    """Route ``argv[0]`` to the matching tool."""
    argv = list(sys.argv[1:] if argv is None else argv)
    handlers = {name: handler for name, _, handler in SUBCOMMANDS}

    if argv and argv[0] in ("-h", "--help"):
        build_root_parser().print_help()
        return 0

    if not argv:
        sys.stderr.write(f"{PROG}: a <command> is required\n\n")
        build_root_parser().print_help(sys.stderr)
        return 2

    if argv[0] not in handlers:
        sys.stderr.write(f"{PROG}: unknown command {argv[0]!r}\n\n")
        build_root_parser().print_help(sys.stderr)
        return 2

    # Repoint sys.argv so the tool's own parser renders
    # `usage: crispytargets <command> ...` and sees only its own arguments.
    sys.argv = [f"{PROG} {argv[0]}"] + argv[1:]
    return handlers[argv[0]]()
