"""Quality report for a species column added to the KY21 dictionary.

Compares the extended dictionary against the file it was built from, measures
coverage for one species, and scores the assignments against a reference
mapping expressed as gene pairs. Writes a metric/value TSV and the standalone
KY21 to species table.

The reference is a two-column TSV with a header: ``gene_id`` (KY21) and the
species label, for example ``KH2013`` or ``ENSEMBL``.

Exit status is non-zero when the backbone is not preserved, which is the
gate the job relies on.
"""

import argparse
import hashlib
import sys

import pandas as pd


def species_columns(species):
    return [f"{species}_primary", f"{species}_secondary"]


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_reference(path, species):
    table = pd.read_csv(path, sep="\t")
    missing = {"gene_id", species} - set(table.columns)
    if missing:
        raise ValueError(f"reference {path} lacks columns: {sorted(missing)}")
    return table[["gene_id", species]].rename(columns={species: "reference"})


def count_base_mismatches(base, extended):
    """Rows differing in any pre-existing column, after sorting by gene_id."""
    left = base.sort_values("gene_id").reset_index(drop=True)
    right = extended.sort_values("gene_id").reset_index(drop=True)
    if len(left) != len(right):
        return -1
    mismatched = 0
    for column in left.columns:
        if column not in right.columns:
            mismatched += len(left)
            continue
        a = left[column].astype(object).fillna("")
        b = right[column].astype(object).fillna("")
        mismatched += int((a != b).sum())
    return mismatched


def any_species_set(extended, columns):
    values = set()
    for column in columns:
        for entry in extended[column].dropna():
            values.update(part for part in str(entry).split("|") if part)
    return values


def reference_agreement(extended, reference, columns):
    """Recall and precision of the primary column against the reference pairs."""
    primary_column, secondary_column = columns
    by_gene = extended.set_index("gene_id")
    recall_hits = 0
    for row in reference.itertuples():
        if row.gene_id not in by_gene.index:
            continue
        primary = by_gene.at[row.gene_id, primary_column]
        secondary = by_gene.at[row.gene_id, secondary_column]
        if row.reference == primary:
            recall_hits += 1
        elif isinstance(secondary, str) and row.reference in secondary.split("|"):
            recall_hits += 1
    reference_genes = set(reference["gene_id"]) & set(by_gene.index)
    reference_gene = reference.drop_duplicates(subset=["gene_id"]).set_index("gene_id")["reference"]
    precision_hits = sum(
        1 for gene in reference_genes
        if by_gene.at[gene, primary_column] == reference_gene.get(gene))
    return recall_hits, len(reference), precision_hits, len(reference_genes)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", required=True, help="Dictionary before the run")
    parser.add_argument("--extended", required=True, help="Dictionary after the run")
    parser.add_argument("--species", required=True, help="Species label, for example KH2013")
    parser.add_argument("--reference", required=True,
                        help="Two-column TSV: gene_id and the species label")
    parser.add_argument("--projection", required=True,
                        help="Standalone KY21 to species TSV")
    parser.add_argument("--report", required=True, help="Metric/value TSV")
    args = parser.parse_args()

    base = pd.read_csv(args.base, sep="\t")
    extended = pd.read_csv(args.extended, sep="\t")
    reference = read_reference(args.reference, args.species)
    columns = species_columns(args.species)
    prefix = args.species.lower()

    missing = [column for column in columns if column not in extended.columns]
    if missing:
        sys.exit(f"extended dictionary lacks columns: {missing}")

    mismatches = count_base_mismatches(base, extended)
    primary = set(extended[columns[0]].dropna())
    any_species = any_species_set(extended, columns)
    primary_rows = int(extended[columns[0]].notna().sum())
    any_rows = int((extended[columns[0]].notna()
                    | extended[columns[1]].notna()).sum())
    recall_hits, recall_total, precision_hits, precision_total = reference_agreement(
        extended, reference, columns)

    metrics = [
        ("rows_base", len(base)),
        ("rows_extended", len(extended)),
        ("base_columns", len(base.columns)),
        ("base_column_mismatches", mismatches),
        (f"{prefix}_primary_rows", primary_rows),
        (f"{prefix}_primary_row_coverage", round(primary_rows / len(extended), 4)),
        (f"{prefix}_any_rows", any_rows),
        (f"{prefix}_any_row_coverage", round(any_rows / len(extended), 4)),
        (f"{prefix}_distinct_primary", len(primary)),
        (f"{prefix}_distinct_any", len(any_species)),
        ("reference_pairs", recall_total),
        ("reference_pairs_matched_primary_or_secondary", recall_hits),
        ("reference_pair_recall", round(recall_hits / recall_total, 4)),
        ("reference_genes_compared", precision_total),
        ("reference_primary_agrees", precision_hits),
        ("reference_primary_precision", round(precision_hits / precision_total, 4)),
        ("sha256_base", sha256(args.base)),
        ("sha256_extended", sha256(args.extended)),
        ("sha256_reference", sha256(args.reference)),
    ]

    report = pd.DataFrame(metrics, columns=["metric", "value"])
    report.to_csv(args.report, sep="\t", index=False)

    projection = extended[["gene_id"] + columns]
    projection.to_csv(args.projection, sep="\t", index=False)

    for metric, value in metrics:
        print(f"{metric}\t{value}")

    if mismatches != 0 or len(extended) != len(base):
        sys.exit("backbone not preserved; see the report")
    return 0


if __name__ == "__main__":
    sys.exit(main())
