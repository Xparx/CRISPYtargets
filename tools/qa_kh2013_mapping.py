"""Quality report for the KH2013 columns added to the KY21 dictionary.

Compares the extended dictionary against the file it was built from, measures
KH2013 coverage, and scores the assignments against GHOST's own KH-to-KY21
mapping. Writes a metric/value TSV and the standalone KY21 to KH2013 table.

Exit status is non-zero when the backbone is not preserved, which is the
gate the job relies on.
"""

import argparse
import gzip
import hashlib
import sys

import pandas as pd

KH_COLUMNS = ["KH2013_primary", "KH2013_secondary"]


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_ghost_pairs(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    rows = []
    with opener(path, "rt") as handle:
        for line in handle:
            parts = line.split()
            if len(parts) >= 2 and parts[0].startswith("KH."):
                rows.append((parts[0], parts[1]))
    return pd.DataFrame(rows, columns=["KH2013", "gene_id"])


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


def primary_set(extended):
    return set(extended["KH2013_primary"].dropna())


def any_kh_set(extended):
    values = set()
    for column in KH_COLUMNS:
        for entry in extended[column].dropna():
            values.update(part for part in str(entry).split("|") if part)
    return values


def ghost_agreement(extended, ghost):
    """Recall and precision of the primary column against GHOST's pairs."""
    by_gene = extended.set_index("gene_id")
    recall_hits = 0
    for row in ghost.itertuples():
        if row.gene_id not in by_gene.index:
            continue
        primary = by_gene.at[row.gene_id, "KH2013_primary"]
        secondary = by_gene.at[row.gene_id, "KH2013_secondary"]
        if row.KH2013 == primary:
            recall_hits += 1
        elif isinstance(secondary, str) and row.KH2013 in secondary.split("|"):
            recall_hits += 1
    ghost_genes = set(ghost["gene_id"]) & set(by_gene.index)
    ghost_kh = ghost.drop_duplicates(subset=["gene_id"]).set_index("gene_id")["KH2013"]
    precision_hits = sum(
        1 for gene in ghost_genes
        if by_gene.at[gene, "KH2013_primary"] == ghost_kh.get(gene))
    return recall_hits, len(ghost), precision_hits, len(ghost_genes)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", required=True, help="Dictionary before the run")
    parser.add_argument("--extended", required=True, help="Dictionary after the run")
    parser.add_argument("--ghost", required=True, help="GHOST KHtoKY21 file")
    parser.add_argument("--projection", required=True, help="Standalone KY21 to KH2013 TSV")
    parser.add_argument("--report", required=True, help="Metric/value TSV")
    args = parser.parse_args()

    base = pd.read_csv(args.base, sep="\t")
    extended = pd.read_csv(args.extended, sep="\t")
    ghost = read_ghost_pairs(args.ghost)

    missing = [column for column in KH_COLUMNS if column not in extended.columns]
    if missing:
        sys.exit(f"extended dictionary lacks columns: {missing}")

    mismatches = count_base_mismatches(base, extended)
    primary = primary_set(extended)
    any_kh = any_kh_set(extended)
    primary_rows = int(extended["KH2013_primary"].notna().sum())
    any_rows = int((extended["KH2013_primary"].notna()
                    | extended["KH2013_secondary"].notna()).sum())
    recall_hits, recall_total, precision_hits, precision_total = ghost_agreement(
        extended, ghost)

    metrics = [
        ("rows_base", len(base)),
        ("rows_extended", len(extended)),
        ("base_columns", len(base.columns)),
        ("base_column_mismatches", mismatches),
        ("kh2013_primary_rows", primary_rows),
        ("kh2013_primary_row_coverage", round(primary_rows / len(extended), 4)),
        ("kh2013_any_rows", any_rows),
        ("kh2013_any_row_coverage", round(any_rows / len(extended), 4)),
        ("kh2013_distinct_primary", len(primary)),
        ("kh2013_distinct_any", len(any_kh)),
        ("ghost_pairs", recall_total),
        ("ghost_pairs_matched_primary_or_secondary", recall_hits),
        ("ghost_pair_recall", round(recall_hits / recall_total, 4)),
        ("ghost_genes_compared", precision_total),
        ("ghost_primary_agrees", precision_hits),
        ("ghost_primary_precision", round(precision_hits / precision_total, 4)),
        ("sha256_base", sha256(args.base)),
        ("sha256_extended", sha256(args.extended)),
        ("sha256_ghost", sha256(args.ghost)),
    ]

    report = pd.DataFrame(metrics, columns=["metric", "value"])
    report.to_csv(args.report, sep="\t", index=False)

    projection = extended[["gene_id"] + KH_COLUMNS]
    projection.to_csv(args.projection, sep="\t", index=False)

    for metric, value in metrics:
        print(f"{metric}\t{value}")

    if mismatches != 0 or len(extended) != len(base):
        sys.exit("backbone not preserved; see the report")
    return 0


if __name__ == "__main__":
    sys.exit(main())
