"""Secondary-hit hygiene in the ortholog mapper.

Distinct BLAST subject IDs can strip to the same accession (for example
``sp|P1|A_HUMAN`` and ``tr|P1|A_HUMAN``). The secondary column must not
repeat the primary entry in that case.
"""

import pandas as pd

from crispytargets.mapping.mapper import create_mapping

GFF3 = (
    "##gff-version 3\n"
    "chr1\ttest\tgene\t1\t1000\t.\t+\t.\tID=G1\n"
    "chr1\ttest\tmRNA\t1\t1000\t.\t+\t.\tID=T1;Parent=G1\n"
)

FIELDS = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
          "qstart", "qend", "sstart", "send", "evalue", "bitscore"]


def write_inputs(tmp_path, rows):
    gff = tmp_path / "mini.gff3"
    gff.write_text(GFF3)
    blast = tmp_path / "hits.tsv"
    blast.write_text("".join("\t".join(str(v) for v in row) + "\n" for row in rows))
    return gff, blast


def hit(sseqid, bitscore, evalue=1e-80):
    return ["T1", sseqid, 90, 100, 0, 0, 1, 100, 1, 100, evalue, bitscore]


def test_secondary_does_not_repeat_the_primary_accession(tmp_path):
    gff, blast = write_inputs(tmp_path, [
        hit("sp|P1|A_HUMAN", 200),
        hit("tr|P1|A_HUMAN", 180),
        hit("sp|P2|B_HUMAN", 150),
    ])
    result = create_mapping(gff, blast, "Human", tmp_path / "out.tsv")
    row = result.set_index("gene_id").loc["G1"]
    assert row["Human_primary"] == "A_HUMAN"
    assert row["Human_secondary"] == "B_HUMAN"


def test_single_subject_leaves_secondary_empty(tmp_path):
    gff, blast = write_inputs(tmp_path, [hit("sp|P1|A_HUMAN", 200)])
    result = create_mapping(gff, blast, "Human", tmp_path / "out.tsv")
    row = result.set_index("gene_id").loc["G1"]
    assert row["Human_primary"] == "A_HUMAN"
    assert row["Human_secondary"] == ""
