"""Subject-accession collapsing in the ortholog mapper.

Ensembl proteomes name BLAST subjects by protein (``ENSCINP``) while the
mapping is wanted at gene level (``ENSCING``). A subject map given to the
mapper must collapse isoforms before primary and secondary are chosen, so an
alternative isoform of the winning gene never lands in the secondary column.
"""

import pandas as pd

from crispytargets.mapping.mapper import create_mapping

GFF3 = (
    "##gff-version 3\n"
    "chr1\ttest\tgene\t1\t1000\t.\t+\t.\tID=G1\n"
    "chr1\ttest\tmRNA\t1\t1000\t.\t+\t.\tID=T1;Parent=G1\n"
)


def hit(sseqid, bitscore, evalue=1e-80):
    return ["T1", sseqid, 90, 100, 0, 0, 1, 100, 1, 100, evalue, bitscore]


def write_inputs(tmp_path, rows, subject_map=None):
    gff = tmp_path / "mini.gff3"
    gff.write_text(GFF3)
    blast = tmp_path / "hits.tsv"
    blast.write_text("".join("\t".join(str(v) for v in row) + "\n" for row in rows))
    map_path = None
    if subject_map:
        map_path = tmp_path / "subject_map.tsv"
        pd.DataFrame(subject_map, columns=["protein_id", "gene_id"]).to_csv(
            map_path, sep="\t", index=False)
    return gff, blast, map_path


def test_isoforms_collapse_to_one_gene(tmp_path):
    gff, blast, subject_map = write_inputs(tmp_path, [
        hit("ENSCINP00000000001", 200),
        hit("ENSCINP00000000002", 180),
        hit("ENSCINP00000000003", 150),
    ], subject_map=[
        ("ENSCINP00000000001", "ENSCING00000000001"),
        ("ENSCINP00000000002", "ENSCING00000000001"),
        ("ENSCINP00000000003", "ENSCING00000000002"),
    ])
    result = create_mapping(gff, blast, "ENSEMBL", tmp_path / "out.tsv",
                            subject_map_path=subject_map)
    row = result.set_index("gene_id").loc["G1"]
    assert row["ENSEMBL_primary"] == "ENSCING00000000001"
    assert row["ENSEMBL_secondary"] == "ENSCING00000000002"


def test_without_a_map_the_accessions_pass_through(tmp_path):
    gff, blast, _ = write_inputs(tmp_path, [hit("ENSCINP00000000001", 200)])
    result = create_mapping(gff, blast, "ENSEMBL", tmp_path / "out.tsv")
    row = result.set_index("gene_id").loc["G1"]
    assert row["ENSEMBL_primary"] == "ENSCINP00000000001"
