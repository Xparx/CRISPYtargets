import pandas as pd
from Bio.Seq import Seq

# IUPAC ambiguity mapping
IUPAC_MAP = {
    'R': 'AG', 'Y': 'CT', 'S': 'GC', 'W': 'AT', 'K': 'GT', 'M': 'AC',
    'B': 'CGT', 'D': 'AGT', 'H': 'ACT', 'V': 'GAT', 'N': 'ATGC'
}

IUPAC_COMPLEMENT = {
    'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
    'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W',
    'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B',
    'H': 'D', 'D': 'H', 'N': 'N'
}

def complement_iupac(char):
    """Returns the complement of an IUPAC ambiguity code."""
    return IUPAC_COMPLEMENT.get(char.upper(), char)


# IUPAC pair mapping for biallelic SNPs
_IUPAC_PAIRS = {
    frozenset(['A', 'C']): 'M',
    frozenset(['A', 'G']): 'R',
    frozenset(['A', 'T']): 'W',
    frozenset(['C', 'G']): 'S',
    frozenset(['C', 'T']): 'Y',
    frozenset(['G', 'T']): 'K',
}

def snp_to_iupac(name):
    """Convert SNP name (chrom:pos:ref/alt) to IUPAC ambiguity code."""
    try:
        ref_alt = name.split(':')[2]
        ref, alt = ref_alt.split('/')
        return _IUPAC_PAIRS.get(frozenset([ref.upper(), alt.upper()]), 'N')
    except Exception:
        return 'N'

ENTROPY_WEIGHTS = {
    'M': 1.0, 'R': 1.0, 'W': 1.0, 'S': 1.0, 'Y': 1.0, 'K': 1.0,
    'V': -1/3, 'H': -1/3, 'D': -1/3, 'B': -1/3,
    'N': 2.0
}

def calculate_entropy(sequence):
    """Calculates sequence entropy using a fast lookup table."""
    return sum(ENTROPY_WEIGHTS.get(char, 0.0) for char in sequence)

def design_oligos(target_seq, strand, target_len=20, pam_len=3):
    """
    Generates fwd and rev oligos.

    Args:
        target_seq (str): The sequence from FlashFry (usually target + PAM).
        strand (str): The strand of the target.
        target_len (int): Expected length of the target sequence (default 20).
        pam_len (int): Expected length of the PAM (default 3).
    """
    expected_len = target_len + pam_len

    if len(target_seq) != expected_len:
        raise ValueError(
            f"Invalid target sequence length: expected {expected_len}bp "
            f"(target {target_len} + PAM {pam_len}), but got {len(target_seq)}bp. "
            f"Sequence: {target_seq}"
        )

    # 1. Remove the PAM (last pam_len characters)
    # 2. Remove the first character and prepend 'G' for T7 promoter efficiency
    # target_seq[1 : expected_len - pam_len] extracts the 19bp core
    core = "G" + target_seq[1 : expected_len - pam_len]
    fwd = "agat" + core

    # Reverse complement the fwd oligo for the rev oligo
    rev_core = str(Seq(core).reverse_complement())
    rev = "aaac" + rev_core

    return fwd, rev

def parse_flashfry_contig(contig_str):
    """Parses 'GeneID|Feature|Chr1:100-123|+' into components."""
    try:
        parts = contig_str.split('|')
        gene_id = parts[0]
        feature = parts[1]
        coords = parts[2]
        strand = parts[3]

        chrom, range_part = coords.split(':')
        start, end = range_part.split('-')
        return gene_id, feature, chrom, int(start), int(end), strand
    except:
        return None, None, None, None, None, None


