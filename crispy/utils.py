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

def calculate_entropy(sequence):
    """Calculates sequence entropy based on IUPAC ambiguity codes."""
    entropy = 0.0
    for char in sequence:
        if char in ['M', 'R', 'W', 'S', 'Y', 'K']:
            entropy += 1.0
        elif char in ['V', 'H', 'D', 'B']:
            entropy += -1.0 * (1/3) # Simplified approximation of R logic
        elif char == 'N':
            entropy += 2.0
    return entropy

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


def convert_output_to_igv_bed(input_file):
    from pathlib import Path

    input_file = Path(input_file)
    # Load the final results
    df = pd.read_csv(input_file, sep='\t')

    # BED format requires: Chromosome, Start, End, Name, Score, Strand
    # We will use 'gene_id' as the name and 'Doench2014OnTarget' as the score
    # bed_df = df[df['over_cutoff']][['Chromosome', 'Start', 'End', 'feature', 'Doench2014OnTarget', 'Strand']]
    df['start'] = df['Start'] + df['start']
    df['stop'] = df['Start'] + df['stop']

    bed_df = df[df['over_cutoff']][['Chromosome', 'start', 'stop', 'feature', 'Doench2014OnTarget', 'Strand']]
    bed_df.columns = ['chrom', 'start', 'end', 'name', 'score', 'strand']
    bed_df = bed_df.drop_duplicates() # (subset=['chrom', 'start', 'end', 'name'])

    out_file = input_file.parent / (input_file.stem + '.bed')
    # Save as BED
    bed_df.to_csv(out_file, sep='\t', index=False, header=False)
    print(f"Converted {input_file.name} to {out_file.name}")
