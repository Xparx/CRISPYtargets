#!/usr/bin/env python3
"""
lift_snps.py - Map SNP positions between genome assemblies (IMPROVED)

Performs whole-genome alignment to establish coordinate correspondence,
then lifts VCF positions from the old assembly to the new assembly.

Usage:
    python3 lift_snps.py \
        --old-genome ../data/ci3_main.fa \
        --new-genome ../data/HT.Ref.fasta \
        --vcf ../data/Ghost_SNPs_main.vcf.gz \
        --output Ghost_SNPs_lifted.vcf \
        --threads 8
"""

import argparse
import bisect
import gzip
import logging
import os
import re
import subprocess
import sys
from collections import defaultdict

try:
    import pysam
except ImportError:
    sys.exit("Error: pysam required. Install with: pip install pysam")

from ..common.exec import run


# ── genome alignment ─────────────────────────────────────────────────────────

def convert_2bit_to_fa(twobit, fa):
    if os.path.exists(fa):
        logging.info("FASTA already exists: %s", fa)
        return
    run(["twoBitToFa", twobit, fa], f"Converting {twobit} → {fa}")
    pysam.faidx(fa)


def index_genome(fasta):
    """Index a FASTA file with pysam if not already indexed."""
    fai = fasta + ".fai"
    if os.path.exists(fai):
        logging.info("Genome already indexed: %s", fai)
        return
    pysam.faidx(fasta)


def check_vcf_file(vcf_path):
    """Check VCF file integrity and ensure it's indexed."""
    logging.info(f"Checking VCF file: {vcf_path}")
    
    # Check if file exists
    if not os.path.exists(vcf_path):
        logging.error(f"VCF file not found: {vcf_path}")
        sys.exit(1)
    
    # Check if indexed
    tbi_path = vcf_path + ".tbi"
    if not os.path.exists(tbi_path):
        logging.info("VCF index not found, creating...")
        if vcf_path.endswith('.gz'):
            run(["tabix", "-f", "-p", "vcf", vcf_path], "Indexing VCF")
        else:
            logging.warning("VCF is not gzipped. Consider compressing with bgzip.")
    
    # Try to open with pysam
    try:
        vcf = pysam.VariantFile(vcf_path)
        header = vcf.header
        n_samples = len(header.samples)
        n_records = sum(1 for _ in vcf)
        vcf.close()
        logging.info(f"VCF valid: {n_samples} samples, {n_records} records")
        return True
    except Exception as e:
        logging.error(f"VCF file appears corrupted: {e}")
        return False


# ── genome alignment ─────────────────────────────────────────────────────────

def align_genomes(new_fa, old_fa, paf_path, threads, preset="asm5"):
    """Align OLD genome (query) to NEW genome (reference) with minimap2."""
    logging.info(f"Aligning genomes (minimap2 -x {preset}) …")

    # -c is required so minimap2 emits the cg:Z: CIGAR. Without it a record
    # spanning an indel looks like one perfectly collinear block, and every
    # position past the indel would lift to the wrong coordinate.
    result = subprocess.run(
        ["minimap2", "-x", preset, "-c",
         "-t", str(threads), new_fa, old_fa],
        capture_output=True, text=True
    )

    if result.returncode != 0:
        logging.error("minimap2 failed:\n%s", result.stderr)
        sys.exit(1)

    with open(paf_path, "w") as f:
        f.write(result.stdout)

    n_records = sum(1 for _ in open(paf_path))
    logging.info("Alignment written to %s (%d records)", paf_path, n_records)


# ── position mapper ──────────────────────────────────────────────────────────

def blocks_from_cigar(qstart, tname, tstart, cigar):
    """Split a PAF CIGAR into collinear blocks.

    Yields (qstart, qend, tname, tstart, tend) for each run of aligned bases
    (M/=/X). Indels advance one offset without emitting a block, so a block
    boundary appears at every indel. Coordinates are 0-based half-open.
    """
    q = qstart
    t = tstart
    for length, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
        length = int(length)
        if op in "M=X":
            yield (q, q + length, tname, t, t + length)
            q += length
            t += length
        elif op == "I":        # query-only bases: deletion in the new assembly
            q += length
        elif op in "DN":       # target-only bases: insertion in the new assembly
            t += length
        # S/H (clipping) and P (padding) consume nothing here: PAF qstart/qend
        # already delimit the aligned region, and minimap2 omits soft clips
        # from the PAF CIGAR for that reason.


class AlignmentBlockMapper:
    """
    Stores aligned blocks from a whole-genome PAF and provides O(log n)
    position look-ups.

    Block = (old_start, old_end, new_chrom, new_start, new_end)
    All coordinates are 0-based, half-open.
    """

    def __init__(self):
        self.blocks = defaultdict(list)       # old_chrom → [block, …]
        self.starts = {}                      # old_chrom → [block_start, …] for fast binary search
        self.records_without_cigar = 0        # records that could not be split

    # ── loading ──────────────────────────────────────────────────────────

    def load_from_paf(self, paf_path):
        """Parse PAF file and extract aligned blocks.

        Each record is split into collinear blocks at its indels using the
        cg:Z: CIGAR, so positions downstream of an indel lift correctly.
        """
        with open(paf_path) as f:
            for line in f:
                fields = line.strip().split("\t")
                if len(fields) < 12:
                    continue

                qname  = fields[0]              # old-genome chromosome
                qstart = int(fields[2])         # 0-based, inclusive
                qend   = int(fields[3])         # exclusive
                tname  = fields[5]              # new-genome chromosome
                tstart = int(fields[7])         # 0-based, inclusive
                tend   = int(fields[8])         # exclusive

                cigar = None
                for tag in fields[12:]:
                    if tag.startswith("cg:Z:"):
                        cigar = tag[len("cg:Z:"):]
                        break

                if cigar:
                    self.blocks[qname].extend(
                        blocks_from_cigar(qstart, tname, tstart, cigar))
                else:
                    # No CIGAR to split on. Falling back to the whole record
                    # treats it as collinear, which is only correct when the
                    # alignment has no indels.
                    self.records_without_cigar += 1
                    self.blocks[qname].append(
                        (qstart, qend, tname, tstart, tend))

        if self.records_without_cigar:
            logging.warning(
                "%d PAF records had no cg:Z: CIGAR; treated as collinear blocks. "
                "Positions downstream of indels in those records will be wrong. "
                "Re-run the alignment with minimap2 -c.",
                self.records_without_cigar,
            )

        # sort + precompute
        for c in self.blocks:
            self.blocks[c].sort()
            self.starts[c] = [b[0] for b in self.blocks[c]]

        n_blocks = sum(len(v) for v in self.blocks.values())
        total_mapped = sum(
            sum(b[1] - b[0] for b in blks)
            for blks in self.blocks.values()
        )
        logging.info("Extracted %d blocks across %d old-genome chromosomes",
                      n_blocks, len(self.blocks))
        logging.info("Total mapped bases: %d bp", total_mapped)
        if n_blocks == 0:
            logging.error("No alignment blocks found!")
            sys.exit(1)

    # ── look-up ──────────────────────────────────────────────────────────

    def map_position(self, old_chrom, old_pos_1based):
        """
        Return (new_chrom, new_pos_1based) or (None, None) on failure.
        """
        q = old_pos_1based - 1                  # work in 0-based

        blks = self.blocks.get(old_chrom)
        if not blks:
            return None, None

        # binary search: find rightmost block whose start ≤ q
        starts = self.starts[old_chrom]
        idx = bisect.bisect_right(starts, q) - 1

        if idx >= 0:
            b_old_s, b_old_e, new_c, b_new_s, b_new_e = blks[idx]
            if b_old_s <= q < b_old_e:
                return new_c, b_new_s + (q - b_old_s) + 1   # → 1-based

        return None, None

    # ── reporting ────────────────────────────────────────────────────────

    def summary(self):
        n_blocks = sum(len(v) for v in self.blocks.values())
        total_bp = sum(
            sum(b[1] - b[0] for b in blks)
            for blks in self.blocks.values()
        )
        
        # Coverage per chromosome
        coverage = {}
        for chrom, blks in self.blocks.items():
            total = sum(b[1] - b[0] for b in blks)
            coverage[chrom] = total
            
        return {"chromosomes": len(self.blocks),
                "blocks": n_blocks,
                "mapped_bp": total_bp,
                "coverage": coverage}

    def new_chromosomes(self):
        """Return the set of new-genome chromosomes seen in alignments."""
        return {b[2] for blks in self.blocks.values() for b in blks}


# ── VCF lifting ──────────────────────────────────────────────────────────────

def lift_vcf(vcf_path, mapper, output_path, output_bed=False):
    # Check if input VCF is gzipped
    is_gzipped = vcf_path.endswith('.gz')
    
    # Open input VCF
    old = pysam.VariantFile(vcf_path)

    # build a new header
    hdr = old.header.copy()
    hdr.add_line("##source=lift_snps.py")
    hdr.add_line("##reference=new_genome_assembly")

    # add contig lines for new-genome chromosomes (natural numeric sort)
    for nc in sorted(mapper.new_chromosomes(),
                     key=lambda x: int(''.join(c for c in x if c.isdigit()) or '9999')):
        hdr.add_line(f"##contig=<ID={nc}>")

    # Determine output format
    if output_path.endswith('.gz'):
        out = pysam.VariantFile(output_path, "wz", header=hdr)
        output_is_gzipped = True
    else:
        out = pysam.VariantFile(output_path, "w", header=hdr)
        output_is_gzipped = False
    
    # Initialize BED output if requested
    bed_out = None
    if output_bed:
        bed_path = output_path.replace('.vcf', '.bed')
        if output_path.endswith('.gz'):
            bed_path = bed_path.replace('.gz', '')
        bed_out = open(bed_path, 'w')
        bed_out.write("#chrom\tstart\tend\tname\tscore\tstrand\n")

    # Build chromosome order map for sorted output
    chrom_order = {c: i for i, c in enumerate(
        sorted(mapper.new_chromosomes(),
               key=lambda x: int(''.join(c for c in x if c.isdigit()) or '9999'))
    )}

    # Collect records for sorted output
    lifted_records = []
    bed_records = []

    lifted = failed = total = 0
    failed_chroms = defaultdict(int)

    for rec in old:
        total += 1
        nc, np_ = mapper.map_position(rec.chrom, rec.pos)

        if nc is not None:
            # Write VCF record
            r = hdr.new_record()
            r.chrom = nc
            r.pos  = np_
            r.id   = rec.id
            r.ref  = rec.ref
            r.alts = rec.alts
            if rec.qual is not None:
                r.qual = rec.qual
            
            # Skip empty filter strings
            for f in rec.filter:
                if f:
                    r.filter.add(f)
            
            # Copy INFO fields, updating END if present
            for k, v in rec.info.items():
                if k == 'END':
                    r.info[k] = np_  # Update END to new position
                else:
                    r.info[k] = v
            
            # Collect VCF record for sorted output
            lifted_records.append((chrom_order.get(nc, 9999), np_, r))
            
            # Collect BED record if requested
            if output_bed:
                bed_start = np_ - 1
                bed_end = np_
                
                if rec.id:
                    snp_name = rec.id
                else:
                    alt_str = ','.join(str(a) for a in rec.alts) if rec.alts else '.'
                    snp_name = f"{rec.chrom}:{rec.pos}:{rec.ref}/{alt_str}"
                
                # Use QUAL as score; default to 1000 if missing
                score = min(int(rec.qual), 1000) if rec.qual else 1000

                bed_records.append((chrom_order.get(nc, 9999), bed_start, bed_end,
                                    snp_name, score, nc))
            
            lifted += 1
        else:
            failed += 1
            failed_chroms[rec.chrom] += 1

    old.close()

    # Sort and write VCF records by chromosome order then position
    lifted_records.sort(key=lambda x: (x[0], x[1]))
    logging.info("Writing %d sorted VCF records …", len(lifted_records))
    for _, _, r in lifted_records:
        out.write(r)
    out.close()

    # Sort and write BED records
    if bed_out:
        bed_records.sort(key=lambda x: (x[0], x[1]))
        for _, start, end, name, score, chrom in bed_records:
            bed_out.write(f"{chrom}\t{start}\t{end}\t{name}\t{score}\t.\n")
        bed_out.close()
        logging.info("BED output written to %s", bed_path)

    logging.info("Lifted %d / %d variants  (%d failed)", lifted, total, failed)
    if failed_chroms:
        logging.warning("Failures by old-genome chromosome:")
        for c in sorted(failed_chroms, key=failed_chroms.get, reverse=True):
            logging.warning("  %s : %d", c, failed_chroms[c])
    
    # Index output if gzipped
    if output_is_gzipped:
        tbi = output_path + ".tbi"
        if os.path.exists(tbi):
            os.remove(tbi)
        run(["tabix", "-p", "vcf", output_path], "Indexing output VCF")
    
    return lifted, failed, total


# ── main ─────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--old-genome", required=True,
                    help="Old assembly  (FASTA or .2bit)")
    ap.add_argument("--new-genome", required=True,
                    help="New assembly FASTA (will be indexed if needed)")
    ap.add_argument("--vcf",        required=True,
                    help="Input VCF(.gz) with old-genome coordinates")
    ap.add_argument("--output",     required=True,
                    help="Output directory for lifted VCF and BED files")
    ap.add_argument("--bed",        action="store_true",
                    help="Also output a BED file with lifted positions")
    ap.add_argument("--threads",    type=int, default=8)
    ap.add_argument("--preset",     default="asm5",
                    choices=["asm5", "asm10", "asm20", "map-ont", "map-hifi"],
                    help="Minimap2 preset for alignment (default: asm5)")
    ap.add_argument("--keep-intermediates", action="store_true",
                    help="Do not delete the intermediate PAF")
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
    )

    # Resolve output directory and derive filenames
    outdir = args.output
    os.makedirs(outdir, exist_ok=True)
    logging.info("Output directory: %s", outdir)

    # Derive output filenames from input VCF
    vcf_base = os.path.splitext(os.path.basename(args.vcf))[0]
    if vcf_base.endswith('.gz'):
        vcf_base = vcf_base[:-3]
    output_vcf = os.path.join(outdir, f"{vcf_base}_lifted.vcf.gz")
    output_bed = os.path.join(outdir, f"{vcf_base}_lifted.bed")

    # ── step 0: check VCF file ───────────────────────────────────────────
    check_vcf_file(args.vcf)

    # ── step 1: convert .2bit → FASTA if needed ──────────────────────────
    old_fa = args.old_genome
    if args.old_genome.endswith(".2bit"):
        old_fa = os.path.join(os.path.dirname(args.output) or ".",
                              os.path.basename(args.old_genome).replace(".2bit", ".fa"))
        convert_2bit_to_fa(args.old_genome, old_fa)

    # ── step 2: ensure genomes are indexed ───────────────────────────────
    index_genome(args.new_genome)
    index_genome(old_fa)

    # ── step 3: whole-genome alignment ───────────────────────────────────
    paf = os.path.join(outdir, "alignment.paf")
    align_genomes(args.new_genome, old_fa, paf, args.threads, args.preset)

    # ── step 4: build position mapper ────────────────────────────────────
    mapper = AlignmentBlockMapper()
    mapper.load_from_paf(paf)
    s = mapper.summary()
    logging.info("Mapper summary: %d chromosomes, %d blocks, ~%d bp mapped",
                 s["chromosomes"], s["blocks"], s["mapped_bp"])
    
    # Log coverage per chromosome
    logging.info("Coverage per old-genome chromosome:")
    for chrom in sorted(s["coverage"].keys()):
        logging.info("  %s: %d bp", chrom, s["coverage"][chrom])

    # ── step 4.5: diagnostic ─────────────────────────────────────────────
    logging.info("Old-genome chromosomes in mapper: %s",
                 sorted(mapper.blocks.keys()))

    # ── step 5: lift VCF ─────────────────────────────────────────────────
    lift_vcf(args.vcf, mapper, output_vcf, output_bed=args.bed)

    # ── cleanup ──────────────────────────────────────────────────────────
    if not args.keep_intermediates:
        if os.path.exists(paf):
            os.remove(paf)

    logging.info("Done.")


if __name__ == "__main__":
    main()
