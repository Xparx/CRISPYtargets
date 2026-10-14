import argparse
import os
import sys
import time
import numpy as np
import pandas as pd
import pyranges as pr
import logging
from pyfaidx import Fasta
from .discovery import discover_and_score, extract_sequences
from ..common.config import load_config
from ..common.exec import run_command
from ..common.logging import setup_logging
from .utils import calculate_entropy, design_oligos, parse_flashfry_contig
from .utils import complement_iupac, snp_to_iupac
from .doench2016 import score_guides
import gc


def _strip_namespace(fid):
    """Strip a GFF3 namespace prefix, e.g. 'KY21:KY21.Chr1.709' -> 'KY21.Chr1.709'.

    Lets bare gene IDs in a user-supplied list match namespaced annotations.
    A no-op for IDs without a namespace prefix.
    """
    return fid.split(':', 1)[-1] if fid and ':' in fid else fid


def _norm_ids(series):
    """Normalize a GFF ID-like column for matching.

    Strips any namespace prefix and maps missing/placeholder values (NaN and
    the 'nan'/'None' strings produced by astype(str)) to None, so they never
    match and never leak into the gene expansion.
    """
    return series.astype(str).map(_strip_namespace).replace(
        {"": None, "nan": None, "None": None})


def _parse_alt(name):
    """Parse the concrete ALT allele from a SNP name like 'chr:pos:REF/ALT'.

    Returns the ALT base (uppercased) or None if the name has no REF/ALT
    suffix. Used to build variant sequences for rescoring.
    """
    try:
        ref_alt = str(name).rsplit(':', 1)[-1]
        if '/' not in ref_alt:
            return None
        alt = ref_alt.split('/')[-1].strip().upper()
        return alt if alt in ('A', 'C', 'G', 'T') else None
    except Exception:
        return None


def _score_variant(guide, pam, flank5, flank3, alt_by_pos, rs3_fn=None):
    """Score a variant guide with ALT alleles substituted.

    alt_by_pos maps a position in the 23-mer (guide 0..19, PAM 20..22) to the
    concrete ALT base. Uses the real flanks and real PAM so the variant score is
    directly comparable to the reference score. Returns (rs2, rs3); rs3 is NaN
    when rs3_fn is None.
    """
    g = list(guide)
    p = list(pam)
    for rp, alt in alt_by_pos.items():
        if rp < len(g):
            g[rp] = alt
        else:
            p[rp - len(g)] = alt
    g = "".join(g)
    p = "".join(p)

    rs2 = float(score_guides([g], flanks=[(flank5, flank3)], pams=[p])[0])
    rs3 = np.nan
    if rs3_fn is not None:
        left = flank5[-4:] if flank5 and len(flank5) >= 4 else "N" * 4
        right = flank3[:3] if flank3 and len(flank3) >= 3 else "N" * 3
        rs3 = float(rs3_fn([(left + g + p + right).upper()])[0])
    return rs2, rs3


class CrispyPipeline:
    def __init__(self, args):
        self.args = args
        self.outdir = args.outdir
        self.ff_jar = args.flashfry_jar
        self.db = os.path.join(args.outdir, "htky_db")
        # Built lazily by extract_sequences; used by post_process for gene_id.
        self.gene_map = {}

    def prepare_targets(self):
        logging.info("Preparing target regions...")
        # Load GFF
        gr = pr.read_gff3(self.args.gff)
        required_gff = {'ID', 'Parent', 'Feature'}
        missing = required_gff - set(gr.df.columns)
        if missing:
            raise ValueError(f"Malformed GFF3 in {self.args.gff}. Missing required columns: {missing}")

        # 1. Get Exons
        exons = gr[(gr.df.Feature == 'exon') | (gr.df.Feature == 'CDS')]
        exons_df = exons.df.drop_duplicates(subset=['Chromosome', 'Start', 'End', 'Strand'])
        exons = pr.PyRanges(exons_df)

        # 2. Get Promoters
        tss_features = self.args.tss_features.split(',')
        upstream = self.args.promoter_upstream
        downstream = self.args.promoter_downstream
        logging.info(f"Identifying promoters using features: {tss_features} (Window: -{upstream} to +{downstream})")
        tpts = gr[gr.df.Feature.isin(tss_features)]
        # tpts = gr[gr.feature == 'transcript']

        # A GFF3 slice with no matching features yields an empty PyRanges,
        # whose .df has no columns at all; indexing Strand/Start/End on it
        # raises KeyError. Guard every access to a strand-filtered slice.
        promoter_dfs = []
        if len(tpts) > 0:
            # Positive strand: TSS is start.
            pos_tpts = tpts[tpts.df['Strand'] == '+']
            if len(pos_tpts) > 0:
                pos_df = pos_tpts.df.copy()
                pos_df['Feature'] = 'promoter' # Explicitly label as promoter
                tss_pos = pos_df['Start']
                pos_df['Start'] = tss_pos - upstream
                pos_df['End'] = tss_pos + downstream
                promoter_dfs.append(pos_df)

            # Negative strand: TSS is end.
            neg_tpts = tpts[tpts.df['Strand'] == '-']
            if len(neg_tpts) > 0:
                neg_df = neg_tpts.df.copy()
                neg_df['Feature'] = 'promoter' # Explicitly label as promoter
                tss_neg = neg_df['End']
                neg_df['Start'] = tss_neg - downstream
                neg_df['End'] = tss_neg + upstream
                promoter_dfs.append(neg_df)

        if promoter_dfs:
            promoters_df = pd.concat(promoter_dfs)
        else:
            logging.warning(
                "No TSS features matched %s on either strand; promoter targets skipped. "
                "Check --tss-features against your GFF3's feature labels.",
                tss_features,
            )
            promoters_df = pd.DataFrame(columns=exons.df.columns)
        promoters = pr.PyRanges(promoters_df)

        # Combine exons and promoters by concatenating their DataFrames
        targets_df = pd.concat([exons.df, promoters.df])

        # Add flanking sequence so FlashFry has context for Doench scoring
        # Doench needs ~17bp of context on each side of the guide
        flank = max(self.args.promoter_upstream, self.args.promoter_downstream, 100)
        targets_df['Start'] = (targets_df['Start'] - flank).clip(lower=0)
        targets_df['End'] = targets_df['End'] + flank
        targets = pr.PyRanges(targets_df)

        # Optional: Filter by provided gene list
        if self.args.genes:
            logging.info(f"Filtering by gene list: {self.args.genes}")
            # Load seed IDs, normalized so bare IDs match namespaced GFF IDs.
            seed_ids = {_strip_namespace(s) for s in
                        pd.read_csv(self.args.genes)['Gene.ID'].astype(str)}

            # We need the full GFF to resolve the hierarchy
            gff_df = gr.df

            # Normalize GFF match keys (strip any 'NAMESPACE:' prefix) so a bare
            # user ID like KY21.Chr1.709 matches KY21:KY21.Chr1.709. No-op for
            # non-namespaced annotations.
            id_norm = _norm_ids(gff_df['ID']) if 'ID' in gff_df.columns else None
            parent_norm = _norm_ids(gff_df['Parent']) if 'Parent' in gff_df.columns else None

            # 1. Initialize expanded_ids with seed_ids.
            # This ensures that if a seed is already a Parent or ID, the loop finds children.
            expanded_ids = set(seed_ids)

            # Also map any features that match seed IDs in Name/name columns to their IDs
            for col in ['Name', 'name']:
                if col in gff_df.columns and id_norm is not None:
                    name_norm = _norm_ids(gff_df[col])
                    expanded_ids.update(id_norm[name_norm.isin(seed_ids)].dropna().tolist())

            # 2. Recursive Expansion: Find all children (Gene -> mRNA -> Exon/CDS)
            while True:
                prev_count = len(expanded_ids)
                if id_norm is not None and parent_norm is not None:
                    # Find features whose Parent is in our current expanded set
                    expanded_ids.update(id_norm[parent_norm.isin(expanded_ids)].dropna().tolist())

                if len(expanded_ids) == prev_count:
                    break

            logging.info(f"Expanded {len(seed_ids)} seed IDs to {len(expanded_ids)} total related features.")

            # 3. Filter targets: Keep if ID or Parent is in the expanded set
            tgt_id = _norm_ids(targets.df['ID']) if 'ID' in targets.df.columns else pd.Series([None] * len(targets.df))
            tgt_parent = _norm_ids(targets.df['Parent']) if 'Parent' in targets.df.columns else pd.Series([None] * len(targets.df))

            combined_mask = tgt_id.isin(expanded_ids) | tgt_parent.isin(expanded_ids)
            targets = targets[combined_mask]

        else:
            logging.info("No gene list provided. Processing all genes in GFF...")

        # Optional Filter BED
        if self.args.filter_bed and os.path.exists(self.args.filter_bed):
            logging.info(f"Applying optional filter: {self.args.filter_bed}")
            filter_gr = pr.read_bed(self.args.filter_bed)
            targets = targets.intersect(filter_gr)

        # --- Build Gene Lineage Map ---
        # This resolves the hierarchy: CDS -> mRNA -> Gene
        self.gene_map = {}
        gff_df = gr.df
        
        # Create quick lookup maps for parents and feature types
        # We use astype(str) to avoid issues with NaN/float types
        parent_lookup = {}
        for fid, parent in zip(gff_df['ID'].astype(str), gff_df['Parent'].astype(str)):
            if parent != 'nan':
                parent_lookup[fid] = parent.split(',')
        feature_lookup = dict(zip(gff_df['ID'].astype(str), gff_df['Feature'].astype(str)))

        def resolve_to_gene(fid):
            """Finds the closest ancestor labeled as 'gene' using BFS."""
            queue = [fid]
            visited = {fid}
            while queue:
                current = queue.pop(0)
                if feature_lookup.get(current) == 'gene':
                    return current
                parents = parent_lookup.get(current, [])
                for p in parents:
                    if p not in visited:
                        visited.add(p)
                        queue.append(p)
            return fid

        # Map every ID in the GFF to its root gene
        for fid in gff_df['ID'].dropna().astype(str):
            self.gene_map[fid] = resolve_to_gene(fid)
        # -----------------------------

        targets.to_bed(os.path.join(self.outdir, "targets.bed"))
        return targets

    def extract_sequences(self, targets_gr):
        return extract_sequences(targets_gr, self.outdir, self.args.ref,
                                 self.gene_map)

    def run_flashfry(self, targets_gr):
        return discover_and_score(self.ff_jar, self.db, self.args.ref,
                                  self.outdir, targets_gr, self.args,
                                  self.gene_map)

    def post_process(self, scored_file):
        logging.info("Post-processing and scoring...")
        df = pd.read_csv(scored_file, sep='\t')
        if 'contig' not in df.columns:
            raise ValueError(f"FlashFry output {scored_file} is missing the 'contig' column. Ensure the --positionOutput flag was used during discovery.")

        # Parse coordinates for all targets (needed for sorting and accessibility)
        parsed_coords = df['contig'].apply(parse_flashfry_contig)
        # Unpack 6 values: gene_id, feature, Chromosome, Start, End, Strand
        df[['gene_id', 'feature', 'Chromosome', 'Start', 'End', 'Strand']] = pd.DataFrame(parsed_coords.tolist(), index=df.index)
        target_gr = pr.PyRanges(df[['Chromosome', 'Start', 'End', 'Strand']])

        # 1. Accessibility (Optional)
        df['accessibility'] = 0
        if self.args.accessome and os.path.exists(self.args.accessome):
            logging.info(f"Calculating accessibility using {self.args.accessome}...")
            acc_gr = pr.read_bed(self.args.accessome)
            if acc_gr.df.shape[1] < 3:
                raise ValueError(f"Accessome BED {self.args.accessome} is malformed. Expected at least 3 columns (Chrom, Start, End).")

            # Preserve original index to map back to df
            target_gr_indexed = target_gr.assign(orig_idx=target_gr.df.index)
            joined = target_gr_indexed.join(acc_gr)

            # Calculate total overlap length (bp) per target
            # joined.df['End'] and joined.df['Start'] are the intersection coordinates
            joined.df['overlap_bp'] = joined.df['End'] - joined.df['Start']
            overlap_sum = joined.df.groupby('orig_idx')['overlap_bp'].sum()

            df['accessibility'] = overlap_sum.reindex(df.index, fill_value=0)
            del joined # Explicitly remove the large joined dataframe
            gc.collect() # Force memory release
        else:
            logging.info("No accessome file provided. Skipping accessibility scoring.")

        # Find sequence column early (needed for SNP annotation too)
        expected_len = self.args.target_len + self.args.pam_len
        seq_col = None
        for col in df.columns:
            first_val = df[col].dropna().iloc[0] if not df[col].dropna().empty else None
            if isinstance(first_val, str) and len(first_val) == expected_len:
                seq_col = col
                break
        if seq_col is None:
            raise KeyError(f"Could not find a sequence column with length {expected_len}. Columns: {list(df.columns)}")
        logging.info(f"Using column '{seq_col}' as the target sequence source.")

        # 1.5 Doench 2016 Rule Set 2 on-target scoring
        logging.info("Computing Doench 2016 on-target scores...")
        guides_20bp = df[seq_col].str[:self.args.target_len].tolist()
        pams = df[seq_col].str[self.args.target_len:self.args.target_len + self.args.pam_len].tolist()

        # Extract flanking sequence from reference for accurate features
        genome = Fasta(self.args.ref)

        def revcomp(seq):
            comp = {'A':'T','T':'A','C':'G','G':'C','N':'N'}
            return ''.join(comp.get(c, 'N') for c in reversed(seq.upper()))

        def get_flank(chrom, guide_start, strand):
            # Returns the genomic flanks around the guide+PAM in guide
            # orientation: (5bp 5' flank, 3bp 3' flank). Consumers slice to the
            # Azimuth/RS3 convention of 4bp 5' + 3bp 3'.
            if strand == '-':
                # On the - strand the 23-mer at fwd[guide_start, guide_start+23)
                # is read as the reverse complement, so the PAM sits at the
                # LOW-coordinate end: PAM is fwd[g, g+3) and the protospacer is
                # fwd[g+3, g+23). The guide-orientation 5' flank is therefore the
                # bases AFTER the protospacer on the forward strand.
                left_start = guide_start + self.args.target_len + self.args.pam_len
                left = str(genome[chrom][left_start:left_start + 5].seq) if left_start + 5 <= len(genome[chrom]) else None
                # 3' flank is BEFORE the PAM on the forward strand.
                right_end = guide_start
                right = str(genome[chrom][right_end - 3:right_end].seq) if right_end >= 3 else None
                return (revcomp(left) if left else None,
                        revcomp(right) if right else None)
            # Forward strand
            left = str(genome[chrom][guide_start - 5:guide_start].seq) if guide_start >= 5 else None
            right_start = guide_start + self.args.target_len + self.args.pam_len
            right = str(genome[chrom][right_start:right_start + 3].seq) if right_start + 3 <= len(genome[chrom]) else None
            return left, right

        guide_starts = (df['Start'].astype(int) + df['start'].astype(int)).values
        strands = df['Strand'].values
        chroms = df['Chromosome'].values
        flanks = [get_flank(c, g, s) for c, g, s in zip(chroms, guide_starts, strands)]

        df['Doench2016OnTarget'] = score_guides(guides_20bp, flanks=flanks, pams=pams)
        logging.info("Doench 2016 scores computed for %d guides", len(df))

        # Rule Set 3 scoring (optional)
        if getattr(self.args, 'rs3', False):
            try:
                from .rs3score import score_guides as rs3_score_guides
                # rs3 needs the Azimuth 30mer layout:
                #   4bp 5' flank + 20bp guide + 3bp PAM + 3bp 3' flank,
                # guide at offset 4..23, PAM at 24..26. (CRISPOR feeds
                # s[26:56] of its 100bp window, which is this layout.)
                rs3_seqs = []
                for i in range(len(df)):
                    flank5 = flanks[i][0] if flanks[i][0] else "NNNN"
                    flank3 = flanks[i][1] if flanks[i][1] else "NNN"
                    left = flank5[-4:] if len(flank5) >= 4 else flank5.rjust(4, "N")
                    right = flank3[:3] if len(flank3) >= 3 else flank3.ljust(3, "N")
                    s30 = (left + guides_20bp[i] + pams[i] + right).upper()
                    rs3_seqs.append(s30)
                df['RuleSet3'] = rs3_score_guides(rs3_seqs)
                logging.info("Rule Set 3 scores computed for %d guides", len(df))
            except ImportError:
                logging.warning("rs3score module or sglearn/lightgbm not available, skipping RS3 scores")
            except Exception as e:
                logging.warning("Rule Set 3 scoring failed: %s", e)

        # 2. SNP Entropy and Annotation
        df['entropy'] = 0.0
        df['snp_overlaps'] = '.'
        df['num_snps'] = 0
        if self.args.snps:
            if not os.path.exists(self.args.snps):
                logging.error(f"SNP file not found: {self.args.snps}")
                sys.exit(1)
            logging.info(f"Calculating SNP entropy using {self.args.snps}...")
            snp_gr = pr.read_bed(self.args.snps)
            if snp_gr.df.shape[1] < 3:
                raise ValueError(f"SNP BED {self.args.snps} is malformed. Expected at least 3 columns (Chrom, Start, End).")

            # Build overlap table using per-chromosome binary search
            # (avoids OOM from pandas merge cross-product)
            target_df = target_gr.df.copy()
            target_df['orig_idx'] = target_df.index

            snp_df = snp_gr.df.copy()
            snp_name_col = 'Name' if 'Name' in snp_df.columns else 'name'

            records = []
            for chrom in target_df['Chromosome'].unique():
                tgt_chr = target_df[target_df['Chromosome'] == chrom]
                snp_chr = snp_df[snp_df['Chromosome'] == chrom]
                if tgt_chr.empty or snp_chr.empty:
                    continue

                tgt_sorted = tgt_chr.sort_values('Start')
                tgt_starts = tgt_sorted['Start'].values
                tgt_ends = tgt_sorted['End'].values
                tgt_orig = tgt_sorted['orig_idx'].values
                tgt_strand = tgt_sorted['Strand'].values

                snp_positions = snp_chr['Start'].values
                snp_names = snp_chr[snp_name_col].values

                right_idx = np.searchsorted(tgt_starts, snp_positions, side='right') - 1

                for i in range(len(snp_positions)):
                    idx = right_idx[i]
                    while idx >= 0 and tgt_starts[idx] <= snp_positions[i]:
                        if tgt_ends[idx] > snp_positions[i]:
                            records.append((
                                tgt_orig[idx], snp_positions[i], snp_names[i],
                                tgt_strand[idx], tgt_starts[idx], tgt_ends[idx],
                            ))
                        idx -= 1

            if records:
                j_df = pd.DataFrame(records, columns=[
                    'orig_idx', 'snp_pos', 'Name', 'Strand', 't_start', 't_end'
                ])
            else:
                j_df = pd.DataFrame(columns=[
                    'orig_idx', 'snp_pos', 'Name', 'Strand', 't_start', 't_end'
                ])

            # 1. Calculate relative positions
            j_df['rel_pos'] = np.where(
                j_df['Strand'] == '+',
                j_df['snp_pos'] - j_df['t_start'],
                j_df['t_end'] - j_df['snp_pos'] - 1
            )

            # 2. Convert SNP alleles to IUPAC code for entropy calculation
            j_df['snp_code'] = j_df['Name'].apply(snp_to_iupac)
            # Keep original for annotation display
            j_df['snp_code_orig'] = j_df['snp_code']
            # Concrete ALT allele (ref/alt parsed from the SNP name) for variant
            # rescoring; complemented for minus-strand guides so it is in the
            # guide's own sequence orientation.
            j_df['alt_allele'] = j_df['Name'].apply(_parse_alt)
            minus_mask = j_df['Strand'] == '-'
            j_df.loc[minus_mask, 'snp_code'] = j_df.loc[minus_mask, 'snp_code'].apply(complement_iupac)
            j_df.loc[minus_mask, 'alt_allele'] = j_df.loc[minus_mask, 'alt_allele'].apply(complement_iupac)

            # 3. GROUPED MUTATION: Only loop over targets that have at least one SNP
            grouped = j_df.groupby('orig_idx')
            entropy_results = {}
            snp_annotations = {}
            snp_variant_info = {}

            for idx, group in grouped:
                # Get original sequence
                seq_list = list(df.at[idx, seq_col])
                seq_len = len(seq_list)
                
                overlaps = []
                variants = []
                for row in group.itertuples():
                    if 0 <= row.rel_pos < seq_len:
                        seq_list[row.rel_pos] = row.snp_code
                        overlaps.append(f"{row.rel_pos}:{row.Name}>{row.snp_code_orig}")
                        if row.alt_allele:
                            variants.append((int(row.rel_pos), row.alt_allele, str(row.Name)))
                
                entropy_results[idx] = calculate_entropy(seq_list)
                snp_annotations[idx] = '|'.join(overlaps) if overlaps else '.'
                snp_variant_info[idx] = variants

            # Map results back to the main dataframe
            df['entropy'] = df.index.map(entropy_results).fillna(0.0)
            df['snp_overlaps'] = df.index.map(snp_annotations).fillna('.')
            df['num_snps'] = df['snp_overlaps'].apply(lambda x: 0 if x == '.' else x.count('|') + 1)

            # Re-score guides that overlap SNPs (if requested): compute RS2/RS3 on
            # variant sequences (concrete ALT allele, real flanks, real PAM) and
            # report per-guide avg/worst single-SNP and all-SNP scores, deltas vs
            # reference, and the min distance from a SNP to the PAM and cut site.
            if getattr(self.args, 'rescore_snps', False) and snp_variant_info:
                rs3_fn = None
                if 'RuleSet3' in df.columns:
                    try:
                        from .rs3score import score_guides as rs3_fn
                    except ImportError:
                        logging.warning("rs3score unavailable; SNP rescoring RS2 only")

                res = {}
                for idx, variants in snp_variant_info.items():
                    if not variants:
                        continue
                    ref_rs2 = float(df.at[idx, 'Doench2016OnTarget'])
                    ref_rs3 = float(df.at[idx, 'RuleSet3']) if 'RuleSet3' in df.columns else np.nan
                    g = guides_20bp[idx]
                    p = pams[idx]
                    f5, f3 = flanks[idx]

                    per = [_score_variant(g, p, f5, f3, {rp: alt}, rs3_fn) for rp, alt, _ in variants]
                    rs2s = [x[0] for x in per]
                    rs3s = [x[1] for x in per]
                    all_rs2, all_rs3 = _score_variant(g, p, f5, f3, {rp: alt for rp, alt, _ in variants}, rs3_fn)

                    row = {
                        'Doench2016OnTarget_avg1snp': float(np.mean(rs2s)),
                        'Doench2016OnTarget_worst1snp': float(np.min(rs2s)),
                        'Doench2016OnTarget_allsnps': all_rs2,
                        'Doench2016_delta_avg': float(np.mean(rs2s)) - ref_rs2,
                        'Doench2016_delta_worst': float(np.min(rs2s)) - ref_rs2,
                        'Doench2016_delta_all': all_rs2 - ref_rs2,
                        # distance from SNP to PAM boundary (23-mer pos 20): >0 in
                        # guide, <=0 inside PAM. Cut site is between guide 17/18.
                        'snp_pam_dist': int(min(20 - rp for rp, _, _ in variants)),
                        'snp_cut_dist': int(min(abs(rp - 17) for rp, _, _ in variants)),
                    }
                    if rs3_fn is not None and not np.isnan(ref_rs3):
                        row.update({
                            'RuleSet3_avg1snp': float(np.mean(rs3s)),
                            'RuleSet3_worst1snp': float(np.min(rs3s)),
                            'RuleSet3_allsnps': all_rs3,
                            'RuleSet3_delta_avg': float(np.mean(rs3s)) - ref_rs3,
                            'RuleSet3_delta_worst': float(np.min(rs3s)) - ref_rs3,
                            'RuleSet3_delta_all': all_rs3 - ref_rs3,
                        })
                    res[idx] = row

                if res:
                    res_df = pd.DataFrame.from_dict(res, orient='index')
                    for col in res_df.columns:
                        df[col] = df.index.map(res_df[col])
                    logging.info("SNP-variant rescoring done for %d guides (RS2%s)",
                                 len(res), " + RS3" if rs3_fn is not None else "")
            elif getattr(self.args, 'rescore_snps', False):
                logging.info("No SNPs found to re-score")
            else:
                logging.info("SNP rescoring skipped (use --rescore-snps to enable)")

        # 3. Oligo Design
        oligos = df.apply(
            lambda row: design_oligos(
                row[seq_col],
                row['Strand'],
                target_len=self.args.target_len,
                pam_len=self.args.pam_len
            ),
            axis=1
        )
        df[['fwd', 'rev']] = pd.DataFrame(oligos.tolist(), index=df.index)

        # 3.5 Deduplication
        # A guide is redundant if it targets the same gene with the same sequence.
        # We keep the first occurrence (which, after sorting, would be the best).
        logging.info(f"Deduplicating targets... Original count: {len(df)}")
        df = df.drop_duplicates(subset=['gene_id', 'target', 'context'])
        logging.info(f"Deduplication complete. Unique guides: {len(df)}")

        # 4. Final Ranking (using Doench 2016, aligned with CRISPOR)
        df['over_cutoff'] = df['Doench2016OnTarget'] > self.args.cutoff

        # Determine sorting priority
        # Base priority: Cutoff (Desc) -> Doench 2016 Score (Desc)
        sort_cols = ['over_cutoff', 'Doench2016OnTarget']
        asc_vals = [False, False]

        # If RS3 is available, add as secondary sort
        if 'RuleSet3' in df.columns:
            sort_cols.append('RuleSet3')
            asc_vals.append(False)

        if self.args.accessome and os.path.exists(self.args.accessome):
            # Priority: Accessibility (Desc)
            sort_cols.append('accessibility')
            asc_vals.append(False)
        else:
            # Priority: Position in gene (Asc)
            sort_cols.append('start')
            asc_vals.append(True)

        # Final priority: Entropy (Asc)
        sort_cols.append('entropy')
        asc_vals.append(True)

        df.sort_values(by=sort_cols, ascending=asc_vals, inplace=True)

        # Return full sorted list for selection
        final = df

        # Strict mode: keep only guides that pass the Doench 2016 cutoff, so the
        # output holds just the usable candidates. Ranking is done above, so the
        # kept rows are already in final order.
        if getattr(self.args, 'strict', False):
            kept = int(final['over_cutoff'].sum())
            dropped = len(final) - kept
            final = final[final['over_cutoff']]
            logging.info("Strict mode: kept %d guides over cutoff %.2f, dropped %d",
                         kept, self.args.cutoff, dropped)

        # Column order: our own columns first, grouped so related scores sit next to
        # each other (identity, coordinates, reference scores, guide sequence, SNP
        # annotation, variant scores). The default off-target metrics (Hsu2013,
        # otCount) live in the reference block; any further FlashFry metric the
        # user requested through --metrics trails at the end, as pure pass-through.
        _own_order = [
            'gene_id', 'feature', 'over_cutoff',
            'Chromosome', 'Start', 'End', 'Strand', 'start', 'stop', 'contig',
            'Doench2016OnTarget', 'RuleSet3', 'Hsu2013', 'otCount', 'accessibility',
            'target', 'context', 'orientation', 'overflow', 'fwd', 'rev',
            'num_snps', 'snp_overlaps', 'entropy', 'snp_pam_dist', 'snp_cut_dist',
            'Doench2016OnTarget_avg1snp', 'Doench2016OnTarget_worst1snp', 'Doench2016OnTarget_allsnps',
            'Doench2016_delta_avg', 'Doench2016_delta_worst', 'Doench2016_delta_all',
            'RuleSet3_avg1snp', 'RuleSet3_worst1snp', 'RuleSet3_allsnps',
            'RuleSet3_delta_avg', 'RuleSet3_delta_worst', 'RuleSet3_delta_all',
        ]
        ordered = [c for c in _own_order if c in final.columns]
        ordered += [c for c in final.columns if c not in ordered]
        final = final[ordered]

        final.to_csv(os.path.join(self.outdir, "final_sgrnas.tsv"), sep='\t', index=False)
        logging.info("Done. Results saved to %s (%d rows)",
                     os.path.join(self.outdir, "final_sgrnas.tsv"), len(final))

        # Genome browser BED output
        if getattr(self.args, 'browser_bed', False):
            self._write_browser_bed(final)

    def _write_browser_bed(self, df):
        """Write a BED6 file for genome browser (IGV, UCSC).

        Colors:
          green  = passes cutoff, no SNPs
          orange = passes cutoff, has SNPs
          red    = fails cutoff
        """
        bed_path = os.path.join(self.outdir, "final_sgrnas.bed")
        # Sort by chromosome then absolute position (numeric)
        bed_df = df.copy()
        bed_df['abs_start'] = bed_df['Start'].astype(int) + bed_df['start'].astype(int)
        bed_df['abs_end'] = bed_df['Start'].astype(int) + bed_df['stop'].astype(int)
        bed_df = bed_df.sort_values(['Chromosome', 'abs_start']).reset_index(drop=True)

        with open(bed_path, 'w') as f:
            f.write("track name=CrispyGuides type=bed visibility=2 itemRgb=On\n")
            for _, row in bed_df.iterrows():
                chrom = row['Chromosome']
                abs_start = int(row['Start']) + int(row['start'])
                abs_end = int(row['Start']) + int(row['stop'])
                strand = row.get('Strand', '.')
                gene = row.get('gene_id', '.')
                seq = row.get('target', '.')
                n_snps = int(row.get('num_snps', 0))
                score_raw = row.get('RuleSet3', row.get('Doench2016OnTarget', 0))
                if pd.notna(score_raw):
                    if 'RuleSet3' in df.columns:
                        score_val = int(score_raw * 100)
                    else:
                        score_val = int(score_raw * 1000)
                else:
                    score_val = 0
                over = row.get('over_cutoff', False)

                if not over:
                    rgb = "200,0,0"       # red — fails cutoff
                elif n_snps > 0:
                    rgb = "255,140,0"     # orange — passes but has SNPs
                else:
                    rgb = "0,180,0"       # green — passes, no SNPs

                snp_tag = f" [SNPs:{n_snps}]" if n_snps > 0 else ""
                name = f"{gene}|{seq}{snp_tag}"
                f.write(f"{chrom}\t{abs_start}\t{abs_end}\t{name}\t{score_val}\t{strand}\t{rgb}\n")
        logging.info("Genome browser BED written to %s", bed_path)


def resolve_flashfry_jar(explicit):
    """Resolve the FlashFry JAR: --flashfry-jar, then $CRISPY_FLASHFRY_JAR.

    The JAR is an untracked external binary and is not shipped inside the
    package, so there is no file-relative default to fall back on.
    """
    if explicit:
        return explicit
    env_jar = os.environ.get("CRISPY_FLASHFRY_JAR")
    if env_jar:
        return env_jar
    sys.exit(
        "ERROR: FlashFry JAR not specified. Pass --flashfry-jar or set the "
        "CRISPY_FLASHFRY_JAR environment variable. Download it from "
        "https://github.com/mckennalab/FlashFry/releases"
    )

def main():
    """CLI entry point: parse args, build the pipeline, and run it."""
    # Pre-parse to check for config file
    pre_parser = argparse.ArgumentParser(add_help=False)
    pre_parser.add_argument("--config", type=str, help="Path to YAML config file")
    pre_args, remaining_args = pre_parser.parse_known_args()

    config = load_config(pre_args.config) if pre_args.config else {}

    parser = argparse.ArgumentParser(description="Crispy CRISPR Target Pipeline")
    parser.add_argument("--config", type=str, help="Path to YAML config file")
    parser.add_argument("--genes", required=False, help="CSV of target genes names (optional)")
    parser.add_argument("--gff", default=config.get('gff'), required=not config.get('gff'), help="Reference GFF file")
    parser.add_argument("--ref", default=config.get('ref'), required=not config.get('ref'), help="Reference genome FASTA")
    parser.add_argument("--accessome", required=False, default=config.get('accessome'), help="accessome.bed file (optional)")
    parser.add_argument("--snps", default=config.get('snps'), help="Optional SNP bed file")
    parser.add_argument("--filter-bed", default=config.get('filter_bed'), help="Optional region filter bed file")
    parser.add_argument("--tss-features", type=str, default=config.get('tss_features', "transcript,mRNA,noncoding_transcript,lncRNA,sRNA,tRNA,rRNA"), help="Comma-separated list of GFF3 features to use as anchors for the TSS; promoter regions are then calculated upstream of these anchors (default: transcript,mRNA,noncoding_transcript,lncRNA,sRNA,tRNA,rRNA)")
    parser.add_argument("--promoter-upstream", type=int, default=config.get('promoter_upstream', 1107), help="Distance upstream of TSS to include in promoter region (default: 1107)")
    parser.add_argument("--promoter-downstream", type=int, default=config.get('promoter_downstream', 107), help="Distance downstream of TSS to include in promoter region (default: 107)")
    parser.add_argument("--cutoff", type=float, default=config.get('cutoff', 0.65), help="Doench score cutoff")
    parser.add_argument("--target-len", type=int, default=config.get('target_len', 20), help="Expected target sequence length (default 20)")
    parser.add_argument("--pam-len", type=int, default=config.get('pam_len', 3), help="Expected PAM length (default 3)")
    parser.add_argument("--mem", type=str, default=config.get('mem', "8g"), help="Memory limit for FlashFry (e.g., 4g, 8g, 16g)")
    parser.add_argument("--metrics", type=str, default=config.get('metrics', "hsu2013"), help="FlashFry scoring metrics, comma-separated (default: hsu2013, the off-target specificity score; add on-target metrics such as doench2014ontarget to request them, they land as trailing columns)")
    parser.add_argument("--rescore-snps", action="store_true", help="Re-score SNP-overlapping guides with RS2/RS3 variant sequences (avg/worst single-SNP and all-SNP scores plus deltas)")
    parser.add_argument("--strict", action="store_true",
                        help="Write only guides over --cutoff to final_sgrnas.tsv (default: keep all, ranked)")
    parser.add_argument("--rs3", action="store_true",
                        help="Compute Rule Set 3 scores (requires sglearn and lightgbm)")
    parser.add_argument("--browser-bed", action="store_true",
                        help="Output a BED file for genome browser visualization")
    parser.add_argument("--outdir", type=str, default="./output",
                        help="Output directory (default: ./output)")
    parser.add_argument("--flashfry-jar", type=str, default=None,
                        help="Path to FlashFry JAR (default: $CRISPY_FLASHFRY_JAR)")
    parser.add_argument("--verbose", action="store_true", help="Show FlashFry Java output")

    args = parser.parse_args()
    args.flashfry_jar = resolve_flashfry_jar(args.flashfry_jar)

    os.makedirs(args.outdir, exist_ok=True)
    # Per-run log lives in the outdir alongside the other outputs.
    setup_logging(log_name=os.path.join(args.outdir, "design.log"))
    pipeline = CrispyPipeline(args)
    t_start = time.perf_counter()

    targets = pipeline.prepare_targets()
    if targets is None or len(targets) == 0:
        logging.error("No target regions were found in the GFF. Check your feature labels (exon, CDS, mRNA, transcript).")
        sys.exit(1)
    t_prepared = time.perf_counter()

    scored_file = pipeline.run_flashfry(targets)
    t_scored = time.perf_counter()

    pipeline.post_process(scored_file)
    t_done = time.perf_counter()

    logging.info("Timing: prepare targets %.1fs | FlashFry index/discover/score %.1fs | post-process %.1fs | total %.1fs",
                 t_prepared - t_start, t_scored - t_prepared, t_done - t_scored, t_done - t_start)


if __name__ == "__main__":
    main()
