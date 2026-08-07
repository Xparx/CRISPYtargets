import argparse
import subprocess
import os
import sys
import numpy as np
import pandas as pd
import pyranges as pr
import logging
import yaml
from pyfaidx import Fasta
from utils import calculate_entropy, design_oligos, parse_flashfry_contig
from utils import complement_iupac, snp_to_iupac
from doench2016 import score_guides
import gc


def run_command(cmd, verbose=False):
    logging.info(f"Executing: {cmd}")
    if verbose:
        subprocess.run(cmd, shell=True, check=True)
    else:
        subprocess.run(cmd, shell=True, check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


class CrispyPipeline:
    def __init__(self, args):
        self.args = args
        self.outdir = args.outdir
        self.ff_jar = args.flashfry_jar
        self.db = os.path.join(args.outdir, "htky_db")

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

        # Positive strand: TSS is start.
        pos_tpts = tpts[tpts.df['Strand'] == '+']
        pos_df = pos_tpts.df.copy()
        pos_df['Feature'] = 'promoter' # Explicitly label as promoter
        tss_pos = pos_df['Start']
        pos_df['Start'] = tss_pos - upstream
        pos_df['End'] = tss_pos + downstream
        pos_promoters = pr.PyRanges(pos_df)

        # Negative strand: TSS is end.
        neg_tpts = tpts[tpts.df['Strand'] == '-']
        neg_df = neg_tpts.df.copy()
        neg_df['Feature'] = 'promoter' # Explicitly label as promoter
        tss_neg = neg_df['End']
        neg_df['Start'] = tss_neg - downstream
        neg_df['End'] = tss_neg + upstream
        neg_promoters = pr.PyRanges(neg_df)

        # Combine positive and negative promoters by concatenating their DataFrames
        promoters_df = pd.concat([pos_promoters.df, neg_promoters.df])
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
            # Load seed IDs from CSV
            seed_ids = set(pd.read_csv(self.args.genes)['Gene.ID'].astype(str).tolist())

            # We need the full GFF to resolve the hierarchy
            gff_df = gr.df

            # 1. Initialize expanded_ids with seed_ids.
            # This ensures that if a seed is already a Parent or ID, the loop finds children.
            expanded_ids = set(seed_ids)

            # Also map any features that match seed IDs in Name/name columns to their IDs
            for col in ['Name', 'name']:
                if col in gff_df.columns:
                    matches = gff_df[gff_df[col].astype(str).isin(seed_ids)]
                    if 'ID' in matches.columns:
                        expanded_ids.update(matches['ID'].dropna().astype(str).tolist())

            # 2. Recursive Expansion: Find all children (Gene -> mRNA -> Exon/CDS)
            while True:
                prev_count = len(expanded_ids)
                if 'Parent' in gff_df.columns:
                    # Find features whose Parent is in our current expanded set
                    children = gff_df[gff_df['Parent'].astype(str).isin(expanded_ids)]
                    if 'ID' in children.columns:
                        expanded_ids.update(children['ID'].dropna().astype(str).tolist())

                if len(expanded_ids) == prev_count:
                    break

            logging.info(f"Expanded {len(seed_ids)} seed IDs to {len(expanded_ids)} total related features.")

            # 3. Filter targets: Keep if ID or Parent is in the expanded set
            # Use vectorized .isin() for performance
            id_col = targets.df['ID'].astype(str) if 'ID' in targets.df.columns else pd.Series([''] * len(targets.df))
            parent_col = targets.df['Parent'].astype(str) if 'Parent' in targets.df.columns else pd.Series([''] * len(targets.df))
            
            combined_mask = id_col.isin(expanded_ids) | parent_col.isin(expanded_ids)
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
        logging.info("Extracting sequences for FlashFry...")
        fasta_out = os.path.join(self.outdir, "targets.fasta")
        genome = Fasta(self.args.ref)

        with open(fasta_out, 'w') as f_out:
            for row in targets_gr.df.itertuples():
                # 1. Identify the feature's own ID or its immediate parent
                fid = getattr(row, 'ID', None)
                if fid is None or str(fid) == 'nan':
                    fid = getattr(row, 'Parent', None)
                
                # 2. Use the gene_map to resolve this ID to the top-level Gene ID
                gene_id = self.gene_map.get(str(fid), "unknown")

                seq = genome[row.Chromosome][row.Start:row.End].seq
                # Include row.Feature in the header to track origin (CDS, exon, promoter)
                f_out.write(f">{gene_id}|{row.Feature}|{row.Chromosome}:{row.Start}-{row.End}|{row.Strand}\n{seq}\n")

        return fasta_out

    def run_flashfry(self, targets_gr):
        logging.info("Running FlashFry...")
        # Index if DB doesn't exist
        if not os.path.exists(self.db):
            # Create a temporary directory for FlashFry indexing
            os.makedirs(os.path.join(self.outdir, "build"), exist_ok=True)
            run_command(f"java -Xmx4g -jar {self.ff_jar} index --database {self.db} --reference {self.args.ref} --enzyme spcas9ngg --tmpLocation {os.path.join(self.outdir, 'build')}",
                        verbose=True)
        else:
            logging.info("FlashFry database already exists: %s", self.db)

        # Extract sequences from PyRanges to FASTA
        targets_fasta = self.extract_sequences(targets_gr)

        # Discover
        run_command(f"java -Xmx{self.args.mem} -jar {self.ff_jar} discover --database {self.db} --fasta {targets_fasta} --positionOutput --output {os.path.join(self.outdir, 'targets.output')}",
                    verbose=getattr(self.args, 'verbose', False))

        # Score
        run_command(f"java -Xmx{self.args.mem} -jar {self.ff_jar} score --input {os.path.join(self.outdir, 'targets.output')} --output {os.path.join(self.outdir, 'targets.output.scored')} --database {self.db} --scoringMetrics={self.args.metrics}",
                    verbose=getattr(self.args, 'verbose', False))

        return os.path.join(self.outdir, "targets.output.scored")

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
            if strand == '-':
                # Reverse strand: upstream on reverse = AFTER guide on forward
                left_start = guide_start + self.args.target_len
                left = str(genome[chrom][left_start:left_start + 5].seq) if left_start + 5 <= len(genome[chrom]) else None
                # Downstream on reverse = BEFORE guide on forward
                right_end = guide_start
                right = str(genome[chrom][right_end - 2:right_end].seq) if right_end >= 2 else None
                return (revcomp(left) if left else None,
                        revcomp(right) if right else None)
            # Forward strand
            left = str(genome[chrom][guide_start - 5:guide_start].seq) if guide_start >= 5 else None
            right_start = guide_start + self.args.target_len + self.args.pam_len
            right = str(genome[chrom][right_start:right_start + 2].seq) if right_start + 2 <= len(genome[chrom]) else None
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
                from rs3score import score_guides as rs3_score_guides
                # rs3 needs 30mers: 5bp flank + 20bp guide + 3bp PAM + 2bp flank
                rs3_seqs = []
                for i in range(len(df)):
                    flank5 = flanks[i][0] if flanks[i][0] else "NNNNN"
                    flank3 = flanks[i][1] if flanks[i][1] else "NN"
                    left = flank5[-5:] if len(flank5) >= 5 else flank5.rjust(5, "N")
                    right = flank3[:2] if len(flank3) >= 2 else flank3.ljust(2, "N")
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
            # Complement for minus-strand guides (needed for sequence substitution)
            minus_mask = j_df['Strand'] == '-'
            j_df.loc[minus_mask, 'snp_code'] = j_df.loc[minus_mask, 'snp_code'].apply(complement_iupac)

            # 3. GROUPED MUTATION: Only loop over targets that have at least one SNP
            grouped = j_df.groupby('orig_idx')
            entropy_results = {}
            snp_annotations = {}
            snp_sequences = {}

            for idx, group in grouped:
                # Get original sequence
                seq_list = list(df.at[idx, seq_col])
                seq_len = len(seq_list)
                
                overlaps = []
                for row in group.itertuples():
                    if 0 <= row.rel_pos < seq_len:
                        seq_list[row.rel_pos] = row.snp_code
                        overlaps.append(f"{row.rel_pos}:{row.Name}>{row.snp_code_orig}")
                
                entropy_results[idx] = calculate_entropy(seq_list)
                snp_annotations[idx] = '|'.join(overlaps) if overlaps else '.'
                snp_sequences[idx] = seq_list

            # Map results back to the main dataframe
            df['entropy'] = df.index.map(entropy_results).fillna(0.0)
            df['snp_overlaps'] = df.index.map(snp_annotations).fillna('.')
            df['num_snps'] = df['snp_overlaps'].apply(lambda x: 0 if x == '.' else x.count('|') + 1)

            # Re-score with SNPs incorporated (if requested)
            if getattr(self.args, 'rescore_snps', False) and snp_sequences:
                expected_len = self.args.target_len + self.args.pam_len
                snp_fasta = os.path.join(self.outdir, "targets_snp.fasta")
                with open(snp_fasta, 'w') as f:
                    for idx, seq_list in snp_sequences.items():
                        row = df.loc[idx]
                        f.write(f">{row.gene_id}|{row.feature}|{row.Chromosome}:{row.Start}-{row.End}|{row.Strand}\n{''.join(seq_list)}\n")

                logging.info("Running FlashFry on SNP-modified sequences...")
                run_command(f"java -Xmx{self.args.mem} -jar {self.ff_jar} discover --database {self.db} --fasta {snp_fasta} --positionOutput --output {os.path.join(self.outdir, 'targets_snp.output')}",
                            verbose=getattr(self.args, 'verbose', False))
                run_command(f"java -Xmx{self.args.mem} -jar {self.ff_jar} score --input {os.path.join(self.outdir, 'targets_snp.output')} --output {os.path.join(self.outdir, 'targets_snp.output.scored')} --database {self.db} --scoringMetrics={self.args.metrics}",
                            verbose=getattr(self.args, 'verbose', False))

                snp_df = pd.read_csv(os.path.join(self.outdir, "targets_snp.output.scored"), sep='\t')
                snp_seq_col = None
                for col in snp_df.columns:
                    first_val = snp_df[col].dropna().iloc[0] if not snp_df[col].dropna().empty else None
                    if isinstance(first_val, str) and len(first_val) == expected_len:
                        snp_seq_col = col
                        break
                if snp_seq_col:
                    # Only merge back guides that directly overlap SNPs
                    # (other guides in the same target region are unaffected)
                    snps_contigs = set(df.loc[j_df['orig_idx'].unique(), 'contig'])
                    snp_df_filtered = snp_df[snp_df['contig'].isin(snps_contigs)]

                    score_cols = [c for c in snp_df_filtered.columns
                                  if c in df.columns
                                  and c not in ['chrom', 'start', 'end', 'stop',
                                                'strand', 'contig']]
                    rename_map = {c: f'{c}_snp' for c in score_cols}
                    snp_to_merge = snp_df_filtered[['contig'] + score_cols].rename(columns=rename_map)
                    df = df.merge(snp_to_merge, on='contig', how='left')
                    logging.info("Re-scored %d / %d guides that directly overlap SNPs",
                                 len(snp_df_filtered), len(snp_df))
                    # Clear _snp scores for guides that don't overlap any SNPs
                    snp_score_cols = [c for c in df.columns if c.endswith('_snp')]
                    df.loc[df['num_snps'] == 0, snp_score_cols] = np.nan
                else:
                    logging.warning("Could not find matching sequence column in SNP output")
            elif getattr(self.args, 'rescore_snps', False):
                logging.info("No SNPs found to re-score")
            else:
                logging.info("SNP re-scoring skipped (use --rescore-snps to enable)")

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

        # Move gene_id to first column
        cols = list(final.columns)
        if 'gene_id' in cols:
            cols.remove('gene_id')
            cols.insert(0, 'gene_id')
            final = final[cols]

        final.to_csv(os.path.join(self.outdir, "final_sgrnas.tsv"), sep='\t', index=False)
        logging.info("Done. Results saved to %s", os.path.join(self.outdir, "final_sgrnas.tsv"))

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


def setup_logging(log_name="crispy.log"):
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(logging.Formatter('%(levelname)s: %(message)s'))
    
    file_handler = logging.FileHandler(log_name)
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
    
    logger.addHandler(console_handler)
    logger.addHandler(file_handler)

def load_config(config_path):
    try:
        with open(config_path, 'r') as f:
            return yaml.safe_load(f) or {}
    except Exception as e:
        # We use print here because logging isn't set up yet
        print(f"ERROR: Failed to load config file {config_path}: {e}")
        return {}

if __name__ == "__main__":
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
    parser.add_argument("--metrics", type=str, default=config.get('metrics', "doench2014ontarget"), help="FlashFry scoring metrics, comma-separated (default: doench2014ontarget)")
    parser.add_argument("--rescore-snps", action="store_true", help="Re-run FlashFry scoring with SNPs incorporated into target sequences")
    parser.add_argument("--rs3", action="store_true",
                        help="Compute Rule Set 3 scores (requires sglearn and lightgbm)")
    parser.add_argument("--browser-bed", action="store_true",
                        help="Output a BED file for genome browser visualization")
    parser.add_argument("--outdir", type=str, default="./output",
                        help="Output directory (default: ./output)")
    parser.add_argument("--flashfry-jar", type=str,
                        default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "FlashFry-assembly-1.15.jar"),
                        help="Path to FlashFry JAR (default: crispy/FlashFry-assembly-1.15.jar)")
    parser.add_argument("--verbose", action="store_true", help="Show FlashFry Java output")

    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    setup_logging(log_name=os.path.join(args.outdir, "crispy.log"))
    pipeline = CrispyPipeline(args)
    # bed_file = pipeline.prepare_targets()
    targets = pipeline.prepare_targets()
    if targets is None or len(targets) == 0:
        logging.error("No target regions were found in the GFF. Check your feature labels (exon, CDS, mRNA, transcript).")
        exit(1)
    scored_file = pipeline.run_flashfry(targets)
    pipeline.post_process(scored_file)
