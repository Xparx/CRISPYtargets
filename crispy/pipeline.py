import argparse
import subprocess
import os
import numpy as np
import pandas as pd
import pyranges as pr
import logging
import yaml
from pyfaidx import Fasta
from utils import calculate_entropy, design_oligos, parse_flashfry_contig
from utils import complement_iupac
import gc


def run_command(cmd):
    print(f"Executing: {cmd}")
    subprocess.run(cmd, shell=True, check=True)


class CrispyPipeline:
    def __init__(self, args):
        self.args = args
        self.ff_jar = "FlashFry-assembly-1.15.jar"
        self.db = "htky_db"

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

        targets.to_bed("targets.bed")
        return targets

    def extract_sequences(self, targets_gr):
        logging.info("Extracting sequences for FlashFry...")
        fasta_out = "targets.fasta"
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
            os.makedirs("build", exist_ok=True)
            run_command(f"java -Xmx4g -jar {self.ff_jar} index --database {self.db} --reference {self.args.ref} --enzyme spcas9ngg --tmpLocation ./build")

        # Extract sequences from PyRanges to FASTA
        targets_fasta = self.extract_sequences(targets_gr)

        # Discover
        run_command(f"java -Xmx{self.args.mem} -jar {self.ff_jar} discover --database {self.db} --fasta {targets_fasta} --positionOutput --output targets.output")

        # Score
        # Use the metrics specified in the arguments (defaults to doench2014ontarget)
        run_command(f"java -Xmx{self.args.mem} -jar {self.ff_jar} score --input targets.output --output targets.output.scored --database {self.db} --scoringMetrics={self.args.metrics}")

        return "targets.output.scored"

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

        # 2. SNP Entropy
        df['entropy'] = 0.0
        if self.args.snps and os.path.exists(self.args.snps):
            logging.info(f"Calculating SNP entropy using {self.args.snps}...")
            snp_gr = pr.read_bed(self.args.snps)
            if snp_gr.df.shape[1] < 3:
                raise ValueError(f"SNP BED {self.args.snps} is malformed. Expected at least 3 columns (Chrom, Start, End).")

            # Add original index and coordinates to target_gr to track back to df
            target_gr_indexed = target_gr.assign(
                orig_idx=target_gr.df.index,
                t_start=target_gr.start,
                t_end=target_gr.end
            )

            # Join targets with SNPs to get all overlapping pairs
            joined = target_gr_indexed.join(snp_gr)
            j_df = joined.df

            # 1. VECTORIZED: Calculate relative positions for all overlaps at once
            j_df['rel_pos'] = np.where(
                j_df['Strand'] == '+',
                j_df['Start_2'] - j_df['t_start'],
                j_df['t_end'] - j_df['Start_2'] - 1
            )

            # 2. VECTORIZED: Handle IUPAC complements for minus strand
            j_df['snp_code'] = j_df['Name']
            minus_mask = j_df['Strand'] == '-'
            j_df.loc[minus_mask, 'snp_code'] = j_df.loc[minus_mask, 'Name'].apply(complement_iupac)

            # 3. GROUPED MUTATION: Only loop over targets that have at least one SNP
            grouped = j_df.groupby('orig_idx')
            entropy_results = {}

            for idx, group in grouped:
                # Get original sequence
                seq_list = list(df.at[idx, 'sequence'])
                seq_len = len(seq_list)
                
                for row in group.itertuples():
                    if 0 <= row.rel_pos < seq_len:
                        seq_list[row.rel_pos] = row.snp_code
                
                entropy_results[idx] = calculate_entropy(seq_list)

            # Map results back to the main dataframe
            df['entropy'] = df.index.map(entropy_results).fillna(0.0)

        # 3. Oligo Design
        # Dynamically find the sequence column (it's the one containing the 23bp target+PAM)
        seq_col = None
        expected_len = self.args.target_len + self.args.pam_len
        for col in df.columns:
            # Check if the first non-null value in the column is a string of the correct length
            first_val = df[col].dropna().iloc[0] if not df[col].dropna().empty else None
            if isinstance(first_val, str) and len(first_val) == expected_len:
                seq_col = col
                break

        if seq_col is None:
            raise KeyError(f"Could not find a sequence column with length {expected_len} in the FlashFry output. Columns found: {list(df.columns)}")

        logging.info(f"Using column '{seq_col}' as the target sequence source.")

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

        # 4. Final Ranking
        df['over_cutoff'] = df['Doench2014OnTarget'] > self.args.cutoff

        # Determine sorting priority
        # Base priority: Cutoff (Desc) -> Doench Score (Desc)
        sort_cols = ['over_cutoff', 'Doench2014OnTarget']
        asc_vals = [False, False]

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

        final.to_csv("final_sgrnas.tsv", sep='\t', index=False)
        logging.info("Done. Results saved to final_sgrnas.tsv")


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

    args = parser.parse_args()

    setup_logging()
    pipeline = CrispyPipeline(args)
    # bed_file = pipeline.prepare_targets()
    targets = pipeline.prepare_targets()
    if targets is None or len(targets) == 0:
        logging.error("No target regions were found in the GFF. Check your feature labels (exon, CDS, mRNA, transcript).")
        exit(1)
    scored_file = pipeline.run_flashfry(targets)
    pipeline.post_process(scored_file)
