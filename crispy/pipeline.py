import argparse
import subprocess
import os
import pandas as pd
import pyranges as pr
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
        print("Preparing target regions...")
        # Load GFF
        gr = pr.read_gff3(self.args.gff)

        # 1. Get Exons
        exons = gr[(gr.df.Feature == 'exon') | (gr.df.Feature == 'CDS')]
        exons_df = exons.df.drop_duplicates(subset=['Chromosome', 'Start', 'End', 'Strand'])
        exons = pr.PyRanges(exons_df)

        # 2. Get Promoters (TSS - 1107 to TSS + 107)
        # Note: In GFF3, the top-level feature is often 'mRNA', 'transcript', or 'noncoding_transcript'
        tpts = gr[gr.df.Feature.isin(['transcript', 'mRNA', 'noncoding_transcript'])]
        # tpts = gr[gr.feature == 'transcript']

        # Positive strand: TSS is start. Range: [start - 1107, start + 107]
        pos_tpts = tpts[tpts.df['Strand'] == '+']
        pos_df = pos_tpts.df.copy()
        pos_df['Feature'] = 'promoter' # Explicitly label as promoter
        tss_pos = pos_df['Start']
        pos_df['Start'] = tss_pos - 1107
        pos_df['End'] = tss_pos + 107
        pos_promoters = pr.PyRanges(pos_df)

        # Negative strand: TSS is end. Range: [end - 107, end + 1107]
        neg_tpts = tpts[tpts.df['Strand'] == '-']
        neg_df = neg_tpts.df.copy()
        neg_df['Feature'] = 'promoter' # Explicitly label as promoter
        tss_neg = neg_df['End']
        neg_df['Start'] = tss_neg - 107
        neg_df['End'] = tss_neg + 1107
        neg_promoters = pr.PyRanges(neg_df)

        # Combine positive and negative promoters by concatenating their DataFrames
        promoters_df = pd.concat([pos_promoters.df, neg_promoters.df])
        promoters = pr.PyRanges(promoters_df)

        # Combine exons and promoters by concatenating their DataFrames
        targets_df = pd.concat([exons.df, promoters.df])
        targets = pr.PyRanges(targets_df)

        # Optional: Filter by provided gene list
        if self.args.genes:
            print(f"Filtering by gene list: {self.args.genes}")
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

            print(f"Expanded {len(seed_ids)} seed IDs to {len(expanded_ids)} total related features.")

            # 3. Filter targets: Keep if ID or Parent is in the expanded set
            def is_related(row):
                # Convert to string to ensure exact match against the set
                row_id = str(getattr(row, 'ID', ''))
                row_parent = str(getattr(row, 'Parent', ''))
                return (row_id in expanded_ids) or (row_parent in expanded_ids)

            targets = targets[targets.df.apply(is_related, axis=1)]

        else:
            print("No gene list provided. Processing all genes in GFF...")

        # Optional Filter BED
        if self.args.filter_bed and os.path.exists(self.args.filter_bed):
            print(f"Applying optional filter: {self.args.filter_bed}")
            filter_gr = pr.read_bed(self.args.filter_bed)
            targets = targets.intersect(filter_gr)

        # --- Build Gene Lineage Map ---
        # This resolves the hierarchy: CDS -> mRNA -> Gene
        self.gene_map = {}
        gff_df = gr.df
        
        # Create quick lookup maps for parents and feature types
        # We use astype(str) to avoid issues with NaN/float types
        parent_lookup = dict(zip(gff_df['ID'].astype(str), gff_df['Parent'].astype(str)))
        feature_lookup = dict(zip(gff_df['ID'].astype(str), gff_df['Feature'].astype(str)))

        def resolve_to_gene(fid):
            """Recursively find the top-level gene ID for a given feature ID."""
            curr = fid
            visited = set() # Prevent infinite loops in malformed GFFs
            while curr in parent_lookup and feature_lookup.get(curr) != 'gene' and curr not in visited:
                visited.add(curr)
                curr = parent_lookup[curr]
            return curr

        # Map every ID in the GFF to its root gene
        for fid in gff_df['ID'].dropna().astype(str):
            self.gene_map[fid] = resolve_to_gene(fid)
        # -----------------------------

        targets.to_bed("targets.bed")
        return targets

    def extract_sequences(self, targets_gr):
        print("Extracting sequences for FlashFry...")
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
        print("Running FlashFry...")
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
        print("Post-processing and scoring...")
        df = pd.read_csv(scored_file, sep='\t')

        # Parse coordinates for all targets (needed for sorting and accessibility)
        parsed_coords = df['contig'].apply(parse_flashfry_contig)
        # Unpack 6 values: gene_id, feature, Chromosome, Start, End, Strand
        df[['gene_id', 'feature', 'Chromosome', 'Start', 'End', 'Strand']] = pd.DataFrame(parsed_coords.tolist(), index=df.index)
        target_gr = pr.PyRanges(df[['Chromosome', 'Start', 'End', 'Strand']])

        # 1. Accessibility (Optional)
        df['accessibility'] = 0
        if self.args.accessome and os.path.exists(self.args.accessome):
            print(f"Calculating accessibility using {self.args.accessome}...")
            acc_gr = pr.read_bed(self.args.accessome)

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
            print("No accessome file provided. Skipping accessibility scoring.")

        # 2. SNP Entropy
        df['entropy'] = 0.0
        if self.args.snps and os.path.exists(self.args.snps):
            print(f"Calculating SNP entropy using {self.args.snps}...")
            snp_gr = pr.read_bed(self.args.snps)

            # Add original index and coordinates to target_gr to track back to df
            target_gr_indexed = target_gr.assign(
                orig_idx=target_gr.df.index,
                t_start=target_gr.start,
                t_end=target_gr.end
            )

            # Join targets with SNPs to get all overlapping pairs
            joined = target_gr_indexed.join(snp_gr)

            # Dictionary to store mutated sequences: {index: sequence_list}
            mutated_seqs = {}

            for row in joined.df.itertuples():
                idx = row.orig_idx
                # Initialize with original sequence if not already mutated
                if idx not in mutated_seqs:
                    mutated_seqs[idx] = list(df.at[idx, 'sequence'])

                # Calculate relative position in the FlashFry sequence
                # row.Start_2 is the SNP start, row.t_start/t_end are target boundaries
                snp_start = row.Start_2
                t_start = row.t_start
                t_end = row.t_end
                strand = row.Strand

                if strand == '+':
                    rel_pos = snp_start - t_start
                else:
                    # For minus strand, the sequence is RC, so position is mirrored
                    rel_pos = t_end - snp_start - 1

                if 0 <= rel_pos < len(mutated_seqs[idx]):
                    snp_code = str(row.Name)
                    if strand == '-':
                        snp_code = complement_iupac(snp_code)
                    mutated_seqs[idx][rel_pos] = snp_code

            # Calculate entropy on the mutated sequences
            for idx, seq_list in mutated_seqs.items():
                df.at[idx, 'entropy'] = calculate_entropy("".join(seq_list))

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

        print(f"Using column '{seq_col}' as the target sequence source.")

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
        print("Done. Results saved to final_sgrnas.tsv")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Crispy CRISPR Target Pipeline")
    parser.add_argument("--genes", required=False, help="CSV of target genes names (optional)")
    parser.add_argument("--gff", required=True, help="Reference GFF file")
    parser.add_argument("--ref", required=True, help="Reference genome FASTA")
    parser.add_argument("--accessome", required=False, help="accessome.bed file (optional)")
    parser.add_argument("--snps", help="Optional SNP bed file")
    parser.add_argument("--filter-bed", help="Optional region filter bed file")
    parser.add_argument("--cutoff", type=float, default=0.65, help="Doench score cutoff")
    parser.add_argument("--target-len", type=int, default=20, help="Expected target sequence length (default 20)")
    parser.add_argument("--pam-len", type=int, default=3, help="Expected PAM length (default 3)")
    parser.add_argument("--mem", type=str, default="8g", help="Memory limit for FlashFry (e.g., 4g, 8g, 16g)")
    parser.add_argument("--metrics", type=str, default="doench2014ontarget", help="FlashFry scoring metrics, comma-separated (default: doench2014ontarget)")

    args = parser.parse_args()

    pipeline = CrispyPipeline(args)
    # bed_file = pipeline.prepare_targets()
    targets = pipeline.prepare_targets()
    if targets is None or len(targets) == 0:
        print("ERROR: No target regions were found in the GFF. Check your feature labels (exon, CDS, mRNA, transcript).")
        exit(1)
    scored_file = pipeline.run_flashfry(targets)
    pipeline.post_process(scored_file)
