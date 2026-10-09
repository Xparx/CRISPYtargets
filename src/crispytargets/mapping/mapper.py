import pandas as pd
import urllib.parse
import logging

class GFF3Resolver:
    def __init__(self, gff3_path):
        self.id_to_parent = {}
        self.id_to_type = {}
        self._parse_gff3(gff3_path)

    def _normalize(self, fid):
        """Removes 'KY21:' prefix to match BLAST qseqids."""
        return fid.split(':', 1)[-1] if fid and ':' in fid else fid

    def _parse_gff3(self, path):
        with open(path, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                cols = line.strip().split('\t')
                if len(cols) < 9: continue

                feat_type = cols[2]
                attrs_str = cols[8]
                attrs = {}
                for item in attrs_str.split(';'):
                    if '=' in item:
                        k, v = item.split('=', 1)
                        attrs[k] = urllib.parse.unquote(v)

                fid = self._normalize(attrs.get('ID'))
                parent_raw = attrs.get('Parent')

                if fid:
                    self.id_to_type[fid] = feat_type
                    if parent_raw:
                        # Store all parents as a list of normalized IDs
                        self.id_to_parent[fid] = [self._normalize(p) for p in parent_raw.split(',')]

    def resolve_to_gene(self, fid):
        """Finds the closest ancestor labeled as 'gene' using BFS."""
        queue = [fid]
        visited = {fid}
        
        while queue:
            current = queue.pop(0)
            
            # Check if the current node is a gene
            if self.id_to_type.get(current) == 'gene':
                return current
            
            # Add parents to the queue to continue climbing
            parents = self.id_to_parent.get(current, [])
            for p in parents:
                if p not in visited:
                    visited.add(p)
                    queue.append(p)
                    
        return fid # Fallback to original ID if no gene ancestor is found

def create_mapping(gff3_path, blast_path, species_label, output_path, evalue_threshold=1e-50, include_top=False):
    resolver = GFF3Resolver(gff3_path)
    cols = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
            'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']

    def process_species(path):
        try:
            df = pd.read_csv(path, sep='\t', names=cols)
            if df.shape[1] != 12:
                raise ValueError(f"Invalid BLAST format in {path}. Expected 12 columns (outfmt 6), found {df.shape[1]}.")
        except FileNotFoundError:
            logging.exception(f"File {path} not found.")
            return pd.DataFrame()

        # 1. Resolve all hits to the root gene_id
        df['gene_id'] = df['qseqid'].apply(resolver.resolve_to_gene)
        
        # 2. Filter by E-value
        is_significant = df['evalue'] <= evalue_threshold
        has_significant = is_significant.groupby(df['gene_id']).transform('any')
        
        # We sort by bitscore descending to ensure the best hit is first
        df = df.sort_values(by=['gene_id', 'bitscore'], ascending=[True, False])
        
        # Identify top hits per gene
        is_top_hit = ~df.duplicated(subset=['gene_id'], keep='first')
        mask = (df['evalue'] <= evalue_threshold) | (~has_significant & is_top_hit & include_top)
        filtered_df = df[mask].copy()

        if filtered_df.empty:
            return pd.DataFrame(columns=['gene_id', 'primary', 'secondary'])

        # 3. Optimal Hit Logic: Group by gene_id to find Primary and Secondary
        results = []
        for gene_id, group in filtered_df.groupby('gene_id'):
            # The first row is the best hit (highest bitscore)
            top_hit = group.iloc[0]['sseqid'].split('|')[-1] if isinstance(group.iloc[0]['sseqid'], str) else group.iloc[0]['sseqid']
            
            # All other unique genes are secondary
            others = group['sseqid'].unique()
            others = [s.split('|')[-1] if isinstance(s, str) else s for s in others]
            others = [o for o in others if o != top_hit]
            
            results.append({
                'gene_id': gene_id,
                'primary': top_hit,
                'secondary': '|'.join(others)
            })

        return pd.DataFrame(results)

    # --- Incremental Backbone Logic ---
    # 1. Load existing output or initialize from GFF3
    try:
        backbone = pd.read_csv(output_path, sep='\t')
        logging.info(f"Loaded existing mapping file: {output_path}")
    except (FileNotFoundError, pd.errors.EmptyDataError):
        logging.info("No existing mapping file found. Initializing backbone from GFF3.")
        # Get all unique genes from the GFF3 to ensure complete coverage
        all_genes = sorted(list(set(resolver.id_to_type.keys())))
        # Only keep those that are actually labeled as 'gene'
        all_genes = [g for g in all_genes if resolver.id_to_type.get(g) == 'gene']
        backbone = pd.DataFrame({'gene_id': all_genes})

    # 2. Process the new species
    species_df = process_species(blast_path)
    if species_df.empty:
        logging.warning(f"No significant hits found for {species_label}. Backbone unchanged.")
        return backbone

    # 3. Rename columns to include species label
    species_df = species_df.rename(columns={
        'primary': f"{species_label}_primary",
        'secondary': f"{species_label}_secondary"
    })

    # 4. Collision Handling: Remove existing columns for this species if they exist
    cols_to_drop = [f"{species_label}_primary", f"{species_label}_secondary"]
    backbone = backbone.drop(columns=[c for c in cols_to_drop if c in backbone.columns])

    # 5. Left Join new species data onto the backbone
    final_df = pd.merge(backbone, species_df, on='gene_id', how='left')

    # Save the updated backbone
    final_df.to_csv(output_path, sep='\t', index=False)
    logging.info(f"Updated {output_path} with {species_label} mappings.")

    return final_df
