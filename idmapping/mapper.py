import pandas as pd
import urllib.parse

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
                parent = self._normalize(attrs.get('Parent'))

                if fid:
                    self.id_to_type[fid] = feat_type
                    if parent:
                        # Handle multiple parents by taking the first one
                        self.id_to_parent[fid] = parent.split(',')[0]

    def resolve_to_gene(self, fid):
        """Climbs the hierarchy to find the top-level Gene ID."""
        current = fid
        visited = set()
        while current in self.id_to_parent:
            if current in visited: break
            visited.add(current)
            parent = self.id_to_parent[current]
            if self.id_to_type.get(parent) == 'gene':
                return parent
            current = parent
        return current if self.id_to_type.get(current) == 'gene' else fid

def create_mapping(gff3_path, human_path, mouse_path, output_path, evalue_threshold=1e-50, include_top=False):
    resolver = GFF3Resolver(gff3_path)
    cols = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
            'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']

    def process_species(path, species_label):
        try:
            df = pd.read_csv(path, sep='\t', names=cols)
        except FileNotFoundError:
            print(f"Warning: File {path} not found. Skipping {species_label}.")
            return pd.DataFrame()

        # Group by query ID to apply the conditional threshold logic
        final_hits = []
        for qid, group in df.groupby('qseqid'):
            # 1. Try to get all hits meeting the threshold
            significant = group[group['evalue'] <= evalue_threshold]

            if not significant.empty:
                final_hits.append(significant)
            elif include_top:
                # 2. Fallback: Keep only the top hit (first row) if no significant hits exist
                final_hits.append(group.head(1))

        if not final_hits:
            return pd.DataFrame()

        res_df = pd.concat(final_hits)
        res_df['gene_id'] = res_df['qseqid'].apply(resolver.resolve_to_gene)
        # Add this line:
        res_df['feature_type'] = res_df['qseqid'].apply(lambda x: resolver.id_to_type.get(x, 'unknown'))

        res_df['species_gene_symbol'] = res_df['sseqid'].apply(lambda x: x.split('|')[-1] if isinstance(x, str) else x)
        res_df['species'] = species_label

        # Update the return columns:
        return res_df[['qseqid', 'feature_type', 'gene_id', 'species_gene_symbol', 'species', 'evalue', 'bitscore']]

    human_df = process_species(human_path, 'Human')
    mouse_df = process_species(mouse_path, 'Mouse')

    final_df = pd.concat([human_df, mouse_df], ignore_index=True)
    if not final_df.empty:
        final_df = final_df.sort_values(by=['gene_id', 'evalue'])
        final_df.to_csv(output_path, sep='\t', index=False)

    return final_df
