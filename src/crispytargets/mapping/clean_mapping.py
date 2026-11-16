import pandas as pd
import argparse
import sys

def clean_locus_id(id_str):
    """
    Collapses a Ciona isoform ID to a Gene/Locus ID.
    Example: KY.Chr1.1.v1.nonSL13-1 -> KY.Chr1.1
    """
    if pd.isna(id_str) or id_str == "":
        return id_str
    
    # Handle pipe-separated lists (common in 'secondary' columns)
    ids = str(id_str).split('|')
    cleaned_ids = []
    
    for cid in ids:
        parts = cid.split('.')
        if len(parts) >= 3:
            # Keep only Species, Chromosome, and Locus (first 3 parts)
            cleaned_ids.append(".".join(parts[:3]))
        else:
            cleaned_ids.append(cid)
            
    # Remove duplicates while preserving original order
    seen = set()
    unique_cleaned = [x for x in cleaned_ids if not (x in seen or seen.add(x))]
    
    return "|".join(unique_cleaned)

def drop_primary(secondary, primary):
    """Remove the primary value from a pipe-separated secondary cell.

    Collapsing isoform IDs to locus IDs can turn a list of alternative
    isoforms into the primary locus itself, so the secondary column would
    repeat the primary. Returns an empty string when nothing else remains.
    """
    if pd.isna(secondary) or pd.isna(primary):
        return secondary
    keep = [x for x in str(secondary).split('|') if x and x != primary]
    return "|".join(keep)


def main():
    parser = argparse.ArgumentParser(description="Clean Ciona isoform IDs in a mapping file.")
    parser.add_argument("-i", "--input", required=True, help="Path to the input mapping TSV file")
    parser.add_argument("-o", "--output", required=True, help="Path to save the cleaned TSV file")
    parser.add_argument("-c", "--columns", required=True, help="Comma-separated list of columns to clean (e.g., 'KY19_primary,KY19_secondary')")
    parser.add_argument("-s", "--sep", default='\t', help="Separator used in the file (default: tab)")

    args = parser.parse_args()

    try:
        # Load the mapping file
        df = pd.read_csv(args.input, sep=args.sep)
        cols_to_clean = args.columns.split(',')

        for col in cols_to_clean:
            if col in df.columns:
                print(f"Cleaning column: {col}...")
                df[col] = df[col].apply(clean_locus_id)
            else:
                print(f"Warning: Column '{col}' not found in file. Skipping.")

        for col in cols_to_clean:
            if col.endswith('_secondary') and col in df.columns:
                primary_col = col[:-len('_secondary')] + '_primary'
                if primary_col in df.columns:
                    df[col] = [drop_primary(s, p) for s, p in
                               zip(df[col], df[primary_col])]

        # Save the cleaned mapping
        df.to_csv(args.output, sep=args.sep, index=False)
        print(f"Successfully saved cleaned mapping to: {args.output}")

    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
