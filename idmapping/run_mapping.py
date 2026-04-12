import argparse
from mapper import create_mapping

def main():
    parser = argparse.ArgumentParser(description="Map GFF3 IDs to Human/Mouse homologs via BLAST.")
    parser.add_argument("--gff", required=True, help="Path to GFF3 file")
    parser.add_argument("--human", required=True, help="Path to human BLAST results")
    parser.add_argument("--mouse", required=True, help="Path to mouse BLAST results")
    parser.add_argument("--output", required=True, help="Output TSV path")
    parser.add_argument("--evalue", type=float, default=1e-50, help="E-value threshold (default: 1e-50)")
    parser.add_argument("--include-top", action="store_true", help="Include the top hit if no hits meet the evalue threshold")

    args = parser.parse_args()

    print(f"Running mapping: E-value={args.evalue}, IncludeTop={args.include_top}")
    create_mapping(
        gff3_path=args.gff,
        human_path=args.human,
        mouse_path=args.mouse,
        output_path=args.output,
        evalue_threshold=args.evalue,
        include_top=args.include_top
    )
    print("Done.")

if __name__ == "__main__":
    main()
