import argparse
import logging
import yaml
from mapper import create_mapping

def load_config(config_path):
    try:
        with open(config_path, 'r') as f:
            return yaml.safe_load(f) or {}
    except Exception as e:
        logging.error(f"Failed to load config file {config_path}: {e}")
        return {}

def setup_logging(log_name="idmapping.log"):
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

def main():
    # Pre-parse to check for config file
    pre_parser = argparse.ArgumentParser(add_help=False)
    pre_parser.add_argument("--config", type=str, help="Path to YAML config file")
    pre_args, remaining_args = pre_parser.parse_known_args()

    config = load_config(pre_args.config) if pre_args.config else {}

    parser = argparse.ArgumentParser(description="Incremental Orthology Mapping Engine")
    parser.add_argument("--config", type=str, help="Path to YAML config file")
    parser.add_argument("--gff", default=config.get('gff'), required=not config.get('gff'), help="Path to GFF3 file")
    parser.add_argument("--blast-file", required=True, help="Path to the BLAST results TSV for the species")
    parser.add_argument("--species", required=True, help="Label for the species (e.g., 'Human', 'Mouse')")
    parser.add_argument("--output", default=config.get('output'), required=not config.get('output'), help="Output TSV path (will be updated incrementally)")
    parser.add_argument("--evalue", type=float, default=config.get('evalue', 1e-50), help="E-value threshold (default: 1e-50)")
    parser.add_argument("--include-top", action="store_true", default=config.get('include_top', False), help="Include the top hit if no hits meet the evalue threshold")

    args = parser.parse_args()

    setup_logging()
    logging.info(f"Running incremental mapping for {args.species}: E-value={args.evalue}")
    create_mapping(
        gff3_path=args.gff,
        blast_path=args.blast_file,
        species_label=args.species,
        output_path=args.output,
        evalue_threshold=args.evalue,
        include_top=args.include_top
    )
    logging.info("Done.")

if __name__ == "__main__":
    main()
