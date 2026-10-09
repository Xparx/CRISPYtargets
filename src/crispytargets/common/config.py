"""YAML configuration loading for tools that accept a --config file."""

import logging

import yaml


def load_config(config_path):
    """Read a YAML config file, returning {} on any failure.

    Logging is used rather than print: by the time this runs, the CLI has
    already called common.setup_logging. Failures fall back to an empty
    config so the caller's argparse defaults still apply.
    """
    try:
        with open(config_path, 'r') as f:
            return yaml.safe_load(f) or {}
    except Exception as e:
        logging.error(f"Failed to load config file {config_path}: {e}")
        return {}