"""Logging setup shared by every CLI tool.

Console at INFO, file at DEBUG, matching the layout both tools always had.
The root logger is configured, and handlers are only added once, so calling
this twice (or from a library context) does not duplicate output lines.
"""

import logging


def setup_logging(log_name="crispytargets.log"):
    """Configure the root logger with a console and a file handler."""
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)

    if getattr(logger, "_crispytargets_configured", False):
        return logger

    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(logging.Formatter('%(levelname)s: %(message)s'))

    file_handler = logging.FileHandler(log_name)
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(
        logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))

    logger.addHandler(console_handler)
    logger.addHandler(file_handler)
    logger._crispytargets_configured = True
    return logger