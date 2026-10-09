"""External tool invocation.

Two shapes, because the tools genuinely differ:

- :func:`run_command` for shell strings (the FlashFry java invocations the
  design pipeline builds as one string). Exit on failure, optional silence.
- :func:`run` for argv lists (tabix, samtools, twoBitToFa in liftover),
  capturing stdout/stderr and logging both on failure.

Both raise or exit loudly rather than returning a silently wrong result.
"""

import logging
import subprocess
import sys


def run_command(cmd, verbose=False):
    """Run a shell command, failing loudly.

    Args:
        cmd: shell string to execute.
        verbose: when False, discard the child's stdout/stderr.
    """
    logging.info(f"Executing: {cmd}")
    if verbose:
        subprocess.run(cmd, shell=True, check=True)
    else:
        subprocess.run(cmd, shell=True, check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def run(cmd, desc="", check=True):
    """Run an argv list, logging and optionally exiting on failure.

    Args:
        cmd: argv list, no shell interpolation.
        desc: human-readable description logged at INFO.
        check: when True, exit(1) on failure; when False, return the
            CompletedProcess (or CalledProcessError) for the caller.
    """
    logging.info(desc)
    logging.debug("  cmd: %s", " ".join(cmd))
    try:
        return subprocess.run(cmd, capture_output=True, text=True, check=check)
    except subprocess.CalledProcessError as e:
        logging.error("Command failed: %s", " ".join(cmd))
        logging.error("STDOUT: %s", e.stdout)
        logging.error("STDERR: %s", e.stderr)
        if check:
            sys.exit(1)
        return e