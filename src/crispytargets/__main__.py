"""Allow ``python -m crispytargets <command> ...`` without installing."""

import sys

from .cli import main

if __name__ == "__main__":
    sys.exit(main())
