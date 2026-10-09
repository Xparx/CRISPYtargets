"""Test configuration.

crispytargets lives under src/, so tests import it as an installed-style
package. Adding src/ to sys.path keeps the suite runnable without
installing the package first.
"""

import sys
from pathlib import Path

SRC = Path(__file__).resolve().parent.parent / "src"

if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))
