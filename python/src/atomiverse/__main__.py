"""``python -m atomiverse`` entry point.

Phase 02 ships a stub: the full ``install`` / ``status`` CLI lands in
phase 03. The module must import and run without the C++ extension so a
pure-Python wheel can still guide the user toward native install.
"""

from __future__ import annotations

import sys
from typing import List, Optional

from atomiverse import __version__
from atomiverse._binding import extension_available


def main(argv: Optional[List[str]] = None) -> int:
    """Print install-oriented guidance and exit 0.

    @param argv: Optional CLI args (unused in the phase-02 stub).
    @return Process exit code (always 0 for the stub).
    """
    _ = argv  # full argparse lands in phase 03
    print(f"atomiverse {__version__}")
    if extension_available():
        print("Native C++ extension (_atomiverse): available")
        print("Full CLI (install / status) arrives in a later release.")
    else:
        print("Native C++ extension (_atomiverse): not installed")
        print()
        print("This is a pure-Python install. To build the compute backend:")
        print("  python -m atomiverse install")
        print()
        print("Until install is available, build from a source checkout:")
        print("  cmake -B build -DATV_BUILD_PYTHON=ON && cmake --build build")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
