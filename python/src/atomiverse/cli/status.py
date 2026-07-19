"""``python -m atomiverse status`` — enabled modules + define fingerprint."""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any, Dict, Optional

from atomiverse.cli.install import state_path


def load_install_state(path: Optional[Path] = None) -> Optional[Dict[str, Any]]:
    """Load install-state JSON written by a previous successful install.

    @param path: Override state file path (tests).
    @return Parsed dict, or None if absent/unreadable.
    """
    p = path if path is not None else state_path()
    if not p.is_file():
        return None
    try:
        return json.loads(p.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None


def format_status(state: Optional[Dict[str, Any]]) -> str:
    """Render a human-readable status report.

    @param state: Install-state dict or None.
    @return Multi-line status text.
    """
    lines = ["atomiverse status"]
    if state is None:
        lines.append("  modules: (none — no native install recorded)")
        lines.append("  fingerprint: (none)")
        lines.append(
            "  hint: run `python -m atomiverse install --md` "
            "(or --config) to build the C++ extension"
        )
        return "\n".join(lines)

    modules = state.get("modules") or []
    if isinstance(modules, list) and modules:
        lines.append("  modules: " + ", ".join(modules))
    else:
        lines.append("  modules: (none)")

    fp = state.get("fingerprint") or ""
    lines.append(f"  fingerprint: {fp if fp else '(none)'}")

    defines = state.get("defines") or {}
    if isinstance(defines, dict) and defines:
        lines.append("  defines:")
        for key in sorted(defines):
            lines.append(f"    {key}={defines[key]}")
    return "\n".join(lines)


def run_status(*, path: Optional[Path] = None, stream=None) -> int:
    """Print status and return 0.

    @param path: Optional state-file override.
    @param stream: Output stream (default stdout).
    @return Always 0.
    """
    out = stream if stream is not None else sys.stdout
    state = load_install_state(path)
    print(format_status(state), file=out)
    return 0
