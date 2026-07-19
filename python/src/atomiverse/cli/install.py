"""Install orchestration: scikit-build via ``pip install`` only.

Never shells out to a hand-written ``cmake -B`` path. The sole compile
entry is::

    pip install --force-reinstall --no-deps \\
        -C wheel.cmake=true \\
        -C cmake.define.ATV_BUILD_PYTHON=ON \\
        -C cmake.define.<KEY>=<VALUE> … \\
        <source>
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Mapping, Optional

from atomiverse.cli.config import (
    ConfigError,
    InstallConfig,
    define_fingerprint,
    enabled_modules_from_defines,
    resolve_defines,
)
from atomiverse.cli.source import SourceSpec


def state_path() -> Path:
    """Path of the install-state JSON next to the ``atomiverse`` package."""
    import atomiverse

    return Path(atomiverse.__file__).resolve().parent / "_install_state.json"


def write_install_state(defines: Mapping[str, str]) -> Path:
    """Persist enabled modules + define fingerprint after a successful install."""
    path = state_path()
    payload = {
        "modules": enabled_modules_from_defines(defines),
        "defines": dict(sorted(defines.items())),
        "fingerprint": define_fingerprint(defines),
    }
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return path


def build_pip_command(
    source: SourceSpec,
    defines: Mapping[str, str],
    *,
    python: Optional[str] = None,
) -> List[str]:
    """Build the exact ``pip install`` argv for a native rebuild.

    Always includes ``ATV_BUILD_PYTHON=ON`` and re-enables
    ``wheel.cmake`` (pure-wheel default is false).

    @param source: Resolved source location.
    @param defines: CMake define map (must already contain ATV_BUILD_PYTHON).
    @param python: Interpreter for ``-m pip`` (default: ``sys.executable``).
    @return Full argv list.
    """
    if defines.get("ATV_BUILD_PYTHON", "").upper() not in ("ON", "TRUE", "1"):
        raise ConfigError("install requires ATV_BUILD_PYTHON=ON in cmake.define")

    exe = python or sys.executable
    cmd: List[str] = [
        exe,
        "-m",
        "pip",
        "install",
        "--force-reinstall",
        "--no-deps",
        "-C",
        "wheel.cmake=true",
    ]
    # Stable order for reproducibility / tests.
    for key in sorted(defines):
        cmd.extend(["-C", f"cmake.define.{key}={defines[key]}"])
    cmd.append(source.location)
    return cmd


def format_defines(defines: Mapping[str, str]) -> str:
    """Pretty-print defines one per line."""
    return "\n".join(f"  {k}={defines[k]}" for k in sorted(defines))


def run_install(
    cfg: InstallConfig,
    source: SourceSpec,
    *,
    dry_run: bool = False,
    python: Optional[str] = None,
) -> int:
    """Resolve defines and invoke pip (or print them on dry-run).

    @param cfg: Parsed install configuration.
    @param source: Source tree / tarball.
    @param dry_run: When True, print defines + pip command and exit 0
        without calling pip or writing state.
    @param python: Interpreter override for pip.
    @return Process exit code.
    """
    defines = resolve_defines(cfg)
    # Belt-and-suspenders: always force PYTHON on.
    defines["ATV_BUILD_PYTHON"] = "ON"

    cmd = build_pip_command(source, defines, python=python)

    print("atomiverse install")
    print(f"  version: {source.version}")
    print(f"  source:  {source.location}")
    print("  modules: " + ", ".join(enabled_modules_from_defines(defines)))
    print("  cmake.define:")
    print(format_defines(defines))
    print("  fingerprint: " + define_fingerprint(defines))

    if dry_run:
        print("  dry-run: pip not invoked")
        print("  pip command:")
        print("    " + " ".join(cmd))
        return 0

    print("  running pip …")
    proc = subprocess.run(cmd, check=False)
    if proc.returncode != 0:
        print(
            f"pip install failed with exit code {proc.returncode}",
            file=sys.stderr,
        )
        return proc.returncode

    try:
        write_install_state(defines)
    except OSError as exc:  # pragma: no cover - filesystem edge
        print(f"warning: could not write install state: {exc}", file=sys.stderr)

    print("install complete")
    return 0


def dry_run_defines(cfg: InstallConfig) -> Dict[str, str]:
    """Resolve defines for tests / regression without invoking pip."""
    defines = resolve_defines(cfg)
    defines["ATV_BUILD_PYTHON"] = "ON"
    return defines
