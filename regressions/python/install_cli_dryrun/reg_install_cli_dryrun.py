#!/usr/bin/env python3
"""Public CLI dry-run regression (python-on-demand-install-03).

Out-of-library consumer: pure sources only, no C++ compile, no real pip.
Exit non-zero on failure.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
SRC = REPO / "python" / "src"
PYTHON_PKG = REPO / "python"


def _run(args: list[str]) -> subprocess.CompletedProcess[str]:
    import os

    env = os.environ.copy()
    # Prefer pure sources; keep any existing PYTHONPATH after.
    prev = env.get("PYTHONPATH", "")
    env["PYTHONPATH"] = str(SRC) + (os.pathsep + prev if prev else "")
    return subprocess.run(
        [sys.executable, "-m", "atomiverse", *args],
        capture_output=True,
        text=True,
        check=False,
        cwd=str(REPO),
        env=env,
    )


def main() -> int:
    if not (SRC / "atomiverse").is_dir():
        print(f"FAIL: pure sources not found at {SRC / 'atomiverse'}")
        return 1

    print("== install_cli_dryrun ==")
    print(f"  python: {sys.executable}")
    print(f"  src:    {SRC}")

    # 1) help / bare
    proc = _run([])
    if proc.returncode != 0:
        print("FAIL: bare invocation non-zero")
        print(proc.stdout, proc.stderr)
        return 1
    if "install" not in (proc.stdout + proc.stderr).lower():
        print("FAIL: bare output missing install guidance")
        print(proc.stdout)
        return 1
    print("ok: bare invocation mentions install")

    proc = _run(["--help"])
    if proc.returncode != 0:
        print("FAIL: --help non-zero")
        return 1
    print("ok: --help")

    # 2) mutual exclusion
    proc = _run(["install", "--config", "x.toml", "--md"])
    if proc.returncode == 0:
        print("FAIL: expected mutual exclusion failure")
        return 1
    blob = (proc.stdout + proc.stderr).lower()
    if "mutual" not in blob and "exclusion" not in blob:
        print("FAIL: mutual exclusion message missing")
        print(proc.stderr)
        return 1
    print("ok: mutual exclusion")

    # 3) empty install
    proc = _run(["install"])
    if proc.returncode == 0:
        print("FAIL: empty install should fail")
        return 1
    if "empty" not in (proc.stdout + proc.stderr).lower():
        print("FAIL: empty install message missing")
        print(proc.stderr)
        return 1
    print("ok: empty install rejected")

    # 4) --hf reject
    proc = _run(["install", "--hf"])
    if proc.returncode == 0:
        print("FAIL: --hf should be rejected")
        return 1
    blob = (proc.stdout + proc.stderr).lower()
    if "hf" not in blob and "scf" not in blob:
        print("FAIL: hf reject message missing")
        print(proc.stderr)
        return 1
    print("ok: --hf rejected")

    # 5) dry-run --md (closure + PYTHON=ON)
    proc = _run(
        [
            "install",
            "--md",
            "--src",
            str(PYTHON_PKG),
            "--dry-run",
        ]
    )
    out = proc.stdout + proc.stderr
    if proc.returncode != 0:
        print("FAIL: dry-run --md failed")
        print(out)
        return 1
    for needle in (
        "ATV_BUILD_MD=ON",
        "ATV_BUILD_POTENTIALS=ON",
        "ATV_BUILD_NEIGHBOR=ON",
        "ATV_BUILD_PYTHON=ON",
    ):
        if needle not in out:
            print(f"FAIL: dry-run output missing {needle}")
            print(out)
            return 1
    if "dry-run" not in out.lower():
        print("FAIL: dry-run marker missing")
        return 1
    print("ok: dry-run --md defines (closure + PYTHON=ON)")

    # 6) status
    proc = _run(["status"])
    if proc.returncode != 0:
        print("FAIL: status non-zero")
        print(proc.stdout, proc.stderr)
        return 1
    st = (proc.stdout + proc.stderr).lower()
    if "module" not in st or "fingerprint" not in st:
        print("FAIL: status missing modules/fingerprint")
        print(proc.stdout)
        return 1
    print("ok: status")

    print("REGRESSION PASSED: python/install_cli_dryrun")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
