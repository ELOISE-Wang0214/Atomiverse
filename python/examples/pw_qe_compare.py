#!/usr/bin/env python3
"""Compare Atomiverse plane-wave DFT against Quantum ESPRESSO for Si diamond.

Requires:
    - atomiverse (with PW enabled): ``pip install -e .``
    - Quantum ESPRESSO ``pw.x`` on PATH
    - Si pseudopotential (Si.upf) in pseudopotentials/ directory

Generates comparison plots:
    1. Total energy vs plane-wave cutoff (ecut)
    2. Total energy vs k-point mesh
    3. Total energy vs lattice constant (EOS)
    4. SCF convergence history (for selected parameters)

Usage::

    python examples/pw_qe_compare.py
    python examples/pw_qe_compare.py --skip-qe       # Atomiverse only
    python examples/pw_qe_compare.py --skip-atomiverse  # QE only
"""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

import molrs
import numpy as np

# ── Paths ──

SCRIPT_DIR = Path(__file__).resolve().parent
PSEUDO_DIR = SCRIPT_DIR / "pseudopotentials"
QE_BIN = "/tmp/qe-install/bin/pw.x"

# ── Si diamond structure builder ──


def build_si_diamond(a_lat: float = 5.43) -> "molrs.Frame":
    """Build the Si diamond primitive FCC cell (2 atoms) as a molrs.Frame."""
    a4 = a_lat / 4.0
    frame = molrs.Frame()
    atoms = molrs.Block()
    atoms.insert("x", np.array([0.0, a4]))
    atoms.insert("y", np.array([0.0, a4]))
    atoms.insert("z", np.array([0.0, a4]))
    atoms.insert("type", np.array([14, 14], dtype=np.int32))
    atoms.insert("mass", np.full(2, 28.0855))
    frame["atoms"] = atoms

    h = a_lat / 2.0
    frame.box = molrs.Box(
        np.array(
            [
                [0.0, h, h],
                [h, 0.0, h],
                [h, h, 0.0],
            ]
        )
    )
    return frame


# ── QE helpers ──

QE_INPUT_TEMPLATE = """&control
    calculation  = 'scf'
    prefix       = 'si'
    pseudo_dir   = '{pseudo_dir}'
    outdir       = '{outdir}'
    verbosity    = 'low'
/
&system
    ibrav        = 2
    celldm(1)    = {alat_bohr}
    nat          = 2
    ntyp         = 1
    ecutwfc      = {ecut_ry}
    occupations  = 'smearing'
    smearing     = 'gaussian'
    degauss      = 0.01
/
&electrons
    conv_thr     = {conv_thr}
    mixing_beta  = 0.4
/
ATOMIC_SPECIES
    Si   28.0855  Si.upf
ATOMIC_POSITIONS crystal
    Si   0.000000  0.000000  0.000000
    Si   0.250000  0.250000  0.250000
K_POINTS automatic
    {kmesh} {kmesh} {kmesh}  0 0 0
"""


def write_qe_input(
    path: str, *, a_lat: float, ecut_ry: float, kmesh: int, conv_thr: float = 1e-6
) -> str:
    """Write a QE input file and return the outdir path."""
    alat_bohr = a_lat / 0.52917721092  # Angstrom → Bohr
    outdir = str(Path(path).parent / "qe_out")

    content = QE_INPUT_TEMPLATE.format(
        pseudo_dir=PSEUDO_DIR,
        outdir=outdir,
        alat_bohr=f"{alat_bohr:.8f}",
        ecut_ry=f"{ecut_ry:.4f}",
        conv_thr=f"{conv_thr:.1e}",
        kmesh=kmesh,
    )
    Path(path).write_text(content)
    return outdir


def run_qe(input_path: str, outdir: str) -> dict | None:
    """Run pw.x and return parsed results, or None on failure."""
    try:
        result = subprocess.run(
            [QE_BIN, "-i", input_path],
            capture_output=True,
            text=True,
            timeout=300,
            env={**os.environ, "TMPDIR": outdir},
        )
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return None

    output = result.stdout + result.stderr

    # Parse total energy
    energy = None
    converged = False
    iterations = 0

    # Look for "!    total energy"
    for line in output.splitlines():
        m = re.search(r"!\s+total energy\s+=\s+([\d\.\-]+)\s+Ry", line)
        if m:
            energy = float(m.group(1))

    # Count SCF iterations
    for line in output.splitlines():
        if "estimated scf accuracy" in line:
            iterations += 1

    # Check convergence
    converged = "convergence has been achieved" in output

    if energy is None:
        return None

    return {
        "energy": energy,
        "converged": converged,
        "iterations": iterations,
    }


# ── Atomiverse helpers ──


def run_atomiverse(
    a_lat: float,
    ecut_ry: float,
    kmesh: int,
    tol: float = 1e-6,
    max_iter: int = 80,
    upf_path: str = "",
) -> dict:
    """Run Atomiverse pw_scf and return results."""
    import atomiverse as av

    frame = build_si_diamond(a_lat)
    r = av.pw_scf(frame, upf_path, ecut_ry, kmesh, tol, max_iter)
    return {
        "energy": r.energy,
        "converged": r.converged,
        "iterations": r.iterations,
    }


# ── Parameter sweep ──


def sweep_ecut(
    skip_qe: bool = False, skip_av: bool = False, upf: str = ""
) -> list[dict]:
    """Energy vs plane-wave cutoff."""
    results = []
    ecuts = [10, 15, 20, 25, 30, 35, 40, 50, 60]

    for ecut in ecuts:
        point = {"param": "ecut", "value": ecut}

        if not skip_av:
            av = run_atomiverse(5.43, ecut, 4, max_iter=80, upf_path=upf)
            point["atomiverse"] = av
            print(
                f"  AV  ecut={ecut:3.0f} Ry  E={av['energy']:15.6f}  iters={av['iterations']}"
            )

        if not skip_qe:
            with tempfile.TemporaryDirectory() as tmpdir:
                inp = os.path.join(tmpdir, "si.in")
                outdir = write_qe_input(inp, a_lat=5.43, ecut_ry=ecut, kmesh=4)
                os.makedirs(outdir, exist_ok=True)
                qe = run_qe(inp, outdir)
                point["qe"] = qe
                if qe:
                    print(
                        f"  QE  ecut={ecut:3.0f} Ry  E={qe['energy']:15.6f}  iters={qe['iterations']}"
                    )

        results.append(point)

    return results


def sweep_kmesh(
    skip_qe: bool = False, skip_av: bool = False, upf: str = ""
) -> list[dict]:
    """Energy vs k-point mesh."""
    results = []
    kmeshes = [1, 2, 3, 4, 5, 6, 8]

    for k in kmeshes:
        point = {"param": "kmesh", "value": k}

        if not skip_av:
            av = run_atomiverse(5.43, 30, k, max_iter=80, upf_path=upf)
            point["atomiverse"] = av
            print(f"  AV  kmesh={k}  E={av['energy']:15.6f}  iters={av['iterations']}")

        if not skip_qe:
            with tempfile.TemporaryDirectory() as tmpdir:
                inp = os.path.join(tmpdir, "si.in")
                outdir = write_qe_input(inp, a_lat=5.43, ecut_ry=30, kmesh=k)
                os.makedirs(outdir, exist_ok=True)
                qe = run_qe(inp, outdir)
                point["qe"] = qe
                if qe:
                    print(
                        f"  QE  kmesh={k}  E={qe['energy']:15.6f}  iters={qe['iterations']}"
                    )

        results.append(point)

    return results


def sweep_alat(
    skip_qe: bool = False, skip_av: bool = False, upf: str = ""
) -> list[dict]:
    """Energy vs lattice constant (EOS scan)."""
    results = []
    alats = [5.20, 5.30, 5.35, 5.40, 5.43, 5.45, 5.50, 5.55, 5.60, 5.70]

    for a in alats:
        point = {"param": "alat", "value": a}

        if not skip_av:
            av = run_atomiverse(a, 30, 4, max_iter=80, upf_path=upf)
            point["atomiverse"] = av
            print(
                f"  AV  alat={a:.2f} A  E={av['energy']:15.6f}  iters={av['iterations']}"
            )

        if not skip_qe:
            with tempfile.TemporaryDirectory() as tmpdir:
                inp = os.path.join(tmpdir, "si.in")
                outdir = write_qe_input(inp, a_lat=a, ecut_ry=30, kmesh=4)
                os.makedirs(outdir, exist_ok=True)
                qe = run_qe(inp, outdir)
                point["qe"] = qe
                if qe:
                    print(
                        f"  QE  alat={a:.2f} A  E={qe['energy']:15.6f}  iters={qe['iterations']}"
                    )

        results.append(point)

    return results


# ── Plotting ──


def plot_comparison(results: dict, save_path: str = ""):
    """Generate multi-panel comparison plot."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle("Atomiverse vs Quantum ESPRESSO — Si Diamond", fontsize=14)

    # Panel 1: Energy vs ecut
    ax = axes[0, 0]
    _plot_compare(
        ax,
        results.get("ecut", []),
        "ecut",
        "ecut [Ry]",
        "E [Ry]",
        "Total energy vs PW cutoff",
    )

    # Panel 2: Energy vs k-mesh
    ax = axes[0, 1]
    _plot_compare(
        ax,
        results.get("kmesh", []),
        "kmesh",
        "k-mesh",
        "E [Ry]",
        "Total energy vs k-point mesh",
    )

    # Panel 3: Energy vs alat (EOS)
    ax = axes[1, 0]
    _plot_compare(
        ax, results.get("alat", []), "alat", "a [Å]", "E [Ry]", "Equation of state"
    )

    # Panel 4: Energy difference
    ax = axes[1, 1]
    _plot_difference(
        ax, results.get("ecut", []), "ecut", "ecut [Ry]", "ΔE [eV] | E_AV − E_QE |"
    )

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"\nSaved: {save_path}")
    else:
        plt.show()


def _plot_compare(ax, data, param, xlabel, ylabel, title):
    """Plot AV and QE results for one parameter sweep."""
    if not data:
        ax.set_title(f"{title}\n(no data)")
        return

    x_vals = [p["value"] for p in data]
    av_vals = [p.get("atomiverse", {}).get("energy") for p in data]
    qe_vals = [p.get("qe", {}).get("energy") for p in data]

    ax.plot(x_vals, av_vals, "o-", label="Atomiverse", markersize=6)
    ax.plot(x_vals, qe_vals, "s--", label="Quantum ESPRESSO", markersize=6)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend()
    ax.grid(True, alpha=0.3)


def _plot_difference(ax, data, param, xlabel, ylabel):
    """Plot absolute energy difference."""
    if not data:
        ax.set_title("Energy difference (no data)")
        return

    x_vals = [p["value"] for p in data]
    diffs = []
    for p in data:
        av_e = p.get("atomiverse", {}).get("energy")
        qe_e = p.get("qe", {}).get("energy")
        if av_e is not None and qe_e is not None:
            diffs.append(abs(av_e - qe_e) * 13.6057)  # Ry → eV
        else:
            diffs.append(None)

    ax.semilogy(x_vals, diffs, "o-", markersize=6)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title("Energy difference")
    ax.grid(True, alpha=0.3)


# ── Main ──


def main():
    parser = argparse.ArgumentParser(
        description="Compare Atomiverse PW-DFT vs Quantum ESPRESSO for Si diamond"
    )
    parser.add_argument("--skip-qe", action="store_true", help="Skip QE calculations")
    parser.add_argument(
        "--skip-av",
        "--skip-atomiverse",
        action="store_true",
        dest="skip_av",
        help="Skip Atomiverse calculations",
    )
    parser.add_argument(
        "--sweep",
        choices=["ecut", "kmesh", "alat", "all"],
        default="all",
        help="Which sweep to run",
    )
    parser.add_argument("--upf", default="", help="Path to Si UPF pseudopotential")
    parser.add_argument("--output", default="", help="Save plot to file")
    args = parser.parse_args()

    # Auto-detect UPF
    upf = args.upf
    if not upf:
        default_upf = PSEUDO_DIR / "Si.upf"
        if default_upf.exists():
            upf = str(default_upf)
            print(f"Using pseudopotential: {upf}")
        else:
            print("No UPF found — using toy local potential (results not physical)")
            print(f"Place Si.upf in {PSEUDO_DIR}/ for meaningful results")

    if args.skip_qe:
        print("Skipping Quantum ESPRESSO calculations\n")
    else:
        if not Path(QE_BIN).exists():
            print(f"QE pw.x not found at {QE_BIN}")
            print("Install QE or use --skip-qe")
            sys.exit(1)
        print(f"Using QE at: {QE_BIN}\n")

    results = {}

    if args.sweep in ("ecut", "all"):
        print("=" * 60)
        print("Sweep: Energy vs plane-wave cutoff")
        print("=" * 60)
        results["ecut"] = sweep_ecut(args.skip_qe, args.skip_av, upf)

    if args.sweep in ("kmesh", "all"):
        print("\n" + "=" * 60)
        print("Sweep: Energy vs k-point mesh")
        print("=" * 60)
        results["kmesh"] = sweep_kmesh(args.skip_qe, args.skip_av, upf)

    if args.sweep in ("alat", "all"):
        print("\n" + "=" * 60)
        print("Sweep: Energy vs lattice constant (EOS)")
        print("=" * 60)
        results["alat"] = sweep_alat(args.skip_qe, args.skip_av, upf)

    # Save results
    out_file = args.output or "pw_qe_results.json"
    with open(out_file, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to: {out_file}")

    # Plot
    plot_path = (
        args.output.replace(".json", ".png") if args.output else "pw_qe_compare.png"
    )
    plot_comparison(results, save_path=plot_path)


if __name__ == "__main__":
    main()
