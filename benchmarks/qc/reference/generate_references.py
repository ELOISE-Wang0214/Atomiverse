#!/usr/bin/env python3
"""Generate PySCF reference data for the QC benchmark suite.

For each xyz file under ``benchmarks/qc/molecules/``, runs:
    1. RKS SCF at the target level of theory on the starting geometry
    2. geometric-based geometry optimization to convergence
    3. One more SCF on the optimized geometry to lock in the reference

Outputs ``<stem>.json`` under ``benchmarks/qc/reference/`` with the
starting geometry, optimized geometry, final energy, and the iteration
counts — numbers the C++ bench (``bench_qc_suite.cpp``) compares against.

Usage:
    python benchmarks/qc/reference/generate_references.py                 # all
    python benchmarks/qc/reference/generate_references.py h2o nh3         # subset

Ground truth is PySCF 2.12.x. The level-of-theory default matches
Phase 1.2 acceptance: PBE / def2-SVP / cart = True (matches atomiverse
Cartesian-d convention).
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

try:
    from pyscf import gto, scf, dft
    from pyscf.geomopt.geometric_solver import optimize as geom_optimize
except ImportError:
    sys.stderr.write("pyscf + geometric required to regenerate references\n")
    raise

REPO_ROOT = Path(__file__).resolve().parents[3]
MOL_DIR = REPO_ROOT / "benchmarks" / "qc" / "molecules"
REF_DIR = REPO_ROOT / "benchmarks" / "qc" / "reference"

Z_FROM_SYMBOL = {
    "H": 1, "He": 2, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7,
    "O": 8, "F": 9, "Ne": 10,
}


def read_xyz(path: Path) -> list[tuple[str, float, float, float]]:
    """Parse a standard XYZ file (atom count, comment, then one line per atom)."""
    text = path.read_text().splitlines()
    n = int(text[0].strip())
    atoms = []
    for line in text[2 : 2 + n]:
        parts = line.split()
        atoms.append((parts[0], float(parts[1]), float(parts[2]), float(parts[3])))
    return atoms


def atoms_to_pyscf_str(atoms):
    return "; ".join(f"{sym} {x} {y} {z}" for sym, x, y, z in atoms)


def pyscf_atoms(mol) -> list[tuple[str, float, float, float]]:
    """Extract optimized atom list from a PySCF Mole in Angstrom."""
    out = []
    coords = mol.atom_coords(unit="Angstrom")
    for i, sym in enumerate(mol.elements):
        x, y, z = coords[i]
        out.append((sym, float(x), float(y), float(z)))
    return out


def _make_mf(mol, method: str, xc: str, tol: float):
    if method.lower() == "rhf":
        mf = scf.RHF(mol)
    elif method.lower() == "rks":
        mf = dft.RKS(mol)
        mf.xc = xc
    else:
        raise ValueError(f"unknown method {method!r} (expected rhf / rks)")
    mf.conv_tol = tol
    return mf


def run_reference(stem: str, basis: str, method: str, xc: str,
                  scf_tol: float = 1e-10) -> dict:
    atoms = read_xyz(MOL_DIR / f"{stem}.xyz")
    mol = gto.M(
        atom=atoms_to_pyscf_str(atoms),
        basis=basis,
        cart=True,  # match atomiverse Cartesian-d convention
        verbose=0,
    )

    mf = _make_mf(mol, method, xc, scf_tol)
    e_start = mf.kernel()
    scf_iters_start = mf.cycles

    t0 = time.perf_counter()
    mol_opt = geom_optimize(mf, maxsteps=200, convergence_energy=1e-6,
                            convergence_grms=3e-4, convergence_gmax=4.5e-4)
    opt_wall = time.perf_counter() - t0

    mf_opt = _make_mf(mol_opt, method, xc, scf_tol)
    e_opt = mf_opt.kernel()

    return {
        "stem": stem,
        "basis": basis,
        "method": method,
        "xc": xc if method.lower() == "rks" else None,
        "cart": True,
        "starting_geometry_angstrom": [
            {"symbol": s, "x": x, "y": y, "z": z} for s, x, y, z in atoms
        ],
        "optimized_geometry_angstrom": [
            {"symbol": s, "x": x, "y": y, "z": z} for s, x, y, z in pyscf_atoms(mol_opt)
        ],
        "energy_hartree_starting": float(e_start),
        "energy_hartree_optimized": float(e_opt),
        "scf_iterations_starting": int(scf_iters_start),
        "opt_wall_seconds": float(opt_wall),
        "pyscf_version": __import__("pyscf").__version__,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("stems", nargs="*",
                        help="Molecule stems to regenerate (default: all xyz files)")
    parser.add_argument("--basis", default="sto-3g",
                        help="PySCF basis name (default: sto-3g)")
    parser.add_argument("--method", default="rhf", choices=["rhf", "rks"],
                        help="rhf or rks (default rhf — matches the atomiverse "
                             "gradient path that exists today)")
    parser.add_argument("--xc", default="pbe",
                        help="XC functional when method=rks (ignored otherwise)")
    args = parser.parse_args()

    if not args.stems:
        args.stems = sorted(p.stem for p in MOL_DIR.glob("*.xyz"))

    REF_DIR.mkdir(parents=True, exist_ok=True)
    label = args.method if args.method == "rhf" else f"{args.method}/{args.xc}"
    for stem in args.stems:
        print(f"[{stem}] generating {label}/{args.basis} reference...")
        t0 = time.perf_counter()
        ref = run_reference(stem, basis=args.basis, method=args.method, xc=args.xc)
        out = REF_DIR / f"{stem}.json"
        out.write_text(json.dumps(ref, indent=2) + "\n")
        print(f"[{stem}] wrote {out.relative_to(REPO_ROOT)}  "
              f"(E_opt={ref['energy_hartree_optimized']:.10f} Ha, "
              f"{time.perf_counter() - t0:.1f} s)")


if __name__ == "__main__":
    main()
