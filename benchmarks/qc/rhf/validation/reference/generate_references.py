#!/usr/bin/env python3
"""Generate PySCF reference data for the QC benchmark suite.

For each xyz file under ``benchmarks/qc/rhf/validation/molecules/``, runs:
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
except ImportError:
    sys.stderr.write("pyscf required to regenerate references\n")
    raise

# geometric is only needed by the geometry-optimization (non-hybrid) path;
# imported lazily there so the hybrid single-point mode runs without it.
geom_optimize = None

# Multiplicity overrides for the radical benchmark entries. Closed-shell
# molecules (spin=0) are handled by the default path; the keys here carry
# 2S for PySCF's `gto.M(spin=...)` convention.
OPEN_SHELL = {
    "o2_triplet":  {"spin": 2},  # triplet  -> S=1
    "no_doublet":  {"spin": 1},  # doublet  -> S=1/2
    "ch3_radical": {"spin": 1},  # doublet
}

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
    m = method.lower()
    open_shell = mol.spin > 0
    if m == "rhf":
        mf = scf.UHF(mol) if open_shell else scf.RHF(mol)
    elif m == "rks":
        mf = dft.UKS(mol) if open_shell else dft.RKS(mol)
        mf.xc = xc
    elif m == "uhf":
        mf = scf.UHF(mol)
    elif m == "uks":
        mf = dft.UKS(mol)
        mf.xc = xc
    else:
        raise ValueError(f"unknown method {method!r} (expected rhf/rks/uhf/uks)")
    mf.conv_tol = tol
    return mf


def run_hybrid_singlepoint(stem: str, xc: str, basis: str = "6-31g*",
                           scf_tol: float = 1e-10) -> dict:
    """Run a single-point hybrid RKS SCF and return an energy/orbital block.

    Runs ``pyscf.dft.RKS(mol, xc=...)`` at the requested global-hybrid
    functional (``b3lyp`` or ``pbe0``) on the STARTING geometry of
    ``<stem>.xyz``. Records the converged total energy and the HOMO /
    LUMO orbital energies. No geometry optimization is performed -- the
    hybrid single-point references pin the C++ KS SCF at a fixed
    geometry (spec phase3-hybrid-dft-02-hybrid-singlepoint).

    Args:
        stem: Molecule stem (xyz file name without extension).
        xc: Hybrid functional name passed to PySCF ('b3lyp' or 'pbe0').
        basis: PySCF basis name (default '6-31g*').
        scf_tol: PySCF SCF energy convergence tolerance.

    Returns:
        A dict with the total energy [Hartree], HOMO/LUMO orbital
        energies [Hartree], basis, functional and pyscf_version.
    """
    atoms = read_xyz(MOL_DIR / f"{stem}.xyz")
    mol = gto.M(
        atom=atoms_to_pyscf_str(atoms),
        basis=basis,
        cart=True,  # match atomiverse Cartesian-d convention
        verbose=0,
    )
    mf = dft.RKS(mol)
    mf.xc = xc
    mf.conv_tol = scf_tol
    e_tot = mf.kernel()

    n_occ = mol.nelectron // 2
    mo_e = mf.mo_energy
    homo = float(mo_e[n_occ - 1])
    lumo = float(mo_e[n_occ])

    return {
        "basis": basis,
        "xc": xc,
        "energy_hartree": float(e_tot),
        "homo_hartree": homo,
        "lumo_hartree": lumo,
        "pyscf_version": __import__("pyscf").__version__,
    }


def run_reference(stem: str, basis: str, method: str, xc: str,
                  scf_tol: float = 1e-10) -> dict:
    atoms = read_xyz(MOL_DIR / f"{stem}.xyz")
    extras = OPEN_SHELL.get(stem, {})
    mol = gto.M(
        atom=atoms_to_pyscf_str(atoms),
        basis=basis,
        cart=True,  # match atomiverse Cartesian-d convention
        verbose=0,
        **extras,
    )

    mf = _make_mf(mol, method, xc, scf_tol)
    e_start = mf.kernel()
    scf_iters_start = mf.cycles

    from pyscf.geomopt.geometric_solver import optimize as geom_optimize

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
        "xc": xc if method.lower() in {"rks", "uks"} else None,
        "spin": int(mol.spin),
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


def run_hybrid_mode(stems: list[str], basis: str) -> None:
    """Generate b3lyp + pbe0 single-point blocks and merge into the JSONs.

    For each stem, runs the hybrid single-point reference for both
    b3lyp and pbe0 and writes the results into ``<stem>.json`` under
    keyed slots ``b3lyp`` / ``pbe0``. The existing top-level PBE
    reference data is preserved untouched -- only the two hybrid keys
    are added or replaced.
    """
    REF_DIR.mkdir(parents=True, exist_ok=True)
    for stem in stems:
        out = REF_DIR / f"{stem}.json"
        if out.exists():
            doc = json.loads(out.read_text())
        else:
            doc = {"stem": stem}
        for xc in ("b3lyp", "pbe0"):
            print(f"[{stem}] generating {xc}/{basis} single-point reference...")
            t0 = time.perf_counter()
            block = run_hybrid_singlepoint(stem, xc=xc, basis=basis)
            doc[xc] = block
            print(f"[{stem}]   {xc}: E={block['energy_hartree']:.10f} Ha "
                  f"HOMO={block['homo_hartree']:.6f} LUMO={block['lumo_hartree']:.6f} "
                  f"({time.perf_counter() - t0:.1f} s)")
        out.write_text(json.dumps(doc, indent=2) + "\n")
        print(f"[{stem}] wrote {out.relative_to(REPO_ROOT)}")


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
    parser.add_argument("--hybrid", action="store_true",
                        help="Hybrid single-point mode: emit b3lyp/pbe0 blocks "
                             "into the per-molecule JSONs (default basis 6-31g*).")
    parser.add_argument("--thermo", action="store_true",
                        help="Thermochemistry mode: emit PySCF thermo() "
                             "references for the vibrational-analysis link "
                             "(spec phase4-vibrations-03). Delegates to "
                             "thermo/generate_thermo_references.py.")
    args = parser.parse_args()

    if args.thermo:
        # Harmonic-thermochemistry references run on fixed geometries
        # (matching the link-2 frequency parity test), so they live in a
        # dedicated generator rather than the xyz-iterating path above.
        from thermo import generate_thermo_references
        generate_thermo_references.main()
        return

    if not args.stems:
        args.stems = sorted(p.stem for p in MOL_DIR.glob("*.xyz"))

    if args.hybrid:
        basis = args.basis if args.basis != "sto-3g" else "6-31g*"
        run_hybrid_mode(args.stems, basis=basis)
        return

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
