"""Generator for PySCF UKS-LSDA reference energies.

Produces the `pyscf_ref` constexpr block in test_uks_solve.cpp.
Regenerate any time the test battery changes.

Requirements: PySCF >= 2.13 installed via `pip install pyscf`.

Usage:
    python tests/cpu/scf/pyscf_reference.py

Reference configuration (matches Atomiverse's UKS solver defaults):
    dft.UKS(mol), xc="lda,vwn", atom_grid=(50, 110), prune=None,
    conv_tol=1e-12, basis="sto-3g"
"""
from __future__ import annotations

from pyscf import gto, dft

CONV_TOL = 1e-12
ATOM_GRID = (50, 110)


def run_uks(atom: str, basis: str, spin: int) -> float:
    mol = gto.M(atom=atom, basis=basis, spin=spin, charge=0)
    mf = dft.UKS(mol)
    mf.xc = "lda,vwn"  # Slater exchange + VWN5 correlation
    mf.grids.atom_grid = {s: ATOM_GRID for s in set(a[0] for a in mol.atom)}
    mf.grids.prune = None
    mf.conv_tol = CONV_TOL
    mf.max_cycle = 200
    return float(mf.kernel())


def main() -> None:
    cases = [
        ("H atom  (doublet)", "H 0 0 0",  "sto-3g", 1),
        ("He atom (singlet)", "He 0 0 0", "sto-3g", 0),
        ("Li atom (doublet)", "Li 0 0 0", "sto-3g", 1),
        ("H2 R=0.74 A (sng)", "H 0 0 0; H 0.74 0 0", "sto-3g", 0),
    ]
    print("# Single-system references")
    for label, atom, basis, spin in cases:
        E = run_uks(atom, basis, spin)
        print(f"# {label}: {E:.16e} Ha")

    print()
    print("# H2 PES scan (singlet, sto-3g, UKS lda,vwn)")
    R_grid = [0.5, 0.6, 0.74, 0.85, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5]
    for R in R_grid:
        atom = f"H 0 0 0; H {R} 0 0"
        E = run_uks(atom, "sto-3g", 0)
        print(f"# R={R:.2f} A : {E:.16e} Ha")


if __name__ == "__main__":
    main()
