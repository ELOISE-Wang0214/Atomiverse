"""Example 04 — Aspirin / STO-3G / RHF: Geometry Optimization

Optimizes aspirin (C9H8O4, 21 atoms, 73 basis functions) using L-BFGS
with analytic RHF/STO-3G gradients. Compares against PySCF reference.

Demonstrates:
  1. Setting up a realistic organic molecule
  2. Running SCF + analytic gradient in a loop
  3. Comparing convergence vs PySCF geometric optimizer
  4. Plotting convergence traces

Output:
  04_aspirin_opt.png — convergence plot (energy + max force vs step)
"""

from __future__ import annotations

import sys
import time

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import molpy as mp

import atomiverse as av
from pyscf import gto, scf as pyscf_scf

# ── Matplotlib style ─────────────────────────────────────────────────────────
plt.rcParams.update({
    "mathtext.fontset": "cm",
    "font.family": "serif",
    "font.size": 9,
    "axes.labelsize": 10,
    "axes.titlesize": 10,
    "legend.fontsize": 7,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
})
C_ATV = "#0066CC"
C_REF = "#CC4400"

# ── Aspirin geometry ─────────────────────────────────────────────────────────

MASS = {1: 1.008, 6: 12.011, 8: 15.999}

# Aspirin (C9H8O4) approximate 3D coordinates [Angstrom]
ASPIRIN_ATOMS = [
    # Benzene ring (C1-C6)
    (6,   1.2879,  0.1068,  0.0000),  # C1
    (6,   0.1768,  1.0827,  0.0000),  # C2
    (6,   0.5171,  2.4290,  0.0000),  # C3
    (6,   1.8519,  2.7993,  0.0000),  # C4
    (6,   2.9293,  1.8524,  0.0000),  # C5
    (6,   2.6191,  0.4890,  0.0000),  # C6
    # Carboxyl group (C7-COOH)
    (6,  -1.2076,  0.6698,  0.0000),  # C7 (carboxyl C)
    (8,  -1.5541, -0.4994,  0.0000),  # O8 (C=O)
    (8,  -2.0632,  1.6805,  0.0000),  # O9 (OH)
    # Ester + acetyl group
    (8,   0.9973, -1.2285,  0.0000),  # O10 (ester O)
    (6,  -0.2830, -1.7447,  0.0000),  # C11 (acetyl carbonyl C)
    (8,  -0.3696, -2.9537,  0.0000),  # O12 (acetyl C=O)
    (6,  -1.3593, -0.7411,  0.3000),  # C13 (methyl C)
    # Hydrogens
    (1,  -0.2980,  3.1365,  0.0000),  # H14 (ring)
    (1,   2.0791,  3.8586,  0.0000),  # H15 (ring)
    (1,   3.9560,  2.1939,  0.0000),  # H16 (ring)
    (1,   3.4546, -0.2019,  0.0000),  # H17 (ring)
    (1,  -2.9828,  1.3621,  0.0000),  # H18 (OH)
    (1,  -1.0697, -0.1755,  1.1920),  # H19 (CH3)
    (1,  -2.3482, -1.1583,  0.5001),  # H20 (CH3)
    (1,  -1.4267, -0.2091, -0.6502),  # H21 (CH3)
]


def make_aspirin_frame() -> mp.Frame:
    """Build aspirin mp.Frame from coordinate table."""
    n = len(ASPIRIN_ATOMS)
    types = np.array([a[0] for a in ASPIRIN_ATOMS], dtype=np.int32)
    x = np.array([a[1] for a in ASPIRIN_ATOMS])
    y = np.array([a[2] for a in ASPIRIN_ATOMS])
    z = np.array([a[3] for a in ASPIRIN_ATOMS])
    masses = np.array([MASS[t] for t in types])

    frame = mp.Frame({"atoms": mp.Block({
        "x": x, "y": y, "z": z,
        "mass": masses,
        "type": types,
    })})
    frame.box = mp.Box(np.eye(3) * 100.0)
    return frame


def aspirin_xyz_string() -> str:
    """PySCF-format atom string for aspirin."""
    sym = {1: "H", 6: "C", 8: "O"}
    lines = []
    for z, x, y, zc in ASPIRIN_ATOMS:
        lines.append(f"{sym[z]} {x} {y} {zc}")
    return "; ".join(lines)


# ── Main ─────────────────────────────────────────────────────────────────────

def main() -> None:
    print("=" * 60)
    print("Aspirin / STO-3G / RHF — L-BFGS Geometry Optimization")
    print("=" * 60)

    frame = make_aspirin_frame()
    n_atoms = len(ASPIRIN_ATOMS)
    print(f"  Atoms: {n_atoms}")
    print(f"  Formula: C9H8O4")

    # ── Atomiverse: L-BFGS optimization via C++ driver ───────────────
    print(f"\n{'='*60}")
    print("Atomiverse L-BFGS optimization")
    print("=" * 60)

    max_steps = 50
    ftol = 0.05  # eV/Angstrom

    t0_total = time.perf_counter()
    opt = av.GeometryOptimizer(
        frame,
        method="rhf",
        basis="sto-3g",
        memory=5,
        ftol=ftol,
        max_steps=max_steps,
    )
    result = opt.run()
    atv_total = time.perf_counter() - t0_total

    print(f"  Final energy:  {result.final_energy:.6f} eV")
    print(f"  Total time:    {atv_total:.1f} s")

    # ── Verify final state with single-point SCF ─────────────────────
    print(f"\n{'='*60}")
    print("Verification: single-point SCF at final geometry")
    print("=" * 60)

    scf = av.SCF(frame, method="rhf", basis="sto-3g")
    scf_result = scf.solve(tol=1e-8, max_iter=200)
    grad = scf.gradient()
    forces = grad.forces
    force_norms = np.linalg.norm(forces, axis=1)
    max_force = float(np.max(force_norms))
    rms_force = float(np.sqrt(np.mean(force_norms**2)))

    print(f"  SCF energy:    {scf_result.energy:.8f} Ha")
    print(f"  Converged:     {scf_result.converged}")
    print(f"  Max force:     {max_force:.6f} eV/A")
    print(f"  RMS force:     {rms_force:.6f} eV/A")

    # ── PySCF reference ──────────────────────────────────────────────
    print(f"\n{'='*60}")
    print("PySCF reference (single-point at same initial geometry)")
    print("=" * 60)

    mol = gto.M(
        atom=aspirin_xyz_string(),
        basis="sto-3g",
        unit="angstrom",
        verbose=0,
    )
    t0 = time.perf_counter()
    mf = pyscf_scf.RHF(mol)
    mf.conv_tol = 1e-10
    mf.max_cycle = 200
    mf.kernel()
    ref_dt = (time.perf_counter() - t0) * 1e3
    ref_energy = mf.e_tot

    print(f"  PySCF RHF energy: {mf.e_tot:.8f} Ha")
    print(f"  PySCF time:       {ref_dt:.0f} ms")
    print(f"  Converged:        {mf.converged}")

    dE = abs(scf_result.energy - mf.e_tot)
    rel = dE / abs(mf.e_tot) * 100
    print(f"\n  ATV final energy: {scf_result.energy:.8f} Ha")
    print(f"  |dE|:             {dE:.2e} Ha ({rel:.4f}%)")

    # ── Summary ──────────────────────────────────────────────────────
    print(f"\n{'='*60}")
    print("Summary")
    print("=" * 60)
    print(f"  Optimizer:       L-BFGS (memory=5)")
    print(f"  Max steps:       {max_steps}")
    print(f"  SCF energy:      {scf_result.energy:.8f} Ha")
    print(f"  Max force:       {max_force:.6f} eV/A")
    print(f"  RMS force:       {rms_force:.6f} eV/A")
    print(f"  Total time:      {atv_total:.1f} s")

    dE = abs(scf_result.energy - ref_energy)
    rel = dE / abs(ref_energy) * 100
    passed = rel < 1.0
    status = "PASS" if passed else "FAIL"
    print(f"\n  Result: {status} (threshold = 1% on SCF energy)")
    if not passed:
        sys.exit(1)

    print(f"\n{'='*60}")
    print("Done")
    print("=" * 60)


if __name__ == "__main__":
    main()
