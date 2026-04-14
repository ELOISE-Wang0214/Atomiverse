"""Single-point energy timing: Atomiverse vs PySCF (aspirin / STO-3G / RHF)."""

import time
import numpy as np
import molpy as mp
import atomiverse as av
from pyscf import gto, scf as pyscf_scf

ASPIRIN_ATOMS = [
    (6,   1.2879,  0.1068,  0.0000),
    (6,   0.1768,  1.0827,  0.0000),
    (6,   0.5171,  2.4290,  0.0000),
    (6,   1.8519,  2.7993,  0.0000),
    (6,   2.9293,  1.8524,  0.0000),
    (6,   2.6191,  0.4890,  0.0000),
    (6,  -1.2076,  0.6698,  0.0000),
    (8,  -1.5541, -0.4994,  0.0000),
    (8,  -2.0632,  1.6805,  0.0000),
    (8,   0.9973, -1.2285,  0.0000),
    (6,  -0.2830, -1.7447,  0.0000),
    (8,  -0.3696, -2.9537,  0.0000),
    (6,  -1.3593, -0.7411,  0.3000),
    (1,  -0.2980,  3.1365,  0.0000),
    (1,   2.0791,  3.8586,  0.0000),
    (1,   3.9560,  2.1939,  0.0000),
    (1,   3.4546, -0.2019,  0.0000),
    (1,  -2.9828,  1.3621,  0.0000),
    (1,  -1.0697, -0.1755,  1.1920),
    (1,  -2.3482, -1.1583,  0.5001),
    (1,  -1.4267, -0.2091, -0.6502),
]
MASS = {1: 1.008, 6: 12.011, 8: 15.999}
SYMBOL = {1: "H", 6: "C", 8: "O"}

types  = np.array([a[0] for a in ASPIRIN_ATOMS], dtype=np.int32)
coords = np.array([[a[1], a[2], a[3]] for a in ASPIRIN_ATOMS], dtype=float)


def run_atomiverse() -> tuple[float, float]:
    masses = np.array([MASS[int(z)] for z in types], dtype=float)
    frame = mp.Frame({"atoms": mp.Block({
        "x": coords[:, 0].copy(), "y": coords[:, 1].copy(), "z": coords[:, 2].copy(),
        "mass": masses, "type": types.copy(),
    })})
    frame.box = mp.Box(np.eye(3) * 100.0)
    t0 = time.perf_counter()
    scf = av.SCF(frame, method="rhf", basis="sto-3g")
    result = scf.solve(tol=1e-8, max_iter=200)
    elapsed = time.perf_counter() - t0
    return float(result.energy), elapsed


def run_pyscf() -> tuple[float, float]:
    geom = [(SYMBOL[int(z)], (float(x), float(y), float(zc)))
            for z, (x, y, zc) in zip(types, coords)]
    mol = gto.M(atom=geom, basis="sto-3g", unit="angstrom", verbose=0)
    mf = pyscf_scf.RHF(mol)
    mf.conv_tol = 1e-8
    mf.max_cycle = 200
    t0 = time.perf_counter()
    energy = float(mf.kernel())
    elapsed = time.perf_counter() - t0
    return energy, elapsed


HARTREE_TO_EV = 27.211386245988

print("Aspirin / STO-3G / RHF  —  single-point energy benchmark")
print("=" * 56)

print("Running PySCF ...")
e_pyscf, t_pyscf = run_pyscf()
print(f"  Energy : {e_pyscf:.10f} Ha  ({e_pyscf * HARTREE_TO_EV:.6f} eV)")
print(f"  Time   : {t_pyscf:.3f} s")

print("Running Atomiverse ...")
e_atv, t_atv = run_atomiverse()
print(f"  Energy : {e_atv / HARTREE_TO_EV:.10f} Ha  ({e_atv:.6f} eV)")
print(f"  Time   : {t_atv:.3f} s")

print("=" * 56)
print(f"Speedup (Atomiverse / PySCF): {t_atv / t_pyscf:.2f}x")
print(f"|dE|: {abs(e_atv / HARTREE_TO_EV - e_pyscf):.3e} Ha")
