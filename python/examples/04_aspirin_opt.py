"""Example 04 — Aspirin / STO-3G / RHF: Geometry-Optimization Benchmark

Optimize the same mildly distorted aspirin geometry with:
  1. Atomiverse L-BFGS + analytic RHF/STO-3G forces
  2. PySCF RHF + analytic gradients + SciPy L-BFGS-B

The script reports:
  - optimization wall time
  - iteration count
  - final RHF energy
  - final max/RMS force
  - final geometry agreement between the two optimizers
"""

from __future__ import annotations

import argparse
import time
from dataclasses import dataclass
import numpy as np
import molpy as mp

import atomiverse as av
from pyscf import gto, scf as pyscf_scf
from scipy import optimize


HARTREE_TO_EV = 27.211386245988
BOHR_TO_ANGSTROM = 0.529177210903
HARTREE_PER_BOHR_TO_EV_PER_ANG = HARTREE_TO_EV / BOHR_TO_ANGSTROM

MASS = {1: 1.008, 6: 12.011, 8: 15.999}
SYMBOL = {1: "H", 6: "C", 8: "O"}

# Aspirin (C9H8O4) approximate 3D coordinates [Angstrom]
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


@dataclass(frozen=True)
class OptimizationSummary:
    name: str
    elapsed_s: float
    iterations: int
    converged: bool
    final_energy_ha: float
    final_max_force_eva: float
    final_rms_force_eva: float
    final_coords: np.ndarray


def reference_types_coords() -> tuple[np.ndarray, np.ndarray]:
    types = np.array([atom[0] for atom in ASPIRIN_ATOMS], dtype=np.int32)
    coords = np.array([[atom[1], atom[2], atom[3]] for atom in ASPIRIN_ATOMS], dtype=float)
    return types, coords


def make_perturbed_coords() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    types, coords = reference_types_coords()
    n_atoms = coords.shape[0]

    perturb = np.zeros_like(coords)
    idx = np.arange(n_atoms, dtype=float)

    # Global gentle puckering and asymmetry.
    perturb[:, 0] += 0.015 * np.cos(0.7 * idx)
    perturb[:, 1] += 0.020 * np.sin(0.9 * idx)
    perturb[:, 2] += 0.040 * np.sin(0.6 * idx + 0.3)

    # Targeted local distortions around the ester / carboxyl region.
    perturb[7]  += np.array([-0.08,  0.10,  0.04])  # O8
    perturb[8]  += np.array([ 0.05, -0.12, -0.03])  # O9
    perturb[10] += np.array([ 0.07,  0.09,  0.02])  # C11
    perturb[12] += np.array([-0.10,  0.05,  0.18])  # C13
    perturb[17] += np.array([-0.14, -0.02,  0.05])  # H18
    perturb[18] += np.array([ 0.06, -0.08,  0.10])  # H19
    perturb[19] += np.array([-0.04,  0.07, -0.06])  # H20
    perturb[20] += np.array([ 0.03, -0.05,  0.08])  # H21

    return types, coords, coords + perturb


def make_frame(types: np.ndarray, coords: np.ndarray) -> mp.Frame:
    masses = np.array([MASS[int(z)] for z in types], dtype=float)
    frame = mp.Frame({"atoms": mp.Block({
        "x": coords[:, 0].copy(),
        "y": coords[:, 1].copy(),
        "z": coords[:, 2].copy(),
        "mass": masses,
        "type": types.copy(),
    })})
    frame.box = mp.Box(np.eye(3) * 100.0)
    return frame


def frame_coords(frame: mp.Frame) -> np.ndarray:
    return np.column_stack((
        np.asarray(frame.atoms.x, dtype=float),
        np.asarray(frame.atoms.y, dtype=float),
        np.asarray(frame.atoms.z, dtype=float),
    ))


def coords_to_pyscf_geom(types: np.ndarray, coords: np.ndarray) -> list[tuple[str, tuple[float, float, float]]]:
    return [
        (SYMBOL[int(z)], (float(x), float(y), float(zc)))
        for z, (x, y, zc) in zip(types, coords, strict=True)
    ]


def force_metrics(forces: np.ndarray) -> tuple[float, float]:
    norms = np.linalg.norm(forces, axis=1)
    return float(np.max(norms)), float(np.sqrt(np.mean(norms**2)))


def displacement_metrics(ref_coords: np.ndarray, coords: np.ndarray) -> tuple[float, float]:
    disp = np.linalg.norm(coords - ref_coords, axis=1)
    return float(np.max(disp)), float(np.sqrt(np.mean(disp**2)))


def compare_geometries(coords_a: np.ndarray, coords_b: np.ndarray) -> tuple[float, float]:
    delta = np.linalg.norm(coords_a - coords_b, axis=1)
    return float(np.max(delta)), float(np.sqrt(np.mean(delta**2)))


def evaluate_atomiverse(types: np.ndarray, coords: np.ndarray) -> tuple[float, np.ndarray]:
    frame = make_frame(types, coords)
    scf = av.SCF(frame, method="rhf", basis="sto-3g")
    result = scf.solve(tol=1e-8, max_iter=200)
    grad = scf.gradient()
    return float(result.energy), np.asarray(grad.forces, dtype=float)


def run_atomiverse(types: np.ndarray, coords0: np.ndarray, max_steps: int, ftol: float, memory: int) -> OptimizationSummary:
    frame = make_frame(types, coords0)

    t0 = time.perf_counter()
    opt = av.GeometryOptimizer(
        frame,
        method="rhf",
        basis="sto-3g",
        memory=memory,
        ftol=ftol,
        max_steps=max_steps,
    )
    opt_result = opt.run()
    elapsed = time.perf_counter() - t0

    final_coords = frame_coords(frame)
    final_energy_ha, final_forces = evaluate_atomiverse(types, final_coords)
    max_force, rms_force = force_metrics(final_forces)

    return OptimizationSummary(
        name="Atomiverse",
        elapsed_s=elapsed,
        iterations=int(opt_result.n_steps),
        converged=max_force <= ftol,
        final_energy_ha=final_energy_ha,
        final_max_force_eva=max_force,
        final_rms_force_eva=rms_force,
        final_coords=final_coords,
    )


def run_pyscf(types: np.ndarray, coords0: np.ndarray, max_steps: int, ftol: float, scf_tol: float) -> OptimizationSummary:
    mol0 = gto.M(
        atom=coords_to_pyscf_geom(types, coords0),
        basis="sto-3g",
        unit="angstrom",
        verbose=0,
    )
    mf0 = pyscf_scf.RHF(mol0)
    mf0.conv_tol = scf_tol
    mf0.max_cycle = 200
    scanner = mf0.nuc_grad_method().as_scanner()

    cache: dict[str, object] = {}

    def evaluate(x_flat: np.ndarray) -> tuple[float, np.ndarray]:
        x = np.asarray(x_flat, dtype=float)
        x_cached = cache.get("x")
        if x_cached is not None and np.array_equal(x, x_cached):
            energy = cache["energy_eva"]
            grad = cache["grad_eva"]
            return float(energy), np.asarray(grad, dtype=float)

        coords = x.reshape(-1, 3)
        geom = coords_to_pyscf_geom(types, coords)
        energy_ha, grad_hb = scanner(geom)
        energy_eva = float(energy_ha) * HARTREE_TO_EV
        grad_eva = np.asarray(grad_hb, dtype=float) * HARTREE_PER_BOHR_TO_EV_PER_ANG

        cache["x"] = x.copy()
        cache["energy_eva"] = energy_eva
        cache["grad_eva"] = grad_eva.copy()
        return energy_eva, grad_eva

    def objective(x_flat: np.ndarray) -> tuple[float, np.ndarray]:
        energy_eva, grad_eva = evaluate(x_flat)
        return energy_eva, grad_eva.reshape(-1)

    x0 = coords0.reshape(-1).copy()
    t0 = time.perf_counter()
    result = optimize.minimize(
        objective,
        x0,
        method="L-BFGS-B",
        jac=True,
        options={
            "maxiter": max_steps,
            "gtol": ftol,
            "maxls": 20,
        },
    )
    elapsed = time.perf_counter() - t0

    final_coords = np.asarray(result.x, dtype=float).reshape(-1, 3)

    mol_final = gto.M(
        atom=coords_to_pyscf_geom(types, final_coords),
        basis="sto-3g",
        unit="angstrom",
        verbose=0,
    )
    mf_final = pyscf_scf.RHF(mol_final)
    mf_final.conv_tol = 1e-10
    mf_final.max_cycle = 200
    final_energy_ha = float(mf_final.kernel())
    final_grad_hb = np.asarray(mf_final.nuc_grad_method().kernel(), dtype=float)
    final_forces = -final_grad_hb * HARTREE_PER_BOHR_TO_EV_PER_ANG
    max_force, rms_force = force_metrics(final_forces)

    return OptimizationSummary(
        name="PySCF",
        elapsed_s=elapsed,
        iterations=int(result.nit),
        converged=bool(result.success) and max_force <= ftol,
        final_energy_ha=final_energy_ha,
        final_max_force_eva=max_force,
        final_rms_force_eva=rms_force,
        final_coords=final_coords,
    )


def print_summary(ref_coords: np.ndarray, perturbed_coords: np.ndarray,
                  atv: OptimizationSummary, pyscf: OptimizationSummary) -> None:
    init_max_disp, init_rms_disp = displacement_metrics(ref_coords, perturbed_coords)
    atv_max_disp, atv_rms_disp = displacement_metrics(ref_coords, atv.final_coords)
    pyscf_max_disp, pyscf_rms_disp = displacement_metrics(ref_coords, pyscf.final_coords)
    cross_max, cross_rms = compare_geometries(atv.final_coords, pyscf.final_coords)
    speedup = atv.elapsed_s / pyscf.elapsed_s if pyscf.elapsed_s > 0 else float("inf")
    dE = abs(atv.final_energy_ha - pyscf.final_energy_ha)

    print("=" * 72)
    print("Aspirin / STO-3G / RHF geometry optimization benchmark")
    print("=" * 72)
    print("Initial geometry")
    print(f"  Atomiverse backend:      {ATOMIVERSE_BACKEND}")
    print(f"  Atoms:                  {len(ref_coords)}")
    print(f"  Initial max displacement: {init_max_disp:.4f} A")
    print(f"  Initial RMS displacement: {init_rms_disp:.4f} A")

    print(f"\n{atv.name}")
    print(f"  Time:                   {atv.elapsed_s:.3f} s")
    print(f"  Steps used:             {atv.iterations}")
    print(f"  Final energy:           {atv.final_energy_ha:.10f} Ha")
    print(f"  Final max force:        {atv.final_max_force_eva:.6f} eV/A")
    print(f"  Final RMS force:        {atv.final_rms_force_eva:.6f} eV/A")
    print(f"  Final max disp vs ref:  {atv_max_disp:.4f} A")
    print(f"  Final RMS disp vs ref:  {atv_rms_disp:.4f} A")
    print(f"  Converged (ftol):       {atv.converged}")

    print(f"\n{pyscf.name}")
    print(f"  Time:                   {pyscf.elapsed_s:.3f} s")
    print(f"  Steps used:             {pyscf.iterations}")
    print(f"  Final energy:           {pyscf.final_energy_ha:.10f} Ha")
    print(f"  Final max force:        {pyscf.final_max_force_eva:.6f} eV/A")
    print(f"  Final RMS force:        {pyscf.final_rms_force_eva:.6f} eV/A")
    print(f"  Final max disp vs ref:  {pyscf_max_disp:.4f} A")
    print(f"  Final RMS disp vs ref:  {pyscf_rms_disp:.4f} A")
    print(f"  Converged (ftol):       {pyscf.converged}")

    print("\nComparison")
    print(f"  Atomiverse / PySCF time ratio: {speedup:.2f}x")
    print(f"  Final |dE|:                   {dE:.3e} Ha")
    print(f"  Final max atom delta:         {cross_max:.4f} A")
    print(f"  Final RMS atom delta:         {cross_rms:.4f} A")
    print("=" * 72)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-steps", type=int, default=50, help="Maximum optimization steps for both engines")
    parser.add_argument("--ftol", type=float, default=0.05, help="Force tolerance in eV/Angstrom")
    parser.add_argument("--memory", type=int, default=5, help="Atomiverse L-BFGS history size")
    parser.add_argument("--pyscf-scf-tol", type=float, default=1e-8, help="PySCF SCF convergence tolerance during optimization")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    types, ref_coords, perturbed_coords = make_perturbed_coords()

    atv = run_atomiverse(types, perturbed_coords, args.max_steps, args.ftol, args.memory)
    pyscf = run_pyscf(types, perturbed_coords, args.max_steps, args.ftol, args.pyscf_scf_tol)

    print_summary(ref_coords, perturbed_coords, atv, pyscf)


if __name__ == "__main__":
    main()
