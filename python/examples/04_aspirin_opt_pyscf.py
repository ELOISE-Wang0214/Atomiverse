"""PySCF-only aspirin geometry-optimization benchmark.

Optimizes the same mildly distorted aspirin geometry used by the
Atomiverse/PySCF comparison script, but runs only the PySCF path.
"""

from __future__ import annotations

import argparse
import time
from dataclasses import dataclass

import numpy as np
from pyscf import gto, lib as pyscf_lib, scf as pyscf_scf
from scipy import optimize

HARTREE_TO_EV = 27.211386245988
BOHR_TO_ANGSTROM = 0.529177210903
HARTREE_PER_BOHR_TO_EV_PER_ANG = HARTREE_TO_EV / BOHR_TO_ANGSTROM

SYMBOL = {1: "H", 6: "C", 8: "O"}

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
class PyscfOptSummary:
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
    idx = np.arange(n_atoms, dtype=float)

    perturb = np.zeros_like(coords)
    perturb[:, 0] += 0.015 * np.cos(0.7 * idx)
    perturb[:, 1] += 0.020 * np.sin(0.9 * idx)
    perturb[:, 2] += 0.040 * np.sin(0.6 * idx + 0.3)

    perturb[7] += np.array([-0.08, 0.10, 0.04])
    perturb[8] += np.array([0.05, -0.12, -0.03])
    perturb[10] += np.array([0.07, 0.09, 0.02])
    perturb[12] += np.array([-0.10, 0.05, 0.18])
    perturb[17] += np.array([-0.14, -0.02, 0.05])
    perturb[18] += np.array([0.06, -0.08, 0.10])
    perturb[19] += np.array([-0.04, 0.07, -0.06])
    perturb[20] += np.array([0.03, -0.05, 0.08])

    return types, coords, coords + perturb


def coords_to_geom(types: np.ndarray, coords: np.ndarray) -> list[tuple[str, tuple[float, float, float]]]:
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


def run_pyscf(types: np.ndarray, coords0: np.ndarray, max_steps: int, ftol: float, scf_tol: float) -> PyscfOptSummary:
    mol0 = gto.M(
        atom=coords_to_geom(types, coords0),
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
            return float(cache["energy_eva"]), np.asarray(cache["grad_eva"], dtype=float)

        coords = x.reshape(-1, 3)
        energy_ha, grad_hb = scanner(coords_to_geom(types, coords))
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
        atom=coords_to_geom(types, final_coords),
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

    return PyscfOptSummary(
        elapsed_s=elapsed,
        iterations=int(result.nit),
        converged=bool(result.success) and max_force <= ftol,
        final_energy_ha=final_energy_ha,
        final_max_force_eva=max_force,
        final_rms_force_eva=rms_force,
        final_coords=final_coords,
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-steps", type=int, default=50, help="Maximum optimization steps")
    parser.add_argument("--ftol", type=float, default=0.05, help="Force tolerance in eV/Angstrom")
    parser.add_argument("--pyscf-scf-tol", type=float, default=1e-8, help="SCF convergence tolerance during optimization")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    types, ref_coords, perturbed_coords = make_perturbed_coords()
    init_max_disp, init_rms_disp = displacement_metrics(ref_coords, perturbed_coords)
    summary = run_pyscf(types, perturbed_coords, args.max_steps, args.ftol, args.pyscf_scf_tol)
    final_max_disp, final_rms_disp = displacement_metrics(ref_coords, summary.final_coords)

    print("=" * 72)
    print("PySCF aspirin / STO-3G / RHF geometry optimization benchmark")
    print("=" * 72)
    print(f"PySCF threads:            {pyscf_lib.num_threads()}")
    print(f"Atoms:                    {len(ref_coords)}")
    print(f"Initial max displacement: {init_max_disp:.4f} A")
    print(f"Initial RMS displacement: {init_rms_disp:.4f} A")
    print(f"Time:                     {summary.elapsed_s:.3f} s")
    print(f"Steps used:               {summary.iterations}")
    print(f"Final energy:             {summary.final_energy_ha:.10f} Ha")
    print(f"Final max force:          {summary.final_max_force_eva:.6f} eV/A")
    print(f"Final RMS force:          {summary.final_rms_force_eva:.6f} eV/A")
    print(f"Final max disp vs ref:    {final_max_disp:.4f} A")
    print(f"Final RMS disp vs ref:    {final_rms_disp:.4f} A")
    print(f"Converged (ftol):         {summary.converged}")
    print("=" * 72)


if __name__ == "__main__":
    main()
