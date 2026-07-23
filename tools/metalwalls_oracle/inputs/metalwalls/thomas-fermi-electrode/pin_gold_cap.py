#!/usr/bin/env python3
"""Offline pin: gold capacitor toy with Thomas–Fermi hardness (MetalWalls model).

Mirrors Atomiverse cpu::md::ConstPotBO matrix algebra:
  Gaussian electrode elastance + 2D xy lattice images + TF diagonal
  A_ii += 4*pi*ke * l_tf**2 * site_density
  + electroneutrality projector.

Never launched from ctest. See regressions/md/thomas-fermi-electrode/REFERENCE.md.
"""
from __future__ import annotations

import math

import numpy as np

KE = 14.399645  # units::coulomb [eV·Å/e²]
SQRT2 = math.sqrt(2.0)
SQRT2_OVER_PI = math.sqrt(2.0 / math.pi)

ETA = 1.0 / 0.56
A0 = 4.07
D_SITE = 4.0 / (A0**3)
A_PLANE = A0 / 2.0
L_GAP = 20.0
LX = LY = A0
R = 40


def pair_a(r: float, eta: float) -> float:
    if r < 1e-14:
        return KE * eta * SQRT2_OVER_PI
    etaij = eta / math.sqrt(2.0)
    return KE * math.erf(etaij * r) / r


def gold_cap_atoms():
    z_surf_bot = 5.0 + 2 * A_PLANE
    z_surf_top = z_surf_bot + L_GAP
    z_bots = [z_surf_bot - 2 * A_PLANE, z_surf_bot - A_PLANE, z_surf_bot]
    z_tops = [z_surf_top, z_surf_top + A_PLANE, z_surf_top + 2 * A_PLANE]
    atoms = [(0.0, 0.0, z, 1) for z in z_bots] + [(0.0, 0.0, z, 2) for z in z_tops]
    return atoms


def build_a(atoms, l_tf: float) -> np.ndarray:
    m = len(atoms)
    a_mat = np.zeros((m, m))
    for i in range(m):
        a_mat[i, i] = KE * ETA * SQRT2_OVER_PI
        xi, yi, zi, _ = atoms[i]
        for j in range(i + 1, m):
            xj, yj, zj, _ = atoms[j]
            aij = 0.0
            for nx in range(-R, R + 1):
                for ny in range(-R, R + 1):
                    dx = xi - xj + nx * LX
                    dy = yi - yj + ny * LY
                    dz = zi - zj
                    r = math.sqrt(dx * dx + dy * dy + dz * dz)
                    if r < 1e-14:
                        continue
                    aij += pair_a(r, ETA)
            a_mat[i, j] = a_mat[j, i] = aij
        self_img = 0.0
        for nx in range(-R, R + 1):
            for ny in range(-R, R + 1):
                if nx == 0 and ny == 0:
                    continue
                r = math.sqrt((nx * LX) ** 2 + (ny * LY) ** 2)
                self_img += KE * math.erf((ETA / SQRT2) * r) / r
        a_mat[i, i] += self_img
    if l_tf > 0.0:
        hardness = 4.0 * math.pi * KE * l_tf * l_tf * D_SITE
        for i in range(m):
            a_mat[i, i] += hardness
    return a_mat


def solve_neutral(a_mat: np.ndarray, psi: np.ndarray) -> np.ndarray:
    inv = np.linalg.inv(a_mat)
    e = np.ones(a_mat.shape[0])
    u = inv @ e
    v = inv @ psi
    return v - (e @ v) / (e @ u) * u


def main() -> None:
    atoms = gold_cap_atoms()
    m = len(atoms)
    psi = np.array([0.0] * 3 + [1.0] * 3)
    print(f"# thomas-fermi-electrode pin  eta={ETA:.12f} d={D_SITE:.12f} R={R}")
    print(f"# atoms z: {[a[2] for a in atoms]}")
    for l_tf in (0.0, 2.0):
        a_mat = build_a(atoms, l_tf)
        q = solve_neutral(a_mat, psi)
        q_top = float(q[3:].sum())
        q_bot = float(q[:3].sum())
        c = abs(q_top) / 1.0
        print(
            f"l_tf={l_tf:.1f}  Q_bot={q_bot:+.12e}  Q_top={q_top:+.12e}  C={c:.12e}  |sumQ|={abs(q.sum()):.3e}"
        )


if __name__ == "__main__":
    main()
