"""Example 01 — H2 / STO-3G / RHF: Atomiverse vs PySCF

H2 dissociation PES with Restricted Hartree-Fock.
Part of a series: 01 (RHF) → 02 (UHF) → 03 (UKS-LSDA).

Output:
  01_h2_rhf.png          — PES + orbital energies + timing benchmark
  01_h2_density_atv.cube — electron density at equilibrium (Atomiverse)
  01_h2_density_ref.cube — electron density at equilibrium (PySCF)

Reference: Szabo & Ostlund Table 3.7 -> E = -1.1175 Ha (STO-3G, R=0.74 A)
"""

from __future__ import annotations

import sys
import time

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import molpy as mp

import atomiverse as av
from pyscf import gto, scf
from pyscf.tools import cubegen

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


# ── Geometry helpers ─────────────────────────────────────────────────────────

def make_h2_frame(r: float) -> mp.Frame:
    frame = mp.Frame({"atoms": mp.Block({
        "x": np.array([0.0, r]),
        "y": np.zeros(2),
        "z": np.zeros(2),
        "mass": np.array([1.008, 1.008]),
        "type": np.array([1, 1], dtype=np.int32),
    })})
    frame.box = mp.Box(np.eye(3) * 100.0)
    return frame


def make_pyscf_h2(r: float) -> gto.Mole:
    return gto.M(atom=f"H 0 0 0; H {r} 0 0", basis="sto-3g",
                 unit="angstrom", verbose=0)


# ── Cube file ────────────────────────────────────────────────────────────────

def write_density_cube(P, mol, filename, label):
    cubegen.density(mol, filename, P.astype(np.float64), nx=60, ny=60, nz=60)
    print(f"  {filename}: density cube written [{label}]")


# ── Plot ─────────────────────────────────────────────────────────────────────

def plot(scan: dict, outfile: str = "01_h2_rhf.png") -> None:
    fig = plt.figure(figsize=(14, 4.2))
    gs = gridspec.GridSpec(1, 3, figure=fig, wspace=0.38)
    R = np.asarray(scan["R"])

    # ── Panel 1: PES ──
    ax1 = fig.add_subplot(gs[0])
    ax1.plot(R, scan["atv_energy"], "-o", color=C_ATV, lw=1.8, ms=3.5,
             label="Atomiverse")
    ax1.plot(R, scan["ref_energy"], "--s", color=C_REF, lw=1.8, ms=3,
             label="PySCF")
    ax1.set_xlabel(r"$R_{\mathrm{H\text{-}H}}$ [\AA]")
    ax1.set_ylabel(r"$E_{\mathrm{tot}}$ [Ha]")
    ax1.set_title(r"$\mathrm{H_2}$ RHF Dissociation")
    ax1.legend()

    # ── Panel 2: orbital energies vs R ──
    ax2 = fig.add_subplot(gs[1])
    atv_homo = np.asarray(scan["atv_homo"])
    atv_lumo = np.asarray(scan["atv_lumo"])
    ref_homo = np.asarray(scan["ref_homo"])
    ref_lumo = np.asarray(scan["ref_lumo"])
    ax2.plot(R, atv_homo, "-o", color=C_ATV, lw=1.8, ms=3.5,
             label=r"ATV $\varepsilon_{\mathrm{HOMO}}$")
    ax2.plot(R, atv_lumo, "-^", color=C_ATV, lw=1.2, ms=3.5, alpha=0.6,
             label=r"ATV $\varepsilon_{\mathrm{LUMO}}$")
    ax2.plot(R, ref_homo, "--s", color=C_REF, lw=1.8, ms=3,
             label=r"PySCF $\varepsilon_{\mathrm{HOMO}}$")
    ax2.plot(R, ref_lumo, "--v", color=C_REF, lw=1.2, ms=3, alpha=0.6,
             label=r"PySCF $\varepsilon_{\mathrm{LUMO}}$")
    ax2.set_xlabel(r"$R_{\mathrm{H\text{-}H}}$ [\AA]")
    ax2.set_ylabel(r"$\varepsilon_i$ [Ha]")
    ax2.set_title(r"Orbital Energies")
    ax2.legend()

    # ── Panel 3: timing ──
    ax3 = fig.add_subplot(gs[2])
    atv_ms = np.asarray(scan["atv_time_ms"])
    ref_ms = np.asarray(scan["ref_time_ms"])
    ax3.plot(R, atv_ms, "-o", color=C_ATV, lw=1.8, ms=3.5,
             label="Atomiverse")
    ax3.plot(R, ref_ms, "--s", color=C_REF, lw=1.8, ms=3,
             label="PySCF")
    ax3.set_xlabel(r"$R_{\mathrm{H\text{-}H}}$ [\AA]")
    ax3.set_ylabel(r"Wall time [ms]")
    ax3.set_title(r"RHF Solve Time")
    ax3.legend()

    mean_speedup = np.mean(ref_ms / np.where(atv_ms > 0, atv_ms, 1))
    ax3.text(0.97, 0.95,
             rf"mean $\times${mean_speedup:.1f}",
             transform=ax3.transAxes, ha="right", va="top",
             fontsize=8, color=C_ATV,
             bbox=dict(boxstyle="round,pad=0.3", fc="white",
                       ec=C_ATV, alpha=0.8))

    dE = np.abs(np.asarray(scan["atv_energy"]) - np.asarray(scan["ref_energy"]))
    fig.suptitle(
        rf"$\mathrm{{H_2}}$ / STO-3G / RHF  —  "
        rf"{len(R)} points, max$|\Delta E|={np.max(dE):.2e}$ Ha",
        fontsize=10,
    )
    fig.savefig(outfile, dpi=150, bbox_inches="tight")
    print(f"\nSaved {outfile}")
    plt.close(fig)


# ── Main ─────────────────────────────────────────────────────────────────────

def main() -> None:
    print("=" * 60)
    print("H2 / STO-3G / RHF  —  Atomiverse vs PySCF")
    print("=" * 60)

    R_vals = np.arange(0.5, 4.01, 0.25)
    scan: dict[str, list] = {
        "R": R_vals,
        "atv_energy": [], "ref_energy": [],
        "atv_homo": [], "atv_lumo": [],
        "ref_homo": [], "ref_lumo": [],
        "atv_time_ms": [], "ref_time_ms": [],
    }

    header = (f"  {'R':>5s}  {'ATV_RHF':>12s}  {'PySCF_RHF':>12s}  "
              f"{'|dE|':>10s}  {'ATV ms':>8s}  {'PySCF ms':>8s}")
    print(header)
    print("  " + "-" * (len(header) - 2))

    for r in R_vals:
        frame = make_h2_frame(r)

        # Atomiverse (timed)
        t0 = time.perf_counter()
        rhf = av.RHF(frame, basis="sto-3g")
        result = rhf.solve(tol=1e-10, max_iter=100)
        atv_dt = (time.perf_counter() - t0) * 1e3

        # PySCF (timed)
        mol = make_pyscf_h2(r)
        t0 = time.perf_counter()
        mf = scf.RHF(mol)
        mf.conv_tol = 1e-12
        mf.kernel()
        ref_dt = (time.perf_counter() - t0) * 1e3

        eig_atv = result.eigenvalues
        eig_ref = mf.mo_energy

        scan["atv_energy"].append(result.energy)
        scan["ref_energy"].append(float(mf.e_tot))
        scan["atv_homo"].append(float(eig_atv[0]))
        scan["atv_lumo"].append(float(eig_atv[1]))
        scan["ref_homo"].append(float(eig_ref[0]))
        scan["ref_lumo"].append(float(eig_ref[1]))
        scan["atv_time_ms"].append(atv_dt)
        scan["ref_time_ms"].append(ref_dt)

        dE = abs(result.energy - mf.e_tot)
        print(f"  {r:5.2f}  {result.energy:12.8f}  {mf.e_tot:12.8f}  "
              f"{dE:10.2e}  {atv_dt:8.1f}  {ref_dt:8.1f}")

    # ── Summary ──
    max_dE = max(abs(a - b) for a, b in
                 zip(scan["atv_energy"], scan["ref_energy"]))
    rel_err = max_dE / abs(min(scan["ref_energy"])) * 100

    atv_total = sum(scan["atv_time_ms"])
    ref_total = sum(scan["ref_time_ms"])

    print(f"\n{'='*60}")
    print("Summary")
    print("=" * 60)
    print(f"  max|dE| = {max_dE:.2e} Ha  ({rel_err:.4f}%)")
    print(f"  Timing: ATV {atv_total:.0f} ms, PySCF {ref_total:.0f} ms", end="")
    if atv_total > 0:
        print(f"  ({ref_total / atv_total:.1f}x)")
    else:
        print()

    THRESHOLD = 1.0
    passed = rel_err < THRESHOLD
    status = "PASS" if passed else "FAIL"
    print(f"  Result: {status}  (threshold {THRESHOLD}%)")

    # ── Cube file at equilibrium ──
    print(f"\n{'='*60}")
    print("Electron density cube (R=0.74 A)")
    print("=" * 60)
    frame_eq = make_h2_frame(0.74)
    P_atv = av.RHF(frame_eq).solve(tol=1e-10).density_matrix
    mol_eq = make_pyscf_h2(0.74)
    mf_eq = scf.RHF(mol_eq)
    mf_eq.conv_tol = 1e-12
    mf_eq.kernel()
    P_ref = mf_eq.make_rdm1()

    write_density_cube(P_atv, mol_eq, "01_h2_density_atv.cube", "Atomiverse")
    write_density_cube(P_ref, mol_eq, "01_h2_density_ref.cube", "PySCF")

    # ── Plot ──
    plot(scan)

    if not passed:
        sys.exit(1)

    print(f"\n{'='*60}")
    print("Done")
    print("=" * 60)


if __name__ == "__main__":
    main()
