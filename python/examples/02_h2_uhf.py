"""Example 02 — H2 / STO-3G / UHF: Atomiverse vs PySCF

Follows the C++ tutorial examples/03_uhf.cpp.

Three cases:
  1. H atom (n_alpha=1, n_beta=0) — only solvable with UHF
  2. H2 at equilibrium (0.74 A) — RHF vs UHF (should agree)
  3. H2 PES scan (0.5 .. 4.0 A) — UHF breaks spin symmetry at stretched R

Output:
  02_h2_uhf.png          — PES + spin contamination + timing benchmark
  02_h2_density_atv.cube — UHF density at R=2.5 A (Atomiverse)
  02_h2_density_ref.cube — UHF density at R=2.5 A (PySCF)

Reference: exact H atom E = -0.5 Ha; STO-3G E = -0.4666 Ha
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
C_UHF = "#888888"


# ── Geometry helpers ─────────────────────────────────────────────────────────

def make_atom_frame(z: int) -> mp.Frame:
    frame = mp.Frame({"atoms": mp.Block({
        "x": np.zeros(1), "y": np.zeros(1), "z": np.zeros(1),
        "mass": np.array([1.008]),
        "type": np.array([z], dtype=np.int32),
    })})
    frame.box = mp.Box(np.eye(3) * 100.0)
    return frame


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


def run_pyscf_uhf(mol: gto.Mole, break_symmetry: bool = False) -> scf.UHF:
    """Run PySCF UHF with optional symmetry breaking initial guess."""
    mf = scf.UHF(mol)
    mf.conv_tol = 1e-12
    mf.max_cycle = 200
    if break_symmetry:
        nao = mol.nao
        dm_a = np.zeros((nao, nao))
        dm_b = np.zeros((nao, nao))
        dm_a[0, 0] = 1.0
        dm_b[nao - 1, nao - 1] = 1.0
        mf.kernel(dm0=np.array([dm_a, dm_b]))
    else:
        mf.kernel()
    return mf


# ── Cube file ────────────────────────────────────────────────────────────────

def write_density_cube(P, mol, filename, label):
    cubegen.density(mol, filename, P.astype(np.float64), nx=60, ny=60, nz=60)
    print(f"  {filename}: density cube written [{label}]")


# ── Plot ─────────────────────────────────────────────────────────────────────

def plot(scan: dict, outfile: str = "02_h2_uhf.png") -> None:
    fig = plt.figure(figsize=(14, 4.2))
    gs = gridspec.GridSpec(1, 3, figure=fig, wspace=0.38)
    R = np.asarray(scan["R"])

    # ── Panel 1: PES — RHF vs UHF (both codes) ──
    ax1 = fig.add_subplot(gs[0])
    ax1.plot(R, scan["atv_rhf"], "-", color=C_UHF, lw=1.2, label="ATV RHF")
    ax1.plot(R, scan["atv_uhf"], "-o", color=C_ATV, lw=1.8, ms=3.5,
             label="ATV UHF")
    ax1.plot(R, scan["ref_rhf"], "--", color=C_UHF, lw=1.2, alpha=0.6,
             label="PySCF RHF")
    ax1.plot(R, scan["ref_uhf"], "--s", color=C_REF, lw=1.8, ms=3,
             label="PySCF UHF")
    exact_limit = 2.0 * (-0.4666)
    ax1.axhline(exact_limit, color="gray", ls=":", lw=0.8,
                label=rf"$2E_{{\mathrm{{H}}}}^{{\mathrm{{STO\text{{-}}3G}}}}={exact_limit:.4f}$")
    ax1.set_xlabel(r"$R_{\mathrm{H\text{-}H}}$ [\AA]")
    ax1.set_ylabel(r"$E_{\mathrm{tot}}$ [Ha]")
    ax1.set_title(r"$\mathrm{H_2}$ Dissociation: RHF vs UHF")
    ax1.legend()

    # ── Panel 2: spin contamination ──
    ax2 = fig.add_subplot(gs[1])
    ax2.plot(R, scan["atv_s2c"], "-o", color=C_ATV, lw=1.8, ms=3.5,
             label="Atomiverse")
    ax2.plot(R, scan["ref_s2c"], "--s", color=C_REF, lw=1.8, ms=3,
             label="PySCF")
    ax2.axhline(0, color="gray", ls=":", lw=0.8)
    ax2.axhline(1, color="red", ls=":", lw=0.8, alpha=0.4)
    ax2.set_xlabel(r"$R_{\mathrm{H\text{-}H}}$ [\AA]")
    ax2.set_ylabel(r"$\langle\hat{S}^2\rangle - S(S+1)$")
    ax2.set_title(r"UHF Spin Contamination")
    ax2.legend()

    # ── Panel 3: timing benchmark ──
    ax3 = fig.add_subplot(gs[2])
    atv_ms = np.asarray(scan["atv_time_ms"])
    ref_ms = np.asarray(scan["ref_time_ms"])
    ax3.plot(R, atv_ms, "-o", color=C_ATV, lw=1.8, ms=3.5,
             label="Atomiverse")
    ax3.plot(R, ref_ms, "--s", color=C_REF, lw=1.8, ms=3,
             label="PySCF")
    ax3.set_xlabel(r"$R_{\mathrm{H\text{-}H}}$ [\AA]")
    ax3.set_ylabel(r"Wall time [ms]")
    ax3.set_title(r"UHF Solve Time")
    ax3.legend()

    mean_speedup = np.mean(ref_ms / np.where(atv_ms > 0, atv_ms, 1))
    ax3.text(0.97, 0.95,
             rf"mean $\times${mean_speedup:.1f}",
             transform=ax3.transAxes, ha="right", va="top",
             fontsize=8, color=C_ATV,
             bbox=dict(boxstyle="round,pad=0.3", fc="white",
                       ec=C_ATV, alpha=0.8))

    dE_uhf = np.abs(np.array(scan["atv_uhf"]) - np.array(scan["ref_uhf"]))
    max_dE = np.max(dE_uhf)
    fig.suptitle(
        rf"$\mathrm{{H_2}}$ / STO-3G / UHF  —  "
        rf"{len(R)} points, max$|\Delta E|={max_dE:.2e}$ Ha",
        fontsize=10,
    )
    fig.savefig(outfile, dpi=150, bbox_inches="tight")
    print(f"\nSaved {outfile}")
    plt.close(fig)


# ── Main ─────────────────────────────────────────────────────────────────────

def main() -> None:

    # =====================================================================
    # Case 1: H atom
    # =====================================================================
    print("=" * 60)
    print("Case 1: H atom / STO-3G / UHF")
    print("=" * 60)

    uhf_h = av.UHF(make_atom_frame(1), multiplicity=2).solve(tol=1e-10, max_iter=50)

    mol_h = gto.M(atom="H 0 0 0", basis="sto-3g", spin=1, verbose=0)
    mf_h = run_pyscf_uhf(mol_h)
    ss_h, _ = mf_h.spin_square()
    ref_s2c_h = ss_h - 0.75

    dE_h = abs(uhf_h.energy - mf_h.e_tot)
    print(f"  Atomiverse: E = {uhf_h.energy:.8f} Ha, "
          f"<S^2>_c = {uhf_h.spin_contamination:.6f}")
    print(f"  PySCF:      E = {mf_h.e_tot:.8f} Ha, "
          f"<S^2>_c = {ref_s2c_h:.6f}")
    print(f"  |dE| = {dE_h:.2e} Ha  (STO-3G reference: -0.4666 Ha)")

    # =====================================================================
    # Case 2: H2 at equilibrium — RHF vs UHF
    # =====================================================================
    print(f"\n{'='*60}")
    print("Case 2: H2 (0.74 Angstrom) / STO-3G / RHF vs UHF")
    print("=" * 60)

    frame_eq = make_h2_frame(0.74)
    rhf_eq = av.RHF(frame_eq).solve(tol=1e-10, max_iter=100)
    uhf_eq = av.UHF(frame_eq, multiplicity=1).solve(tol=1e-10, max_iter=200)

    mol_eq = make_pyscf_h2(0.74)
    mf_rhf_eq = scf.RHF(mol_eq)
    mf_rhf_eq.conv_tol = 1e-12
    mf_rhf_eq.kernel()
    mf_uhf_eq = run_pyscf_uhf(mol_eq)
    ss_eq, _ = mf_uhf_eq.spin_square()

    print(f"  Atomiverse RHF: E = {rhf_eq.energy:.8f} Ha")
    print(f"  Atomiverse UHF: E = {uhf_eq.energy:.8f} Ha, "
          f"<S^2>_c = {uhf_eq.spin_contamination:.6f}")
    print(f"  PySCF RHF:      E = {mf_rhf_eq.e_tot:.8f} Ha")
    print(f"  PySCF UHF:      E = {mf_uhf_eq.e_tot:.8f} Ha, "
          f"<S^2>_c = {ss_eq:.6f}")
    dE_eq = abs(uhf_eq.energy - rhf_eq.energy)
    print(f"  |E_UHF - E_RHF| = {dE_eq:.2e} Ha (should be ~0 at equilibrium)")

    # =====================================================================
    # Case 3: H2 PES scan — RHF vs UHF + timing benchmark
    # =====================================================================
    print(f"\n{'='*60}")
    print("Case 3: H2 PES scan / STO-3G / RHF vs UHF")
    print("=" * 60)

    R_vals = np.arange(0.5, 4.01, 0.25)
    scan: dict[str, list] = {
        "R": R_vals,
        "atv_rhf": [], "atv_uhf": [], "atv_s2c": [],
        "ref_rhf": [], "ref_uhf": [], "ref_s2c": [],
        "atv_time_ms": [], "ref_time_ms": [],
    }

    header = (f"  {'R':>5s}  {'ATV_UHF':>12s}  {'PySCF_UHF':>12s}  "
              f"{'dE':>10s}  {'ATV ms':>8s}  {'PySCF ms':>8s}  "
              f"{'S2c_ATV':>8s}")
    print(header)
    print("  " + "-" * (len(header) - 2))

    for r in R_vals:
        frame = make_h2_frame(r)

        # Atomiverse (timed: integrals + solve together)
        t0 = time.perf_counter()
        rhf_r = av.RHF(frame).solve(tol=1e-10, max_iter=100)
        uhf_r = av.UHF(frame, multiplicity=1).solve(tol=1e-10, max_iter=400)
        atv_dt = (time.perf_counter() - t0) * 1e3

        # PySCF (timed)
        mol = make_pyscf_h2(r)
        t0 = time.perf_counter()
        mf_r = scf.RHF(mol)
        mf_r.conv_tol = 1e-12
        mf_r.kernel()
        mf_u = run_pyscf_uhf(mol, break_symmetry=True)
        ref_dt = (time.perf_counter() - t0) * 1e3

        ss, _ = mf_u.spin_square()
        ref_s2c = ss - 0.0

        scan["atv_rhf"].append(rhf_r.energy)
        scan["atv_uhf"].append(uhf_r.energy)
        scan["atv_s2c"].append(uhf_r.spin_contamination)
        scan["ref_rhf"].append(float(mf_r.e_tot))
        scan["ref_uhf"].append(float(mf_u.e_tot))
        scan["ref_s2c"].append(float(ref_s2c))
        scan["atv_time_ms"].append(atv_dt)
        scan["ref_time_ms"].append(ref_dt)

        dE = abs(scan["atv_uhf"][-1] - scan["ref_uhf"][-1])
        print(f"  {r:5.2f}  {scan['atv_uhf'][-1]:12.8f}  "
              f"{scan['ref_uhf'][-1]:12.8f}  {dE:10.2e}  "
              f"{atv_dt:8.1f}  {ref_dt:8.1f}  "
              f"{scan['atv_s2c'][-1]:8.4f}")

    # ── Summary ──
    max_dE_uhf = max(abs(a - b)
                     for a, b in zip(scan["atv_uhf"], scan["ref_uhf"]))
    max_dE_rhf = max(abs(a - b)
                     for a, b in zip(scan["atv_rhf"], scan["ref_rhf"]))

    print(f"\n  Exact dissociation limit: 2 * E_H = {2*(-0.4666):.4f} Ha")
    print(f"\n{'='*60}")
    print("Agreement (Atomiverse vs PySCF)")
    print("=" * 60)
    print(f"  max|dE| RHF = {max_dE_rhf:.2e} Ha")
    print(f"  max|dE| UHF = {max_dE_uhf:.2e} Ha")

    atv_total = sum(scan["atv_time_ms"])
    ref_total = sum(scan["ref_time_ms"])
    print(f"\n  Timing (full scan):")
    print(f"    Atomiverse: {atv_total:.0f} ms")
    print(f"    PySCF:      {ref_total:.0f} ms")
    if atv_total > 0:
        print(f"    Speedup:    {ref_total / atv_total:.1f}x")

    rel_err_uhf = max_dE_uhf / abs(min(scan["ref_uhf"])) * 100
    THRESHOLD = 1.0  # 1%
    passed = rel_err_uhf < THRESHOLD
    status = "PASS" if passed else "FAIL"
    print(f"\n  Result: {status}  (max relative error {rel_err_uhf:.4f}% "
          f"vs {THRESHOLD}% threshold)")

    # ── Cube files at stretched R ──
    print(f"\n{'='*60}")
    print("Electron density cube files (R=2.50 Angstrom)")
    print("=" * 60)
    r_cube = 2.50
    uhf_c = av.UHF(make_h2_frame(r_cube), multiplicity=1) \
              .solve(tol=1e-10, max_iter=400)
    P_atv = uhf_c.density_matrix_alpha + uhf_c.density_matrix_beta

    mol_c = make_pyscf_h2(r_cube)
    mf_c = run_pyscf_uhf(mol_c, break_symmetry=True)
    P_ref = mf_c.make_rdm1()[0] + mf_c.make_rdm1()[1]

    write_density_cube(P_atv, mol_c, "02_h2_density_atv.cube", "Atomiverse UHF")
    write_density_cube(P_ref, mol_c, "02_h2_density_ref.cube", "PySCF UHF")

    # ── Plot ──
    plot(scan)

    if not passed:
        sys.exit(1)

    print(f"\n{'='*60}")
    print("Done")
    print("=" * 60)


if __name__ == "__main__":
    main()
