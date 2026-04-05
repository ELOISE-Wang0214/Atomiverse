"""Example 02 — H2 / STO-3G / UHF: Atomiverse vs PySCF

Follows the C++ tutorial examples/tutorial/03_uhf.cpp.

Three cases:
  1. H atom (n_alpha=1, n_beta=0) — only solvable with UHF
  2. H2 at equilibrium (0.74 A) — RHF vs UHF (should agree)
  3. H2 PES scan (0.5 .. 4.0 A) — UHF breaks spin symmetry at stretched R

Output:
  02_h2_uhf.png          — PES + spin contamination + convergence traces
  02_h2_density_atv.cube — UHF density at R=2.5 A (Atomiverse)
  02_h2_density_ref.cube — UHF density at R=2.5 A (PySCF)

Reference: exact H atom E = -0.5 Ha; STO-3G E = -0.4666 Ha
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from utils.timing import timed

import atomiverse as av
from pyscf import gto, scf
from pyscf.tools import cubegen


# ── Geometry helpers ─────────────────────────────────────────────────────────

def _detect_dtype() -> np.dtype:
    f = av.CpuFrame(1)
    return f.get_x().dtype


DTYPE = None  # lazy-init


def get_dtype() -> np.dtype:
    global DTYPE
    if DTYPE is None:
        DTYPE = _detect_dtype()
    return DTYPE


def make_atom_frame(z: int) -> av.CpuFrame:
    dtype = get_dtype()
    frame = av.CpuFrame(1)
    frame.set_positions(
        np.zeros(1, dtype=dtype),
        np.zeros(1, dtype=dtype),
        np.zeros(1, dtype=dtype),
    )
    frame.set_mass(np.array([1.008], dtype=dtype))
    frame.set_type(np.array([z], dtype=np.int32))
    frame.set_box(np.eye(3, dtype=dtype) * 100.0)
    return frame


def make_h2_frame(r: float) -> av.CpuFrame:
    dtype = get_dtype()
    frame = av.CpuFrame(2)
    frame.set_positions(
        np.array([0.0, r], dtype=dtype),
        np.zeros(2, dtype=dtype),
        np.zeros(2, dtype=dtype),
    )
    frame.set_mass(np.array([1.008, 1.008], dtype=dtype))
    frame.set_type(np.array([1, 1], dtype=np.int32))
    frame.set_box(np.eye(3, dtype=dtype) * 100.0)
    return frame


def make_pyscf_h2(r: float) -> gto.Mole:
    return gto.M(atom=f"H 0 0 0; H {r} 0 0", basis="sto-3g",
                 unit="angstrom", verbose=0)


def run_pyscf_uhf(mol: gto.Mole, break_symmetry: bool = False) -> scf.UHF:
    """Run PySCF UHF with optional symmetry breaking initial guess."""
    trace: list[dict] = []

    def callback(envs: dict) -> None:
        trace.append({
            "iter": envs["cycle"] + 1,
            "energy": envs["e_tot"],
            "residual": envs.get("norm_ddm", float("nan")),
        })

    mf = scf.UHF(mol)
    mf.conv_tol = 1e-12
    mf.max_cycle = 200
    mf.callback = callback

    if break_symmetry:
        # Localize initial guess: alpha on atom A, beta on atom B
        nao = mol.nao
        dm_a = np.zeros((nao, nao))
        dm_b = np.zeros((nao, nao))
        dm_a[0, 0] = 1.0
        dm_b[nao - 1, nao - 1] = 1.0
        mf.kernel(dm0=np.array([dm_a, dm_b]))
    else:
        mf.kernel()

    mf._trace = trace
    return mf


# ── Cube file ────────────────────────────────────────────────────────────────

def write_density_cube(
    P: np.ndarray,
    mol: gto.Mole,
    filename: str,
    label: str,
) -> None:
    cubegen.density(mol, filename, P.astype(np.float64),
                    nx=60, ny=60, nz=60)
    print(f"  {filename}: density cube written [{label}]")


# ── Plot ─────────────────────────────────────────────────────────────────────

def plot(
    scan: dict,
    outfile: str = "02_h2_uhf.png",
) -> None:
    fig = plt.figure(figsize=(14, 4))
    gs = gridspec.GridSpec(1, 3, figure=fig, wspace=0.4)
    R = scan["R"]

    # Panel 1: PES — RHF vs UHF (both codes)
    ax1 = fig.add_subplot(gs[0])
    ax1.plot(R, scan["atv_rhf"], "-", color="#0066CC", lw=1.5,
             label="ATV RHF")
    ax1.plot(R, scan["atv_uhf"], "-o", color="#0066CC", lw=2, ms=3,
             label="ATV UHF")
    ax1.plot(R, scan["ref_rhf"], "--", color="#CC4400", lw=1.5,
             label="PySCF RHF")
    ax1.plot(R, scan["ref_uhf"], "--s", color="#CC4400", lw=2, ms=3,
             label="PySCF UHF")
    exact_limit = 2.0 * (-0.4666)
    ax1.axhline(exact_limit, color="gray", ls=":", lw=1,
                label=f"2*E_H = {exact_limit:.4f}")
    ax1.set_xlabel("R(H-H) [Angstrom]")
    ax1.set_ylabel("Energy [Hartree]")
    ax1.set_title("H2 Dissociation: RHF vs UHF")
    ax1.legend(fontsize=6)

    # Panel 2: spin contamination
    ax2 = fig.add_subplot(gs[1])
    ax2.plot(R, scan["atv_s2c"], "-o", color="#0066CC", ms=3,
             label="Atomiverse")
    ax2.plot(R, scan["ref_s2c"], "--s", color="#CC4400", ms=3,
             label="PySCF")
    ax2.axhline(0, color="gray", ls=":", lw=0.8)
    ax2.axhline(1, color="red", ls=":", lw=0.8, alpha=0.5)
    ax2.set_xlabel("R(H-H) [Angstrom]")
    ax2.set_ylabel("<S^2> - S(S+1)")
    ax2.set_title("UHF Spin Contamination")
    ax2.legend()

    # Panel 3: energy difference |ATV - PySCF|
    ax3 = fig.add_subplot(gs[2])
    dE_uhf = np.abs(np.array(scan["atv_uhf"]) - np.array(scan["ref_uhf"]))
    dE_rhf = np.abs(np.array(scan["atv_rhf"]) - np.array(scan["ref_rhf"]))
    ax3.semilogy(R, dE_rhf, "-", color="gray", lw=1.5, label="|dE| RHF")
    ax3.semilogy(R, dE_uhf, "-o", color="#0066CC", ms=3, label="|dE| UHF")
    ax3.set_xlabel("R(H-H) [Angstrom]")
    ax3.set_ylabel("|dE| [Hartree]")
    ax3.set_title("Agreement (ATV vs PySCF)")
    ax3.legend()

    max_dE = max(dE_uhf)
    fig.suptitle(
        f"H2 / STO-3G / UHF  —  {len(R)} points, max|dE|={max_dE:.2e} Ha",
        fontsize=9,
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

    frame_h = make_atom_frame(1)
    ip_h = av.build_integrals_sto3g(frame_h)
    uhf_h = av.uhf_solve(ip_h, n_alpha=1, n_beta=0, tol=1e-10, max_iter=50)

    # PySCF H atom
    mol_h = gto.M(atom="H 0 0 0", basis="sto-3g", spin=1, verbose=0)
    mf_h = run_pyscf_uhf(mol_h)
    ss_h, _ = mf_h.spin_square()
    # spin_contamination = <S^2> - S(S+1), S=0.5 -> S(S+1)=0.75
    ref_s2c_h = ss_h - 0.75

    dE_h = abs(float(uhf_h.energy) - mf_h.e_tot)
    print(f"  Atomiverse: E = {float(uhf_h.energy):.8f} Ha, "
          f"<S^2>_c = {float(uhf_h.spin_contamination):.6f}")
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
    ip_eq = av.build_integrals_sto3g(frame_eq)
    rhf_eq = av.rhf_solve(ip_eq, tol=1e-10, max_iter=100)
    uhf_eq = av.uhf_solve(ip_eq, n_alpha=1, n_beta=1, tol=1e-10, max_iter=200)

    mol_eq = make_pyscf_h2(0.74)
    mf_rhf_eq = scf.RHF(mol_eq)
    mf_rhf_eq.conv_tol = 1e-12
    mf_rhf_eq.kernel()
    mf_uhf_eq = run_pyscf_uhf(mol_eq)
    ss_eq, _ = mf_uhf_eq.spin_square()

    print(f"  Atomiverse RHF: E = {float(rhf_eq.energy):.8f} Ha")
    print(f"  Atomiverse UHF: E = {float(uhf_eq.energy):.8f} Ha, "
          f"<S^2>_c = {float(uhf_eq.spin_contamination):.6f}")
    print(f"  PySCF RHF:      E = {mf_rhf_eq.e_tot:.8f} Ha")
    print(f"  PySCF UHF:      E = {mf_uhf_eq.e_tot:.8f} Ha, "
          f"<S^2>_c = {ss_eq:.6f}")
    dE_eq = abs(float(uhf_eq.energy) - float(rhf_eq.energy))
    print(f"  |E_UHF - E_RHF| = {dE_eq:.2e} Ha (should be ~0 at equilibrium)")

    # =====================================================================
    # Case 3: H2 PES scan — RHF vs UHF
    # =====================================================================
    print(f"\n{'='*60}")
    print("Case 3: H2 PES scan / STO-3G / RHF vs UHF")
    print("=" * 60)

    R_vals = np.arange(0.5, 4.01, 0.25)
    scan: dict[str, list] = {
        "R": R_vals,
        "atv_rhf": [], "atv_uhf": [], "atv_s2c": [],
        "ref_rhf": [], "ref_uhf": [], "ref_s2c": [],
    }

    header = (f"  {'R':>6s}  {'ATV_RHF':>12s}  {'ATV_UHF':>12s}  "
              f"{'PySCF_UHF':>12s}  {'dE_UHF':>10s}  "
              f"{'ATV_S2c':>8s}  {'PySCF_S2c':>10s}")
    print(header)
    print("  " + "-" * (len(header) - 2))

    for r in R_vals:
        # Atomiverse
        frame = make_h2_frame(r)
        ip = av.build_integrals_sto3g(frame)
        rhf = av.rhf_solve(ip, tol=1e-10, max_iter=100)
        uhf = av.uhf_solve(ip, n_alpha=1, n_beta=1, tol=1e-10, max_iter=400)

        # PySCF (break symmetry to match Atomiverse)
        mol = make_pyscf_h2(r)
        mf_r = scf.RHF(mol)
        mf_r.conv_tol = 1e-12
        mf_r.kernel()
        mf_u = run_pyscf_uhf(mol, break_symmetry=True)
        ss, _ = mf_u.spin_square()
        # For H2 singlet: S=0, S(S+1)=0, so contamination = <S^2>
        ref_s2c = ss - 0.0

        scan["atv_rhf"].append(float(rhf.energy))
        scan["atv_uhf"].append(float(uhf.energy))
        scan["atv_s2c"].append(float(uhf.spin_contamination))
        scan["ref_rhf"].append(float(mf_r.e_tot))
        scan["ref_uhf"].append(float(mf_u.e_tot))
        scan["ref_s2c"].append(float(ref_s2c))

        dE = abs(scan["atv_uhf"][-1] - scan["ref_uhf"][-1])
        print(f"  {r:6.3f}  {scan['atv_rhf'][-1]:12.8f}  "
              f"{scan['atv_uhf'][-1]:12.8f}  {scan['ref_uhf'][-1]:12.8f}  "
              f"{dE:10.2e}  {scan['atv_s2c'][-1]:8.4f}  "
              f"{scan['ref_s2c'][-1]:10.4f}")

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
    frame_c = make_h2_frame(r_cube)
    ip_c = av.build_integrals_sto3g(frame_c)
    uhf_c = av.uhf_solve(ip_c, n_alpha=1, n_beta=1, tol=1e-10, max_iter=400)
    P_atv = np.array(uhf_c.get_P_alpha()) + np.array(uhf_c.get_P_beta())

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
