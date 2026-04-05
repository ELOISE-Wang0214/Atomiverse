"""Example 01 — H2 / STO-3G / RHF: Atomiverse vs PySCF

Follows the C++ tutorial examples/tutorial/02_rhf.cpp.

End-to-end RHF workflow:
  1. Build H2 geometry
  2. Compute STO-3G molecular integrals
  3. Run RHF SCF
  4. Print results and compare with PySCF

Output:
  01_h2_rhf.png          — convergence + orbital energy + timing bar chart
  01_h2_density_atv.cube — electron density from Atomiverse P
  01_h2_density_ref.cube — electron density from PySCF P

Reference: Szabo & Ostlund Table 3.7 -> E = -1.1175 Ha (STO-3G)
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

# ── Geometry ─────────────────────────────────────────────────────────────────
BOND_LENGTH = 0.74  # Angstrom (H-H equilibrium)


# ── Atomiverse ───────────────────────────────────────────────────────────────

def make_h2_frame(bond: float = BOND_LENGTH) -> av.CpuFrame:
    """Build H2 cpu::Frame, matching the C++ tutorial."""
    frame = av.CpuFrame(2)
    dtype = frame.get_x().dtype  # detect Real precision
    frame.set_positions(
        np.array([0.0, bond], dtype=dtype),
        np.zeros(2, dtype=dtype),
        np.zeros(2, dtype=dtype),
    )
    frame.set_mass(np.array([1.008, 1.008], dtype=dtype))
    frame.set_type(np.array([1, 1], dtype=np.int32))
    frame.set_box(np.eye(3, dtype=dtype) * 100.0)
    return frame


@timed(label="Atomiverse integrals")
def atv_integrals(frame: av.CpuFrame) -> av.IntegralPack:
    return av.build_integrals_sto3g(frame)


@timed(label="Atomiverse RHF")
def atv_rhf(ip: av.IntegralPack) -> av.RHFResult:
    return av.rhf_solve(ip, tol=1e-10, max_iter=100)


# ── PySCF ────────────────────────────────────────────────────────────────────

def make_pyscf_mol(bond: float = BOND_LENGTH) -> gto.Mole:
    return gto.M(
        atom=f"H 0 0 0; H {bond} 0 0",
        basis="sto-3g",
        unit="angstrom",
        verbose=0,
    )


@timed(label="PySCF RHF")
def pyscf_rhf(mol: gto.Mole) -> tuple[scf.RHF, list[dict]]:
    """Run PySCF RHF and capture per-iteration convergence trace."""
    trace: list[dict] = []

    def callback(envs: dict) -> None:
        trace.append({
            "iter": envs["cycle"] + 1,
            "energy": envs["e_tot"],
            "residual": envs.get("norm_ddm", float("nan")),
        })

    mf = scf.RHF(mol)
    mf.conv_tol = 1e-12
    mf.callback = callback
    mf.kernel()
    return mf, trace


# ── Cube file ────────────────────────────────────────────────────────────────

def write_density_cube(
    P: np.ndarray,
    mol: gto.Mole,
    filename: str,
    label: str,
) -> None:
    """Write electron density cube file via PySCF cubegen.

    Uses the PySCF Mole for AO evaluation; the density matrix P can come
    from either Atomiverse or PySCF.
    """
    cubegen.density(mol, filename, P.astype(np.float64),
                    nx=60, ny=60, nz=60)
    print(f"  {filename}: density cube written [{label}]")


# ── Plot ─────────────────────────────────────────────────────────────────────

def plot(
    atv_result: av.RHFResult,
    atv_ip: av.IntegralPack,
    mf: scf.RHF,
    pyscf_trace: list[dict],
    outfile: str = "01_h2_rhf.png",
) -> None:
    fig = plt.figure(figsize=(14, 4))
    gs = gridspec.GridSpec(1, 3, figure=fig, wspace=0.4)

    # --- Panel 1: SCF convergence ---
    ax1 = fig.add_subplot(gs[0])
    atv_trace = atv_result.trace
    atv_iters = [t["iter"] for t in atv_trace]
    atv_resid = [t["residual"] for t in atv_trace]
    ax1.semilogy(atv_iters, atv_resid, "-o", color="#0066CC",
                 markersize=5, label="Atomiverse")
    if pyscf_trace:
        ref_iters = [t["iter"] for t in pyscf_trace]
        ref_resid = [t["residual"] for t in pyscf_trace]
        ax1.semilogy(ref_iters, ref_resid, "-s", color="#CC4400",
                     markersize=5, label="PySCF")
    ax1.set_xlabel("SCF iteration")
    ax1.set_ylabel("||dP||_F")
    ax1.set_title("SCF Convergence")
    ax1.legend()

    # --- Panel 2: orbital energies ---
    ax2 = fig.add_subplot(gs[1])
    eig_atv = np.array(atv_result.get_eigenvalues())
    eig_ref = mf.mo_energy[:len(eig_atv)]
    n_mo = len(eig_atv)
    x = np.arange(n_mo)
    w = 0.3
    ax2.bar(x - w / 2, eig_atv, w, label="Atomiverse", color="#0066CC")
    ax2.bar(x + w / 2, eig_ref, w, label="PySCF", color="#CC4400")
    ax2.set_xticks(x)
    ax2.set_xticklabels([f"e{i+1}" for i in range(n_mo)])
    ax2.set_ylabel("Energy [Hartree]")
    ax2.set_title("Orbital Energies")
    ax2.legend()

    # --- Panel 3: timing ---
    ax3 = fig.add_subplot(gs[2])
    labels = ["Integrals", "RHF solve"]
    atv_times = [atv_integrals.elapsed_ms, atv_rhf.elapsed_ms]
    pyscf_times = [0.0, pyscf_rhf.elapsed_ms]  # PySCF bundles integrals in kernel()
    x = np.arange(len(labels))
    ax3.bar(x - w / 2, atv_times, w, label="Atomiverse", color="#0066CC")
    ax3.bar(x + w / 2, pyscf_times, w, label="PySCF", color="#CC4400")
    ax3.set_xticks(x)
    ax3.set_xticklabels(labels)
    ax3.set_ylabel("Time [ms]")
    ax3.set_title("Wall-clock Time")
    ax3.legend()

    e_atv = float(atv_result.energy)
    e_ref = float(mf.e_tot)
    dE = abs(e_atv - e_ref)
    fig.suptitle(
        f"H2 / STO-3G / RHF  —  "
        f"Atomiverse E={e_atv:.8f}  PySCF E={e_ref:.8f}  "
        f"|dE|={dE:.2e} Ha",
        fontsize=9,
    )
    fig.savefig(outfile, dpi=150, bbox_inches="tight")
    print(f"\nSaved {outfile}")
    plt.close(fig)


# ── Main ─────────────────────────────────────────────────────────────────────

def main() -> None:
    print("=" * 60)
    print("H2 / STO-3G / RHF  —  Atomiverse vs PySCF")
    print("=" * 60)
    print(f"Geometry: H-H at {BOND_LENGTH} Angstrom\n")

    # ── 1. Build H2 geometry ──
    frame = make_h2_frame()

    # ── 2. Compute molecular integrals ──
    ip = atv_integrals(frame)
    print(f"Basis functions : {ip.n_basis}")
    print(f"Electrons       : {ip.n_electrons}")
    print(f"Occupied MOs    : {ip.n_occ} (RHF closed-shell)")
    print(f"V_nuc           : {float(ip.V_nuc):.10f} Hartree")

    S = np.array(ip.get_S())
    Hcore = np.array(ip.get_Hcore())
    ERI = np.array(ip.get_ERI())

    print(f"\nOverlap S:\n{S}")
    print(f"\nCore Hamiltonian Hcore = T + V:\n{Hcore}")
    print(f"\nKey ERI (chemist's notation):")
    print(f"  (11|11) = {ERI[0,0,0,0]:.6f}")
    print(f"  (11|22) = {ERI[0,0,1,1]:.6f}")
    print(f"  (12|12) = {ERI[0,1,0,1]:.6f}")

    # ── 3. Run Atomiverse RHF ──
    result = atv_rhf(ip)

    print(f"\n{'='*60}")
    print("Atomiverse RHF Result")
    print(f"{'='*60}")
    print(f"Converged       : {'yes' if result.converged else 'no'}"
          f" ({result.iterations} iterations)")
    print(f"E_electronic    : {float(result.electronic_energy):.10f} Hartree")
    print(f"V_nuclear       : {float(ip.V_nuc):.10f} Hartree")
    print(f"E_total         : {float(result.energy):.10f} Hartree")

    eig = np.array(result.get_eigenvalues())
    print(f"\nOrbital energies:")
    for i, e in enumerate(eig):
        occ = "  (occupied)" if i < ip.n_occ else ""
        print(f"  epsilon[{i}] = {e:12.6f} Hartree{occ}")

    P_atv = np.array(result.get_P())
    tr_PS = np.sum(P_atv * S)
    print(f"\nDensity check: Tr(P*S) = {tr_PS:.6f} (expected: {ip.n_electrons})")
    print(f"Timing: integrals {atv_integrals.elapsed_ms:.2f} ms, "
          f"RHF {atv_rhf.elapsed_ms:.2f} ms")

    # ── 4. PySCF reference ──
    print(f"\n{'='*60}")
    print("PySCF RHF Result")
    print(f"{'='*60}")
    mol = make_pyscf_mol()
    mf, pyscf_trace = pyscf_rhf(mol)
    P_ref = mf.make_rdm1()
    S_ref = mol.intor("int1e_ovlp")

    print(f"Converged       : {'yes' if mf.converged else 'no'}")
    print(f"E_electronic    : {mf.e_tot - mol.energy_nuc():.10f} Hartree")
    print(f"V_nuclear       : {mol.energy_nuc():.10f} Hartree")
    print(f"E_total         : {mf.e_tot:.10f} Hartree")

    print(f"\nOrbital energies:")
    for i, e in enumerate(mf.mo_energy):
        occ = "  (occupied)" if i < mol.nelectron // 2 else ""
        print(f"  epsilon[{i}] = {e:12.6f} Hartree{occ}")

    tr_PS_ref = np.trace(P_ref @ S_ref)
    print(f"\nDensity check: Tr(P*S) = {tr_PS_ref:.6f} (expected: {mol.nelectron})")
    print(f"Timing: {pyscf_rhf.elapsed_ms:.2f} ms")

    # ── 5. Convergence trace comparison ──
    print(f"\n{'='*60}")
    print("SCF Convergence Trace")
    print(f"{'='*60}")
    atv_trace = result.trace
    n_rows = max(len(atv_trace), len(pyscf_trace))
    print(f"  {'iter':>4s}  {'Atomiverse ||dP||':>20s}  {'PySCF ||dP||':>20s}")
    print(f"  {'----':>4s}  {'--------------------':>20s}  {'--------------------':>20s}")
    for i in range(n_rows):
        atv_val = f"{atv_trace[i]['residual']:.6e}" if i < len(atv_trace) else ""
        ref_val = f"{pyscf_trace[i]['residual']:.6e}" if i < len(pyscf_trace) else ""
        it = i + 1
        print(f"  {it:4d}  {atv_val:>20s}  {ref_val:>20s}")

    # Final residual comparison
    atv_final = atv_trace[-1]["residual"] if atv_trace else float("nan")
    ref_final = pyscf_trace[-1]["residual"] if pyscf_trace else float("nan")
    print(f"\n  Final residual — Atomiverse: {atv_final:.6e}, PySCF: {ref_final:.6e}")

    # ── 6. Agreement ──
    e_atv = float(result.energy)
    e_ref = float(mf.e_tot)
    dE = abs(e_atv - e_ref)
    rel_err = abs(dE / e_ref) * 100

    dP_max = np.max(np.abs(P_atv.astype(np.float64) - P_ref))

    eig_ref = mf.mo_energy[:len(eig)]
    d_eig_max = np.max(np.abs(eig.astype(np.float64) - eig_ref))

    print(f"\n{'='*60}")
    print("Agreement (Atomiverse vs PySCF)")
    print(f"{'='*60}")
    print(f"  |dE_total|       = {dE:.2e} Ha  ({rel_err:.4f}%)")
    print(f"  max|dP|          = {dP_max:.2e}")
    print(f"  max|d_epsilon|   = {d_eig_max:.2e} Ha")
    print(f"  Textbook (S&O)   = -1.1175 Ha")

    THRESHOLD = 0.01  # 1%
    passed = rel_err < THRESHOLD * 100
    status = "PASS" if passed else "FAIL"
    print(f"\n  Result: {status}  (relative error {rel_err:.4f}% vs {THRESHOLD*100}% threshold)")

    # ── 7. Cube files ──
    print(f"\n{'='*60}")
    print("Electron density cube files")
    print(f"{'='*60}")
    write_density_cube(P_atv, mol, "01_h2_density_atv.cube", "Atomiverse")
    write_density_cube(P_ref, mol, "01_h2_density_ref.cube", "PySCF")

    # ── 8. Plot ──
    plot(result, ip, mf, pyscf_trace)

    if not passed:
        sys.exit(1)

    print(f"\n{'='*60}")
    print("Done")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
