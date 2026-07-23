#!/usr/bin/env python3
"""Offline PySCF pin helper for Atomiverse regressions.

NEVER imported by CMake/ctest. Developers / autonomous agents run this,
then paste numbers into REFERENCE.md and hard-coded regression asserts.

Examples:
  python tools/pyscf_oracle/pin_scf.py rhf --geom water --basis sto-3g
  python tools/pyscf_oracle/pin_scf.py uhf --geom o2 --basis sto-3g --spin 2
  python tools/pyscf_oracle/pin_scf.py ump2 --geom o2 --basis sto-3g --spin 2
  python tools/pyscf_oracle/pin_scf.py mp2-fc --geom water --basis 6-31g* --frozen 1
  python tools/pyscf_oracle/pin_scf.py rks --geom water --basis 6-31g* --xc b3lyp
  python tools/pyscf_oracle/pin_scf.py tda --geom h2co --basis 6-31g* --nstates 3
  python tools/pyscf_oracle/pin_scf.py cisd --geom h2co --basis sto-3g --nstates 3
  python tools/pyscf_oracle/pin_scf.py rks-tda --geom h2 --basis sto-3g --xc b3lyp --nstates 3
  python tools/pyscf_oracle/pin_scf.py rks-tddft --geom h2 --basis sto-3g --xc b3lyp --nstates 3
"""
from __future__ import annotations

import argparse
import json
import sys

# CODATA Ha→eV (same constant Atomiverse units::kHartreeInEv uses for gates).
_HA_TO_EV = 27.211386245988


def _mol(geom: str, basis: str, charge: int, spin: int, cart: bool | None):
    from pyscf import gto

    presets = {
        "water": """
O  0.000000  0.000000  0.117300
H  0.000000  0.757200 -0.469200
H  0.000000 -0.757200 -0.469200
""",
        "h2": """
H  0.000000  0.000000  0.000000
H  0.000000  0.000000  0.740000
""",
        "o2": """
O  0.000000  0.000000  0.000000
O  0.000000  0.000000  1.207500
""",
        "h2co": """
C  0.000000  0.000000  0.000000
O  0.000000  0.000000  1.210000
H  0.000000  0.943000 -0.587000
H  0.000000 -0.943000 -0.587000
""",
        "ch4": """
C  0.000000  0.000000  0.000000
H  0.629118  0.629118  0.629118
H -0.629118 -0.629118  0.629118
H  0.629118 -0.629118 -0.629118
H -0.629118  0.629118 -0.629118
""",
    }
    if geom not in presets:
        raise SystemExit(f"unknown geom preset {geom!r}; choose {sorted(presets)}")

    # Pople * bases: Cartesian d is Atomiverse default; pure for def2.
    if cart is None:
        cart = not basis.lower().startswith("def2")

    return gto.M(
        atom=presets[geom],
        basis=basis,
        charge=charge,
        spin=spin,
        unit="Angstrom",
        cart=cart,
        verbose=0,
    )


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "method",
        choices=[
            "rhf",
            "uhf",
            "rohf",
            "rks",
            "uks",
            "mp2",
            "mp2-fc",
            "ump2",
            "tda",
            "tdhf",
            "cisd",
            "rks-tda",
            "rks-tddft",
        ],
    )
    p.add_argument("--geom", default="water")
    p.add_argument("--basis", default="sto-3g")
    p.add_argument("--charge", type=int, default=0)
    p.add_argument("--spin", type=int, default=0, help="2S (PySCF spin)")
    p.add_argument("--xc", default="b3lyp")
    p.add_argument("--nstates", type=int, default=3)
    p.add_argument("--conv-tol", type=float, default=1e-12)
    p.add_argument(
        "--frozen",
        type=int,
        default=None,
        help="Number of lowest occupied MOs to freeze (mp2-fc default: 1 = O 1s on water)",
    )
    p.add_argument("--cart", action=argparse.BooleanOptionalAction, default=None)
    p.add_argument("--grad", action="store_true")
    args = p.parse_args()

    try:
        import pyscf
        from pyscf import scf, dft, mp, tdscf
    except ImportError as e:
        print(json.dumps({"ok": False, "error": f"pyscf missing: {e}"}))
        return 2

    mol = _mol(args.geom, args.basis, args.charge, args.spin, args.cart)
    out: dict = {
        "ok": True,
        "pyscf": pyscf.__version__,
        "method": args.method,
        "geom": args.geom,
        "basis": args.basis,
        "charge": args.charge,
        "spin": args.spin,
        "cart": bool(mol.cart),
        "conv_tol": args.conv_tol,
        "natm": mol.natm,
        "nao": mol.nao,
    }

    m = args.method
    if m == "rhf":
        mf = scf.RHF(mol)
    elif m in ("uhf", "ump2"):
        mf = scf.UHF(mol)
    elif m == "rohf":
        mf = scf.ROHF(mol)
    elif m in ("rks", "rks-tda", "rks-tddft"):
        mf = dft.RKS(mol)
        mf.xc = args.xc
        out["xc"] = args.xc
    elif m == "uks":
        mf = dft.UKS(mol)
        mf.xc = args.xc
        out["xc"] = args.xc
    else:
        # tda / tdhf / mp2 ground state: RHF
        mf = scf.RHF(mol)

    mf.conv_tol = args.conv_tol
    e = mf.kernel()
    out["e_tot_ha"] = float(e)
    out["converged"] = bool(mf.converged)

    if args.grad and m in ("rhf", "uhf", "rohf", "rks", "uks"):
        g = mf.nuc_grad_method().kernel()
        out["grad_ha_bohr"] = g.tolist()

    if m == "mp2":
        m2 = mp.MP2(mf).run()
        out["e_corr_ha"] = float(m2.e_corr)
        out["e_tot_mp2_ha"] = float(m2.e_tot)
    elif m == "mp2-fc":
        # Frozen-core RMP2: freeze lowest N occupied MOs (default 1 = O 1s on H2O).
        n_frozen = 1 if args.frozen is None else int(args.frozen)
        m2 = mp.MP2(mf, frozen=n_frozen).run()
        out["frozen"] = n_frozen
        out["e_hf_ha"] = float(e)
        out["e_corr_ha"] = float(m2.e_corr)
        out["e_tot_mp2_ha"] = float(m2.e_tot)
    elif m == "ump2":
        # Explicit UMP2 on the UHF reference (triplet O2, etc.).
        m2 = mp.UMP2(mf).run()
        out["e_hf_ha"] = float(e)
        out["e_corr_ha"] = float(m2.e_corr)
        out["e_tot_mp2_ha"] = float(m2.e_tot)
        if hasattr(mf, "spin_square"):
            ss, mult = mf.spin_square()
            out["s2"] = float(ss)
            out["multiplicity"] = float(mult)
    elif m in ("tda", "tdhf", "rks-tda", "rks-tddft"):
        # rks-tda: hybrid/pure TDDFT-TDA on the KS reference (exact-exchange
        # kernel for hybrids via PySCF's Casida A matrix).
        # rks-tddft: full RPA TDDFT (Casida A,B) — tdscf.TDDFT, not TDA.
        if m == "tdhf":
            td = tdscf.TDHF(mf)
        elif m == "rks-tddft":
            td = tdscf.TDDFT(mf)
        else:
            td = tdscf.TDA(mf)
        td.nstates = args.nstates
        td.kernel()
        out["excitation_ha"] = [float(x) for x in td.e]
        out["excitation_ev"] = [float(x) * _HA_TO_EV for x in td.e]
        if hasattr(td, "converged"):
            out["td_converged"] = [bool(c) for c in td.converged]
        if hasattr(td, "oscillator_strength"):
            out["oscillator_strength"] = [float(x) for x in td.oscillator_strength()]
    elif m == "cisd":
        # CIS(D): Head-Gordon et al., Chem. Phys. Lett. 219, 21 (1994).
        # PySCF has no product CIS(D); evaluate spin-orbital eqs. (8)/(13)/(14)
        # on pyscf.tdscf.TDA amplitudes + MO ERIs (offline pin only).
        from pyscf import ao2mo
        import numpy as np

        td = tdscf.TDA(mf)
        td.nstates = args.nstates
        td.conv_tol = min(args.conv_tol, 1e-10)
        td.kernel()
        eps = mf.mo_energy
        C = mf.mo_coeff
        nmo = C.shape[1]
        nocc = mol.nelectron // 2
        nvir = nmo - nocc
        eri_mo = ao2mo.restore(1, ao2mo.full(mf._eri, C), nmo)
        occ_so = [(i, 0) for i in range(nocc)] + [(i, 1) for i in range(nocc)]
        vir_so = [(a, 0) for a in range(nocc, nmo)] + [
            (a, 1) for a in range(nocc, nmo)
        ]
        nso_o, nso_v = len(occ_so), len(vir_so)

        def so_eri(p, q, r, s):
            (P, sp), (Q, sq), (R, sr), (S, ss) = p, q, r, s
            val = 0.0
            if sp == sr and sq == ss:
                val += eri_mo[P, R, Q, S]
            if sp == ss and sq == sr:
                val -= eri_mo[P, S, Q, R]
            return val

        a_amp = np.zeros((nso_o, nso_o, nso_v, nso_v))
        for ii, i in enumerate(occ_so):
            for jj, j in enumerate(occ_so):
                for aa, a in enumerate(vir_so):
                    for bb, b in enumerate(vir_so):
                        denom = (
                            eps[a[0]] + eps[b[0]] - eps[i[0]] - eps[j[0]]
                        )
                        a_amp[ii, jj, aa, bb] = -so_eri(a, b, i, j) / denom

        roots = []
        for s in range(args.nstates):
            X = td.xy[s][0]
            bvec = np.zeros((nso_o, nso_v))
            for i in range(nocc):
                for a in range(nvir):
                    bvec[i, a] = X[i, a]
                    bvec[nocc + i, nvir + a] = X[i, a]
            omega = float(td.e[s])
            u = np.zeros((nso_o, nso_o, nso_v, nso_v))
            for ii, i in enumerate(occ_so):
                for jj, j in enumerate(occ_so):
                    for aa, a in enumerate(vir_so):
                        for bb, bv in enumerate(vir_so):
                            ssum = 0.0
                            for cc, c in enumerate(vir_so):
                                ssum += so_eri(a, bv, c, j) * bvec[ii, cc]
                                ssum -= so_eri(a, bv, c, i) * bvec[jj, cc]
                            for kk, k in enumerate(occ_so):
                                ssum += so_eri(k, a, i, j) * bvec[kk, bb]
                                ssum -= so_eri(k, bv, i, j) * bvec[kk, aa]
                            u[ii, jj, aa, bb] = ssum
            V = 0.0
            for ii, i in enumerate(occ_so):
                for jj, j in enumerate(occ_so):
                    for aa, a in enumerate(vir_so):
                        for bb, bv in enumerate(vir_so):
                            denom = (
                                eps[a[0]]
                                + eps[bv[0]]
                                - eps[i[0]]
                                - eps[j[0]]
                                - omega
                            )
                            V += u[ii, jj, aa, bb] ** 2 / denom
            V *= -0.25
            v = np.zeros((nso_o, nso_v))
            for ii in range(nso_o):
                for aa in range(nso_v):
                    ssum = 0.0
                    for jj in range(nso_o):
                        for kk in range(nso_o):
                            for bb in range(nso_v):
                                for cc in range(nso_v):
                                    jklbc = so_eri(
                                        occ_so[jj],
                                        occ_so[kk],
                                        vir_so[bb],
                                        vir_so[cc],
                                    )
                                    ssum += jklbc * (
                                        bvec[ii, bb] * a_amp[jj, kk, cc, aa]
                                        + bvec[jj, aa] * a_amp[ii, kk, cc, bb]
                                        + 2.0
                                        * bvec[jj, bb]
                                        * a_amp[ii, kk, aa, cc]
                                    )
                    v[ii, aa] = 0.5 * ssum
            U = float(np.sum(bvec * v))
            roots.append(
                {
                    "cis_ha": omega,
                    "cis_ev": omega * _HA_TO_EV,
                    "U_ha": U,
                    "V_ha": V,
                    "corr_ha": U + V,
                    "cisd_ha": omega + U + V,
                    "cisd_ev": (omega + U + V) * _HA_TO_EV,
                }
            )
        out["roots"] = roots
        out["excitation_cis_ev"] = [r["cis_ev"] for r in roots]
        out["excitation_cisd_ev"] = [r["cisd_ev"] for r in roots]
        if hasattr(td, "converged"):
            out["td_converged"] = [bool(c) for c in td.converged]

    print(json.dumps(out, indent=2, sort_keys=True))
    return 0 if out.get("converged", True) else 1


if __name__ == "__main__":
    sys.exit(main())
