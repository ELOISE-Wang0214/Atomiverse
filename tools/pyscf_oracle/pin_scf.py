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
  python tools/pyscf_oracle/pin_scf.py rks-tda --geom h2 --basis sto-3g --xc b3lyp --nstates 3
  python tools/pyscf_oracle/pin_scf.py rks-tddft --geom h2 --basis sto-3g --xc b3lyp --nstates 3
  python tools/pyscf_oracle/pin_scf.py ghf --geom h3 --basis sto-3g --spin 1 --grad
"""
from __future__ import annotations

import argparse
import json
import sys

# CODATA Ha→eV (same constant Atomiverse units::kHartreeInEv uses for gates).
_HA_TO_EV = 27.211386245988


def _mol(geom: str, basis: str, charge: int, spin: int, cart: bool | None):
    from pyscf import gto
    import numpy as np

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
        # Equilateral H3, side 1.0 Å (GHF spin-frustration geometry).
        "h3": f"""
H  0.000000  0.000000  0.000000
H  1.000000  0.000000  0.000000
H  0.500000  {0.5 * np.sqrt(3):.12f}  0.000000
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
            "ghf",
            "rks",
            "uks",
            "mp2",
            "mp2-fc",
            "ump2",
            "tda",
            "tdhf",
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
    elif m == "ghf":
        mf = scf.GHF(mol)
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
    # GHF on frustrated H3 needs a noncollinear Néel guess (collinear traps).
    if m == "ghf" and args.geom == "h3":
        import numpy as np

        n = mol.nao
        dm0 = np.zeros((2 * n, 2 * n))
        dm0[:n, :n] = np.eye(n)

        def _rot(dm, theta):
            G = np.zeros((2 * n, 2 * n))
            for mu in range(n):
                th = theta * mu
                c, s = np.cos(th / 2), np.sin(th / 2)
                G[mu, mu] = c
                G[mu, n + mu] = -s
                G[n + mu, mu] = s
                G[n + mu, n + mu] = c
            return G @ dm @ G.T

        e = mf.kernel(dm0=_rot(dm0, 2.0 * np.pi / 3.0))
    else:
        e = mf.kernel()
    out["e_tot_ha"] = float(e)
    out["converged"] = bool(mf.converged)

    if args.grad and m in ("rhf", "uhf", "rohf", "rks", "uks"):
        g = mf.nuc_grad_method().kernel()
        out["grad_ha_bohr"] = g.tolist()
    elif args.grad and m == "ghf":
        # PySCF has no scf.GHF.nuc_grad_method; pin analytic dE/dR at the
        # converged spinor density with central-FD dAO integrals (Ha/Bohr).
        import numpy as np

        n = mol.nao
        dm = mf.make_rdm1()
        Paa, Pab = dm[:n, :n], dm[:n, n:]
        Pba, Pbb = dm[n:, :n], dm[n:, n:]
        Pc = Paa + Pbb
        mo, e_mo, occ = mf.mo_coeff, mf.mo_energy, mf.mo_occ
        W = np.zeros_like(dm)
        for i in range(mo.shape[1]):
            if occ[i] < 0.5:
                continue
            W += e_mo[i] * np.outer(mo[:, i], mo[:, i]) * occ[i]
        Wc = W[:n, :n] + W[n:, n:]
        h = 1e-5
        coords0 = mol.atom_coords(unit="Bohr").copy()
        grad = np.zeros((mol.natm, 3))
        for A in range(mol.natm):
            for d in range(3):

                def _mats(sign, _A=A, _d=d):
                    coords = coords0.copy()
                    coords[_A, _d] += sign * h
                    m2 = mol.copy()
                    m2.set_geom_(coords, unit="Bohr")
                    m2.build(0, 0)
                    S = m2.intor("int1e_ovlp")
                    H = m2.intor("int1e_kin") + m2.intor("int1e_nuc")
                    eri = m2.intor("int2e", aosym="s1").reshape(n, n, n, n)
                    return S, H, eri, m2.energy_nuc()

                Sp, Hp, ep, Vp = _mats(+1)
                Sm, Hm, em, Vm = _mats(-1)
                dS = (Sp - Sm) / (2 * h)
                dH = (Hp - Hm) / (2 * h)
                deri = (ep - em) / (2 * h)
                dVn = (Vp - Vm) / (2 * h)
                gA = float(np.sum(Pc * dH) - np.sum(Wc * dS) + dVn)
                gA += 0.5 * float(np.einsum("ij,kl,ijkl->", Pc, Pc, deri))
                for D, E in ((Paa, Paa), (Pbb, Pbb), (Pab, Pba), (Pba, Pab)):
                    gA += -0.5 * float(np.einsum("jk,li,ijkl->", D, E, deri))
                grad[A, d] = gA
        out["grad_ha_bohr"] = grad.tolist()
        out["grad_note"] = (
            "analytic GHF at fixed dm with FD dAO (PySCF has no GHF nuc_grad)"
        )

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

    print(json.dumps(out, indent=2, sort_keys=True))
    return 0 if out.get("converged", True) else 1


if __name__ == "__main__":
    sys.exit(main())
