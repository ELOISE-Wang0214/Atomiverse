#!/usr/bin/env python3
"""Offline PySCF pin helper for Atomiverse regressions.

NEVER imported by CMake/ctest. Developers / autonomous agents run this,
then paste numbers into REFERENCE.md and hard-coded regression asserts.

Examples:
  python tools/pyscf_oracle/pin_scf.py rhf --geom water --basis sto-3g
  python tools/pyscf_oracle/pin_scf.py uhf --geom o2 --basis sto-3g --spin 2
  python tools/pyscf_oracle/pin_scf.py ump2 --geom o2 --basis sto-3g --spin 2
  python tools/pyscf_oracle/pin_scf.py romp2 --geom o2 --basis sto-3g --spin 2
  python tools/pyscf_oracle/pin_scf.py mp2-fc --geom water --basis 6-31g* --frozen 1
  python tools/pyscf_oracle/pin_scf.py scs-mp2 --geom water --basis sto-3g
  python tools/pyscf_oracle/pin_scf.py sos-mp2 --geom water --basis sto-3g
  python tools/pyscf_oracle/pin_scf.py mp2 --geom water --basis sto-3g --grad
  python tools/pyscf_oracle/pin_scf.py mp3 --geom h2 --basis sto-3g
  python tools/pyscf_oracle/pin_scf.py ccsd --geom h2 --basis sto-3g
  python tools/pyscf_oracle/pin_scf.py b2plyp --geom water --basis sto-3g
  python tools/pyscf_oracle/pin_scf.py rks --geom water --basis 6-31g* --xc b3lyp
  python tools/pyscf_oracle/pin_scf.py roks --geom o2 --basis sto-3g --spin 2 --xc LDA,VWN
  python tools/pyscf_oracle/pin_scf.py tda --geom h2co --basis 6-31g* --nstates 3
  python tools/pyscf_oracle/pin_scf.py cisd --geom h2co --basis sto-3g --nstates 3
  python tools/pyscf_oracle/pin_scf.py adc2 --geom h2co --basis sto-3g --nstates 3
  python tools/pyscf_oracle/pin_scf.py rks-tda --geom h2 --basis sto-3g --xc b3lyp --nstates 3
  python tools/pyscf_oracle/pin_scf.py rks-tddft --geom h2 --basis sto-3g --xc b3lyp --nstates 3
  python tools/pyscf_oracle/pin_scf.py uks-tda --geom o2 --basis sto-3g --spin 2 --xc LDA,VWN --nstates 5
  python tools/pyscf_oracle/pin_scf.py sf-tda --geom o2 --basis sto-3g --spin 2 --nstates 5
  python tools/pyscf_oracle/pin_scf.py d3bj --geom water --basis sto-3g --xc pbe
  python tools/pyscf_oracle/pin_scf.py g0w0 --geom h2 --basis sto-3g
  python tools/pyscf_oracle/pin_scf.py delta-scf --geom n2 --basis sto-3g
  python tools/pyscf_oracle/pin_scf.py rhf --geom water --basis 6-31g* --grad --pcm
"""
from __future__ import annotations

import argparse
import json
import sys

import numpy as np

# CODATA Ha→eV (same constant Atomiverse units::kHartreeInEv uses for gates).
_HA_TO_EV = 27.211386245988


def _as_float(x) -> float:
    """Coerce scalar / 0-d / length-1 array to float (PySCF dispersion returns)."""
    if isinstance(x, dict):
        if "energy" in x:
            return _as_float(x["energy"])
        raise TypeError(f"dict without 'energy': {x!r}")
    arr = np.asarray(x)
    if arr.ndim == 0:
        return float(arr)
    return float(arr.ravel()[0])


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
        "n2": """
N  0.000000  0.000000  0.000000
N  0.000000  0.000000  1.097700
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
        choices=["rhf", "uhf", "rohf", "rks", "uks", "roks", "mp2", "mp2-fc", "mp3", "scs-mp2", "sos-mp2", "ump2", "romp2", "tda", "tdhf", "cisd", "adc2", "rks-tda", "rks-tddft", "uks-tda", "sf-tda", "d3bj", "g0w0", "ccsd", "b2plyp", "delta-scf"],
    )
    p.add_argument("--geom", default="water")
    p.add_argument("--basis", default="sto-3g")
    p.add_argument("--charge", type=int, default=0)
    p.add_argument("--spin", type=int, default=0, help="2S (PySCF spin)")
    p.add_argument("--xc", default="b3lyp")
    p.add_argument("--nstates", type=int, default=3)
    p.add_argument(
        "--from-mo",
        type=int,
        default=None,
        help="Occupied MO index to depopulate for delta-scf (default: HOMO)",
    )
    p.add_argument(
        "--to-mo",
        type=int,
        default=None,
        help="Virtual MO index to populate for delta-scf (default: LUMO)",
    )
    p.add_argument("--conv-tol", type=float, default=1e-12)
    p.add_argument(
        "--frozen",
        type=int,
        default=None,
        help="Number of lowest occupied MOs to freeze (mp2-fc default: 1 = O 1s on water)",
    )
    p.add_argument("--cart", action=argparse.BooleanOptionalAction, default=None)
    p.add_argument("--grad", action="store_true")
    p.add_argument(
        "--pcm",
        action="store_true",
        help="Attach PySCF solvent.PCM (C-PCM) before SCF / gradient",
    )
    p.add_argument("--pcm-eps", type=float, default=78.4, help="PCM dielectric ε")
    p.add_argument(
        "--pcm-scale",
        type=float,
        default=1.2,
        help="PCM vdw radii scale (Bondi × scale)",
    )
    p.add_argument(
        "--pcm-lebedev-order",
        type=int,
        default=29,
        help="Lebedev angular order for PCM cavity",
    )
    p.add_argument(
        "--pcm-method",
        default="C-PCM",
        help="PySCF solvent method string (default C-PCM)",
    )
    args = p.parse_args()

    try:
        import pyscf
        from pyscf import scf, dft, mp, tdscf, cc as pyscf_cc
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
    if m == "d3bj":
        # Geometry-only DFT-D3(BJ). Requires pyscf-dispersion (simple-dftd3).
        # E_disp is independent of SCF; still report optional e_tot with/without.
        try:
            from pyscf.dispersion import dftd3 as dftd3_mod
        except ImportError as e:
            print(
                json.dumps(
                    {
                        "ok": False,
                        "error": f"pyscf-dispersion missing (need dftd3): {e}",
                    }
                )
            )
            return 2
        xc = args.xc if args.xc else "pbe"
        out["xc"] = xc
        model = dftd3_mod.DFTD3Dispersion(mol, xc=xc, version="d3bj")
        ed = _as_float(model.get_dispersion())
        out["e_disp_ha"] = ed
        # Optional KS totals for cross-check (same geometry / basis / xc).
        mf0 = dft.RKS(mol)
        mf0.xc = xc
        mf0.conv_tol = args.conv_tol
        e0 = float(mf0.kernel())
        mf1 = dft.RKS(mol)
        mf1.xc = xc
        mf1.conv_tol = args.conv_tol
        mf1.disp = "d3bj"
        e1 = float(mf1.kernel())
        out["e_tot_no_disp_ha"] = e0
        out["e_tot_with_disp_ha"] = e1
        out["e_tot_ha"] = e1
        out["converged"] = bool(mf1.converged)
        print(json.dumps(out, indent=2, sort_keys=True))
        return 0

    if m == "delta-scf":
        # UHF ΔSCF with IMOM occupation locking (PySCF scf.addons.mom_occ_).
        # Default promotion: α HOMO → LUMO on a closed-shell RHF reference.
        import numpy as np
        from pyscf.scf import addons

        mf0 = scf.RHF(mol)
        mf0.conv_tol = args.conv_tol
        e0 = float(mf0.kernel())
        nmo = int(mf0.mo_coeff.shape[1])
        nocc = int(mol.nelectron // 2)
        from_mo = int(args.from_mo) if args.from_mo is not None else nocc - 1
        to_mo = int(args.to_mo) if args.to_mo is not None else nocc
        if from_mo < 0 or from_mo >= nocc or to_mo < nocc or to_mo >= nmo:
            print(
                json.dumps(
                    {
                        "ok": False,
                        "error": f"bad promotion from={from_mo} to={to_mo} "
                        f"(nocc={nocc}, nmo={nmo})",
                    }
                )
            )
            return 2
        mo0 = np.asarray(mf0.mo_coeff)
        occ_a = np.zeros(nmo)
        occ_b = np.zeros(nmo)
        occ_a[:nocc] = 1.0
        occ_b[:nocc] = 1.0
        occ_a[from_mo] = 0.0
        occ_a[to_mo] = 1.0

        def _rdm1(C, occ):
            return np.einsum("pi,i,qi->pq", C, occ, C)

        dm = (_rdm1(mo0, occ_a), _rdm1(mo0, occ_b))
        mfu = scf.UHF(mol)
        mfu.conv_tol = args.conv_tol
        mfu.max_cycle = 200
        mfu = addons.mom_occ_(mfu, (mo0, mo0), (occ_a.copy(), occ_b.copy()))
        e1 = float(mfu.kernel(dm))
        de_ha = e1 - e0
        out["e_ground_ha"] = e0
        out["e_excited_ha"] = e1
        out["delta_e_ha"] = de_ha
        out["delta_e_ev"] = de_ha * _HA_TO_EV
        out["from_mo"] = from_mo
        out["to_mo"] = to_mo
        out["nocc"] = nocc
        out["nmo"] = nmo
        out["converged"] = bool(mfu.converged)
        if hasattr(mfu, "spin_square"):
            ss, mult = mfu.spin_square()
            out["s2"] = float(ss)
            out["multiplicity"] = float(mult)
        out["oracle"] = "pyscf_scf_addons_mom_occ_imom"
        print(json.dumps(out, indent=2, sort_keys=True))
        return 0

    if m == "g0w0":
        # Spectral (frequency-exact) G0W0 with dRPA screening — PySCF
        # gw.gw_exact.GWExact twin of Atomiverse RHF::g0w0(). HF density via
        # RKS(xc='hf') (gw_exact asserts a DFT mean-field object).
        # Local import keeps `np` bound even though later branches also
        # `import numpy as np` (which makes np function-local in Python).
        import numpy as np
        from pyscf.gw import gw_exact
        from pyscf import tdscf as tdscf_mod

        mf = dft.RKS(mol)
        mf.xc = "hf"
        mf.conv_tol = args.conv_tol
        e = float(mf.kernel())
        nocc = int(mol.nelectron // 2)
        nvir = int(mf.mo_coeff.shape[1] - nocc)
        td = tdscf_mod.dRPA(mf)
        td.nstates = max(1, nocc * nvir)
        td.verbose = 0
        td.kernel()
        ge = gw_exact.GWExact(mf, tdmf=td)
        ge.verbose = 0
        ge.kernel()
        qp = np.asarray(ge.mo_energy, dtype=float)
        homo = nocc - 1
        lumo = nocc
        gap_ha = float(qp[lumo] - qp[homo])
        out["e_tot_ha"] = e
        out["converged"] = bool(mf.converged)
        out["eps_hf_ha"] = np.asarray(mf.mo_energy, dtype=float).tolist()
        out["qp_energy_ha"] = qp.tolist()
        out["qp_homo_ha"] = float(qp[homo])
        out["qp_lumo_ha"] = float(qp[lumo])
        out["gap_ha"] = gap_ha
        out["gap_ev"] = gap_ha * _HA_TO_EV
        out["ip_ev"] = float(-qp[homo] * _HA_TO_EV)
        out["ea_ev"] = float(qp[lumo] * _HA_TO_EV)
        out["drpa_omega_ha"] = np.asarray(td.e, dtype=float).tolist()
        # AC companion (default gw.GW on RHF) for near-twin comparison only.
        try:
            from pyscf import gw as gw_mod

            mf_rhf = scf.RHF(mol)
            mf_rhf.conv_tol = args.conv_tol
            mf_rhf.kernel()
            gac = gw_mod.GW(mf_rhf)
            gac.verbose = 0
            gac.kernel()
            qac = np.asarray(gac.mo_energy, dtype=float)
            out["qp_ac_ha"] = qac.tolist()
            out["gap_ac_ev"] = float((qac[lumo] - qac[homo]) * _HA_TO_EV)
        except Exception as ac_err:  # pragma: no cover — optional twin
            out["ac_error"] = str(ac_err)
        print(json.dumps(out, indent=2, sort_keys=True))
        return 0 if out.get("converged", True) else 1

    if m == "rhf":
        mf = scf.RHF(mol)
    elif m in ("uhf", "ump2", "sf-tda"):
        mf = scf.UHF(mol)
    elif m in ("rohf", "romp2"):
        mf = scf.ROHF(mol)
    elif m in ("rks", "rks-tda", "rks-tddft"):
        mf = dft.RKS(mol)
        mf.xc = args.xc
        out["xc"] = args.xc
    elif m == "b2plyp":
        # Grimme B2PLYP double hybrid (JCP 124, 034108 (2006)):
        #   a_x = 0.53 HF, (1-a_x) B88, (1-a_c) LYP, a_c = 0.27 PT2/MP2.
        # Hosted libxc often lacks the named B2PLYP composite; compose the
        # hybrid-KS part explicitly and scale MP2 by a_c offline.
        mf = dft.RKS(mol)
        mf.xc = "0.53*HF + 0.47*GGA_X_B88 + 0.73*GGA_C_LYP"
        out["xc"] = mf.xc
        out["a_x"] = 0.53
        out["a_c"] = 0.27
    elif m in ("uks", "uks-tda"):
        mf = dft.UKS(mol)
        mf.xc = args.xc
        out["xc"] = args.xc
    elif m == "roks":
        mf = dft.ROKS(mol)
        mf.xc = args.xc
        out["xc"] = args.xc
    else:
        # tda / tdhf / mp2 / mp3 / ccsd ground state: RHF
        mf = scf.RHF(mol)


    if args.pcm:
        if m not in ("rhf", "uhf", "rohf", "rks", "uks", "roks"):
            print(
                json.dumps(
                    {
                        "error": f"--pcm not supported for method {m!r}",
                    },
                    indent=2,
                )
            )
            return 2
        mf = mf.PCM()
        mf.with_solvent.method = args.pcm_method
        mf.with_solvent.eps = args.pcm_eps
        mf.with_solvent.vdw_scale = args.pcm_scale
        mf.with_solvent.lebedev_order = args.pcm_lebedev_order
        out["pcm"] = {
            "method": args.pcm_method,
            "eps": args.pcm_eps,
            "vdw_scale": args.pcm_scale,
            "lebedev_order": args.pcm_lebedev_order,
        }

    mf.conv_tol = args.conv_tol
    e = mf.kernel()
    out["e_tot_ha"] = float(e)
    out["converged"] = bool(mf.converged)
    if args.pcm and hasattr(mf, "with_solvent"):
        out["e_solvent_ha"] = float(getattr(mf.with_solvent, "e", 0.0) or 0.0)

    if m in ("rohf", "roks", "uks", "uks-tda", "uhf", "ump2", "romp2", "sf-tda") and hasattr(
        mf, "spin_square"
    ):
        ss, mult = mf.spin_square()
        out["s2"] = float(ss)
        out["multiplicity"] = float(mult)

    if args.grad and m in ("rhf", "uhf", "rohf", "rks", "uks", "roks"):
        g = mf.nuc_grad_method().kernel()
        out["grad_ha_bohr"] = g.tolist()

    if m == "mp2":
        m2 = mp.MP2(mf).run()
        out["e_corr_ha"] = float(m2.e_corr)
        out["e_tot_mp2_ha"] = float(m2.e_tot)
        out["e_corr_ss_ha"] = float(m2.e_corr_ss)
        out["e_corr_os_ha"] = float(m2.e_corr_os)
        if args.grad:
            g = m2.nuc_grad_method().kernel()
            out["grad_ha_bohr"] = g.tolist()
    elif m == "mp2-fc":
        # Frozen-core RMP2: freeze lowest N occupied MOs (default 1 = O 1s on H2O).
        n_frozen = 1 if args.frozen is None else int(args.frozen)
        m2 = mp.MP2(mf, frozen=n_frozen).run()
        out["frozen"] = n_frozen
        out["e_hf_ha"] = float(e)
        out["e_corr_ha"] = float(m2.e_corr)
        out["e_tot_mp2_ha"] = float(m2.e_tot)
    elif m == "scs-mp2":
        # Spin-component-scaled MP2 (Grimme, JCP 118, 9095 (2003)):
        # E_corr^SCS = (1/3) E_SS + 1.2 E_OS on the RHF-MP2 spin channels.
        m2 = mp.MP2(mf).run()
        out["e_hf_ha"] = float(e)
        out["e_corr_mp2_ha"] = float(m2.e_corr)
        out["e_corr_ss_ha"] = float(m2.e_corr_ss)
        out["e_corr_os_ha"] = float(m2.e_corr_os)
        out["e_corr_ha"] = float(m2.emp2_scs)
        out["e_tot_mp2_ha"] = float(m2.e_tot_scs)
        out["c_os"] = 1.2
        out["c_ss"] = 1.0 / 3.0
    elif m == "sos-mp2":
        # Scaled-opposite-spin MP2 (Jung et al., JCP 121, 9793 (2004)):
        # E_corr^SOS = 1.3 E_OS + 0 * E_SS (opposite-spin only).
        # PySCF has no emp2_sos property; scale the RMP2 OS channel offline.
        m2 = mp.MP2(mf).run()
        c_os, c_ss = 1.3, 0.0
        e_corr_sos = c_os * float(m2.e_corr_os) + c_ss * float(m2.e_corr_ss)
        out["e_hf_ha"] = float(e)
        out["e_corr_mp2_ha"] = float(m2.e_corr)
        out["e_corr_ss_ha"] = float(m2.e_corr_ss)
        out["e_corr_os_ha"] = float(m2.e_corr_os)
        out["e_corr_ha"] = e_corr_sos
        out["e_tot_mp2_ha"] = float(e) + e_corr_sos
        out["c_os"] = c_os
        out["c_ss"] = c_ss
    elif m == "mp3":
        # Closed-shell RMP3. PySCF 2.x has no mp.MP3 — rebuild Szabo/Ostlund
        # (Psi4NumPy MP3.py) spatial diagrams on PySCF RHF MOs + ao2mo ERIs.
        # E_corr = E^(2) + E^(3); E^(2) must match mp.MP2.
        from pyscf import ao2mo
        import numpy as np

        m2 = mp.MP2(mf).run()
        C = np.asarray(mf.mo_coeff)
        eps = np.asarray(mf.mo_energy)
        nmo = C.shape[1]
        ndocc = mol.nelectron // 2
        eri_chem = ao2mo.restore(1, ao2mo.full(mol, C), nmo)
        # Physicist <pq|rs> = (pr|qs)
        MO = np.einsum("prqs->pqrs", eri_chem)
        o = slice(0, ndocc)
        v = slice(ndocc, nmo)
        eocc = eps[:ndocc]
        evirt = eps[ndocc:]
        den = (
            eocc.reshape(-1, 1, 1, 1)
            + eocc.reshape(-1, 1, 1)
            - evirt.reshape(-1, 1)
            - evirt
        )
        if np.any(den >= 0.0):
            print(
                json.dumps(
                    {
                        "ok": False,
                        "error": "mp3: non-negative orbital-energy denominator",
                    }
                )
            )
            return 2
        eps_inv = 1.0 / den
        e2 = 2.0 * np.einsum(
            "abrs,rsab,abrs", MO[o, o, v, v], MO[v, v, o, o], eps_inv
        )
        e2 -= np.einsum(
            "abrs,rsba,abrs", MO[o, o, v, v], MO[v, v, o, o], eps_inv
        )
        # 12 third-order diagrams (Psi4NumPy / Szabo Table 6.2).
        e3 = 2.0 * np.einsum(
            "abru,ruts,tsab,abru,abts",
            MO[o, o, v, v],
            MO[v, v, v, v],
            MO[v, v, o, o],
            eps_inv,
            eps_inv,
        )
        e3 += 2.0 * np.einsum(
            "adrs,cbad,rscb,adrs,cbrs",
            MO[o, o, v, v],
            MO[o, o, o, o],
            MO[v, v, o, o],
            eps_inv,
            eps_inv,
        )
        e3 += -4.0 * np.einsum(
            "acrt,rbsc,stab,acrt,abst",
            MO[o, o, v, v],
            MO[v, o, v, o],
            MO[v, v, o, o],
            eps_inv,
            eps_inv,
        )
        e3 += -4.0 * np.einsum(
            "bcrt,rasb,stac,bcrt,acst",
            MO[o, o, v, v],
            MO[v, o, v, o],
            MO[v, v, o, o],
            eps_inv,
            eps_inv,
        )
        e3 += 8.0 * np.einsum(
            "acrt,btsc,rsab,acrt,abrs",
            MO[o, o, v, v],
            MO[o, v, v, o],
            MO[v, v, o, o],
            eps_inv,
            eps_inv,
        )
        e3 += 2.0 * np.einsum(
            "cbrt,atsc,rsab,cbrt,abrs",
            MO[o, o, v, v],
            MO[o, v, v, o],
            MO[v, v, o, o],
            eps_inv,
            eps_inv,
        )
        e3 += -1.0 * np.einsum(
            "acrs,dbac,srdb,acrs,dbrs",
            MO[o, o, v, v],
            MO[o, o, o, o],
            MO[v, v, o, o],
            eps_inv,
            eps_inv,
        )
        e3 += -1.0 * np.einsum(
            "abrt,trus,usab,abtr,abus",
            MO[o, o, v, v],
            MO[v, v, v, v],
            MO[v, v, o, o],
            eps_inv,
            eps_inv,
        )
        e3 += 2.0 * np.einsum(
            "bcrt,arbs,tsac,cbrt,acst",
            MO[o, o, v, v],
            MO[o, v, o, v],
            MO[v, v, o, o],
            eps_inv,
            eps_inv,
        )
        e3 += 2.0 * np.einsum(
            "cbrt,rasb,stac,cbrt,acst",
            MO[o, o, v, v],
            MO[v, o, v, o],
            MO[v, v, o, o],
            eps_inv,
            eps_inv,
        )
        e3 += -4.0 * np.einsum(
            "abrs,scat,rtbc,abrs,cbrt",
            MO[o, o, v, v],
            MO[v, o, o, v],
            MO[v, v, o, o],
            eps_inv,
            eps_inv,
        )
        e3 += -4.0 * np.einsum(
            "bcrt,atsc,rsab,bctr,abrs",
            MO[o, o, v, v],
            MO[o, v, v, o],
            MO[v, v, o, o],
            eps_inv,
            eps_inv,
        )
        e_corr = float(e2 + e3)
        out["e_hf_ha"] = float(e)
        out["e_corr_mp2_ha"] = float(e2)
        out["e_corr_mp3_ha"] = float(e3)
        out["e_corr_ha"] = e_corr
        out["e_tot_mp3_ha"] = float(e) + e_corr
        out["e_corr_mp2_pyscf_ha"] = float(m2.e_corr)
        out["oracle"] = "szabo_psi4numpy_on_pyscf_rhf_ao2mo"
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
    elif m == "romp2":
        # ROHF reference → PySCF converts to UHF view + UMP2 (semi-canonical
        # spin-Fock diagonal denominators). No separate mp.ROMP2 class.
        m2 = mp.UMP2(mf).run()
        out["e_hf_ha"] = float(e)
        out["e_corr_ha"] = float(m2.e_corr)
        out["e_tot_mp2_ha"] = float(m2.e_tot)
        out["e_corr_ss_ha"] = float(m2.e_corr_ss)
        out["e_corr_os_ha"] = float(m2.e_corr_os)
        if hasattr(mf, "spin_square"):
            ss, mult = mf.spin_square()
            out["s2"] = float(ss)
            out["multiplicity"] = float(mult)
    elif m in ("tda", "tdhf", "rks-tda", "rks-tddft", "uks-tda"):
        # rks-tda: hybrid/pure TDDFT-TDA on the KS reference (exact-exchange
        # kernel for hybrids via PySCF's Casida A matrix).
        # rks-tddft: full RPA TDDFT (Casida A,B) — tdscf.TDDFT, not TDA.
        # uks-tda: open-shell spin-conserving TDA on a UKS reference
        # (PySCF tdscf.TDA on dft.UKS; no closed-shell singlet/triplet split).
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
    elif m == "sf-tda":
        # Collinear SF-CIS / SF-TDA (Shao, Head-Gordon, Krylov JCP 118, 4807
        # (2003)). PySCF mainline has no product SF-TDA; build the Ms-decreasing
        # singles block with PySCF UHF MOs + ao2mo (offline pin / fair twin).
        # A_ia,jb = δ_ij δ_ab (ε_a^min - ε_i^maj) - (ij|ab).
        from pyscf import ao2mo
        import numpy as np

        mo_a, mo_b = mf.mo_coeff
        ea, eb = mf.mo_energy
        occa, occb = mf.mo_occ
        oa = np.where(occa > 0)[0]
        va = np.where(occa == 0)[0]
        ob = np.where(occb > 0)[0]
        vb = np.where(occb == 0)[0]
        alpha_to_beta = len(oa) >= len(ob)
        if alpha_to_beta:
            o_idx, v_idx = oa, vb
            e_o, e_v = ea, eb
            Co, Cv = mo_a[:, oa], mo_b[:, vb]
            out["sf_channel"] = "alpha_to_beta"
        else:
            o_idx, v_idx = ob, va
            e_o, e_v = eb, ea
            Co, Cv = mo_b[:, ob], mo_a[:, va]
            out["sf_channel"] = "beta_to_alpha"
        no, nv = len(o_idx), len(v_idx)
        nov = no * nv
        nstates = min(args.nstates, nov)
        ijab = ao2mo.general(mol, (Co, Co, Cv, Cv), compact=False).reshape(
            no, no, nv, nv
        )
        A = np.zeros((nov, nov))
        for i in range(no):
            for a in range(nv):
                p = i * nv + a
                A[p, p] = e_v[v_idx[a]] - e_o[o_idx[i]]
                for j in range(no):
                    for b in range(nv):
                        q = j * nv + b
                        A[p, q] -= ijab[i, j, a, b]
        w = np.linalg.eigvalsh(A)
        out["n_ov"] = int(nov)
        out["excitation_ha"] = [float(x) for x in w[:nstates]]
        out["excitation_ev"] = [float(x) * _HA_TO_EV for x in w[:nstates]]
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
    elif m == "adc2":
        # Strict ADC(2) electronic excitation energies (singlet, EE).
        # Schirmer PRA 26, 2395 (1982); Trofimov & Schirmer JPB 28, 2299 (1995).
        from pyscf import adc as pyscf_adc

        myadc = pyscf_adc.ADC(mf)
        myadc.method = "adc(2)"
        myadc.method_type = "ee"
        myadc.conv_tol = min(args.conv_tol, 1e-10)
        e_adc, _v, _p, _x = myadc.kernel(nroots=args.nstates)
        out["e_corr_mp2_ha"] = float(myadc.e_corr)
        out["excitation_ha"] = [float(x) for x in e_adc]
        out["excitation_ev"] = [float(x) * _HA_TO_EV for x in e_adc]
        out["s1_ha"] = float(e_adc[0])
        out["s1_ev"] = float(e_adc[0]) * _HA_TO_EV
    elif m == "ccsd":
        # Closed-shell CCSD on the RHF reference (PySCF cc.CCSD).
        # Local import: other branches also re-bind `np`, so use a fresh name.
        import numpy as _np

        mycc = pyscf_cc.CCSD(mf)
        mycc.conv_tol = min(args.conv_tol, 1e-10)
        e_corr, t1, t2 = mycc.kernel()
        out["e_hf_ha"] = float(e)
        out["e_corr_ha"] = float(e_corr)
        out["e_tot_ccsd_ha"] = float(mycc.e_tot)
        out["e_tot_ha"] = float(mycc.e_tot)
        out["ccsd_converged"] = bool(mycc.converged)
        t1a = _np.asarray(t1).ravel()
        t2a = _np.asarray(t2).ravel()
        out["t1_norm"] = float(_np.sqrt(t1a.dot(t1a)))
        out["t2_norm"] = float(_np.sqrt(t2a.dot(t2a)))
    elif m == "b2plyp":
        # Double-hybrid total energy: E_KS + a_c * E_MP2 on the KS orbitals.
        m2 = mp.MP2(mf).run()
        a_c = 0.27
        e_corr_pt2 = a_c * float(m2.e_corr)
        out["e_ks_ha"] = float(e)
        out["e_corr_mp2_ha"] = float(m2.e_corr)
        out["e_corr_ss_ha"] = float(m2.e_corr_ss)
        out["e_corr_os_ha"] = float(m2.e_corr_os)
        out["e_corr_pt2_ha"] = e_corr_pt2
        out["e_corr_ha"] = e_corr_pt2
        out["e_tot_ha"] = float(e) + e_corr_pt2
        out["a_c"] = a_c
        out["a_x"] = 0.53

    print(json.dumps(out, indent=2, sort_keys=True))
    return 0 if out.get("converged", True) else 1


if __name__ == "__main__":
    sys.exit(main())
