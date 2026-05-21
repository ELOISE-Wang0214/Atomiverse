#!/usr/bin/env python3
"""Generate PySCF harmonic_analysis reference frequencies.

Runs a genuine analytic-Hessian vibrational analysis for H2O, NH3, and
CH4 at RHF / 6-31G and emits two artefacts that the C++ parity test
(``test_frequencies_vs_pyscf.cpp``) consumes:

  * ``pyscf_frequencies.json``   — provenance + per-molecule data
  * ``pyscf_frequencies_data.h`` — the same data as a C++ header (the
    test build has no JSON library, so the numbers are embedded as a
    generated header; the JSON is the human-readable provenance copy).

The geometries here are the EXACT geometries the C++ test uses, so the
numerical Hessian and the analytic Hessian are evaluated at the same
nuclear configuration and the only gap is numerical-vs-analytic.

Usage:
    python tests/cpu/vib/reference/generate_pyscf_frequencies.py
"""

import json
import os
import textwrap

import numpy
import pyscf
from pyscf import gto, scf
from pyscf.hessian import thermo

# ── Molecule definitions ────────────────────────────────────────────────
# Each geometry is given in Angstrom and MUST match the C++ test verbatim.
MOLECULES = {
    "h2o": {
        "atoms": [
            ("O", (0.0000000000,  0.0000000000,  0.0000000000)),
            ("H", (0.0000000000,  0.7572000000,  0.5865000000)),
            ("H", (0.0000000000, -0.7572000000,  0.5865000000)),
        ],
    },
    "nh3": {
        "atoms": [
            ("N", (0.0000000000,  0.0000000000,  0.0000000000)),
            ("H", (0.0000000000,  0.9377000000, -0.3816000000)),
            ("H", (0.8121000000, -0.4688000000, -0.3816000000)),
            ("H", (-0.8121000000, -0.4688000000, -0.3816000000)),
        ],
    },
    "ch4": {
        "atoms": [
            ("C", (0.0000000000,  0.0000000000,  0.0000000000)),
            ("H", (0.6276000000,  0.6276000000,  0.6276000000)),
            ("H", (-0.6276000000, -0.6276000000,  0.6276000000)),
            ("H", (-0.6276000000,  0.6276000000, -0.6276000000)),
            ("H", (0.6276000000, -0.6276000000, -0.6276000000)),
        ],
    },
}

BASIS = "6-31g"
METHOD = "RHF"


ANG2BOHR = 1.8897261245650618


def build_mol(symbols, coords_ang):
    mol = gto.Mole()
    mol.atom = [(symbols[i], tuple(coords_ang[i]))
                for i in range(len(symbols))]
    mol.unit = "Angstrom"
    mol.basis = BASIS
    mol.verbose = 0
    mol.build()
    return mol


def rhf_energy_grad(symbols, coords_ang):
    """RHF energy [Hartree] and nuclear gradient [Hartree/Angstrom]."""
    mol = build_mol(symbols, coords_ang)
    mf = scf.RHF(mol)
    mf.conv_tol = 1e-12
    mf.conv_tol_grad = 1e-9
    e = mf.kernel()
    g_bohr = mf.nuc_grad_method().kernel()  # Hartree/Bohr
    g_ang = g_bohr * ANG2BOHR               # Hartree/Angstrom
    return e, g_ang, mf, mol


def optimize_geometry(symbols, coords_ang):
    """Plain BFGS geometry optimization to a tight gradient norm.

    A tight stationary point is required: the C++ numerical Hessian and
    the PySCF analytic Hessian must be evaluated at the SAME stationary
    geometry, otherwise the translation/rotation projection mixes
    gradient contamination into the low modes.
    """
    x = numpy.asarray(coords_ang, dtype=float).reshape(-1)
    n = x.size
    H = numpy.eye(n)  # inverse-Hessian estimate.
    e, g, _, _ = rhf_energy_grad(symbols, x.reshape(-1, 3))
    g = g.reshape(-1)
    for it in range(200):
        gmax = numpy.abs(g).max()
        if gmax < 1e-6:  # Hartree/Angstrom — tight stationary point.
            break
        p = -H.dot(g)
        # Backtracking line search on the energy.
        alpha = 1.0
        for _ in range(30):
            x_new = x + alpha * p
            e_new, g_new, _, _ = rhf_energy_grad(
                symbols, x_new.reshape(-1, 3))
            if e_new <= e + 1e-4 * alpha * g.dot(p):
                break
            alpha *= 0.5
        g_new = g_new.reshape(-1)
        s = x_new - x
        y = g_new - g
        sy = s.dot(y)
        if sy > 1e-12:  # BFGS inverse-Hessian update.
            rho = 1.0 / sy
            I = numpy.eye(n)
            H = ((I - rho * numpy.outer(s, y)).dot(H)
                 .dot(I - rho * numpy.outer(y, s))
                 + rho * numpy.outer(s, s))
        x, e, g = x_new, e_new, g_new
    return x.reshape(-1, 3), e, numpy.abs(g).max()


def analyze(name, spec):
    """Optimize the geometry, then run the analytic harmonic analysis."""
    symbols = [a[0] for a in spec["atoms"]]
    coords0 = [list(a[1]) for a in spec["atoms"]]

    opt_coords, e_opt, gmax = optimize_geometry(symbols, coords0)
    print(f"    {name}: optimized, |grad|_max = {gmax:.2e} Ha/Ang")

    mol = build_mol(symbols, opt_coords)
    mf = scf.RHF(mol)
    mf.conv_tol = 1e-12
    mf.conv_tol_grad = 1e-9
    mf.kernel()

    hobj = mf.Hessian()
    hess = hobj.kernel()  # analytic Hessian, [natm,natm,3,3], Hartree/Bohr^2
    res = thermo.harmonic_analysis(mol, hess, imaginary_freq=False)
    freq = numpy.asarray(res["freq_wavenumber"]).real
    freq = numpy.sort(freq)
    return {
        "method": METHOD,
        "basis": BASIS,
        "natoms": mol.natm,
        "geometry_angstrom": [
            [symbols[i], [float(c) for c in opt_coords[i]]]
            for i in range(len(symbols))
        ],
        "n_modes": int(freq.size),
        "frequencies_cm1": [float(x) for x in freq],
        "scf_energy_hartree": float(mf.e_tot),
        "grad_max_hartree_per_ang": float(gmax),
    }


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    data = {
        "provenance": {
            "pyscf_version": pyscf.__version__,
            "numpy_version": numpy.__version__,
            "method": METHOD,
            "basis": BASIS,
            "hessian": "analytic (mf.Hessian().kernel())",
            "harmonic_analysis": "pyscf.hessian.thermo.harmonic_analysis",
            "note": "Reference frequencies are from a genuine PySCF run.",
        },
        "molecules": {},
    }
    for name, spec in MOLECULES.items():
        print(f"  analyzing {name} ...")
        data["molecules"][name] = analyze(name, spec)
        freqs = data["molecules"][name]["frequencies_cm1"]
        print(f"    {name}: {len(freqs)} modes -> "
              f"{[round(f, 1) for f in freqs]}")

    json_path = os.path.join(here, "pyscf_frequencies.json")
    with open(json_path, "w") as fh:
        json.dump(data, fh, indent=2)
    print(f"wrote {json_path}")

    # Emit the C++ header consumed by the parity test.
    header_path = os.path.join(here, "pyscf_frequencies_data.h")
    write_header(header_path, data)
    print(f"wrote {header_path}")


def write_header(path, data):
    lines = []
    lines.append("/// @file pyscf_frequencies_data.h")
    lines.append("/// @brief GENERATED — PySCF harmonic_analysis reference "
                 "frequencies.")
    lines.append("///")
    lines.append("/// Generated by "
                 "tests/cpu/vib/reference/generate_pyscf_frequencies.py.")
    lines.append(f"/// PySCF {data['provenance']['pyscf_version']}, "
                 f"{data['provenance']['method']} / "
                 f"{data['provenance']['basis']}, analytic Hessian.")
    lines.append("/// DO NOT EDIT BY HAND — regenerate via the script.")
    lines.append("#pragma once")
    lines.append("")
    lines.append("#include <vector>")
    lines.append("")
    lines.append("namespace pyscf_ref {")
    lines.append("")
    lines.append("/// @brief One atom in a reference geometry.")
    lines.append("struct RefAtom {")
    lines.append("    int Z;                       ///< Atomic number.")
    lines.append("    double x, y, z;              ///< Coordinates [Ang].")
    lines.append("};")
    lines.append("")
    lines.append("/// @brief One molecule: optimized geometry + reference "
                 "frequencies.")
    lines.append("struct RefMolecule {")
    lines.append("    const char* name;            ///< Molecule label.")
    lines.append("    int natoms;                  ///< Atom count.")
    lines.append("    std::vector<RefAtom> geom;   ///< Optimized geometry.")
    lines.append("    std::vector<double> freq;    ///< Wavenumbers [cm^-1].")
    lines.append("};")
    lines.append("")
    lines.append(f"/// @brief PySCF version that produced these references.")
    lines.append(f'inline const char* kPyscfVersion = '
                 f'"{data["provenance"]["pyscf_version"]}";')
    lines.append("")
    lines.append("/// @brief Reference frequencies for all parity molecules.")
    lines.append("inline const std::vector<RefMolecule>& reference() {")
    lines.append("    static const std::vector<RefMolecule> kRef = {")
    sym2z = {"H": 1, "C": 6, "N": 7, "O": 8, "F": 9}
    for name, m in data["molecules"].items():
        freqs = ", ".join(f"{f:.10f}" for f in m["frequencies_cm1"])
        atoms = ", ".join(
            f"{{{sym2z[sym]}, {c[0]:.12f}, {c[1]:.12f}, {c[2]:.12f}}}"
            for sym, c in m["geometry_angstrom"])
        lines.append(f'        {{"{name}", {m["natoms"]},')
        lines.append(f'         {{{atoms}}},')
        lines.append(f'         {{{freqs}}}}},')
    lines.append("    };")
    lines.append("    return kRef;")
    lines.append("}")
    lines.append("")
    lines.append("} // namespace pyscf_ref")
    lines.append("")
    with open(path, "w") as fh:
        fh.write("\n".join(lines))


if __name__ == "__main__":
    main()
