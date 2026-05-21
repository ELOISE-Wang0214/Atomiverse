#!/usr/bin/env python3
"""Generate PySCF ``thermo()`` reference data for the vibrational-analysis
thermochemistry link (spec phase4-vibrations-03).

For each of H2O / NH3 / CH4 this runs a genuine PySCF analytic-Hessian
vibrational analysis followed by ``pyscf.hessian.thermo.thermo()`` at
T = 298.15 K, P = 101325 Pa, and records:

  * the zero-point energy (ZPE)
  * the thermal Gibbs correction ``G_tot - E0``
  * the per-component entropies ``S_trans`` / ``S_rot`` / ``S_vib``
  * the per-component / total ``thermo()`` dict entries needed by the
    C++ parity test (``test_thermochem_vs_pyscf.cpp``)

The geometries are the EXACT optimized geometries that the link-2
frequency parity test uses (``tests/cpu/vib/reference/
pyscf_frequencies.json``), so the C++ ``cpu::vib::Thermochem`` and PySCF
``thermo()`` operate on the same nuclear configuration and the same
harmonic frequencies — the only gap is the thermochemistry port itself.

The symmetry numbers are pinned per molecule (H2O = 2, NH3 = 3,
CH4 = 12) and passed verbatim to ``thermo()``; sigma auto-detection is
out of scope for this link.

Output:
  * ``pyscf_thermo_references.json`` — provenance + per-molecule data.
  * ``tests/cpu/vib/reference/pyscf_thermo_data.h`` — the same data as a
    generated C++ header (the test build has no JSON library, so the
    numbers are embedded; the JSON is the human-readable provenance copy).

Usage:
    python benchmarks/qc/reference/thermo/generate_thermo_references.py
"""

from __future__ import annotations

import json
import os

import numpy
import pyscf
from pyscf import gto, scf
from pyscf.hessian import thermo

# Optimized geometries (Angstrom) — copied verbatim from
# tests/cpu/vib/reference/pyscf_frequencies.json so the thermochemistry
# port is validated at the SAME stationary geometry as the link-2
# frequency parity test.
MOLECULES = {
    "h2o": {
        "symmetry_number": 2,
        "atoms": [
            ("O", (-0.000000000000, -0.000000000000, 0.034902683452)),
            ("H", (-0.000000000000, 0.785166496369, 0.569048658274)),
            ("H", (0.000000000000, -0.785166496369, 0.569048658274)),
        ],
    },
    "nh3": {
        "symmetry_number": 3,
        "atoms": [
            ("N", (0.000000000000, 0.000023233161, -0.138019413462)),
            ("H", (-0.000000000000, 0.971477016737, -0.335582217196)),
            ("H", (0.841301488611, -0.485700124949, -0.335599184671)),
            ("H", (-0.841301488611, -0.485700124949, -0.335599184671)),
        ],
    },
    "ch4": {
        "symmetry_number": 12,
        "atoms": [
            ("C", (-0.000000000000, 0.000000000000, 0.000000000000)),
            ("H", (0.624754213345, 0.624754213345, 0.624754213345)),
            ("H", (-0.624754213345, -0.624754213345, 0.624754213345)),
            ("H", (-0.624754213345, 0.624754213345, -0.624754213345)),
            ("H", (0.624754213345, -0.624754213345, -0.624754213345)),
        ],
    },
}

BASIS = "6-31g"
METHOD = "RHF"
TEMPERATURE = 298.15   # K
PRESSURE = 101325.0    # Pa (1 atm)


def build_mol(symbols, coords_ang):
    """Build a PySCF Mole at the given Angstrom geometry."""
    mol = gto.Mole()
    mol.atom = [(symbols[i], tuple(coords_ang[i]))
                for i in range(len(symbols))]
    mol.unit = "Angstrom"
    mol.basis = BASIS
    mol.verbose = 0
    mol.build()
    return mol


def analyze(name, spec):
    """Run analytic Hessian + ``thermo()`` for one molecule."""
    symbols = [a[0] for a in spec["atoms"]]
    coords = [list(a[1]) for a in spec["atoms"]]
    sigma = spec["symmetry_number"]

    mol = build_mol(symbols, coords)
    mf = scf.RHF(mol)
    mf.conv_tol = 1e-12
    mf.conv_tol_grad = 1e-9
    e0 = mf.kernel()

    hobj = mf.Hessian()
    hess = hobj.kernel()  # analytic Hessian [natm,natm,3,3], Hartree/Bohr^2
    harm = thermo.harmonic_analysis(mol, hess, imaginary_freq=False)
    freq_au = numpy.asarray(harm["freq_au"]).real

    info = thermo.thermo(mf, harm["freq_au"], TEMPERATURE, PRESSURE)

    # thermo() returns dict entries as (value, unit) tuples.
    zpe = float(info["ZPE"][0])
    e_0k = float(info["E_0K"][0])
    e_tot = float(info["E_tot"][0])
    h_tot = float(info["H_tot"][0])
    g_tot = float(info["G_tot"][0])
    s_tot = float(info["S_tot"][0])
    s_elec = float(info["S_elec"][0])
    s_trans = float(info["S_trans"][0])
    s_rot = float(info["S_rot"][0])
    s_vib = float(info["S_vib"][0])

    return {
        "method": METHOD,
        "basis": BASIS,
        "natoms": mol.natm,
        "symmetry_number": sigma,
        "multiplicity": int(mol.spin) + 1,
        "temperature_K": TEMPERATURE,
        "pressure_Pa": PRESSURE,
        "geometry_angstrom": [
            [symbols[i], [float(c) for c in coords[i]]]
            for i in range(len(symbols))
        ],
        "scf_energy_hartree": float(e0),
        "n_modes": int(freq_au.size),
        "freq_au": [float(x) for x in freq_au],
        # Reference thermochemistry — all in Hartree / Hartree-per-K.
        "zpe_hartree": zpe,
        "e_0k_hartree": e_0k,
        "e_tot_hartree": e_tot,
        "h_tot_hartree": h_tot,
        "g_tot_hartree": g_tot,
        "s_tot_hartree_per_K": s_tot,
        "s_elec_hartree_per_K": s_elec,
        "s_trans_hartree_per_K": s_trans,
        "s_rot_hartree_per_K": s_rot,
        "s_vib_hartree_per_K": s_vib,
        # The acceptance compares this electronic-energy-independent
        # quantity: the thermal correction to Gibbs.
        "g_thermal_hartree": g_tot - float(e0),
        "zpe_check_hartree": e_0k - float(e0),  # must equal ZPE.
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
            "thermo": "pyscf.hessian.thermo.thermo",
            "temperature_K": TEMPERATURE,
            "pressure_Pa": PRESSURE,
            "note": ("Reference thermochemistry from a genuine PySCF "
                     "thermo() run. Geometries match "
                     "tests/cpu/vib/reference/pyscf_frequencies.json."),
        },
        "molecules": {},
    }
    for name, spec in MOLECULES.items():
        print(f"  analyzing {name} (sigma={spec['symmetry_number']}) ...")
        block = analyze(name, spec)
        data["molecules"][name] = block
        print(f"    {name}: ZPE={block['zpe_hartree']:.8f} Ha  "
              f"G-E0={block['g_thermal_hartree']:.8f} Ha  "
              f"S_tot={block['s_tot_hartree_per_K']:.10f} Ha/K")

    out = os.path.join(here, "pyscf_thermo_references.json")
    with open(out, "w") as fh:
        json.dump(data, fh, indent=2)
    print(f"wrote {out}")

    # Emit the C++ header consumed by the parity test.
    header = os.path.normpath(os.path.join(
        here, "..", "..", "..", "..",
        "tests", "cpu", "vib", "reference", "pyscf_thermo_data.h"))
    write_header(header, data)
    print(f"wrote {header}")


def write_header(path, data):
    """Emit the generated C++ reference header for the parity test."""
    sym2z = {"H": 1, "C": 6, "N": 7, "O": 8, "F": 9}
    lines = []
    lines.append("/// @file pyscf_thermo_data.h")
    lines.append("/// @brief GENERATED — PySCF thermo() reference "
                 "thermochemistry.")
    lines.append("///")
    lines.append("/// Generated by benchmarks/qc/reference/thermo/"
                 "generate_thermo_references.py.")
    prov = data["provenance"]
    lines.append(f"/// PySCF {prov['pyscf_version']}, {prov['method']} / "
                 f"{prov['basis']}, analytic Hessian + thermo().")
    lines.append(f"/// T = {prov['temperature_K']} K, "
                 f"P = {prov['pressure_Pa']} Pa.")
    lines.append("/// DO NOT EDIT BY HAND — regenerate via the script.")
    lines.append("#pragma once")
    lines.append("")
    lines.append("#include <vector>")
    lines.append("")
    lines.append("namespace pyscf_thermo_ref {")
    lines.append("")
    lines.append("/// @brief One atom of a reference geometry.")
    lines.append("struct RefAtom {")
    lines.append("    int Z;            ///< Atomic number.")
    lines.append("    double x, y, z;   ///< Coordinates [Angstrom].")
    lines.append("};")
    lines.append("")
    lines.append("/// @brief One molecule: geometry + PySCF thermo() "
                 "reference values.")
    lines.append("struct RefThermo {")
    lines.append("    const char* name;            ///< Molecule label.")
    lines.append("    int natoms;                  ///< Atom count.")
    lines.append("    int multiplicity;            ///< Spin multiplicity.")
    lines.append("    int symmetry_number;         ///< Rotational sigma.")
    lines.append("    double temperature_K;        ///< Temperature [K].")
    lines.append("    double pressure_Pa;          ///< Pressure [Pa].")
    lines.append("    std::vector<RefAtom> geom;   ///< Geometry [Ang].")
    lines.append("    double zpe_hartree;          ///< ZPE [Hartree].")
    lines.append("    double g_thermal_hartree;    ///< G_tot - E0 [Hartree].")
    lines.append("    double s_trans;              ///< S_trans [Eh/K].")
    lines.append("    double s_rot;                ///< S_rot [Eh/K].")
    lines.append("    double s_vib;                ///< S_vib [Eh/K].")
    lines.append("    double s_tot;                ///< S_tot [Eh/K].")
    lines.append("};")
    lines.append("")
    lines.append("/// @brief PySCF version that produced these references.")
    lines.append(f'inline const char* kPyscfVersion = '
                 f'"{prov["pyscf_version"]}";')
    lines.append("")
    lines.append("/// @brief PySCF thermo() references for all parity "
                 "molecules.")
    lines.append("inline const std::vector<RefThermo>& reference() {")
    lines.append("    static const std::vector<RefThermo> kRef = {")
    for name, m in data["molecules"].items():
        atoms = ", ".join(
            f"{{{sym2z[sym]}, {c[0]:.12f}, {c[1]:.12f}, {c[2]:.12f}}}"
            for sym, c in m["geometry_angstrom"])
        lines.append(f'        {{"{name}", {m["natoms"]}, '
                     f'{m["multiplicity"]}, {m["symmetry_number"]},')
        lines.append(f'         {m["temperature_K"]}, {m["pressure_Pa"]},')
        lines.append(f'         {{{atoms}}},')
        lines.append(f'         {m["zpe_hartree"]:.15e}, '
                     f'{m["g_thermal_hartree"]:.15e},')
        lines.append(f'         {m["s_trans_hartree_per_K"]:.15e}, '
                     f'{m["s_rot_hartree_per_K"]:.15e},')
        lines.append(f'         {m["s_vib_hartree_per_K"]:.15e}, '
                     f'{m["s_tot_hartree_per_K"]:.15e}}},')
    lines.append("    };")
    lines.append("    return kRef;")
    lines.append("}")
    lines.append("")
    lines.append("} // namespace pyscf_thermo_ref")
    lines.append("")
    with open(path, "w") as fh:
        fh.write("\n".join(lines))


if __name__ == "__main__":
    main()
