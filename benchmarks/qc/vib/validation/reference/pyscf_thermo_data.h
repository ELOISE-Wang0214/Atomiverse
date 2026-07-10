/// @file pyscf_thermo_data.h
/// @brief GENERATED — PySCF thermo() reference thermochemistry.
///
/// Pinned from a one-off PySCF run. Generator removed; do not regenerate.
/// PySCF 2.13.0, RHF / 6-31g, analytic Hessian + thermo().
/// T = 298.15 K, P = 101325.0 Pa.
/// DO NOT EDIT BY HAND — regenerate via the script.
#pragma once

#include <vector>

namespace pyscf_thermo_ref {

/// @brief One atom of a reference geometry.
struct RefAtom {
    int Z;          ///< Atomic number.
    double x, y, z; ///< Coordinates [Angstrom].
};

/// @brief One molecule: geometry + PySCF thermo() reference values.
struct RefThermo {
    const char* name;          ///< Molecule label.
    int natoms;                ///< Atom count.
    int multiplicity;          ///< Spin multiplicity.
    int symmetry_number;       ///< Rotational sigma.
    double temperature_K;      ///< Temperature [K].
    double pressure_Pa;        ///< Pressure [Pa].
    std::vector<RefAtom> geom; ///< Geometry [Ang].
    double zpe_hartree;        ///< ZPE [Hartree].
    double g_thermal_hartree;  ///< G_tot - E0 [Hartree].
    double s_trans;            ///< S_trans [Eh/K].
    double s_rot;              ///< S_rot [Eh/K].
    double s_vib;              ///< S_vib [Eh/K].
    double s_tot;              ///< S_tot [Eh/K].
};

/// @brief PySCF version that produced these references.
inline const char* kPyscfVersion = "2.13.0";

/// @brief PySCF thermo() references for all parity molecules.
inline const std::vector<RefThermo>& reference() {
    static const std::vector<RefThermo> kRef = {
        {"h2o",
         3,
         1,
         2,
         298.15,
         101325.0,
         {{8, -0.000000000000, -0.000000000000, 0.034902683452},
          {1, -0.000000000000, 0.785166496369, 0.569048658274},
          {1, 0.000000000000, -0.785166496369, 0.569048658274}},
         2.248565969680725e-02,
         4.911971424306216e-03,
         5.515294252088886e-05,
         1.645601332456206e-05,
         6.807322317448142e-09,
         7.161576316776837e-05},
        {"nh3",
         4,
         1,
         3,
         298.15,
         101325.0,
         {{7, 0.000000000000, 0.000023233161, -0.138019413462},
          {1, -0.000000000000, 0.971477016737, -0.335582217196},
          {1, 0.841301488611, -0.485700124949, -0.335599184671},
          {1, -0.841301488611, -0.485700124949, -0.335599184671}},
         3.638906106977707e-02,
         1.838109336023308e-02,
         5.488612512686279e-05,
         1.799668961583614e-05,
         7.335733605506259e-07,
         7.361638810324955e-05},
        {"ch4",
         5,
         1,
         12,
         298.15,
         101325.0,
         {{6, -0.000000000000, 0.000000000000, 0.000000000000},
          {1, 0.624754213345, 0.624754213345, 0.624754213345},
          {1, -0.624754213345, -0.624754213345, 0.624754213345},
          {1, -0.624754213345, 0.624754213345, -0.624754213345},
          {1, 0.624754213345, -0.624754213345, -0.624754213345}},
         4.792693888592365e-02,
         3.063379764980567e-02,
         5.460223990284688e-05,
         1.605852090887093e-05,
         6.780410914117017e-08,
         7.072856492085897e-05},
    };
    return kRef;
}

} // namespace pyscf_thermo_ref
