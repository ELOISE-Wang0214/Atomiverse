/// @file pyscf_frequencies_data.h
/// @brief GENERATED — PySCF harmonic_analysis reference frequencies.
///
/// Pinned from a one-off PySCF run. Generator removed; do not regenerate.
/// PySCF 2.13.0, RHF / 6-31g, analytic Hessian.
/// DO NOT EDIT BY HAND — regenerate via the script.
#pragma once

#include <vector>

namespace pyscf_ref {

/// @brief One atom in a reference geometry.
struct RefAtom {
    int Z;          ///< Atomic number.
    double x, y, z; ///< Coordinates [Ang].
};

/// @brief One molecule: optimized geometry + reference frequencies.
struct RefMolecule {
    const char* name;          ///< Molecule label.
    int natoms;                ///< Atom count.
    std::vector<RefAtom> geom; ///< Optimized geometry.
    std::vector<double> freq;  ///< Wavenumbers [cm^-1].
};

/// @brief PySCF version that produced these references.
inline const char* kPyscfVersion = "2.13.0";

/// @brief Reference frequencies for all parity molecules.
inline const std::vector<RefMolecule>& reference() {
    static const std::vector<RefMolecule> kRef = {
        {"h2o",
         3,
         {{8, -0.000000000000, -0.000000000000, 0.034902683452},
          {1, -0.000000000000, 0.785166496369, 0.569048658274},
          {1, 0.000000000000, -0.785166496369, 0.569048658274}},
         {1736.8524497254, 3988.1512229854, 4145.0600563142}},
        {"nh3",
         4,
         {{7, 0.000000000000, 0.000023233161, -0.138019413462},
          {1, -0.000000000000, 0.971477016737, -0.335582217196},
          {1, 0.841301488611, -0.485700124949, -0.335599184671},
          {1, -0.841301488611, -0.485700124949, -0.335599184671}},
         {597.3615374423, 1814.4315527013, 1814.4318020332, 3779.4845062762,
          3983.6202587061, 3983.6218741190}},
        {"ch4",
         5,
         {{6, -0.000000000000, 0.000000000000, 0.000000000000},
          {1, 0.624754213345, 0.624754213345, 0.624754213345},
          {1, -0.624754213345, -0.624754213345, 0.624754213345},
          {1, -0.624754213345, 0.624754213345, -0.624754213345},
          {1, 0.624754213345, -0.624754213345, -0.624754213345}},
         {1516.7504826961, 1516.7504830831, 1516.7504838415, 1708.6742156259,
          1708.6742158787, 3181.6142770523, 3296.0934350147, 3296.0934356081,
          3296.0934363777}},
    };
    return kRef;
}

} // namespace pyscf_ref
