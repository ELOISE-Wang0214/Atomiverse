/// @file cpu/scf/grid.h
/// @brief Molecular integration grid for DFT: Lebedev angular + radial + Becke partitioning.
///
/// Atom-centered quadrature grids for numerical XC integration.
///
/// References:
///   - Lebedev-Laikov: Doklady Math. 59, 477 (1999)
///   - Mura-Knowles radial: JCP 104, 9848 (1996)
///   - Becke partitioning: JCP 88, 2547 (1988)
#pragma once

#include <vector>
#include "atomiverse/cpu/scf/native_integrals.h"

namespace ATV {
namespace cpu {
namespace scf {

/// @brief Single quadrature grid point with position and weight.
struct GridPoint {
    double x, y, z;   ///< Position in Bohr.
    double weight;     ///< Quadrature weight (includes radial, angular, and Becke partition).
};

/// @brief Molecular integration grid.
struct MolecularGrid {
    std::vector<GridPoint> points;  ///< All grid points.
    int n_radial;                   ///< Radial points per atom.
    int n_angular;                  ///< Lebedev angular points per atom.
};

/// @brief Build atom-centered molecular grid with Becke partitioning.
///
/// Constructs atom-centered grids using Euler-Maclaurin radial quadrature
/// and Lebedev angular quadrature, combined with Becke fuzzy-cell partitioning
/// for molecular integration.
///
/// @param atoms Array of atom descriptors (Z, position in Bohr).
/// @param natm  Number of atoms.
/// @param n_radial Number of radial quadrature points per atom.
/// @param n_angular Number of Lebedev angular points per atom (6, 14, 26, 50, 110, 194, 302, 590).
/// @return MolecularGrid with all weighted grid points.
MolecularGrid build_molecular_grid(const native::Atom* atoms, int natm,
                                   int n_radial = 50, int n_angular = 110);

} // namespace scf
} // namespace cpu
} // namespace ATV
