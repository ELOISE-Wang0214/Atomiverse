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
/// Constructs atom-centered grids using Mura-Knowles radial quadrature
/// (JCP 104, 9848 (1996)) and tabulated Lebedev-Laikov angular quadrature
/// (Doklady Math. 59, 477 (1999)), combined with Becke fuzzy-cell
/// partitioning (JCP 88, 2547 (1988)) for molecular integration.
///
/// @param atoms Array of atom descriptors (Z, position in Bohr).
/// @param natm  Number of atoms.
/// @param n_radial Number of radial quadrature points per atom.
/// @param n_angular Number of Lebedev angular points per atom; one of
///        {6, 50, 110, 194, 302, 590}. Other values round down to the
///        next supported order.
/// @return MolecularGrid with all weighted grid points.
MolecularGrid build_molecular_grid(const native::Atom* atoms, int natm,
                                   int n_radial = 50, int n_angular = 110);

/// @brief Bragg-Slater atomic radius for element Z, in Bohr.
///
/// Throws std::invalid_argument for Z outside the supported range [1, 86].
/// Values are Slater 1964 (JCP 41, 3199) with noble-gas radii from
/// Clementi-Raimondi-Reinhardt 1967 (JCP 47, 1300) and At from Cordero 2008.
///
/// @param Z Atomic number (1 <= Z <= 86).
/// @return Bragg-Slater radius in Bohr.
double bragg_slater_radius(int Z);

/// @brief A single radial quadrature point.
struct RadialPoint {
    double r;  ///< Radius in Bohr.
    double w;  ///< Quadrature weight, includes r^2 Jacobian (so integrals
               ///< of f(r) over 3-space are sum_i w_i * f(r_i) * 4*pi).
};

/// @brief Mura-Knowles radial quadrature on (0, infinity).
///
/// Reference: M. E. Mura and P. J. Knowles, "Improved radial grids for
/// quadrature in molecular density-functional calculations", J. Chem.
/// Phys. 104, 9848 (1996), DOI 10.1063/1.471749.
///
/// Uses the log3 mapping r = -R * ln(1 - x^3) with x = (i + 0.5)/n
/// (uniform open Newton-Cotes) and weights w_i = R * 3 x_i^2 / (n (1 - x_i^3)) * r_i^2.
/// Matches PySCF's default radial grid scheme.
///
/// @param n Number of radial points.
/// @param R_bragg Element-specific Bragg-Slater radius in Bohr.
/// @return Vector of n radial points and weights.
std::vector<RadialPoint> radial_grid_mura_knowles(int n, double R_bragg);

} // namespace scf
} // namespace cpu
} // namespace ATV
