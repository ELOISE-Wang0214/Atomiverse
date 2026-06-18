/// @file cpu/gpw/internal/one_electron.h
/// @brief Analytic overlap (S) and kinetic-energy (T) integrals over primitive
///        Cartesian Gaussians, for the GPW KS matrix.
///
/// These are the "hard" (analytic) one-electron parts of the GPW Hamiltonian
/// H_KS = T + V_PP + V_soft[ρ]: the kinetic energy and the overlap metric of the
/// generalized eigenproblem H C = S C ε. Implemented with the McMurchie-Davidson
/// Hermite-expansion scheme (arbitrary angular momentum), operating on the shared
/// `GridPrimitive` type so the GPW basis (`GridBasisFunction`) is used directly.
#pragma once

#include "atomiverse/cpu/gpw/internal/collocate.h" // GridPrimitive, GridBasisFunction

namespace ATV {
namespace cpu {
namespace gpw {

/// @brief Overlap integral ⟨a|b⟩ of two primitive Cartesian Gaussians.
double prim_overlap(const GridPrimitive& a, const GridPrimitive& b);

/// @brief Kinetic-energy integral ⟨a|-½∇²|b⟩ of two primitive Gaussians.
double prim_kinetic(const GridPrimitive& a, const GridPrimitive& b);

/// @brief Overlap of two contracted basis functions (sum over primitive pairs).
double overlap(const GridBasisFunction& a, const GridBasisFunction& b);

/// @brief Kinetic energy of two contracted basis functions.
double kinetic(const GridBasisFunction& a, const GridBasisFunction& b);

} // namespace gpw
} // namespace cpu
} // namespace ATV
