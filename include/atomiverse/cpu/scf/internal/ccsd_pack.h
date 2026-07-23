/// @file cpu/scf/internal/ccsd_pack.h
/// @brief Private pack-based closed-shell CCSD entry.
#pragma once

#include "atomiverse/cpu/scf/ccsd.h"
#include "atomiverse/cpu/scf/internal/integrals.h"
#include "atomiverse/cpu/scf/slater_determinant.h"

namespace ATV {
namespace cpu {
namespace scf {

/// @brief Closed-shell CCSD correlation energy from a converged RHF reference.
///
/// Expands the spatial MO ERI tensor to spin-orbitals, builds Gauss–Stanton
/// intermediates, and iterates T₁/T₂ by plain Jacobi updates until
/// max|Δt| < `tol` or `max_iter` is exhausted.
///
/// @param ints     Integral pack with the dense ERI tensor materialized.
/// @param det      Converged closed-shell determinant (MO coeffs, ε, n_e).
/// @param e_scf_ha Converged SCF total energy [Hartree].
/// @param tol      Amplitude residual tolerance (dimensionless max-norm).
/// @param max_iter Jacobi iteration budget (dimensionless).
/// @return CCSD correlation and total energies [Hartree].
/// @throws std::invalid_argument on open-shell / missing ERI / non-aufbau
///         denominators / non-convergence (no silent fallback).
CCSDResult ccsd_energy(const IntegralPack& ints, const SlaterDeterminant& det,
                       double e_scf_ha, double tol = 1e-10, int max_iter = 100);

} // namespace scf
} // namespace cpu
} // namespace ATV
