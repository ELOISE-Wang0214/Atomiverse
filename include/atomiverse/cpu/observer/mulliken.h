/// @file cpu/observer/mulliken.h
/// @brief Mulliken population analysis — pure observer (compute only, no I/O).
///
/// Computes Mulliken gross populations by projecting a density matrix onto
/// the atom-centred basis.  The algorithm is:
///
///   N_A = Tr_A(P S)  =  sum_{mu in A} (PS)_{mu mu}
///
/// Charge: q_A = Z_A - N_A
/// Spin:   sigma_A = Tr_A((P^alpha - P^beta) S)
///
/// The functions are stateless and side-effect free: they take the density
/// matrices and overlap matrix and return per-atom vectors.  Analysis is
/// delegated to molrs via the CXX bridge.
///
/// @see cpu::scf::IntegralPack for the S matrix
#pragma once

#include <vector>

#include <xtensor/xtensor.hpp>

#include "atomiverse/core/real.h"
#include "atomiverse/cpu/frame.h"

namespace ATV {
namespace cpu {

/// @brief Compute Mulliken net charges: q_A = Z_A - Tr_A(P S).
///
/// Delegates to molrs::mulliken_charge via CXX bridge.
///
/// @param frame   CPU Frame supplying type_id[] and n_atoms().
/// @param P_flat  Density matrix as flat row-major array [n_basis * n_basis].
/// @param S_flat  Overlap matrix as flat row-major array [n_basis * n_basis].
/// @param n_basis Number of basis functions.
/// @return Per-atom Mulliken net charges [e], size == frame.n_atoms().
std::vector<Real> mulliken_charge(const cpu::Frame& frame,
                                  const Real* P_flat,
                                  const Real* S_flat,
                                  int n_basis);

/// @brief Compute Mulliken spin densities: sigma_A = Tr_A((P^alpha - P^beta) S).
///
/// @param frame     CPU Frame supplying type_id[] and n_atoms().
/// @param Pa_flat   Alpha density matrix as flat row-major array.
/// @param Pb_flat   Beta  density matrix as flat row-major array.
/// @param S_flat    Overlap matrix as flat row-major array.
/// @param n_basis   Number of basis functions.
/// @return Per-atom Mulliken spin densities [e], size == frame.n_atoms().
std::vector<Real> mulliken_spin(const cpu::Frame& frame,
                                const Real* Pa_flat,
                                const Real* Pb_flat,
                                const Real* S_flat,
                                int n_basis);

// ── Overloads accepting xtensor matrices (for CPU-side callers) ──────────────

/// @brief Convenience overload taking xtensor density and overlap matrices.
std::vector<Real> mulliken_charge(const cpu::Frame& frame,
                                  const xt::xtensor<Real, 2>& P,
                                  const xt::xtensor<Real, 2>& S);

/// @brief Convenience overload taking xtensor alpha/beta density matrices.
std::vector<Real> mulliken_spin(const cpu::Frame& frame,
                                const xt::xtensor<Real, 2>& P_alpha,
                                const xt::xtensor<Real, 2>& P_beta,
                                const xt::xtensor<Real, 2>& S);

} // namespace cpu
} // namespace ATV
