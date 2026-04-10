/// @file cpu/scf/rys_roots.h
/// @brief Rys quadrature roots and weights — public API for the native
///        Rys-quadrature ERI engine.
///
/// Implementation is a namespaced copy of libcint's rys_roots.c (BSD-2,
/// Copyright (c) 2012 Qiming Sun). Only the closed-form rys_root1..5
/// branches are retained; d-shells (Ltot ≤ 8) need nroots ≤ 5.
/// Higher angular momenta would require the iterative Jacobi/Laguerre
/// fallbacks (not ported).
///
/// @see rys_roots.cpp for the ported implementation.
/// @see rys_eri.cpp   for the caller (2D-VRR + HRR + primitive loop).
#pragma once

namespace ATV {
namespace cpu {
namespace scf {
namespace native {

/// @brief Compute Rys quadrature roots u[i] and weights w[i] for the Boys
///        integrand t^{2m} exp(-x t^2).
///
/// @details The Rys roots are the zeros of orthogonal polynomials with respect
///   to the weight exp(-x t^2) on [0, 1]. They convert the Boys integral
///   evaluation into a discrete sum of products over `nroots` samples,
///   which is the algorithmic foundation of the Rys-quadrature ERI scheme
///   used by libcint / PySCF.
///
///   Output `u[i]` stores t_i^2 / (1 - t_i^2) (libcint convention), not the
///   raw Rys root t_i.
///
/// @param nroots Number of roots (1..5 supported; sufficient for Ltot ≤ 8,
///        i.e. integral class (dd|dd) and below).
/// @param x      Boys parameter T = ρ·|P-Q|² (non-negative).
/// @param u      Output array [nroots] — transformed roots t_i²/(1-t_i²).
/// @param w      Output array [nroots] — quadrature weights.
void rys_roots(int nroots, double x, double *u, double *w);

} // namespace native
} // namespace scf
} // namespace cpu
} // namespace ATV
