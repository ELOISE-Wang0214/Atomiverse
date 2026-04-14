/// @file cpu/scf/rys_eri.h
/// @brief Native Rys-quadrature electron repulsion integrals.
///
/// This engine replicates the algorithm used by libcint / PySCF for
/// evaluating two-electron integrals (ab|cd) over contracted Cartesian
/// Gaussians. It is algorithmically equivalent to libcint (same Rys
/// roots, same 2D-VRR, same HRR) and is intended to supersede the
/// McMurchie-Davidson path in native_integrals.cpp for all s/p/d shells.
///
/// Data flow per shell quartet (A B | C D):
///   1. For each primitive-pair quartet (p,q,r,s):
///        a. Compute bra/ket Gaussian-product centres P, Q and the Boys
///           parameter x = ρ·|P-Q|².
///        b. Call native::rys_roots(nroots, x, u, w) — closed-form
///           piecewise polynomial approximation ported from libcint.
///        c. Per root, form the 9 Rys coefficients b00, b10, b01,
///           c00{xyz}, c0p{xyz}.
///        d. Build the 2D Rys intermediate gx[i][k], gy[i][k], gz[i][k]
///           via the standard Obara–Saika-style recurrence along the
///           bra (i) and ket (k) axes — see `CINTg0_2e_2d` in libcint.
///   2. Apply the 2D→4D HRR to shift angular momentum from
///      (li+lj, 0 | lk+ll, 0) → (li, lj | lk, ll). HRR is a pure
///      geometric transform and does not depend on the primitive
///      quartet.
///   3. For each Cartesian component tuple (a, b, c, d), assemble the
///      contracted ERI as
///          ERI(a,b,c,d) = Σ_root gx(...) · gy(...) · gz(...)
///      and scatter into `out` with the primitive coefficient product.
///
/// Public API: one function, `rys_eri_block`. The output buffer layout
/// is identical to `native::eri_block`:
///
///     out[ia*nb*nc*nd + ib*nc*nd + ic*nd + id] += (ia ib | ic id)
///
/// @see rys_roots.h for the ported root/weight generator.
/// @see native_integrals.h for the ShellBlock data structure.
#pragma once

#include "atomiverse/cpu/scf/native_integrals.h"

namespace ATV {
namespace cpu {
namespace scf {
namespace native {

/// @brief Maximum supported angular momentum per shell (inclusive).
/// l = 0 (s), 1 (p), 2 (d). Shell quartets with higher angular momentum
/// need nroots > 5 and the iterative Jacobi / Laguerre fallback paths,
/// which are not currently ported.
constexpr int kRysMaxL = 2;

/// @brief Maximum Rys roots supported (d-d-d-d has Ltot=8 → nroots=5).
constexpr int kRysMaxRoots = 5;

/// @brief Compute the full Cartesian ERI block (AB|CD) using Rys
///        quadrature.
///
/// The output buffer is accumulated, not overwritten — callers can
/// pre-zero once and build many block contributions into the same
/// buffer. Input shells may have any angular momentum up to @ref
/// kRysMaxL and any contraction length supported by `ShellBlock`
/// (alpha/coeff inline arrays).
///
/// The function is algorithmically equivalent to libcint's
/// `cint2e_cart` for point-nucleus, non-range-separated Coulomb.
///
/// @param A,B Shell blocks for the bra (electron 1), in that order.
/// @param C,D Shell blocks for the ket (electron 2), in that order.
/// @param out Pre-zeroed flat buffer of size
///            A.n_cart·B.n_cart·C.n_cart·D.n_cart.
void rys_eri_block(const ShellBlock& A, const ShellBlock& B,
                   const ShellBlock& C, const ShellBlock& D,
                   double* out);

/// @brief Block-level derivative of the ERI w.r.t. one nuclear centre.
///
/// Computes ∂(AB|CD)/∂R_X — all three Cartesian components — into an
/// axis-split output buffer. Uses a single extended 2D Rys grid per
/// primitive quartet (bra or ket direction extended by +1 angular
/// momentum on the deriv side) and a nabla recurrence
///   f[i] = i·g[i-1] - 2α·g[i+1]
/// to read out the derivative from the already-built grid. All three
/// Cartesian axes share the same grid; this is the same algorithmic
/// skeleton used by libcint's `int2e_ip1`.
///
/// @param A,B,C,D Shell blocks in the conventional ERI ordering (AB|CD).
/// @param deriv_center 0 = ∂/∂R_A, 1 = ∂/∂R_B, 2 = ∂/∂R_C, 3 = ∂/∂R_D.
/// @param out_axes Pre-zeroed flat buffer of size
///        3 · A.n_cart · B.n_cart · C.n_cart · D.n_cart, interpreted as
///        `out_axes[axis * (na*nb*nc*nd) + ia*nbcd + ib*ncd + ic*nd + id]`.
///        Output is accumulated, not overwritten (mirrors
///        `rys_eri_block`).
///
/// @note Only Cartesian Gaussians up to l = 2 (s/p/d) are supported.
///       STO-3G uses only l ∈ {0, 1}, so a derivative adds at most one
///       angular step and stays within the `kRysMaxRoots = 5` budget.
void rys_eri_block_deriv(const ShellBlock& A, const ShellBlock& B,
                          const ShellBlock& C, const ShellBlock& D,
                          int deriv_center,
                          double* out_axes);

} // namespace native
} // namespace scf
} // namespace cpu
} // namespace ATV
