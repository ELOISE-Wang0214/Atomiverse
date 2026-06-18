/// @file cpu/scf/scf_linalg.h
/// @brief SCF linear-algebra facade: re-exports the shared dense-LA primitives
///        from `cpu::linalg` into `cpu::scf`, plus the SCF-specific helpers that
///        do not belong in the generic module (broken-symmetry guess, spin
///        contamination).
///
/// The method-agnostic primitives (orthogonalization, generalized eigensolve,
/// density-from-MO, RealMatrix) were lifted to `cpu/linalg/dense.h` so the GPW
/// path can share them; this header keeps existing `cpu::scf` call sites working
/// unchanged. (Replaces the former, vaguely-named `scf_utils.h`.)
#pragma once

#include <cmath>

#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/views/xview.hpp>

#include "atomiverse/cpu/linalg/dense.h"

namespace ATV {
namespace cpu {
namespace scf {

// Shared dense LA, re-exported so `cpu::scf` code uses these names directly.
using linalg::build_density_matrix;
using linalg::build_orthogonalization_matrix;
using linalg::diagonalize_fock;
using linalg::diagonalize_fock_generalized;
using linalg::RealMatrix;
using linalg::RealVector;

/// @brief Apply HOMO-LUMO mixing to break spin symmetry in the initial guess.
///
/// Rotates the HOMO and LUMO columns of C by `theta`, creating an asymmetric
/// starting density that can locate broken-symmetry solutions (e.g. stretched
/// H2). SCF-specific; not part of the generic LA module.
inline void apply_homo_lumo_mixing(RealMatrix& C, int n_occ, int n,
                                   double theta = 0.3490658503988659) {
    if (n_occ <= 0 || n_occ >= n) return;
    const double ca = std::cos(theta);
    const double sa = std::sin(theta);
    const int homo = n_occ - 1;
    const int lumo = n_occ;
    for (int r = 0; r < n; ++r) {
        const double ci = C(r, homo);
        const double cj = C(r, lumo);
        C(r, homo) = ca * ci + sa * cj;
        C(r, lumo) = -sa * ci + ca * cj;
    }
}

/// @brief Spin contamination <S^2> - S(S+1) (zero for pure spin states).
inline double compute_spin_contamination(const RealMatrix& C_alpha,
                                         const RealMatrix& C_beta,
                                         const RealMatrix& S, int n_alpha,
                                         int n_beta) {
    if (n_alpha == 0 || n_beta == 0) return 0.0;
    auto Ca = xt::view(C_alpha, xt::all(), xt::range(0, n_alpha));
    auto Cb = xt::view(C_beta, xt::all(), xt::range(0, n_beta));
    RealMatrix SCb = xt::linalg::dot(S, Cb);
    RealMatrix O = xt::linalg::dot(xt::transpose(Ca), SCb);
    return static_cast<double>(n_beta) - xt::sum(O * O)();
}

} // namespace scf
} // namespace cpu
} // namespace ATV
