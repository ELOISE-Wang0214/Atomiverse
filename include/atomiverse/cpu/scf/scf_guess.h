/// @file cpu/scf/scf_guess.h
/// @brief SCF initial-guess perturbations. Broken-symmetry HOMO-LUMO mixing for
///        locating spin-broken solutions (e.g. stretched H2).
#pragma once

#include <cmath>

#include "atomiverse/cpu/scf/dense_matrix.h"

namespace ATV {
namespace cpu {
namespace scf {

/// @brief Rotate the HOMO and LUMO columns of C by `theta` to seed a broken-
///        symmetry guess (an asymmetric starting density).
/// @param C MO coefficient matrix [n x n], modified in place.
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

} // namespace scf
} // namespace cpu
} // namespace ATV
