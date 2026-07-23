/// @file cpu/planewave/internal/smearing.h
/// @brief Metallic occupation smearing (Fermi–Dirac / Gaussian / Methfessel–Paxton).
///
/// Spin-unpolarized convention: spinless \(f\in[0,1]\) is doubled to a band
/// occupation in \([0,2]\). Width `degauss` is in the same energy unit as the
/// eigenvalues (Hartree inside the PW SCF).
#pragma once

#include "atomiverse/core/precision.h"

#include <utility>
#include <vector>

namespace ATV {
namespace cpu {
namespace pw {
namespace internal {

/// @brief Smearing kind (mirrors public `cpu::pw::Smearing`, without `none`).
enum class SmearKind : int {
    fermi_dirac = 1,
    methfessel_paxton = 2,
    gaussian = 3,
};

/// @brief Spinless occupation \(f(x)\in[0,1]\), \(x=(\varepsilon-\mu)/\sigma\).
/// @param x        Reduced energy (dimensionless).
/// @param kind     Smearing kind.
/// @param mp_order Methfessel–Paxton order \(N\ge 1\) (ignored otherwise).
/// @return Occupation in \([0,1]\).
AccReal smear_f_spinless(AccReal x, SmearKind kind, int mp_order = 1);

/// @brief Closed-shell band occupation in \([0,2]\): \(2\,f(x)\).
AccReal smear_f_band(AccReal eig, AccReal mu, AccReal degauss, SmearKind kind,
                     int mp_order = 1);

/// @brief Entropy weight \(s(x)\) such that \(TS = \sigma\sum w\,s\) per band
///        with spin degeneracy already folded into \(s\) (band form).
///
/// For Fermi–Dirac: \(s = -2[f\ln f+(1-f)\ln(1-f)]\) with spinless \(f\).
/// For Gaussian/MP: QE-style \(w1gauss\) contribution \(\times 2\).
AccReal smear_entropy_band(AccReal eig, AccReal mu, AccReal degauss,
                           SmearKind kind, int mp_order = 1);

/// @brief Find \(\mu\) so \(\sum_i w_i f(\varepsilon_i,\mu)=n_{\mathrm{elec}}\).
///
/// @param evals   Eigenvalues (same unit as degauss).
/// @param weights Per-eigenvalue weight (k-point weight, one entry per eval).
/// @param n_elec  Target electron count (dimensionless).
/// @param degauss Smearing width (>0).
/// @param kind    Smearing kind.
/// @param mp_order MP order.
/// @return Fermi level \(\mu\).
AccReal find_fermi_level(const std::vector<AccReal>& evals,
                         const std::vector<AccReal>& weights, AccReal n_elec,
                         AccReal degauss, SmearKind kind, int mp_order = 1);

/// @brief Occupations, electron sum, and entropy sum for a spectrum at fixed \(\mu\).
/// @return `{occupations [0,2], n_occ, entropy_sum}` where
///         \(E_{\mathrm{free}} = E_{\mathrm{KS}} - \sigma\cdot\mathrm{entropy\_sum}\).
struct SmearFill {
    std::vector<AccReal> occ; ///< Band occupations in [0, 2].
    AccReal n_occ = 0;        ///< \(\sum w f\).
    AccReal entropy = 0;      ///< Dimensionless \(\sum w s\) (band form).
};

SmearFill fill_occupations(const std::vector<AccReal>& evals,
                           const std::vector<AccReal>& weights, AccReal mu,
                           AccReal degauss, SmearKind kind, int mp_order = 1);

} // namespace internal
} // namespace pw
} // namespace cpu
} // namespace ATV
