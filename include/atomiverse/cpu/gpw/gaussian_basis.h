/// @file cpu/gpw/gaussian_basis.h
/// @brief GaussianBasis — contracted Gaussian AO basis, constructed by
///        (Element, name) lookup into the in-tree static tables
///        (basis_data/gpw_tables.h). No file or string parsing anywhere.
#pragma once

#include <cstddef>
#include <string>
#include <string_view>
#include <vector>

#include "atomiverse/molrs.hpp"

namespace ATV {
namespace cpu {
namespace gpw {

/// @brief One contraction "set" of a Gaussian basis (shared exponents).
///
/// A set groups all angular momenta from `lmin` to `lmax` that share the same
/// `exps`. The contraction matrix has one row per exponent and one column per
/// contracted shell, columns ordered by angular momentum then shell index
/// (s-shells first, then p, then d, ...).
struct GaussianBasisSet {
    int n = 0;                  ///< principal quantum number label.
    int lmin = 0;               ///< lowest angular momentum in the set.
    int lmax = 0;               ///< highest angular momentum in the set.
    std::vector<double> exps;   ///< Gaussian exponents (nexp).
    std::vector<int> nshell;    ///< contracted shells per l, indexed [l - lmin].
    int ncoef = 0;              ///< number of contraction columns = sum(nshell).
    std::vector<double> coef;   ///< nexp*ncoef contraction matrix, row-major.

    /// @brief Number of primitive Gaussian exponents in this set.
    int nexp() const { return static_cast<int>(exps.size()); }

    /// @brief Contraction coefficient for exponent `iexp`, column `icol`.
    double coef_at(int iexp, int icol) const {
        return coef[static_cast<std::size_t>(iexp) * ncoef + icol];
    }

    /// @brief Number of contracted shells for angular momentum `l`.
    int nshell_l(int l) const { return nshell[static_cast<std::size_t>(l - lmin)]; }
};

/// @brief A complete contracted Gaussian basis for one element/name.
struct GaussianBasis {
    Element element{};                  ///< chemical element (molrs enum, Z = value).
    std::string name;                   ///< basis name, e.g. "DZVP-MOLOPT-GTH".
    std::vector<GaussianBasisSet> sets; ///< contraction sets.

    GaussianBasis() = default;

    /// @brief Look up a published basis set from the in-tree table.
    ///
    /// Provenance: MOLOPT — VandeVondele & Hutter, JCP 127, 114105 (2007);
    /// GTH-PADE bases — the CP2K BASIS_SET distribution. Throws
    /// std::runtime_error naming the supported entries when (element, name) is
    /// not in the table — no fallback.
    /// @param element Chemical element (e.g. `Element::O`).
    /// @param name    Basis name (e.g. "DZVP-MOLOPT-SR-GTH").
    GaussianBasis(Element element, std::string_view name);
};

} // namespace gpw
} // namespace cpu
} // namespace ATV
