/// @file cpu/gpw/internal/cp2k_basis.h
/// @brief Parser for CP2K Gaussian basis-set data files (GPW spec D1).
///
/// GPW keeps the Kohn-Sham problem in CP2K's contracted Gaussian AO basis, so
/// it must ingest the *exact* exponents and contraction coefficients CP2K uses
/// (DZVP-GTH-PADE for the M1 LDA target, DZVP-MOLOPT-GTH for the H2O tutorial).
/// This reader is the basis-set companion to the GTH potential reader in
/// `cp2k_io.h`; it mirrors that file's style (std-only, throw on error).
///
/// File grammar (one entry, e.g. data/BASIS_SET, data/BASIS_MOLOPT):
///   <Element> <name> [<alias> ...]
///       <nset>
///       <n> <lmin> <lmax> <nexp> <nshell_lmin> ... <nshell_lmax>
///            <exp_1> <coef_1...ncoef>     # nexp rows, ncoef = sum_l nshell_l
///            ...                          #   columns ordered by l then shell
///       <next set ...>
/// '#' lines separate entries.
#pragma once

#include <cstddef>
#include <string>
#include <vector>

namespace ATV {
namespace cpu {
namespace gpw {

/// @brief One contraction "set" of a CP2K Gaussian basis (shared exponents).
///
/// A set groups all angular momenta from `lmin` to `lmax` that share the same
/// `exps`. The contraction matrix has one row per exponent and one column per
/// contracted shell, columns ordered by angular momentum then shell index
/// (s-shells first, then p, then d, ...), exactly as laid out in the file.
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
    std::string element;              ///< chemical symbol, e.g. "O".
    std::string name;                 ///< matched name, e.g. "DZVP-MOLOPT-GTH".
    std::vector<GaussianBasisSet> sets;
};

/// @brief Parse one Gaussian basis entry from a CP2K basis-set file.
///
/// Matches the entry whose header line starts with `element` and lists `name`
/// among its (alias) names. Throws std::runtime_error if the file cannot be
/// opened or no matching entry is found.
///
/// @param path    Path to the CP2K basis file (e.g. BASIS_SET, BASIS_MOLOPT).
/// @param element Chemical symbol (e.g. "Ar", "O", "H").
/// @param name    Basis name/alias (e.g. "DZVP-GTH-PADE", "DZVP-MOLOPT-GTH").
/// @return Parsed GaussianBasis.
GaussianBasis parse_gth_basis(const std::string& path, const std::string& element,
                              const std::string& name);

} // namespace gpw
} // namespace cpu
} // namespace ATV
