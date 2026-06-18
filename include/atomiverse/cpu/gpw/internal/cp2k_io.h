/// @file cpu/gpw/internal/cp2k_io.h
/// @brief Parsers for CP2K's GTH pseudopotential and Gaussian basis data files.
///
/// GPW must consume the *exact* GTH parameters and contracted-GTO coefficients
/// CP2K uses, so that Atomiverse reproduces CP2K Quickstep results (see
/// `.claude/specs/gpw-gaussian-plane-wave.md`, design item D1). This header
/// declares the data structures and the readers for CP2K's `GTH_POTENTIALS` /
/// `POTENTIAL` files (analytic Goedecker-Teter-Hutter form). The Gaussian
/// basis-file reader is added alongside in a later task.
///
/// File grammar (one entry):
///   <Element> <name> [<alias> ...]
///       <n_l0> <n_l1> ...                 # valence occupation per l; sum = Zion
///       <rloc> <nexp> <C_1> ... <C_nexp>  # local part
///       <nchannels>                       # number of nonlocal projector l-channels
///       <r_l> <nproj> <h row 0 ...>       # per channel: radius, #projectors,
///       <h row 1 ...>                     #   symmetric h matrix, upper triangle
///       ...                               #   spread over nproj lines
#pragma once

#include <string>
#include <vector>

namespace ATV {
namespace cpu {
namespace gpw {

/// @brief One nonlocal projector angular-momentum channel of a GTH potential.
struct GthNonlocalChannel {
    int l = 0;             ///< angular momentum (0=s, 1=p, ...), by channel order.
    double r = 0.0;        ///< projector radius r_l [Bohr].
    int nproj = 0;         ///< number of projectors in this channel.
    std::vector<double> h; ///< nproj*nproj symmetric coupling matrix, row-major.

    /// @brief h(i,j) accessor into the row-major dense matrix.
    double hij(int i, int j) const { return h[static_cast<std::size_t>(i) * nproj + j]; }
};

/// @brief Analytic GTH pseudopotential for one element/parameterization.
///
/// Local part: V_loc(r) = -Zion/r erf(r/(sqrt(2) rloc))
///   + exp(-(r/rloc)^2/2) * sum_i C[i] (r/rloc)^(2(i-1)).
/// Nonlocal part: separable Gaussian projectors with couplings `channels[*].h`.
struct GthPotential {
    std::string element;             ///< chemical symbol, e.g. "O".
    std::string name;                ///< matched name, e.g. "GTH-PBE-q6".
    int zion = 0;                    ///< number of valence electrons (= occupation sum).
    double rloc = 0.0;               ///< local Gaussian radius [Bohr].
    std::vector<double> cloc;        ///< local polynomial coefficients C_1..C_nexp.
    std::vector<GthNonlocalChannel> channels; ///< nonlocal projector channels.
};

/// @brief Parse one GTH potential entry from a CP2K `GTH_POTENTIALS`/`POTENTIAL` file.
///
/// Matches the entry whose header line starts with `element` and lists `name`
/// among its (alias) names. Throws std::runtime_error if the file cannot be
/// opened or no matching entry is found.
///
/// @param path    Path to the CP2K potential file.
/// @param element Chemical symbol (e.g. "O", "Ar", "H").
/// @param name    Potential name/alias (e.g. "GTH-PBE-q6", "GTH-PADE-q8").
/// @return Parsed GthPotential.
GthPotential parse_gth_potential(const std::string& path, const std::string& element,
                                 const std::string& name);

} // namespace gpw
} // namespace cpu
} // namespace ATV
