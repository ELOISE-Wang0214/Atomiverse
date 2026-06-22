/// @file cpu/gpw/internal/dftd3_params.h
/// @brief Grimme DFT-D3 reference parameter tables (zero-damping), parsed to
///        match CP2K's `qs_dispersion_*` exactly.
///
/// Three sources feed the D3 model (Grimme, Antony, Ehrlich, Krieg,
/// J. Chem. Phys. 132, 154104 (2010)):
///   1. The tabulated, coordination-number-resolved reference C6 coefficients,
///      read from CP2K's `dftd3.dat` (same file, same encoding as
///      `qs_dispersion_d3::dftd3_c6_param`).
///   2. The element covalent radii `rcov` (Pyykkö & Atsumi, Chem. Eur. J. 15,
///      2009, 188; metals reduced 10%), pre-scaled by k2 = 4/3 and converted to
///      Bohr — exactly the `dispersion_env%rcov = k2*rcov*bohr` of CP2K.
///   3. The `r2r4` values (PBE0/def2-QZVP) post-processed to
///      `sqrt(0.5*r2r4_raw*sqrt(Z))` so that C8 = 3*C6*r2r4_A*r2r4_B, and the
///      packed lower-triangular `r0ab(za,zb)` cutoff-radius table [Bohr].
///
/// All tables cover Z = 1..94 (`kMaxElem`). This header exposes the parsed,
/// ready-to-use parameters; physics lives in `dftd3.h`.
#pragma once

#include <array>
#include <string>
#include <vector>

namespace ATV {
namespace cpu {
namespace gpw {

/// @brief Largest atomic number supported by the D3 tables.
constexpr int kMaxElem = 94;

/// @brief Reference-state entry of the CN-resolved C6 table for one element pair.
struct D3C6Ref {
    double c6 = 0.0;  ///< Reference C6^AB [a.u.] for this (CN_A, CN_B) pair.
    double cn_a = 0.0; ///< Reference coordination number of A.
    double cn_b = 0.0; ///< Reference coordination number of B.
};

/// @brief Parsed DFT-D3 parameter set (everything geometry-independent).
///
/// `rcov`, `r2r4` and `r0ab` are indexed by atomic number Z (1-based; index 0
/// unused). `c6ab[za][zb]` is the list of reference states for that ordered
/// element pair (symmetric: both (za,zb) and (zb,za) are populated).
struct D3Params {
    /// @brief k2-folded covalent radii [Bohr], rcov[Z] = (4/3)*rcov_Pyykko*bohr.
    std::array<double, kMaxElem + 1> rcov{};
    /// @brief Processed r2r4 [a.u.], r2r4[Z] = sqrt(0.5*r2r4_raw*sqrt(Z)).
    std::array<double, kMaxElem + 1> r2r4{};
    /// @brief Cutoff radii r0ab[Za][Zb] [Bohr] (symmetric).
    std::array<std::array<double, kMaxElem + 1>, kMaxElem + 1> r0ab{};
    /// @brief CN-resolved reference C6 states per ordered element pair.
    std::array<std::array<std::vector<D3C6Ref>, kMaxElem + 1>, kMaxElem + 1> c6ab{};
};

/// @brief Build the full D3 parameter set by parsing CP2K's `dftd3.dat`.
///
/// The covalent-radius, r2r4 and r0ab tables are compiled-in (identical to
/// CP2K's `setrcov`/`setr0ab`); only the CN-resolved C6 coefficients are read
/// from file. The file uses CP2K's encoding: a header line `nl nlines`, then
/// `nlines` records of 5 reals `(C6, iat, jat, cnA, cnB)` where `iat`/`jat`
/// encode `Z + 1000*(ref_index-1)` (the 1000-stride form used by all current
/// `dftd3.dat`).
///
/// @param dftd3_dat_path Absolute path to CP2K's `dftd3.dat`.
/// @return Fully populated `D3Params`.
/// @throw std::runtime_error if the file cannot be opened or is malformed.
D3Params load_dftd3_params(const std::string& dftd3_dat_path);

} // namespace gpw
} // namespace cpu
} // namespace ATV
