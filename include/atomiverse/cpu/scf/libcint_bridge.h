/// @file cpu/scf/libcint_bridge.h
/// @brief Minimal STO-3G -> libcint argument bridge for CPU SCF paths.
///
/// Mirrors the basis normalization and atm/bas/env layout used by PySCF's
/// `make_env` / `make_bas_env` helpers for point nuclei and contracted GTOs.
#pragma once

#include <vector>

#include "atomiverse/cpu/frame.h"

namespace ATV {
namespace cpu {
namespace scf {

#ifdef ATV_HAS_LIBCINT

/// @brief Prepared libcint arguments for the current STO-3G molecule.
struct LibcintSto3gContext {
    std::vector<int> atm;            ///< libcint atm table [natm * ATM_SLOTS].
    std::vector<int> bas;            ///< libcint bas table [nshells * BAS_SLOTS].
    std::vector<double> env;         ///< libcint env array (Bohr / normalized coeffs).
    std::vector<int> shell_offsets;  ///< AO offset per shell.
    int natm = 0;                    ///< Number of atoms.
    int nshells = 0;                 ///< Number of libcint shells.
    int nao = 0;                     ///< Number of cartesian AOs (s/p only).
};

/// @brief Build libcint atm/bas/env arrays for the molecule in STO-3G.
///
/// For STO-3G SP shells, the bridge expands one raw shell into two libcint
/// shells: one contracted s shell and one contracted p shell. Because the
/// current basis only uses l <= 1, cartesian and spherical AO dimensions match.
LibcintSto3gContext build_libcint_sto3g_context(const Frame& frame);

#endif

} // namespace scf
} // namespace cpu
} // namespace ATV
