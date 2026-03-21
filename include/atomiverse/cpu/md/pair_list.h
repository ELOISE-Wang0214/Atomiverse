/// @file cpu/md/pair_list.h
/// @brief PairList -- LAMMPS-style CSR neighbor list format for CPU backend.
#pragma once

#include <xtensor/xtensor.hpp>

#include "atomiverse/core/int.h"
#include "atomiverse/core/real.h"
#include "atomiverse/cpu/frame.h"

namespace ATV {
namespace cpu {
namespace md {

/// @brief LAMMPS-style CSR neighbor list (half or full).
///
/// Force components iterate:
/// @code{.cpp}
/// for (Int ii = 0; ii < pairs.inum; ++ii) {
///     Int i = pairs.ilist[ii];
///     Int jnum = pairs.numneigh[i];
///     Int off = pairs.firstneigh[i];
///     for (Int jj = 0; jj < jnum; ++jj) {
///         Int j = pairs.neighbors[off + jj];
///         // compute pair (i,j)
///     }
/// }
/// @endcode
struct PairList {
    Int inum = 0;                          ///< Number of I atoms with neighbor lists.
    xt::xtensor<Int, 1> ilist;             ///< I atom indices [inum].
    xt::xtensor<Int, 1> numneigh;          ///< Number of neighbors per atom [n_atoms].
    xt::xtensor<Int, 1> firstneigh;        ///< Offset into neighbors[] per atom [n_atoms].
    xt::xtensor<Int, 1> neighbors;         ///< Flat neighbor indices [total_pairs].
    Int total_pairs = 0;                   ///< Total entries in neighbors[].

    /// @brief Zero-overhead pair iteration (half list).
    /// @tparam F Callable with signature void(Int i, Int j).
    template <typename F>
    void for_each_pair(F&& func) const {
        for (Int ii = 0; ii < inum; ++ii) {
            Int i = ilist[ii];
            Int jnum_i = numneigh[i];
            Int offset = firstneigh[i];
            for (Int jj = 0; jj < jnum_i; ++jj) {
                Int j = neighbors[offset + jj];
                func(i, j);
            }
        }
    }

    /// @brief OpenMP-parallel pair iteration.
    /// @tparam F Callable with signature void(Int i, Int j).
    template <typename F>
    void for_each_pair_omp(F&& func) const {
        #pragma omp parallel for schedule(static)
        for (Int ii = 0; ii < inum; ++ii) {
            Int i = ilist[ii];
            Int jnum_i = numneigh[i];
            Int offset = firstneigh[i];
            for (Int jj = 0; jj < jnum_i; ++jj) {
                Int j = neighbors[offset + jj];
                func(i, j);
            }
        }
    }
};

} // namespace md
} // namespace cpu
} // namespace ATV
