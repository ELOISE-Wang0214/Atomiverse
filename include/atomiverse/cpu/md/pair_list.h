/// @file cpu/md/pair_list.h
/// @brief PairList -- LAMMPS-style CSR neighbor list format for CPU backend.
#pragma once

#include <xtensor/containers/xtensor.hpp>

#include "atomiverse/cpu/frame.h"

namespace ATV {
namespace cpu {
namespace md {

/// @brief LAMMPS-style CSR neighbor list (half or full).
///
/// Force components iterate:
/// @code{.cpp}
/// for (int ii = 0; ii < pairs.inum; ++ii) {
///     int i = pairs.ilist[ii];
///     int jnum = pairs.numneigh[i];
///     int off = pairs.firstneigh[i];
///     for (int jj = 0; jj < jnum; ++jj) {
///         int j = pairs.neighbors[off + jj];
///         // compute pair (i,j)
///     }
/// }
/// @endcode
struct PairList {
    int inum = 0;                          ///< Number of I atoms with neighbor lists.
    xt::xtensor<int, 1> ilist;             ///< I atom indices [inum].
    xt::xtensor<int, 1> numneigh;          ///< Number of neighbors per atom [n_atoms].
    xt::xtensor<int, 1> firstneigh;        ///< Offset into neighbors[] per atom [n_atoms].
    xt::xtensor<int, 1> neighbors;         ///< Flat neighbor indices [total_pairs].
    int total_pairs = 0;                   ///< Total entries in neighbors[].

    /// @brief Zero-overhead pair iteration (half list).
    /// @tparam F Callable with signature void(int i, int j).
    template <typename F>
    void for_each_pair(F&& func) const {
        for (int ii = 0; ii < inum; ++ii) {
            int i = ilist[ii];
            int jnum_i = numneigh[i];
            int offset = firstneigh[i];
            for (int jj = 0; jj < jnum_i; ++jj) {
                int j = neighbors[offset + jj];
                func(i, j);
            }
        }
    }

    /// @brief OpenMP-parallel pair iteration.
    /// @tparam F Callable with signature void(int i, int j).
    template <typename F>
    void for_each_pair_omp(F&& func) const {
        #pragma omp parallel for schedule(static)
        for (int ii = 0; ii < inum; ++ii) {
            int i = ilist[ii];
            int jnum_i = numneigh[i];
            int offset = firstneigh[i];
            for (int jj = 0; jj < jnum_i; ++jj) {
                int j = neighbors[offset + jj];
                func(i, j);
            }
        }
    }
};

} // namespace md
} // namespace cpu
} // namespace ATV
