/// @file cpu/md/cell_list.h
/// @brief CellList -- LAMMPS-style linked-cell spatial index for CPU backend.
///
/// O(N) neighbor list construction via cell binning + stencil search.
/// Produces LAMMPS-style CSR PairList (half or full).
#pragma once

#include <vector>
#include <xtensor/xtensor.hpp>

#include "atomiverse/core/int.h"
#include "atomiverse/core/real.h"
#include "atomiverse/cpu/frame.h"
#include "atomiverse/cpu/md/pair_list.h"

namespace ATV {
namespace cpu {
namespace md {

/// @brief Cell list spatial index -- LAMMPS bin style.
///
/// Algorithm:
/// 1. Bin atoms into cells of side >= cutoff + skin
/// 2. Precompute stencil of cell offsets
/// 3. Build CSR half-neighbor list via stencil search
/// 4. Track max displacement for rebuild check
class CellList {
public:
    /// @param cutoff Interaction cutoff [Angstrom].
    /// @param skin Verlet skin width [Angstrom]. Negative = auto (0.3 * cutoff).
    explicit CellList(Real cutoff, Real skin = Real(-1));

    /// @brief One-time setup: compute grid dimensions, allocate cell arrays.
    void init(const cpu::Frame& frame, Real cutoff_skin);

    /// @brief Full rebuild: bin atoms, build stencil.
    void build(const cpu::Frame& frame);

    /// @brief Check if any atom moved more than skin/2.
    bool needs_rebuild(const cpu::Frame& frame) const;

    /// @brief Rebuild if needed, else no-op.
    void update(const cpu::Frame& frame);

    /// @brief Produce LAMMPS-style half-neighbor list (CSR). Newton on.
    PairList build_half(const cpu::Frame& frame) const;

    /// @brief Produce full-neighbor list (CSR). Newton off.
    PairList build_full(const cpu::Frame& frame) const;

private:
    Real cutoff_, skin_, cutoff_skin_;

    // Grid
    Int ncx_ = 0, ncy_ = 0, ncz_ = 0;   ///< Cells per axis.
    Real csx_ = 0, csy_ = 0, csz_ = 0;   ///< Cell dimensions.
    Int total_cells_ = 0;

    // Binning
    xt::xtensor<Int, 1> atom_cell_;       ///< Atom -> cell index [n_atoms].
    xt::xtensor<Int, 1> cell_count_;      ///< Atoms per cell [total_cells].
    xt::xtensor<Int, 1> cell_offset_;     ///< Prefix sum [total_cells].
    xt::xtensor<Int, 1> binned_atoms_;    ///< Sorted atom indices [n_atoms].

    // Stencil
    std::vector<Int> half_stencil_;       ///< Cell offsets for half list.
    std::vector<Int> full_stencil_;       ///< Cell offsets for full list.

    // Displacement tracking
    xt::xtensor<Real, 1> x_save_, y_save_, z_save_;

    void compute_grid(const cpu::Frame& frame);
    void bin_atoms(const cpu::Frame& frame);
    void build_stencil();
    void save_positions(const cpu::Frame& frame);
};

} // namespace md
} // namespace cpu
} // namespace ATV
