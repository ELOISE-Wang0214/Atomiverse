/// @file cpu/md/cell_list.h
/// @brief CellList -- LAMMPS-style linked-cell spatial index for CPU backend.
///
/// O(N) neighbor list construction via cell binning + stencil search.
/// Produces LAMMPS-style CSR PairList (half or full).
///
/// **Multi-rank ghost-aware mode** (`init_extended`):
/// the bin grid is widened by `cutoff_skin` on each side along every axis
/// `d` with `dims[d] > 1`, so ghost atoms whose positions lie outside the
/// box (after PBC pre-shift by `HaloExchangeForward`) bin into dedicated
/// halo cells instead of wrapping back onto owner cells. Iteration bounds
/// are explicit: `bin_atoms` walks `[0, n_active)` (= owners + current
/// ghosts), and the half-list outer loop walks `[0, n_active)` so callers
/// can later trim `inum = n_local` for force kernels.
#pragma once

#include <vector>
#include <xtensor/containers/xtensor.hpp>

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
    explicit CellList(double cutoff, double skin = -1.0);

    /// @brief One-time setup (single-rank): grid spans `[0, L)` per axis,
    /// ghost positions are wrapped during binning.
    /// @param frame Frame whose box defines the grid.
    /// @param cutoff_skin Interaction cutoff plus Verlet skin [Angstrom].
    void init(const cpu::Frame& frame, double cutoff_skin);

    /// @brief One-time setup (multi-rank): grid is widened by `cutoff_skin`
    /// on each side along every axis `d` with `dims[d] > 1`. Ghost atoms
    /// whose positions lie outside the box are binned into the halo cells
    /// without PBC wrap-around. Single-rank axes (`dims[d] == 1`) keep the
    /// `[0, L)` range and are wrapped during binning.
    /// @param frame Frame whose box defines the grid.
    /// @param cutoff_skin Interaction cutoff plus Verlet skin [Angstrom].
    /// @param dims MPI Cartesian dimensions per axis (≥1).
    /// @param capacity Maximum number of active atoms (owners + ghost slack);
    ///        per-atom bookkeeping arrays are sized to this.
    void init_extended(const cpu::Frame& frame, double cutoff_skin,
                        const int dims[3], int capacity);

    /// @brief Full rebuild: bin atoms, build stencil.
    /// @param frame Frame to bin.
    /// @param n_active Number of active atoms (owners + ghosts) to iterate.
    ///        `-1` means `frame.atoms.size()` (single-rank legacy).
    void build(const cpu::Frame& frame, int n_active = -1);

    /// @brief Check if any atom moved more than skin/2.
    /// @param frame Frame whose positions are compared to last `build()`.
    /// @param n_active Number of active atoms to scan; `-1` = `frame.atoms.size()`.
    bool needs_rebuild(const cpu::Frame& frame, int n_active = -1) const;

    /// @brief Rebuild if needed, else no-op.
    /// @param frame Frame to potentially rebin.
    /// @param n_active Number of active atoms; `-1` = `frame.atoms.size()`.
    void update(const cpu::Frame& frame, int n_active = -1);

    /// @brief Produce LAMMPS-style half-neighbor list (CSR). Newton on.
    /// @param frame Frame for distance evaluation.
    /// @param n_active Outer-loop / array-size bound; `-1` = `frame.atoms.size()`.
    ///        Caller (e.g. `NeighborList`) typically trims `inum = n_local` after.
    PairList build_half(const cpu::Frame& frame, int n_active = -1) const;

    /// @brief Produce full-neighbor list (CSR). Newton off.
    /// @param frame Frame for distance evaluation.
    /// @param n_active Outer-loop / array-size bound; `-1` = `frame.atoms.size()`.
    PairList build_full(const cpu::Frame& frame, int n_active = -1) const;

private:
    double cutoff_, skin_, cutoff_skin_;

    // Grid
    int ncx_ = 0, ncy_ = 0, ncz_ = 0;   ///< Cells per axis.
    double csx_ = 0, csy_ = 0, csz_ = 0;   ///< Cell dimensions.
    int total_cells_ = 0;
    /// @brief Grid origin per axis [Angstrom]. `0.0` for single-rank axes;
    /// `-cutoff_skin` for axes extended to cover ghost halo.
    double origin_[3] = {0.0, 0.0, 0.0};
    /// @brief Per-axis flag: `true` if grid is extended for multi-rank ghosts
    /// (skip PBC wrap during binning); `false` for single-rank (wrap into `[0, L)`).
    bool axis_extended_[3] = {false, false, false};
    /// @brief Bookkeeping capacity. Per-atom arrays sized to this; ≥ any
    /// `n_active` ever passed.
    int capacity_ = 0;

    // Binning
    xt::xtensor<int, 1> atom_cell_;       ///< Atom -> cell index [capacity].
    xt::xtensor<int, 1> cell_count_;      ///< Atoms per cell [total_cells].
    xt::xtensor<int, 1> cell_offset_;     ///< Prefix sum [total_cells].
    xt::xtensor<int, 1> binned_atoms_;    ///< Sorted atom indices [capacity].

    // Stencil
    std::vector<int> half_stencil_;       ///< Cell offsets for half list.
    std::vector<int> full_stencil_;       ///< Cell offsets for full list.

    // Displacement tracking
    xt::xtensor<double, 1> x_save_, y_save_, z_save_;

    void compute_grid(const cpu::Frame& frame);
    void bin_atoms(const cpu::Frame& frame, int n_active);
    void build_stencil();
    void save_positions(const cpu::Frame& frame, int n_active);
};

} // namespace md
} // namespace cpu
} // namespace ATV
