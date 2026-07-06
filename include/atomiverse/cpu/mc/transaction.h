/// @file cpu/mc/transaction.h
/// @brief cpu::mc::Transaction — per-move changeset for Monte Carlo trial moves.
#pragma once

#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

#include "atomiverse/core/precision.h"

namespace ATV {
namespace cpu {
namespace mc {

/// @brief Pure-data changeset for one Monte Carlo trial move.
///
/// Records which atoms a Mover displaced and their positions [Å] before
/// (`old_*`) and after (`new_*`) the proposal, so the sampling driver can
/// either commit the move (keep the frame) or roll it back on rejection
/// (restore `old_*`). It is **pure data**: it carries no logic that reaches
/// into the neighbor list or the Frame — consumers derive their own views
/// (CLAUDE.md: "hot-path data structures carry only fundamental physical
/// quantities"). It is **never** a member of `cpu::Frame` (Frame-purity); the
/// driver owns it and publishes a non-owning pointer through
/// `cpu::Context::proposal`.
///
/// The seven arrays are parallel and equal-length: element `k` describes the
/// displacement of atom `moved[k]` from `(old_x,old_y,old_z)[k]` to
/// `(new_x,new_y,new_z)[k]`.
struct Transaction {
    std::vector<LocalInt> moved;    ///< Local indices of the displaced atoms.
    std::vector<StorageReal> old_x; ///< Pre-move x [Å], parallel to `moved`.
    std::vector<StorageReal> old_y; ///< Pre-move y [Å], parallel to `moved`.
    std::vector<StorageReal> old_z; ///< Pre-move z [Å], parallel to `moved`.
    std::vector<StorageReal> new_x; ///< Post-move x [Å], parallel to `moved`.
    std::vector<StorageReal> new_y; ///< Post-move y [Å], parallel to `moved`.
    std::vector<StorageReal> new_z; ///< Post-move z [Å], parallel to `moved`.

    /// @brief Construct an empty changeset (the driver owns it; a Mover fills
    ///        it each `propose` phase).
    Transaction() = default;

    /// @brief Construct a changeset from parallel arrays, validating that every
    ///        position array has the same length as `moved`.
    /// @param moved_ Local indices of the displaced atoms.
    /// @param old_x_ Pre-move x [Å]. @param old_y_ Pre-move y [Å].
    /// @param old_z_ Pre-move z [Å]. @param new_x_ Post-move x [Å].
    /// @param new_y_ Post-move y [Å]. @param new_z_ Post-move z [Å].
    /// @throws std::invalid_argument if any array length differs from
    ///         `moved.size()`.
    Transaction(std::vector<LocalInt> moved_, std::vector<StorageReal> old_x_,
                std::vector<StorageReal> old_y_, std::vector<StorageReal> old_z_,
                std::vector<StorageReal> new_x_, std::vector<StorageReal> new_y_,
                std::vector<StorageReal> new_z_)
        : moved(std::move(moved_)), old_x(std::move(old_x_)),
          old_y(std::move(old_y_)), old_z(std::move(old_z_)),
          new_x(std::move(new_x_)), new_y(std::move(new_y_)),
          new_z(std::move(new_z_)) {
        const std::size_t n = moved.size();
        if (old_x.size() != n || old_y.size() != n || old_z.size() != n
            || new_x.size() != n || new_y.size() != n || new_z.size() != n) {
            throw std::invalid_argument(
                "mc::Transaction: every position array must match moved.size()");
        }
    }

    /// @brief Number of displaced atoms recorded in this changeset.
    /// @return `moved.size()`.
    std::size_t size() const noexcept { return moved.size(); }

    /// @brief Reset to an empty changeset, retaining the backing storage so a
    ///        Mover reuses it across moves without per-move reallocation.
    void clear() noexcept {
        moved.clear();
        old_x.clear();
        old_y.clear();
        old_z.clear();
        new_x.clear();
        new_y.clear();
        new_z.clear();
    }
};

} // namespace mc
} // namespace cpu
} // namespace ATV
