/// @file cpu/dynamics/fix_freeze.cpp
/// @brief FixFreeze — zero force / freeze position and velocity by atom type.

#include "atomiverse/cpu/dynamics/fix_freeze.h"

#include <stdexcept>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/cpu/component/owned_range.h"

namespace ATV {
namespace cpu {
namespace md {

FixFreeze::FixFreeze(int atom_type) : type_(atom_type) {
    if (atom_type < 0)
        throw std::invalid_argument("FixFreeze: atom_type must be >= 0");
}

int FixFreeze::n_constraints(const cpu::Context& /*ctx*/,
                             const cpu::Frame& /*frame*/) const {
    // Three Cartesian constraints per frozen atom (global count cached at build).
    return 3 * n_frozen_global_;
}

void FixFreeze::build(const cpu::BuildContext& ctx) {
    const int n = owned_count(ctx);
    const auto type = ctx.frame.span(kAtom_type);
    const auto x = ctx.frame.span(kAtom_x);
    const auto y = ctx.frame.span(kAtom_y);
    const auto z = ctx.frame.span(kAtom_z);

    local_idx_.clear();
    x0_.clear();
    y0_.clear();
    z0_.clear();
    local_idx_.reserve(static_cast<std::size_t>(n));
    x0_.reserve(static_cast<std::size_t>(n));
    y0_.reserve(static_cast<std::size_t>(n));
    z0_.reserve(static_cast<std::size_t>(n));

    for (int i = 0; i < n; ++i) {
        if (type[i] == type_) {
            local_idx_.push_back(i);
            x0_.push_back(x[i]);
            y0_.push_back(y[i]);
            z0_.push_back(z[i]);
        }
    }
    // Single-rank v1: global frozen count = local. MPI Allreduce is a follow-up.
    n_frozen_global_ = static_cast<int>(local_idx_.size());
}

void FixFreeze::launch(const cpu::LaunchContext& ctx) {
    auto& f = ctx.frame;
    StorageReal* __restrict vx = f.span(kAtom_vx);
    StorageReal* __restrict vy = f.span(kAtom_vy);
    StorageReal* __restrict vz = f.span(kAtom_vz);
    StorageReal* __restrict x = f.span(kAtom_x);
    StorageReal* __restrict y = f.span(kAtom_y);
    StorageReal* __restrict z = f.span(kAtom_z);
    AccReal* __restrict fx = f.span(kAtom_fx);
    AccReal* __restrict fy = f.span(kAtom_fy);
    AccReal* __restrict fz = f.span(kAtom_fz);

    if (ctx.at(ATV::md::Phase::post_force)) {
        for (int k = 0; k < static_cast<int>(local_idx_.size()); ++k) {
            const int i = local_idx_[static_cast<std::size_t>(k)];
            fx[i] = AccReal(0);
            fy[i] = AccReal(0);
            fz[i] = AccReal(0);
        }
        return;
    }

    // Integrate halves: restore positions + kill velocities (after VV / NH).
    if (ctx.at(ATV::md::Phase::initial_integrate) ||
        ctx.at(ATV::md::Phase::final_integrate)) {
        for (int k = 0; k < static_cast<int>(local_idx_.size()); ++k) {
            const int i = local_idx_[static_cast<std::size_t>(k)];
            x[i] = x0_[static_cast<std::size_t>(k)];
            y[i] = y0_[static_cast<std::size_t>(k)];
            z[i] = z0_[static_cast<std::size_t>(k)];
            vx[i] = StorageReal(0);
            vy[i] = StorageReal(0);
            vz[i] = StorageReal(0);
        }
    }
}

} // namespace md
} // namespace cpu
} // namespace ATV
