/// @file src/cpu/mc/volume_move.cpp
/// @brief cpu::mc::VolumeMove implementation.

#include "atomiverse/cpu/mc/volume_move.h"

#include <cmath>
#include <random>
#include <stdexcept>

#include "atomiverse/core/precision.h"
#include "atomiverse/cpu/mc/transaction.h"

namespace ATV {
namespace cpu {
namespace mc {

VolumeMove::VolumeMove(double dv_max) : dv_max_(dv_max) {
    if (!(dv_max > 0.0)) {
        throw std::invalid_argument("VolumeMove: dv_max must be > 0 Å³");
    }
}

void VolumeMove::launch(const cpu::LaunchContext& ctx) {
    if (ctx.step < 0)
        return; // energy-only bootstrap pass

    std::mt19937_64* rng = ctx.runtime.mc.rng;
    mc::Transaction* txn = ctx.runtime.mc.proposal;
    if (rng == nullptr || txn == nullptr) {
        throw std::runtime_error(
            "VolumeMove::launch: ctx.mc.rng / ctx.mc.proposal not published — "
            "VolumeMove must run under an MC preset");
    }

    const double v_old = ctx.frame.box.volume();
    std::uniform_real_distribution<double> dv(-dv_max_, dv_max_);
    const double v_new = v_old + dv(*rng);
    if (!(v_new > 0.0)) {
        // A non-positive trial volume has zero Boltzmann weight — record an
        // empty proposal the policy rejects outright.
        txn->begin(Transaction::Verb::None, ctx.frame);
        return;
    }

    const double s = std::cbrt(v_new / v_old);
    txn->begin(Transaction::Verb::Box, ctx.frame);
    txn->record_box(ctx.frame);
    auto& a = ctx.frame.atoms;
    for (int i = 0; i < a.nlocal; ++i) {
        txn->record_touch(ctx.frame, i); // snapshot BEFORE rescaling the row
        a.x[i] = static_cast<StorageReal>(static_cast<double>(a.x[i]) * s);
        a.y[i] = static_cast<StorageReal>(static_cast<double>(a.y[i]) * s);
        a.z[i] = static_cast<StorageReal>(static_cast<double>(a.z[i]) * s);
    }
    ctx.frame.box.scale_isotropic(s);
}

} // namespace mc
} // namespace cpu
} // namespace ATV
