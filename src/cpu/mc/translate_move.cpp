/// @file src/cpu/mc/translate_move.cpp
/// @brief cpu::mc::TranslateMove implementation.

#include "atomiverse/cpu/mc/translate_move.h"

#include <random>
#include <stdexcept>

#include "atomiverse/core/precision.h"
#include "atomiverse/cpu/mc/transaction.h"

namespace ATV {
namespace cpu {
namespace mc {

TranslateMove::TranslateMove(double delta) : delta_(delta) {
    if (!(delta > 0.0)) {
        throw std::invalid_argument("TranslateMove: delta must be > 0 Å");
    }
}

void TranslateMove::launch(const cpu::LaunchContext& ctx) {
    // step < 0 is the driver's energy-only bootstrap pass — propose nothing so
    // the evaluator sees the unperturbed initial configuration.
    if (ctx.step < 0)
        return;

    std::mt19937_64* rng = ctx.runtime.rng;
    mc::Transaction* txn = ctx.runtime.proposal;
    // Both are published by SamplingDriver::run(); a null here means the move
    // was driven outside a sampling driver, which is a usage error.
    if (rng == nullptr || txn == nullptr) {
        throw std::runtime_error(
            "TranslateMove::launch: Context::rng / Context::proposal not "
            "published — TranslateMove must run under a SamplingDriver");
    }

    const int n = ctx.frame.atoms.nlocal;
    std::uniform_int_distribution<int> pick(0, n - 1);
    std::uniform_real_distribution<double> disp(-delta_, delta_);

    const int i = pick(*rng);
    const StorageReal ox = ctx.frame.atoms.x[i];
    const StorageReal oy = ctx.frame.atoms.y[i];
    const StorageReal oz = ctx.frame.atoms.z[i];
    const StorageReal nx = ox + static_cast<StorageReal>(disp(*rng));
    const StorageReal ny = oy + static_cast<StorageReal>(disp(*rng));
    const StorageReal nz = oz + static_cast<StorageReal>(disp(*rng));

    ctx.frame.atoms.x[i] = nx;
    ctx.frame.atoms.y[i] = ny;
    ctx.frame.atoms.z[i] = nz;

    txn->clear();
    txn->moved.push_back(static_cast<LocalInt>(i));
    txn->old_x.push_back(ox);
    txn->old_y.push_back(oy);
    txn->old_z.push_back(oz);
    txn->new_x.push_back(nx);
    txn->new_y.push_back(ny);
    txn->new_z.push_back(nz);
}

} // namespace mc
} // namespace cpu
} // namespace ATV
