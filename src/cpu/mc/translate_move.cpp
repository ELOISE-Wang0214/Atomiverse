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

    std::mt19937_64* rng = ctx.runtime.mc.rng;
    mc::Transaction* txn = ctx.runtime.mc.proposal;
    // Both are published by cpu::mc::MC::run(); a null here means the move was driven
    // outside an MC run, which is a usage error.
    if (rng == nullptr || txn == nullptr) {
        throw std::runtime_error(
            "TranslateMove::launch: ctx.mc.rng / ctx.mc.proposal not "
            "published — TranslateMove must run under cpu::mc::MC");
    }

    const int n = ctx.frame.atoms.nlocal;
    std::uniform_int_distribution<int> pick(0, n - 1);
    std::uniform_real_distribution<double> disp(-delta_, delta_);

    const int i = pick(*rng);
    const StorageReal ox = ctx.frame.atoms.x[i];
    const StorageReal oy = ctx.frame.atoms.y[i];
    const StorageReal oz = ctx.frame.atoms.z[i];
    double nx = static_cast<double>(ox) + disp(*rng);
    double ny = static_cast<double>(oy) + disp(*rng);
    double nz = static_cast<double>(oz) + disp(*rng);

    // Wrap the displaced atom back into the primary cell so coordinates stay
    // bounded over long runs — the COM otherwise random-walks off to infinity
    // (a displacement perpendicular to a pair bond is ΔE=0 and always accepted),
    // which eventually drives an atom outside the box and blows up the cell
    // list. Energy is unaffected: the evaluator uses the minimum image. A
    // degenerate (zero-length) box is left as-is (wrap_axis is a no-op).
    ctx.frame.box.wrap_into(nx, ny, nz);

    const StorageReal wx = static_cast<StorageReal>(nx);
    const StorageReal wy = static_cast<StorageReal>(ny);
    const StorageReal wz = static_cast<StorageReal>(nz);
    ctx.frame.atoms.x[i] = wx;
    ctx.frame.atoms.y[i] = wy;
    ctx.frame.atoms.z[i] = wz;

    txn->clear();
    txn->moved.push_back(static_cast<LocalInt>(i));
    txn->old_x.push_back(ox);
    txn->old_y.push_back(oy);
    txn->old_z.push_back(oz);
    txn->new_x.push_back(wx);
    txn->new_y.push_back(wy);
    txn->new_z.push_back(wz);
}

} // namespace mc
} // namespace cpu
} // namespace ATV
