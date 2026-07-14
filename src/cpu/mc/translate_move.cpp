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
    // Snapshot the row into the patch BEFORE mutating it.
    txn->begin(Transaction::Verb::Displace, ctx.frame);
    txn->record_touch(ctx.frame, i);

    double nx = static_cast<double>(ctx.frame.atoms.x[i]) + disp(*rng);
    double ny = static_cast<double>(ctx.frame.atoms.y[i]) + disp(*rng);
    double nz = static_cast<double>(ctx.frame.atoms.z[i]) + disp(*rng);

    // Wrap the displaced atom back into the primary cell so coordinates stay
    // bounded over long runs — the COM otherwise random-walks off to infinity
    // (a displacement perpendicular to a pair bond is ΔE=0 and always accepted),
    // which eventually drives an atom outside the box and blows up the cell
    // list. Energy is unaffected: the evaluator uses the minimum image. A
    // degenerate (zero-length) box is left as-is (wrap_axis is a no-op).
    ctx.frame.box.wrap_into(nx, ny, nz);

    ctx.frame.atoms.x[i] = static_cast<StorageReal>(nx);
    ctx.frame.atoms.y[i] = static_cast<StorageReal>(ny);
    ctx.frame.atoms.z[i] = static_cast<StorageReal>(nz);
}

} // namespace mc
} // namespace cpu
} // namespace ATV
