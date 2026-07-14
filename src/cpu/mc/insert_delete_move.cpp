/// @file src/cpu/mc/insert_delete_move.cpp
/// @brief cpu::mc::InsertDeleteMove implementation.

#include "atomiverse/cpu/mc/insert_delete_move.h"

#include <random>
#include <stdexcept>

#include "atomiverse/core/precision.h"
#include "atomiverse/cpu/mc/transaction.h"

namespace ATV {
namespace cpu {
namespace mc {

InsertDeleteMove::InsertDeleteMove(int type, double mass_amu)
    : type_(type), mass_amu_(mass_amu) {
    if (type <= 0)
        throw std::invalid_argument("InsertDeleteMove: type must be > 0");
    if (!(mass_amu > 0.0))
        throw std::invalid_argument("InsertDeleteMove: mass must be > 0 amu");
}

void InsertDeleteMove::launch(const cpu::LaunchContext& ctx) {
    if (ctx.step < 0)
        return; // energy-only bootstrap pass

    std::mt19937_64* rng = ctx.runtime.mc.rng;
    mc::Transaction* txn = ctx.runtime.mc.proposal;
    if (rng == nullptr || txn == nullptr) {
        throw std::runtime_error(
            "InsertDeleteMove::launch: ctx.mc.rng / ctx.mc.proposal not "
            "published — InsertDeleteMove must run under an MC preset");
    }

    auto& a = ctx.frame.atoms;
    std::uniform_real_distribution<double> u01(0.0, 1.0);

    if (u01(*rng) < 0.5) {
        // ── Insert at a uniform fractional position ─────────────────────
        txn->begin(Transaction::Verb::Insert, ctx.frame);
        if (a.nlocal + 1 > a.nmax)
            a.grow(a.nlocal + 1);
        const int slot = a.nlocal;
        const double u0 = u01(*rng), u1 = u01(*rng), u2 = u01(*rng);
        const auto& c = ctx.frame.box.cell();
        a.x[slot] = static_cast<StorageReal>(u0 * c(0, 0) + u1 * c(1, 0)
                                             + u2 * c(2, 0));
        a.y[slot] = static_cast<StorageReal>(u0 * c(0, 1) + u1 * c(1, 1)
                                             + u2 * c(2, 1));
        a.z[slot] = static_cast<StorageReal>(u0 * c(0, 2) + u1 * c(1, 2)
                                             + u2 * c(2, 2));
        // Zero everything else — this slot may hold a stale deleted row.
        a.vx[slot] = 0;
        a.vy[slot] = 0;
        a.vz[slot] = 0;
        a.fx[slot] = 0;
        a.fy[slot] = 0;
        a.fz[slot] = 0;
        a.charge[slot] = 0;
        a.image(slot, 0) = 0;
        a.image(slot, 1) = 0;
        a.image(slot, 2) = 0;
        a.num_bond[slot] = 0;
        a.num_angle[slot] = 0;
        a.num_dihedral[slot] = 0;
        a.nspecial(slot, 0) = 0;
        a.nspecial(slot, 1) = 0;
        a.nspecial(slot, 2) = 0;
        a.mass[slot] = static_cast<StorageReal>(mass_amu_);
        a.type[slot] = type_;
        a.tag[slot] = a.acquire_tag();
        ++a.nlocal;
    } else {
        // ── Delete a uniformly chosen atom ──────────────────────────────
        if (a.nlocal == 0) {
            // "A trial deletion from N = 0 is simply rejected" — empty
            // proposal; the policy rolls it back without a Boltzmann draw.
            txn->begin(Transaction::Verb::None, ctx.frame);
            return;
        }
        std::uniform_int_distribution<int> pick(0, a.nlocal - 1);
        const int victim = pick(*rng);
        txn->begin(Transaction::Verb::Delete, ctx.frame);
        txn->record_touch(ctx.frame, victim); // full row, topology included
        a.copy_row(victim, a.nlocal - 1);
        --a.nlocal;
    }
}

} // namespace mc
} // namespace cpu
} // namespace ATV
