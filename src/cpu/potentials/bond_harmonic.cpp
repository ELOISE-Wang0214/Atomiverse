/// @file cpu/potentials/bond_harmonic.cpp
/// @brief CPU BondHarmonic — V = K(r − r₀)² implementation.

#include "atomiverse/cpu/potentials/bond_harmonic.h"
#include "atomiverse/cpu/md/resolve_view.h"

#include <algorithm>
#include <cmath>
#include <cstdint>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace ATV {
namespace cpu {

BondHarmonic::BondHarmonic(std::vector<double> K, std::vector<double> r0)
    : K_(std::move(K)), r0_(std::move(r0)) {}

void BondHarmonic::build(const BuildContext& ctx) {
    int n_threads = 1;
#ifdef _OPENMP
    n_threads = omp_get_max_threads();
#endif
    // Per-thread accumulators are sized to the full backing-store capacity
    // (`nmax`); only the live `[0, n_force)` rows are touched per step (a
    // partner row may be a ghost, so the range must span owners + ghosts).
    scratch_.build(n_threads, ctx.frame.atoms.nmax);
}

void BondHarmonic::compute(const LaunchContext& ctx) {
    auto& f = ctx.frame;
    auto& a = f.atoms;

    // Single source of truth for owned/ghost counts + bounded tag→local
    // resolution. `n_local` = owners to iterate; `n_force` = owners + ghosts,
    // the extent the per-thread force buffers and the reduction touch (NOT
    // `atoms.n_total()`, which is the global atom count in multi-rank).
    // `resolve()` returns -1 for any tag that is absent or maps outside
    // `[0, n_force)` — a stale ghost partner is dropped, not read as garbage.
    const md::ResolveView rv(ctx);
    const int nlocal = rv.n_local();
    const int n_force = rv.n_force();

    // Quick exit if no bonds anywhere on this rank.
    int max_nb = 0;
    for (int i = 0; i < nlocal; ++i)
        max_nb = std::max(max_nb, a.num_bond[i]);
    if (max_nb == 0)
        return;

    double energy_total = 0;
    md::Voigt6 vir;
    std::uint64_t dropped = 0, evaluated = 0;

    const double* __restrict x_ptr = a.x.data();
    const double* __restrict y_ptr = a.y.data();
    const double* __restrict z_ptr = a.z.data();

// Stiff bonds contribute significantly to the stress tensor. The per-pair
// virial accumulates into a thread-private `vir_local` and is combined under
// a critical section (once per thread per step) — an OMP user-defined
// reduction on Voigt6 would heap-allocate per step, violating the
// "no allocation in launch()" invariant.
#pragma omp parallel reduction(+ : energy_total, dropped, evaluated)
    {
        int tid = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num();
#endif
        const md::ThreadForce tf = scratch_.begin_thread(tid, n_force);
        md::Voigt6 vir_local;

        // LAMMPS-style per-atom iteration: each owned atom `i` iterates its
        // own bond list (Newton-bond-ON: only bonds where `i` is the
        // smaller-tag owner). Partner `j` may be an owner or a ghost; ghost
        // force halves ship home via HaloExchangeReverse.
#pragma omp for schedule(static)
        for (int i = 0; i < nlocal; ++i) {
            const int nb = a.num_bond[i];
            for (int m = 0; m < nb; ++m) {
                const int j = rv.resolve(a.bond_atom(i, m));
                if (j < 0) {
                    ++dropped; // partner not resolvable on this rank
                    continue;
                }
                ++evaluated;
                const int type = a.bond_type(i, m);

                double dx = x_ptr[i] - x_ptr[j];
                double dy = y_ptr[i] - y_ptr[j];
                double dz = z_ptr[i] - z_ptr[j];
                f.box.min_image(dx, dy, dz);

                const double r = std::sqrt(dx * dx + dy * dy + dz * dz);
                const double dr = r - r0_[type];

                energy_total += K_[type] * dr * dr;

                const double fpair = -double(2) * K_[type] * dr / r;
                const double fx_val = fpair * dx;
                const double fy_val = fpair * dy;
                const double fz_val = fpair * dz;

                tf.fx[i] += fx_val;
                tf.fx[j] -= fx_val;
                tf.fy[i] += fy_val;
                tf.fy[j] -= fy_val;
                tf.fz[i] += fz_val;
                tf.fz[j] -= fz_val;

                vir_local.add_pair(dx, dy, dz, fx_val, fy_val, fz_val);
            }
        }
#pragma omp critical
        vir += vir_local;
    }

    scratch_.reduce_into(f, n_force);

    f.scalars.energy += energy_total;
    f.scalars.energy_bond += energy_total;
    vir.flush_to(f.scalars);

    drops_.dropped += dropped;
    drops_.evaluated += evaluated;
}

void BondHarmonic::destroy(const DestroyContext& /*ctx*/) {
    drops_.report("BondHarmonic");
    drops_.reset();
    scratch_.clear();
    // Owned-tuple cache is owned by the cpu::Bond base; nothing to clear here.
}

} // namespace cpu
} // namespace ATV
