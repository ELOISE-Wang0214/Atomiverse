/// @file cpu/potentials/bond_harmonic.cpp
/// @brief CPU BondHarmonic — V = K(r − r₀)² implementation.

#include "atomiverse/cpu/potentials/bond_harmonic.h"
#include "atomiverse/core/mat3.h"
#include "atomiverse/cpu/md/mpi/domain.h"

#include <algorithm>
#include <cstdint>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace ATV {
namespace cpu {

BondHarmonic::BondHarmonic(std::vector<double> K, std::vector<double> r0)
    : K_(std::move(K)), r0_(std::move(r0)) {}

void BondHarmonic::build(const BuildContext& ctx) {
    n_threads_ = 1;
#ifdef _OPENMP
    n_threads_ = omp_get_max_threads();
#endif
    // Size per-thread accumulators to the full backing-store capacity
    // (`nmax`), not the owner count (`nlocal`). Force kernels accumulate
    // into ghost rows too (the partner row may be a ghost), so the
    // accumulator must accommodate indices up to `nlocal + nghost`,
    // which is bounded by `nmax`.
    int n = ctx.frame.atoms.nmax;
    fx_thr_.resize(n_threads_);
    fy_thr_.resize(n_threads_);
    fz_thr_.resize(n_threads_);
    for (int t = 0; t < n_threads_; ++t) {
        fx_thr_[t].resize(n, double(0));
        fy_thr_[t].resize(n, double(0));
        fz_thr_[t].resize(n, double(0));
    }
}

void BondHarmonic::compute(const LaunchContext& ctx) {
    auto& f = ctx.frame;
    auto& a = f.atoms;

    // Resolve partner tag → local row index, and pick the right owner-row
    // count: multi-rank uses `dd->n_local` (the authoritative DD-managed
    // count); single-rank uses `atoms.nlocal`. `atoms.nlocal` is left at
    // the Frame-ctor default in multi-rank because syncing it broke the
    // MPI_Alltoallv buffer math in atom_migration (the pack stride was
    // tied to the old `frame.atoms.size()` value). The `dd->tag_to_local`
    // map already covers owners + ghosts post-halo.
    const auto* dd = ctx.runtime.dd;
    const bool multi_rank = (dd != nullptr) && (ctx.runtime.mpi_size > 1);
    const int nlocal = multi_rank ? dd->n_local : a.nlocal;
    auto resolve = [&](std::int64_t tag) -> int {
        if (multi_rank) {
            auto it = dd->tag_to_local.find(tag);
            return (it == dd->tag_to_local.end()) ? -1 : it->second;
        }
        auto it = a.map.find(tag);
        return (it == a.map.end()) ? -1 : it->second;
    };

    // Bond virial accumulators (Voigt xx yy zz xy xz yz). Stiff bonds
    // contribute significantly to the stress tensor; without this the
    // NPT barostat over-expands the box. Matches LAMMPS
    // `bond_harmonic::compute`.
    double v_xx = 0, v_yy = 0, v_zz = 0, v_xy = 0, v_xz = 0, v_yz = 0;

    // Quick exit if no bonds anywhere on this rank.
    int max_nb = 0;
    for (int i = 0; i < nlocal; ++i)
        max_nb = std::max(max_nb, a.num_bond[i]);
    if (max_nb == 0)
        return;

    const int n_atoms = a.n_total();
    double energy_total = 0;

#pragma omp parallel reduction(+ : energy_total) reduction(+ : v_xx)           \
    reduction(+ : v_yy) reduction(+ : v_zz) reduction(+ : v_xy)                \
    reduction(+ : v_xz) reduction(+ : v_yz)
    {
        int tid = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num();
#endif
        auto& fxt = fx_thr_[tid];
        auto& fyt = fy_thr_[tid];
        auto& fzt = fz_thr_[tid];
        std::fill(fxt.begin(), fxt.begin() + n_atoms, double(0));
        std::fill(fyt.begin(), fyt.begin() + n_atoms, double(0));
        std::fill(fzt.begin(), fzt.begin() + n_atoms, double(0));

        const double* __restrict x_ptr = a.x.data();
        const double* __restrict y_ptr = a.y.data();
        const double* __restrict z_ptr = a.z.data();

        // LAMMPS-style per-atom iteration: each owned atom `i` iterates
        // its own bond list (which by Newton-bond-ON convention contains
        // only bonds where `i` is the smaller-tag owner). Partner atom
        // `j` is resolved via `atoms.map` — it may be owned or a ghost.
        // Ghost partners get their force half via the reverse-comm step.
#pragma omp for schedule(static)
        for (int i = 0; i < nlocal; ++i) {
            const int nb = a.num_bond[i];
            for (int m = 0; m < nb; ++m) {
                const std::int64_t partner_tag = a.bond_atom(i, m);
                const int j = resolve(partner_tag);
                if (j < 0)
                    continue; // partner not on this rank — dropped silently
                              // (should not happen if halo cutoff is correct)
                const int type = a.bond_type(i, m);

                const double ri[3] = {x_ptr[i], y_ptr[i], z_ptr[i]};
                const double rj[3] = {x_ptr[j], y_ptr[j], z_ptr[j]};
                double dx = ri[0] - rj[0];
                double dy = ri[1] - rj[1];
                double dz = ri[2] - rj[2];
                f.box.min_image(dx, dy, dz);

                double r = std::sqrt(dx * dx + dy * dy + dz * dz);
                double dr = r - r0_[type];

                energy_total += K_[type] * dr * dr;

                double fpair = -double(2) * K_[type] * dr / r;
                double fx_val = fpair * dx;
                double fy_val = fpair * dy;
                double fz_val = fpair * dz;

                fxt[i] += fx_val;
                fxt[j] -= fx_val;
                fyt[i] += fy_val;
                fyt[j] -= fy_val;
                fzt[i] += fz_val;
                fzt[j] -= fz_val;

                v_xx += dx * fx_val;
                v_yy += dy * fy_val;
                v_zz += dz * fz_val;
                v_xy += dx * fy_val;
                v_xz += dx * fz_val;
                v_yz += dy * fz_val;
            }
        }

#pragma omp for schedule(static)
        for (int i = 0; i < n_atoms; ++i) {
            for (int t = 0; t < n_threads_; ++t) {
                f.atoms.fx[i] += fx_thr_[t][i];
                f.atoms.fy[i] += fy_thr_[t][i];
                f.atoms.fz[i] += fz_thr_[t][i];
            }
        }
    }

    f.scalars.energy += energy_total;
    f.scalars.energy_bond += energy_total;
    f.scalars.stress[0] += v_xx;
    f.scalars.stress[1] += v_yy;
    f.scalars.stress[2] += v_zz;
    f.scalars.stress[3] += v_xy;
    f.scalars.stress[4] += v_xz;
    f.scalars.stress[5] += v_yz;
}

void BondHarmonic::destroy(const DestroyContext& /*ctx*/) {
    fx_thr_.clear();
    fy_thr_.clear();
    fz_thr_.clear();
    // Owned-tuple cache is owned by the cpu::Bond base; nothing to clear here.
}

} // namespace cpu
} // namespace ATV
