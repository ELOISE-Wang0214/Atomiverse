/// @file src/cpu/potentials/coul_long.cpp
/// @brief Real-space screened Coulomb implementation (CPU).
#include "atomiverse/cpu/potentials/coul_long.h"

#include <cmath>
#include <stdexcept>

#include "atomiverse/core/units.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/component/owned_range.h"
#include "atomiverse/cpu/neighbor/pair_list.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef ATV_HAS_MPI
#include "atomiverse/cpu/comm/domain.h" // DomainDecomposition (dd->n_local/n_ghost)
#endif

namespace ATV {
namespace cpu {

// 2/sqrt(pi), the erfc-derivative prefactor (LAMMPS EWALD_F).
static constexpr double kEwaldF = 1.1283791670955126;

CoulLong::CoulLong(double alpha, double cutoff)
    : alpha_(alpha), cutoff_(cutoff), cutoff_sq_(cutoff * cutoff) {
    if (alpha <= 0.0)
        throw std::invalid_argument("CoulLong: alpha must be > 0");
    if (cutoff <= 0.0)
        throw std::invalid_argument("CoulLong: cutoff must be > 0");
}

void CoulLong::build(const BuildContext& ctx) {
    n_threads_ = 1;
#ifdef _OPENMP
    n_threads_ = omp_get_max_threads();
#endif
    const int n = ctx.frame.atoms.nmax;
    fx_thr_.assign(n_threads_, std::vector<double>(n, 0.0));
    fy_thr_.assign(n_threads_, std::vector<double>(n, 0.0));
    fz_thr_.assign(n_threads_, std::vector<double>(n, 0.0));
}

void CoulLong::launch(const LaunchContext& ctx) {
    auto& f = ctx.frame;
    if (!ctx.runtime.pair_list) {
        throw std::runtime_error(
            "CoulLong::launch: runtime.pair_list is null. Add a NeighborList "
            "to the force phase ahead of CoulLong.");
    }
    const auto& pairs = *ctx.runtime.pair_list;
    const auto& box = f.box;

    const auto* dd = ctx.runtime.dd;
    const bool mpi = dd && ctx.runtime.mpi_size > 1;
    const int n_force = mpi ? dd->n_local + dd->n_ghost : f.atoms.size();

    // Ghost charges arrive via the build_halo border message (x,y,z,tag,q), so
    // f.atoms.charge[j] is valid for ghost rows too — no global tag→charge map.

    // Grow per-thread force arenas to cover the ghost region if needed.
    if (static_cast<int>(fx_thr_.empty() ? 0 : fx_thr_[0].size()) < n_force) {
        for (int t = 0; t < n_threads_; ++t) {
            fx_thr_[t].resize(n_force, 0.0);
            fy_thr_[t].resize(n_force, 0.0);
            fz_thr_[t].resize(n_force, 0.0);
        }
    }

    const int* __restrict ilist = pairs.ilist.data();
    const int* __restrict numneigh = pairs.numneigh.data();
    const int* __restrict firstneigh = pairs.firstneigh.data();
    const int* __restrict neighbors = pairs.neighbors.data();
    const bool has_shifts = pairs.has_shifts;
    const double* __restrict shx = has_shifts ? pairs.shift_x.data() : nullptr;
    const double* __restrict shy = has_shifts ? pairs.shift_y.data() : nullptr;
    const double* __restrict shz = has_shifts ? pairs.shift_z.data() : nullptr;

    const double* __restrict x_ptr = f.atoms.x.data();
    const double* __restrict y_ptr = f.atoms.y.data();
    const double* __restrict z_ptr = f.atoms.z.data();
    const double* __restrict q_ptr = f.atoms.charge.data();

    const double alpha = alpha_;
    const double rc_sq = cutoff_sq_;
    const double qqr2e = ATV::kCoulomb;

    double energy_total = 0;
    double v_xx = 0, v_yy = 0, v_zz = 0, v_xy = 0, v_xz = 0, v_yz = 0;

#pragma omp parallel reduction(+ : energy_total) reduction(+ : v_xx)           \
    reduction(+ : v_yy) reduction(+ : v_zz) reduction(+ : v_xy)                \
    reduction(+ : v_xz) reduction(+ : v_yz)
    {
        int tid = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num();
#endif
        double* __restrict fxt = fx_thr_[tid].data();
        double* __restrict fyt = fy_thr_[tid].data();
        double* __restrict fzt = fz_thr_[tid].data();
        std::fill(fx_thr_[tid].begin(), fx_thr_[tid].begin() + n_force, 0.0);
        std::fill(fy_thr_[tid].begin(), fy_thr_[tid].begin() + n_force, 0.0);
        std::fill(fz_thr_[tid].begin(), fz_thr_[tid].begin() + n_force, 0.0);

        const double lx = box(0, 0), ly = box(1, 1), lz = box(2, 2);
        const double inv_lx = 1.0 / lx, inv_ly = 1.0 / ly, inv_lz = 1.0 / lz;

#pragma omp for schedule(static)
        for (int ii = 0; ii < pairs.inum; ++ii) {
            const int i = ilist[ii];
            const double xi = x_ptr[i], yi = y_ptr[i], zi = z_ptr[i];
            const double qi = q_ptr[i];
            if (qi == 0.0)
                continue;
            const int jnum = numneigh[i];
            const int off = firstneigh[i];
            double fi_x = 0, fi_y = 0, fi_z = 0;
            for (int jj = 0; jj < jnum; ++jj) {
                const int j = neighbors[off + jj];
                double dx, dy, dz;
                if (has_shifts) {
                    dx = xi - (x_ptr[j] + shx[off + jj]);
                    dy = yi - (y_ptr[j] + shy[off + jj]);
                    dz = zi - (z_ptr[j] + shz[off + jj]);
                } else {
                    dx = xi - x_ptr[j];
                    dy = yi - y_ptr[j];
                    dz = zi - z_ptr[j];
                    dx -= lx * std::nearbyint(dx * inv_lx);
                    dy -= ly * std::nearbyint(dy * inv_ly);
                    dz -= lz * std::nearbyint(dz * inv_lz);
                }
                const double r2 = dx * dx + dy * dy + dz * dz;
                if (r2 >= rc_sq || r2 == 0.0)
                    continue;
                // Ghost charge is carried in f.atoms.charge by the border halo,
                // so a single read covers both owned and ghost j.
                const double qj = q_ptr[j];
                if (qj == 0.0)
                    continue;

                const double r = std::sqrt(r2);
                const double r2inv = 1.0 / r2;
                const double grij = alpha * r;
                const double expm2 = std::exp(-grij * grij);
                const double erfcv = std::erfc(grij);
                const double prefactor = qqr2e * qi * qj / r;
                const double forcecoul =
                    prefactor * (erfcv + kEwaldF * grij * expm2);
                const double fpair = forcecoul * r2inv;

                const double fxp = dx * fpair;
                const double fyp = dy * fpair;
                const double fzp = dz * fpair;
                fi_x += fxp;
                fi_y += fyp;
                fi_z += fzp;
                fxt[j] -= fxp;
                fyt[j] -= fyp;
                fzt[j] -= fzp;

                v_xx += dx * fxp;
                v_yy += dy * fyp;
                v_zz += dz * fzp;
                v_xy += dx * fyp;
                v_xz += dx * fzp;
                v_yz += dy * fzp;

                energy_total += prefactor * erfcv;
            }
            fxt[i] += fi_x;
            fyt[i] += fi_y;
            fzt[i] += fi_z;
        }

#pragma omp for schedule(static)
        for (int i = 0; i < n_force; ++i) {
            for (int t = 0; t < n_threads_; ++t) {
                f.atoms.fx[i] += fx_thr_[t][i];
                f.atoms.fy[i] += fy_thr_[t][i];
                f.atoms.fz[i] += fz_thr_[t][i];
            }
        }
    }

    f.scalars.energy += energy_total;
    f.scalars.energy_pair += energy_total;
    f.scalars.stress[0] += v_xx;
    f.scalars.stress[1] += v_yy;
    f.scalars.stress[2] += v_zz;
    f.scalars.stress[3] += v_xy;
    f.scalars.stress[4] += v_xz;
    f.scalars.stress[5] += v_yz;
}

void CoulLong::destroy(const DestroyContext&) {
    fx_thr_.clear();
    fy_thr_.clear();
    fz_thr_.clear();
}

} // namespace cpu
} // namespace ATV
