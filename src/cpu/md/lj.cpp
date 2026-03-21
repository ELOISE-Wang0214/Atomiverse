/// @file cpu/md/lj.cpp
/// @brief CPU LJ 12-6 (shifted) with OpenMP per-thread force reduction.

#include "atomiverse/cpu/md/lj.h"
#include "atomiverse/cpu/md/pair_list.h"
#include "atomiverse/core/mat3.h"

#include <algorithm>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace ATV {
namespace cpu {
namespace md {

static inline void minimum_image(Real& dx, Real& dy, Real& dz, const Mat3& box) {
    Real lx = box(0, 0), ly = box(1, 1), lz = box(2, 2);
    if (dx >  lx * Real(0.5)) dx -= lx;
    if (dx < -lx * Real(0.5)) dx += lx;
    if (dy >  ly * Real(0.5)) dy -= ly;
    if (dy < -ly * Real(0.5)) dy += ly;
    if (dz >  lz * Real(0.5)) dz -= lz;
    if (dz < -lz * Real(0.5)) dz += lz;
}

LJ::LJ(Real epsilon, Real sigma, Real cutoff)
    : epsilon_(epsilon), sigma_(sigma), cutoff_(cutoff) {
    sigma6_ = sigma * sigma * sigma * sigma * sigma * sigma;
    cutoff_sq_ = cutoff * cutoff;

    // Compute V_shift = V(r_cutoff)
    Real r2inv = Real(1) / cutoff_sq_;
    Real r6inv = r2inv * r2inv * r2inv * sigma6_;
    v_shift_ = Real(4) * epsilon_ * (r6inv * r6inv - r6inv);
}

void LJ::build(const cpu::BuildContext& ctx) {
    n_threads_ = 1;
#ifdef _OPENMP
    n_threads_ = omp_get_max_threads();
#endif
    Int n = ctx.frame.n_atoms();
    fx_thr_.resize(n_threads_);
    fy_thr_.resize(n_threads_);
    fz_thr_.resize(n_threads_);
    for (int t = 0; t < n_threads_; ++t) {
        fx_thr_[t].resize(n, Real(0));
        fy_thr_[t].resize(n, Real(0));
        fz_thr_[t].resize(n, Real(0));
    }
}

void LJ::launch(const cpu::LaunchContext& ctx) {
    auto& f = ctx.frame;
    auto& pairs = *static_cast<const PairList*>(f.pair_list);
    const auto& box = f.box;
    Int n = f.n_atoms();
    Real energy_total = 0;

    #pragma omp parallel reduction(+:energy_total)
    {
        int tid = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num();
#endif
        auto& fxt = fx_thr_[tid];
        auto& fyt = fy_thr_[tid];
        auto& fzt = fz_thr_[tid];
        std::fill(fxt.begin(), fxt.end(), Real(0));
        std::fill(fyt.begin(), fyt.end(), Real(0));
        std::fill(fzt.begin(), fzt.end(), Real(0));

        #pragma omp for schedule(static)
        for (Int ii = 0; ii < pairs.inum; ++ii) {
            Int i = pairs.ilist[ii];
            Int jnum = pairs.numneigh[i];
            Int off = pairs.firstneigh[i];
            for (Int jj = 0; jj < jnum; ++jj) {
                Int j = pairs.neighbors[off + jj];
                Real dx = f.x[i] - f.x[j];
                Real dy = f.y[i] - f.y[j];
                Real dz = f.z[i] - f.z[j];
                minimum_image(dx, dy, dz, box);
                Real r2 = dx * dx + dy * dy + dz * dz;
                if (r2 >= cutoff_sq_) continue;

                Real r2inv = Real(1) / r2;
                Real r6inv = r2inv * r2inv * r2inv * sigma6_;
                Real fpair = Real(24) * epsilon_ * r2inv * (Real(2) * r6inv * r6inv - r6inv);

                fxt[i] += fpair * dx;  fxt[j] -= fpair * dx;
                fyt[i] += fpair * dy;  fyt[j] -= fpair * dy;
                fzt[i] += fpair * dz;  fzt[j] -= fpair * dz;

                energy_total += Real(4) * epsilon_ * (r6inv * r6inv - r6inv) - v_shift_;
            }
        }

        // Reduce per-thread forces into frame
        #pragma omp for schedule(static)
        for (Int i = 0; i < n; ++i) {
            for (int t = 0; t < n_threads_; ++t) {
                f.fx[i] += fx_thr_[t][i];
                f.fy[i] += fy_thr_[t][i];
                f.fz[i] += fz_thr_[t][i];
            }
        }
    }

    f.energy += energy_total;
}

void LJ::destroy(const cpu::DestroyContext& ctx) {
    fx_thr_.clear();
    fy_thr_.clear();
    fz_thr_.clear();
}

} // namespace md
} // namespace cpu
} // namespace ATV
