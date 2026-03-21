/// @file cpu/bonded/harmonic_bond.cpp
/// @brief CPU harmonic bond stretching implementation.

#include "atomiverse/cpu/bonded/harmonic_bond.h"
#include "atomiverse/core/mat3.h"

#include <algorithm>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace ATV {
namespace cpu {
namespace bonded {

static inline void minimum_image(Real& dx, Real& dy, Real& dz, const Mat3& box) {
    Real lx = box(0, 0), ly = box(1, 1), lz = box(2, 2);
    if (dx >  lx * Real(0.5)) dx -= lx;
    if (dx < -lx * Real(0.5)) dx += lx;
    if (dy >  ly * Real(0.5)) dy -= ly;
    if (dy < -ly * Real(0.5)) dy += ly;
    if (dz >  lz * Real(0.5)) dz -= lz;
    if (dz < -lz * Real(0.5)) dz += lz;
}

HarmonicBond::HarmonicBond(std::vector<Real> K, std::vector<Real> r0)
    : K_(std::move(K)), r0_(std::move(r0)) {}

void HarmonicBond::build(const cpu::BuildContext& ctx) {
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

void HarmonicBond::launch(const cpu::LaunchContext& ctx) {
    auto& f = ctx.frame;
    Int n_bonds = static_cast<Int>(f.bond_i.size());
    if (n_bonds == 0) return;

    Int n_atoms = f.n_atoms();
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
        for (Int b = 0; b < n_bonds; ++b) {
            Int i = f.bond_i[b];
            Int j = f.bond_j[b];
            Int type = f.bond_type[b];

            Real dx = f.x[i] - f.x[j];
            Real dy = f.y[i] - f.y[j];
            Real dz = f.z[i] - f.z[j];
            minimum_image(dx, dy, dz, f.box);

            Real r = std::sqrt(dx * dx + dy * dy + dz * dz);
            Real dr = r - r0_[type];

            // V = K(r - r0)^2, F = -2K(r - r0) * r_hat
            energy_total += K_[type] * dr * dr;

            Real fpair = -Real(2) * K_[type] * dr / r;
            Real fx_val = fpair * dx;
            Real fy_val = fpair * dy;
            Real fz_val = fpair * dz;

            fxt[i] += fx_val;  fxt[j] -= fx_val;
            fyt[i] += fy_val;  fyt[j] -= fy_val;
            fzt[i] += fz_val;  fzt[j] -= fz_val;
        }

        #pragma omp for schedule(static)
        for (Int i = 0; i < n_atoms; ++i) {
            for (int t = 0; t < n_threads_; ++t) {
                f.fx[i] += fx_thr_[t][i];
                f.fy[i] += fy_thr_[t][i];
                f.fz[i] += fz_thr_[t][i];
            }
        }
    }

    f.energy += energy_total;
}

void HarmonicBond::destroy(const cpu::DestroyContext& ctx) {
    fx_thr_.clear();
    fy_thr_.clear();
    fz_thr_.clear();
}

} // namespace bonded
} // namespace cpu
} // namespace ATV
