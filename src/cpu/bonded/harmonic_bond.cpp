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

static inline void minimum_image(double& dx, double& dy, double& dz, const Mat3& box) {
    double lx = box(0, 0), ly = box(1, 1), lz = box(2, 2);
    if (dx >  lx * double(0.5)) dx -= lx;
    if (dx < -lx * double(0.5)) dx += lx;
    if (dy >  ly * double(0.5)) dy -= ly;
    if (dy < -ly * double(0.5)) dy += ly;
    if (dz >  lz * double(0.5)) dz -= lz;
    if (dz < -lz * double(0.5)) dz += lz;
}

HarmonicBond::HarmonicBond(std::vector<double> K, std::vector<double> r0)
    : K_(std::move(K)), r0_(std::move(r0)) {}

void HarmonicBond::build(const cpu::BuildContext& ctx) {
    n_threads_ = 1;
#ifdef _OPENMP
    n_threads_ = omp_get_max_threads();
#endif
    int n = ctx.frame.atoms.size();
    fx_thr_.resize(n_threads_);
    fy_thr_.resize(n_threads_);
    fz_thr_.resize(n_threads_);
    for (int t = 0; t < n_threads_; ++t) {
        fx_thr_[t].resize(n, double(0));
        fy_thr_[t].resize(n, double(0));
        fz_thr_[t].resize(n, double(0));
    }
}

void HarmonicBond::launch(const cpu::LaunchContext& ctx) {
    auto& f = ctx.frame;
    int n_bonds = static_cast<int>(f.bonds.atomi.size());
    if (n_bonds == 0) return;

    int n_atoms = f.atoms.size();
    double energy_total = 0;

    #pragma omp parallel reduction(+:energy_total)
    {
        int tid = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num();
#endif
        auto& fxt = fx_thr_[tid];
        auto& fyt = fy_thr_[tid];
        auto& fzt = fz_thr_[tid];
        std::fill(fxt.begin(), fxt.end(), double(0));
        std::fill(fyt.begin(), fyt.end(), double(0));
        std::fill(fzt.begin(), fzt.end(), double(0));

        #pragma omp for schedule(static)
        for (int b = 0; b < n_bonds; ++b) {
            int i = f.bonds.atomi[b];
            int j = f.bonds.atomj[b];
            int type = f.bonds.type[b];

            double dx = f.atoms.x[i] - f.atoms.x[j];
            double dy = f.atoms.y[i] - f.atoms.y[j];
            double dz = f.atoms.z[i] - f.atoms.z[j];
            minimum_image(dx, dy, dz, f.box);

            double r = std::sqrt(dx * dx + dy * dy + dz * dz);
            double dr = r - r0_[type];

            // V = K(r - r0)^2, F = -2K(r - r0) * r_hat
            energy_total += K_[type] * dr * dr;

            double fpair = -double(2) * K_[type] * dr / r;
            double fx_val = fpair * dx;
            double fy_val = fpair * dy;
            double fz_val = fpair * dz;

            fxt[i] += fx_val;  fxt[j] -= fx_val;
            fyt[i] += fy_val;  fyt[j] -= fy_val;
            fzt[i] += fz_val;  fzt[j] -= fz_val;
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
