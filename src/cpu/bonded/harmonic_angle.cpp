/// @file cpu/bonded/harmonic_angle.cpp
/// @brief CPU harmonic angle bending implementation.

#include "atomiverse/cpu/bonded/harmonic_angle.h"
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

HarmonicAngle::HarmonicAngle(std::vector<Real> K, std::vector<Real> theta0)
    : K_(std::move(K)), theta0_(std::move(theta0)) {}

void HarmonicAngle::build(const cpu::BuildContext& ctx) {
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

void HarmonicAngle::launch(const cpu::LaunchContext& ctx) {
    auto& f = ctx.frame;
    Int n_angles = static_cast<Int>(f.angle_i.size());
    if (n_angles == 0) return;

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
        for (Int a = 0; a < n_angles; ++a) {
            Int i = f.angle_i[a];
            Int j = f.angle_j[a]; // center atom
            Int k = f.angle_k[a];
            Int type = f.angle_type[a];

            // r_ji = r_i - r_j
            Real dx_ji = f.x[i] - f.x[j];
            Real dy_ji = f.y[i] - f.y[j];
            Real dz_ji = f.z[i] - f.z[j];
            minimum_image(dx_ji, dy_ji, dz_ji, f.box);

            // r_jk = r_k - r_j
            Real dx_jk = f.x[k] - f.x[j];
            Real dy_jk = f.y[k] - f.y[j];
            Real dz_jk = f.z[k] - f.z[j];
            minimum_image(dx_jk, dy_jk, dz_jk, f.box);

            Real r_ji = std::sqrt(dx_ji * dx_ji + dy_ji * dy_ji + dz_ji * dz_ji);
            Real r_jk = std::sqrt(dx_jk * dx_jk + dy_jk * dy_jk + dz_jk * dz_jk);

            Real cos_theta = (dx_ji * dx_jk + dy_ji * dy_jk + dz_ji * dz_jk) / (r_ji * r_jk);
            cos_theta = std::max(Real(-1), std::min(Real(1), cos_theta));
            Real theta = std::acos(cos_theta);

            Real dtheta = theta - theta0_[type];
            energy_total += K_[type] * dtheta * dtheta;

            // dV/dtheta = 2 * K * dtheta
            Real dV = Real(2) * K_[type] * dtheta;

            // sin(theta) for force projection
            Real sin_theta = std::sqrt(Real(1) - cos_theta * cos_theta);
            if (sin_theta < Real(1e-10)) continue;

            Real a_coeff = -dV / (r_ji * sin_theta);
            Real c_coeff = -dV / (r_jk * sin_theta);

            // Force on i: a * (r_jk/r_jk - cos*r_ji/r_ji)
            Real fi_x = a_coeff * (dx_jk / r_jk - cos_theta * dx_ji / r_ji);
            Real fi_y = a_coeff * (dy_jk / r_jk - cos_theta * dy_ji / r_ji);
            Real fi_z = a_coeff * (dz_jk / r_jk - cos_theta * dz_ji / r_ji);

            // Force on k: c * (r_ji/r_ji - cos*r_jk/r_jk)
            Real fk_x = c_coeff * (dx_ji / r_ji - cos_theta * dx_jk / r_jk);
            Real fk_y = c_coeff * (dy_ji / r_ji - cos_theta * dy_jk / r_jk);
            Real fk_z = c_coeff * (dz_ji / r_ji - cos_theta * dz_jk / r_jk);

            // Force on j: -(fi + fk) (Newton's 3rd law)
            fxt[i] += fi_x;  fyt[i] += fi_y;  fzt[i] += fi_z;
            fxt[k] += fk_x;  fyt[k] += fk_y;  fzt[k] += fk_z;
            fxt[j] -= fi_x + fk_x;
            fyt[j] -= fi_y + fk_y;
            fzt[j] -= fi_z + fk_z;
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

void HarmonicAngle::destroy(const cpu::DestroyContext& ctx) {
    fx_thr_.clear();
    fy_thr_.clear();
    fz_thr_.clear();
}

} // namespace bonded
} // namespace cpu
} // namespace ATV
