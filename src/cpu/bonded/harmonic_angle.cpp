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

static inline void minimum_image(double& dx, double& dy, double& dz, const Mat3& box) {
    double lx = box(0, 0), ly = box(1, 1), lz = box(2, 2);
    if (dx >  lx * double(0.5)) dx -= lx;
    if (dx < -lx * double(0.5)) dx += lx;
    if (dy >  ly * double(0.5)) dy -= ly;
    if (dy < -ly * double(0.5)) dy += ly;
    if (dz >  lz * double(0.5)) dz -= lz;
    if (dz < -lz * double(0.5)) dz += lz;
}

HarmonicAngle::HarmonicAngle(std::vector<double> K, std::vector<double> theta0)
    : K_(std::move(K)), theta0_(std::move(theta0)) {}

void HarmonicAngle::build(const cpu::BuildContext& ctx) {
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

void HarmonicAngle::launch(const cpu::LaunchContext& ctx) {
    auto& f = ctx.frame;
    int n_angles = static_cast<int>(f.angles.atomi.size());
    if (n_angles == 0) return;

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
        for (int a = 0; a < n_angles; ++a) {
            int i = f.angles.atomi[a];
            int j = f.angles.atomj[a]; // center atom
            int k = f.angles.atomk[a];
            int type = f.angles.type[a];

            // r_ji = r_i - r_j
            double dx_ji = f.atoms.x[i] - f.atoms.x[j];
            double dy_ji = f.atoms.y[i] - f.atoms.y[j];
            double dz_ji = f.atoms.z[i] - f.atoms.z[j];
            minimum_image(dx_ji, dy_ji, dz_ji, f.box);

            // r_jk = r_k - r_j
            double dx_jk = f.atoms.x[k] - f.atoms.x[j];
            double dy_jk = f.atoms.y[k] - f.atoms.y[j];
            double dz_jk = f.atoms.z[k] - f.atoms.z[j];
            minimum_image(dx_jk, dy_jk, dz_jk, f.box);

            double r_ji = std::sqrt(dx_ji * dx_ji + dy_ji * dy_ji + dz_ji * dz_ji);
            double r_jk = std::sqrt(dx_jk * dx_jk + dy_jk * dy_jk + dz_jk * dz_jk);

            double cos_theta = (dx_ji * dx_jk + dy_ji * dy_jk + dz_ji * dz_jk) / (r_ji * r_jk);
            cos_theta = std::max(double(-1), std::min(double(1), cos_theta));
            double theta = std::acos(cos_theta);

            double dtheta = theta - theta0_[type];
            energy_total += K_[type] * dtheta * dtheta;

            // dV/dtheta = 2 * K * dtheta
            double dV = double(2) * K_[type] * dtheta;

            // sin(theta) for force projection
            double sin_theta = std::sqrt(double(1) - cos_theta * cos_theta);
            if (sin_theta < double(1e-10)) continue;

            // F_i = -dV/dr_i = +(dV/dθ)/(r_ji·sin θ) · (r_jk_hat - cos θ·r_ji_hat)
            // (the - in d(cos θ)/dr_i and the - in dθ/d(cos θ) = -1/sin θ
            //  combine into a + here -- the previous version dropped one sign).
            double a_coeff = dV / (r_ji * sin_theta);
            double c_coeff = dV / (r_jk * sin_theta);

            double fi_x = a_coeff * (dx_jk / r_jk - cos_theta * dx_ji / r_ji);
            double fi_y = a_coeff * (dy_jk / r_jk - cos_theta * dy_ji / r_ji);
            double fi_z = a_coeff * (dz_jk / r_jk - cos_theta * dz_ji / r_ji);

            double fk_x = c_coeff * (dx_ji / r_ji - cos_theta * dx_jk / r_jk);
            double fk_y = c_coeff * (dy_ji / r_ji - cos_theta * dy_jk / r_jk);
            double fk_z = c_coeff * (dz_ji / r_ji - cos_theta * dz_jk / r_jk);

            // Force on j: -(fi + fk) (Newton's 3rd law)
            fxt[i] += fi_x;  fyt[i] += fi_y;  fzt[i] += fi_z;
            fxt[k] += fk_x;  fyt[k] += fk_y;  fzt[k] += fk_z;
            fxt[j] -= fi_x + fk_x;
            fyt[j] -= fi_y + fk_y;
            fzt[j] -= fi_z + fk_z;
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

void HarmonicAngle::destroy(const cpu::DestroyContext& ctx) {
    fx_thr_.clear();
    fy_thr_.clear();
    fz_thr_.clear();
}

} // namespace bonded
} // namespace cpu
} // namespace ATV
