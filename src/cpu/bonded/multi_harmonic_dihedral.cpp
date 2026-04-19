/// @file cpu/bonded/multi_harmonic_dihedral.cpp
/// @brief CPU multi-harmonic dihedral implementation.

#include "atomiverse/cpu/bonded/multi_harmonic_dihedral.h"
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

// Cross product
static inline void cross(double ax, double ay, double az,
                          double bx, double by, double bz,
                          double& cx, double& cy, double& cz) {
    cx = ay * bz - az * by;
    cy = az * bx - ax * bz;
    cz = ax * by - ay * bx;
}

MultiHarmonicDihedral::MultiHarmonicDihedral(std::vector<std::vector<double>> Cn)
    : Cn_(std::move(Cn)) {}

void MultiHarmonicDihedral::build(const cpu::BuildContext& ctx) {
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

void MultiHarmonicDihedral::launch(const cpu::LaunchContext& ctx) {
    auto& f = ctx.frame;
    int n_dihedrals = static_cast<int>(f.dihedrals.atomi.size());
    if (n_dihedrals == 0) return;

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
        for (int d = 0; d < n_dihedrals; ++d) {
            int i1 = f.dihedrals.atomi[d];
            int i2 = f.dihedrals.atomj[d];
            int i3 = f.dihedrals.atomk[d];
            int i4 = f.dihedrals.atoml[d];
            int type = f.dihedrals.type[d];

            // Vectors: b1 = r2-r1, b2 = r3-r2, b3 = r4-r3
            double b1x = f.atoms.x[i2] - f.atoms.x[i1], b1y = f.atoms.y[i2] - f.atoms.y[i1], b1z = f.atoms.z[i2] - f.atoms.z[i1];
            double b2x = f.atoms.x[i3] - f.atoms.x[i2], b2y = f.atoms.y[i3] - f.atoms.y[i2], b2z = f.atoms.z[i3] - f.atoms.z[i2];
            double b3x = f.atoms.x[i4] - f.atoms.x[i3], b3y = f.atoms.y[i4] - f.atoms.y[i3], b3z = f.atoms.z[i4] - f.atoms.z[i3];
            minimum_image(b1x, b1y, b1z, f.box);
            minimum_image(b2x, b2y, b2z, f.box);
            minimum_image(b3x, b3y, b3z, f.box);

            // Normal vectors: n1 = b1 x b2, n2 = b2 x b3
            double n1x, n1y, n1z, n2x, n2y, n2z;
            cross(b1x, b1y, b1z, b2x, b2y, b2z, n1x, n1y, n1z);
            cross(b2x, b2y, b2z, b3x, b3y, b3z, n2x, n2y, n2z);

            double n1_sq = n1x * n1x + n1y * n1y + n1z * n1z;
            double n2_sq = n2x * n2x + n2y * n2y + n2z * n2z;
            if (n1_sq < double(1e-20) || n2_sq < double(1e-20)) continue;

            double n1_inv = double(1) / std::sqrt(n1_sq);
            double n2_inv = double(1) / std::sqrt(n2_sq);

            double cos_phi = (n1x * n2x + n1y * n2y + n1z * n2z) * n1_inv * n2_inv;
            cos_phi = std::max(double(-1), std::min(double(1), cos_phi));

            // V = sum(Cn * cos^n(phi))
            // dV/dcos(phi) = sum(n * Cn * cos^(n-1)) -- accumulated via a
            // running cos_pow_minus1 to avoid the /cos_phi singularity at
            // phi = 90°.
            const auto& coeffs = Cn_[type];
            double energy = 0;
            double dV_dcos = 0;
            double cos_pow = 1;          // cos^n at start of iteration
            double cos_pow_minus1 = 0;   // cos^(n-1); 0 for n=0 (term has factor n)
            for (size_t n = 0; n < coeffs.size(); ++n) {
                energy += coeffs[n] * cos_pow;
                if (n > 0)
                    dV_dcos += static_cast<double>(n) * coeffs[n] * cos_pow_minus1;
                cos_pow_minus1 = cos_pow;
                cos_pow *= cos_phi;
            }
            energy_total += energy;

            // sin(phi) via the identity (b1 · (b2 × b3)) * |b2| / (|n1||n2|).
            // b2 × b3 == n2, so b1 · n2 supplies the numerator.
            double b2_len = std::sqrt(b2x * b2x + b2y * b2y + b2z * b2z);
            double sin_phi = (b1x * n2x + b1y * n2y + b1z * n2z) * b2_len
                              * n1_inv * n2_inv;

            // dV/dphi = -sin(phi) * dV/dcos(phi)
            double dV_dphi = -sin_phi * dV_dcos;

            // Force computation via chain rule:
            //   F1 = (dV/dphi) * b2_len / n1_sq * n1
            //   F4 = -(dV/dphi) * b2_len / n2_sq * n2
            // (sign convention chosen so F1 + F2 + F3 + F4 = 0).
            double f1_scale = dV_dphi * b2_len / n1_sq;
            double f4_scale = -dV_dphi * b2_len / n2_sq;

            double f1x = f1_scale * n1x, f1y = f1_scale * n1y, f1z = f1_scale * n1z;
            double f4x = f4_scale * n2x, f4y = f4_scale * n2y, f4z = f4_scale * n2z;

            // Intermediate forces (approximate for i2, i3)
            double b1_dot_b2 = b1x * b2x + b1y * b2y + b1z * b2z;
            double b3_dot_b2 = b3x * b2x + b3y * b2y + b3z * b2z;
            double b2_sq = b2x * b2x + b2y * b2y + b2z * b2z;

            double s1 = b1_dot_b2 / b2_sq;
            double s2 = b3_dot_b2 / b2_sq;

            double f2x = -f1x + s1 * f1x - s2 * f4x;
            double f2y = -f1y + s1 * f1y - s2 * f4y;
            double f2z = -f1z + s1 * f1z - s2 * f4z;
            double f3x = -f4x - s1 * f1x + s2 * f4x;
            double f3y = -f4y - s1 * f1y + s2 * f4y;
            double f3z = -f4z - s1 * f1z + s2 * f4z;

            fxt[i1] += f1x; fyt[i1] += f1y; fzt[i1] += f1z;
            fxt[i2] += f2x; fyt[i2] += f2y; fzt[i2] += f2z;
            fxt[i3] += f3x; fyt[i3] += f3y; fzt[i3] += f3z;
            fxt[i4] += f4x; fyt[i4] += f4y; fzt[i4] += f4z;
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

void MultiHarmonicDihedral::destroy(const cpu::DestroyContext& ctx) {
    fx_thr_.clear();
    fy_thr_.clear();
    fz_thr_.clear();
}

} // namespace bonded
} // namespace cpu
} // namespace ATV
