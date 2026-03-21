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

static inline void minimum_image(Real& dx, Real& dy, Real& dz, const Mat3& box) {
    Real lx = box(0, 0), ly = box(1, 1), lz = box(2, 2);
    if (dx >  lx * Real(0.5)) dx -= lx;
    if (dx < -lx * Real(0.5)) dx += lx;
    if (dy >  ly * Real(0.5)) dy -= ly;
    if (dy < -ly * Real(0.5)) dy += ly;
    if (dz >  lz * Real(0.5)) dz -= lz;
    if (dz < -lz * Real(0.5)) dz += lz;
}

// Cross product
static inline void cross(Real ax, Real ay, Real az,
                          Real bx, Real by, Real bz,
                          Real& cx, Real& cy, Real& cz) {
    cx = ay * bz - az * by;
    cy = az * bx - ax * bz;
    cz = ax * by - ay * bx;
}

MultiHarmonicDihedral::MultiHarmonicDihedral(std::vector<std::vector<Real>> Cn)
    : Cn_(std::move(Cn)) {}

void MultiHarmonicDihedral::build(const cpu::BuildContext& ctx) {
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

void MultiHarmonicDihedral::launch(const cpu::LaunchContext& ctx) {
    auto& f = ctx.frame;
    Int n_dihedrals = static_cast<Int>(f.dihedral_i.size());
    if (n_dihedrals == 0) return;

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
        for (Int d = 0; d < n_dihedrals; ++d) {
            Int i1 = f.dihedral_i[d];
            Int i2 = f.dihedral_j[d];
            Int i3 = f.dihedral_k[d];
            Int i4 = f.dihedral_l[d];
            Int type = f.dihedral_type[d];

            // Vectors: b1 = r2-r1, b2 = r3-r2, b3 = r4-r3
            Real b1x = f.x[i2] - f.x[i1], b1y = f.y[i2] - f.y[i1], b1z = f.z[i2] - f.z[i1];
            Real b2x = f.x[i3] - f.x[i2], b2y = f.y[i3] - f.y[i2], b2z = f.z[i3] - f.z[i2];
            Real b3x = f.x[i4] - f.x[i3], b3y = f.y[i4] - f.y[i3], b3z = f.z[i4] - f.z[i3];
            minimum_image(b1x, b1y, b1z, f.box);
            minimum_image(b2x, b2y, b2z, f.box);
            minimum_image(b3x, b3y, b3z, f.box);

            // Normal vectors: n1 = b1 x b2, n2 = b2 x b3
            Real n1x, n1y, n1z, n2x, n2y, n2z;
            cross(b1x, b1y, b1z, b2x, b2y, b2z, n1x, n1y, n1z);
            cross(b2x, b2y, b2z, b3x, b3y, b3z, n2x, n2y, n2z);

            Real n1_sq = n1x * n1x + n1y * n1y + n1z * n1z;
            Real n2_sq = n2x * n2x + n2y * n2y + n2z * n2z;
            if (n1_sq < Real(1e-20) || n2_sq < Real(1e-20)) continue;

            Real n1_inv = Real(1) / std::sqrt(n1_sq);
            Real n2_inv = Real(1) / std::sqrt(n2_sq);

            Real cos_phi = (n1x * n2x + n1y * n2y + n1z * n2z) * n1_inv * n2_inv;
            cos_phi = std::max(Real(-1), std::min(Real(1), cos_phi));

            // V = sum(Cn * cos^n(phi))
            const auto& coeffs = Cn_[type];
            Real energy = 0;
            Real dV_dcos = 0;
            Real cos_n = 1;  // cos^0 = 1
            for (size_t n = 0; n < coeffs.size(); ++n) {
                energy += coeffs[n] * cos_n;
                if (n > 0)
                    dV_dcos += static_cast<Real>(n) * coeffs[n] * cos_n / cos_phi;
                cos_n *= cos_phi;
            }
            energy_total += energy;

            // Force computation via chain rule (simplified)
            Real b2_len = std::sqrt(b2x * b2x + b2y * b2y + b2z * b2z);

            // F1 = -(dV/dphi) * b2_len / n1_sq * n1
            Real f1_scale = -dV_dcos * b2_len / n1_sq;
            // F4 = (dV/dphi) * b2_len / n2_sq * n2
            Real f4_scale = dV_dcos * b2_len / n2_sq;

            Real f1x = f1_scale * n1x, f1y = f1_scale * n1y, f1z = f1_scale * n1z;
            Real f4x = f4_scale * n2x, f4y = f4_scale * n2y, f4z = f4_scale * n2z;

            // Intermediate forces (approximate for i2, i3)
            Real b1_dot_b2 = b1x * b2x + b1y * b2y + b1z * b2z;
            Real b3_dot_b2 = b3x * b2x + b3y * b2y + b3z * b2z;
            Real b2_sq = b2x * b2x + b2y * b2y + b2z * b2z;

            Real s1 = b1_dot_b2 / b2_sq;
            Real s2 = b3_dot_b2 / b2_sq;

            Real f2x = -f1x + s1 * f1x - s2 * f4x;
            Real f2y = -f1y + s1 * f1y - s2 * f4y;
            Real f2z = -f1z + s1 * f1z - s2 * f4z;
            Real f3x = -f4x - s1 * f1x + s2 * f4x;
            Real f3y = -f4y - s1 * f1y + s2 * f4y;
            Real f3z = -f4z - s1 * f1z + s2 * f4z;

            fxt[i1] += f1x; fyt[i1] += f1y; fzt[i1] += f1z;
            fxt[i2] += f2x; fyt[i2] += f2y; fzt[i2] += f2z;
            fxt[i3] += f3x; fyt[i3] += f3y; fzt[i3] += f3z;
            fxt[i4] += f4x; fyt[i4] += f4y; fzt[i4] += f4z;
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

void MultiHarmonicDihedral::destroy(const cpu::DestroyContext& ctx) {
    fx_thr_.clear();
    fy_thr_.clear();
    fz_thr_.clear();
}

} // namespace bonded
} // namespace cpu
} // namespace ATV
