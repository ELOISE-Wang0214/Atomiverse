/// @file cpu/potentials/dihedral_multi_harmonic.cpp
/// @brief CPU DihedralMultiHarmonic — V = Σₙ Cₙ cosⁿ(φ) implementation.

#include "atomiverse/cpu/potentials/dihedral_multi_harmonic.h"
#include "atomiverse/cpu/potentials/resolve_view.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace ATV {
namespace cpu {

// Cross product
static inline void cross(double ax, double ay, double az, double bx, double by,
                         double bz, double& cx, double& cy, double& cz) {
    cx = ay * bz - az * by;
    cy = az * bx - ax * bz;
    cz = ax * by - ay * bx;
}

DihedralMultiHarmonic::DihedralMultiHarmonic(
    std::vector<std::vector<double>> Cn)
    : Cn_(std::move(Cn)) {}

void DihedralMultiHarmonic::build(const BuildContext& ctx) {
    int n_threads = 1;
#ifdef _OPENMP
    n_threads = omp_get_max_threads();
#endif
    scratch_.build(n_threads, ctx.frame.atoms.nmax);
}

void DihedralMultiHarmonic::compute(const LaunchContext& ctx) {
    auto& f = ctx.frame;
    auto& a = f.atoms;

    const md::ResolveView rv(ctx);
    const int nlocal = rv.n_local();
    const int n_force = rv.n_force();

    // Tuple cache: tag→local hash resolution happens only when row
    // indices may have moved (migration / halo full rebuild), not every
    // step. Layout: {i1, i2, i3, i4, type} per dihedral (i2 = owner).
    if (cache_.stale(ctx)) {
        cache_.flat.clear();
        cache_.dropped_at_build = 0;
        for (int i2 = 0; i2 < nlocal; ++i2) {
            const int nd = a.num_dihedral[i2];
            for (int m = 0; m < nd; ++m) {
                const int i1 = rv.resolve(a.dihedral_atom(i2, 3 * m + 0));
                const int i3 = rv.resolve(a.dihedral_atom(i2, 3 * m + 1));
                const int i4 = rv.resolve(a.dihedral_atom(i2, 3 * m + 2));
                if (i1 < 0 || i3 < 0 || i4 < 0) {
                    ++cache_.dropped_at_build;
                    continue;
                }
                cache_.flat.push_back(i1);
                cache_.flat.push_back(i2);
                cache_.flat.push_back(i3);
                cache_.flat.push_back(i4);
                cache_.flat.push_back(a.dihedral_type(i2, m));
            }
        }
        cache_.mark_built(ctx);
        drops_.dropped += cache_.dropped_at_build;
    }
    const int n_tuples = static_cast<int>(cache_.flat.size()) / 5;
    if (n_tuples == 0)
        return;
    const int* __restrict tup = cache_.flat.data();

    double energy_total = 0;
    // Dihedral virial = (-b1)⊗F1 + b2⊗(F3+F4) + b3⊗F4 (with i2 as the
    // reference atom and Newton's 3rd F2 = -(F1+F3+F4)) — translation/PBC
    // invariant. Each term is one Voigt6::add_pair.
    md::Voigt6 vir;

    const StorageReal* __restrict x_ptr = f.span(kAtom_x);
    const StorageReal* __restrict y_ptr = f.span(kAtom_y);
    const StorageReal* __restrict z_ptr = f.span(kAtom_z);

    // Resolve the frame's force spans ONCE, before the parallel region opens
    // (never inside it — see BondedForceScratch::bind).
    scratch_.bind(f);

// Thread-private virial combined under a critical section (no per-step heap
// alloc, unlike an OMP user-defined reduction on Voigt6).
#pragma omp parallel reduction(+ : energy_total)
    {
        int tid = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num();
#endif
        const md::ThreadForce tf = scratch_.begin_thread(tid, n_force);
        md::Voigt6 vir_local;

        // Hot loop over pre-resolved tuples: atom `i2` (2nd of the
        // i-j-k-l quadruple) owns the dihedral; partners may be ghosts.
#pragma omp for schedule(static)
        for (int t = 0; t < n_tuples; ++t) {
            {
                const int i1 = tup[5 * t + 0];
                const int i2 = tup[5 * t + 1];
                const int i3 = tup[5 * t + 2];
                const int i4 = tup[5 * t + 3];
                const int type = tup[5 * t + 4];

                const double r1[3] = {x_ptr[i1], y_ptr[i1], z_ptr[i1]};
                const double r2[3] = {x_ptr[i2], y_ptr[i2], z_ptr[i2]};
                const double r3[3] = {x_ptr[i3], y_ptr[i3], z_ptr[i3]};
                const double r4[3] = {x_ptr[i4], y_ptr[i4], z_ptr[i4]};

                // Vectors: b1 = r2-r1, b2 = r3-r2, b3 = r4-r3
                double b1x = r2[0] - r1[0], b1y = r2[1] - r1[1],
                       b1z = r2[2] - r1[2];
                double b2x = r3[0] - r2[0], b2y = r3[1] - r2[1],
                       b2z = r3[2] - r2[2];
                double b3x = r4[0] - r3[0], b3y = r4[1] - r3[1],
                       b3z = r4[2] - r3[2];
                f.box.min_image(b1x, b1y, b1z);
                f.box.min_image(b2x, b2y, b2z);
                f.box.min_image(b3x, b3y, b3z);

                // Normal vectors: n1 = b1 x b2, n2 = b2 x b3
                double n1x, n1y, n1z, n2x, n2y, n2z;
                cross(b1x, b1y, b1z, b2x, b2y, b2z, n1x, n1y, n1z);
                cross(b2x, b2y, b2z, b3x, b3y, b3z, n2x, n2y, n2z);

                const double n1_sq = n1x * n1x + n1y * n1y + n1z * n1z;
                const double n2_sq = n2x * n2x + n2y * n2y + n2z * n2z;
                if (n1_sq < double(1e-20) || n2_sq < double(1e-20))
                    continue;

                const double n1_inv = double(1) / std::sqrt(n1_sq);
                const double n2_inv = double(1) / std::sqrt(n2_sq);

                double cos_phi =
                    (n1x * n2x + n1y * n2y + n1z * n2z) * n1_inv * n2_inv;
                cos_phi = std::max(double(-1), std::min(double(1), cos_phi));

                // V = Σ Cn cos^n; dV/dcos = Σ n·Cn·cos^(n-1). Running
                // cos_pow_minus1 avoids dividing by cos(phi) (vanishes at ±90°).
                const auto& coeffs = Cn_[type];
                double energy = 0;
                double dV_dcos = 0;
                double cos_pow = 1;
                double cos_pow_minus1 = 0;
                for (std::size_t n = 0; n < coeffs.size(); ++n) {
                    energy += coeffs[n] * cos_pow;
                    if (n > 0)
                        dV_dcos += static_cast<double>(n) * coeffs[n]
                                   * cos_pow_minus1;
                    cos_pow_minus1 = cos_pow;
                    cos_pow *= cos_phi;
                }
                energy_total += energy;

                // F_i = -dV/dcos · ∂cos_phi/∂r_i (see test_pe_fd_force.cpp).
                const double n1h_x = n1x * n1_inv, n1h_y = n1y * n1_inv,
                             n1h_z = n1z * n1_inv;
                const double n2h_x = n2x * n2_inv, n2h_y = n2y * n2_inv,
                             n2h_z = n2z * n2_inv;

                const double fh1_x = (n2h_x - cos_phi * n1h_x) * n1_inv;
                const double fh1_y = (n2h_y - cos_phi * n1h_y) * n1_inv;
                const double fh1_z = (n2h_z - cos_phi * n1h_z) * n1_inv;

                const double fh2_x = (n1h_x - cos_phi * n2h_x) * n2_inv;
                const double fh2_y = (n1h_y - cos_phi * n2h_y) * n2_inv;
                const double fh2_z = (n1h_z - cos_phi * n2h_z) * n2_inv;

                const double db1_x = b2y * fh1_z - b2z * fh1_y;
                const double db1_y = b2z * fh1_x - b2x * fh1_z;
                const double db1_z = b2x * fh1_y - b2y * fh1_x;

                const double db3_x = fh2_y * b2z - fh2_z * b2y;
                const double db3_y = fh2_z * b2x - fh2_x * b2z;
                const double db3_z = fh2_x * b2y - fh2_y * b2x;

                const double db2_x =
                    (fh1_y * b1z - fh1_z * b1y) + (b3y * fh2_z - b3z * fh2_y);
                const double db2_y =
                    (fh1_z * b1x - fh1_x * b1z) + (b3z * fh2_x - b3x * fh2_z);
                const double db2_z =
                    (fh1_x * b1y - fh1_y * b1x) + (b3x * fh2_y - b3y * fh2_x);

                const double f1x = dV_dcos * db1_x;
                const double f1y = dV_dcos * db1_y;
                const double f1z = dV_dcos * db1_z;

                const double f4x = -dV_dcos * db3_x;
                const double f4y = -dV_dcos * db3_y;
                const double f4z = -dV_dcos * db3_z;

                const double gx = dV_dcos * db2_x;
                const double gy = dV_dcos * db2_y;
                const double gz = dV_dcos * db2_z;

                const double f2x = -f1x + gx;
                const double f2y = -f1y + gy;
                const double f2z = -f1z + gz;

                const double f3x = -gx - f4x;
                const double f3y = -gy - f4y;
                const double f3z = -gz - f4z;

                tf.fx[i1] += f1x;
                tf.fy[i1] += f1y;
                tf.fz[i1] += f1z;
                tf.fx[i2] += f2x;
                tf.fy[i2] += f2y;
                tf.fz[i2] += f2z;
                tf.fx[i3] += f3x;
                tf.fy[i3] += f3y;
                tf.fz[i3] += f3z;
                tf.fx[i4] += f4x;
                tf.fy[i4] += f4y;
                tf.fz[i4] += f4z;

                // Virial: (-b1)⊗F1 + b2⊗(F3+F4) + b3⊗F4.
                vir_local.add_pair(-b1x, -b1y, -b1z, f1x, f1y, f1z);
                vir_local.add_pair(b2x, b2y, b2z, f3x + f4x, f3y + f4y,
                                   f3z + f4z);
                vir_local.add_pair(b3x, b3y, b3z, f4x, f4y, f4z);
            }
        }
#pragma omp critical
        vir += vir_local;
    }

    scratch_.reduce_into(n_force);

    f.scalars.energy += energy_total;
    f.scalars.energy_dihedral += energy_total;
    vir.flush_to(f.scalars);

    drops_.evaluated += static_cast<std::uint64_t>(n_tuples);
}

void DihedralMultiHarmonic::destroy(const DestroyContext& /*ctx*/) {
    drops_.report("DihedralMultiHarmonic");
    drops_.reset();
    cache_.clear();
    scratch_.clear();
    // Owned-tuple cache is owned by the cpu::Dihedral base.
}

} // namespace cpu
} // namespace ATV
