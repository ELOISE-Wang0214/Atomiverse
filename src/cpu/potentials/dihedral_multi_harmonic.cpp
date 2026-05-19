/// @file cpu/potentials/dihedral_multi_harmonic.cpp
/// @brief CPU DihedralMultiHarmonic — V = Σₙ Cₙ cosⁿ(φ) implementation.

#include "atomiverse/cpu/potentials/dihedral_multi_harmonic.h"
#include "atomiverse/core/mat3.h"
#include "atomiverse/cpu/md/mpi/domain.h"

#include <algorithm>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace ATV {
namespace cpu {

// Cross product
static inline void cross(double ax, double ay, double az,
                          double bx, double by, double bz,
                          double& cx, double& cy, double& cz) {
    cx = ay * bz - az * by;
    cy = az * bx - ax * bz;
    cz = ax * by - ay * bx;
}

DihedralMultiHarmonic::DihedralMultiHarmonic(std::vector<std::vector<double>> Cn)
    : Cn_(std::move(Cn)) {}

void DihedralMultiHarmonic::build(const BuildContext& ctx) {
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

void DihedralMultiHarmonic::compute(const LaunchContext& ctx) {
    auto& f = ctx.frame;
    const int n_dihedrals = static_cast<int>(f.dihedrals.atomi.size());
    if (n_dihedrals == 0) return;

    const auto* dd = ctx.runtime.dd;
    const bool multi_rank = (dd != nullptr) && (ctx.runtime.mpi_size > 1);
    int n_atoms = multi_rank ? (dd->n_local + dd->n_ghost) : f.atoms.size();

    // Owner-tag rule for dihedral (i,j,k,l): 2nd atom (j). LAMMPS convention.
    const auto& owned = select_owned_tuples(
        ctx, n_dihedrals,
        [&](int d) -> std::int64_t {
            return static_cast<std::int64_t>(f.dihedrals.atomj[d]);
        });
    const int n_work = static_cast<int>(owned.size());
    if (n_work == 0) return;

    double energy_total = 0;
    // Dihedral virial: rewriting Σ_atoms r_α·F_β with r_i2 as reference
    // (and Newton's 3rd F_i2 = -(F_i1+F_i3+F_i4)) gives
    //   virial = -b1·F1 + b2·(F3+F4) + b3·F4
    // which is translation/PBC-invariant.
    double v_xx = 0, v_yy = 0, v_zz = 0, v_xy = 0, v_xz = 0, v_yz = 0;

    #pragma omp parallel reduction(+:energy_total) \
        reduction(+:v_xx) reduction(+:v_yy) reduction(+:v_zz) \
        reduction(+:v_xy) reduction(+:v_xz) reduction(+:v_yz)
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

        // Pointer hoist: drop xtensor's strided indexing on every access.
        const double* __restrict x_ptr = f.atoms.x.data();
        const double* __restrict y_ptr = f.atoms.y.data();
        const double* __restrict z_ptr = f.atoms.z.data();

        #pragma omp for schedule(static)
        for (int id = 0; id < n_work; ++id) {
            const int d = owned[static_cast<std::size_t>(id)];
            int i1 = f.dihedrals.atomi[d];
            int i2 = f.dihedrals.atomj[d];
            int i3 = f.dihedrals.atomk[d];
            int i4 = f.dihedrals.atoml[d];
            int type = f.dihedrals.type[d];

            // Multi-rank: dihedrals.atom{i,j,k,l} store global tags. Owner
            // = rank holding the 2nd atom (j). Skip if not owned, or if
            // any of the four atoms is beyond this rank's halo skin.
            if (multi_rank) {
                auto it2 = dd->tag_to_local.find(static_cast<std::int64_t>(i2));
                if (it2 == dd->tag_to_local.end()) continue;
                auto it1 = dd->tag_to_local.find(static_cast<std::int64_t>(i1));
                auto it3 = dd->tag_to_local.find(static_cast<std::int64_t>(i3));
                auto it4 = dd->tag_to_local.find(static_cast<std::int64_t>(i4));
                if (it1 == dd->tag_to_local.end() ||
                    it3 == dd->tag_to_local.end() ||
                    it4 == dd->tag_to_local.end()) continue;
                i1 = it1->second;
                i2 = it2->second;
                i3 = it3->second;
                i4 = it4->second;
            }

            // Pack the dihedral tuple's positions into stack-resident triples
            // (B.3): four contiguous loads per atom let the compiler keep all
            // 12 components in registers across the cross-product math below,
            // removing the strided indirection through xtensor's `operator[]`.
            // Bit-identical to the inlined form (same loads, same FP order).
            const double r1[3] = { x_ptr[i1], y_ptr[i1], z_ptr[i1] };
            const double r2[3] = { x_ptr[i2], y_ptr[i2], z_ptr[i2] };
            const double r3[3] = { x_ptr[i3], y_ptr[i3], z_ptr[i3] };
            const double r4[3] = { x_ptr[i4], y_ptr[i4], z_ptr[i4] };

            // Vectors: b1 = r2-r1, b2 = r3-r2, b3 = r4-r3
            double b1x = r2[0] - r1[0], b1y = r2[1] - r1[1], b1z = r2[2] - r1[2];
            double b2x = r3[0] - r2[0], b2y = r3[1] - r2[1], b2z = r3[2] - r2[2];
            double b3x = r4[0] - r3[0], b3y = r4[1] - r3[1], b3z = r4[2] - r3[2];
            f.box.min_image(b1x, b1y, b1z);
            f.box.min_image(b2x, b2y, b2z);
            f.box.min_image(b3x, b3y, b3z);

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

            // V = sum Cn cos^n(phi); dV/dcos = sum n·Cn·cos^(n-1).
            // Use a running cos_pow_minus1 to avoid dividing by cos(phi)
            // (which vanishes at phi = ±90°).
            const auto& coeffs = Cn_[type];
            double energy = 0;
            double dV_dcos = 0;
            double cos_pow = 1;          // cos^n
            double cos_pow_minus1 = 0;   // cos^(n-1); 0 for n=0
            for (size_t n = 0; n < coeffs.size(); ++n) {
                energy += coeffs[n] * cos_pow;
                if (n > 0)
                    dV_dcos += static_cast<double>(n) * coeffs[n] * cos_pow_minus1;
                cos_pow_minus1 = cos_pow;
                cos_pow *= cos_phi;
            }
            energy_total += energy;

            // Force via direct chain rule: F_i = -dV/dcos · ∂cos_phi/∂r_i.
            //
            // Let f_hat_1 = ∂(cos_phi)/∂n1, f_hat_2 = ∂(cos_phi)/∂n2 as 3-vectors:
            //   f_hat_1 = (n2/|n2| - cos_phi·n1/|n1|) / |n1|
            //   f_hat_2 = (n1/|n1| - cos_phi·n2/|n2|) / |n2|
            // Then
            //   ∂cos_phi/∂b1 = b2 × f_hat_1
            //   ∂cos_phi/∂b3 = f_hat_2 × b2
            //   ∂cos_phi/∂b2 = (f_hat_1 × b1) + (b3 × f_hat_2)
            // and finally
            //   F_1 =  dV/dcos · ∂cos/∂b1
            //   F_4 = -dV/dcos · ∂cos/∂b3
            //   F_2 = -F_1 + dV/dcos · ∂cos/∂b2
            //   F_3 = -dV/dcos · ∂cos/∂b2 - F_4
            // (checked against a central-difference gradient -- see
            //  tests/cpu/md/test_pe_fd_force.cpp).
            const double n1h_x = n1x * n1_inv, n1h_y = n1y * n1_inv, n1h_z = n1z * n1_inv;
            const double n2h_x = n2x * n2_inv, n2h_y = n2y * n2_inv, n2h_z = n2z * n2_inv;

            const double fh1_x = (n2h_x - cos_phi * n1h_x) * n1_inv;
            const double fh1_y = (n2h_y - cos_phi * n1h_y) * n1_inv;
            const double fh1_z = (n2h_z - cos_phi * n1h_z) * n1_inv;

            const double fh2_x = (n1h_x - cos_phi * n2h_x) * n2_inv;
            const double fh2_y = (n1h_y - cos_phi * n2h_y) * n2_inv;
            const double fh2_z = (n1h_z - cos_phi * n2h_z) * n2_inv;

            // ∂cos/∂b1 = b2 × fh1
            const double db1_x = b2y * fh1_z - b2z * fh1_y;
            const double db1_y = b2z * fh1_x - b2x * fh1_z;
            const double db1_z = b2x * fh1_y - b2y * fh1_x;

            // ∂cos/∂b3 = fh2 × b2
            const double db3_x = fh2_y * b2z - fh2_z * b2y;
            const double db3_y = fh2_z * b2x - fh2_x * b2z;
            const double db3_z = fh2_x * b2y - fh2_y * b2x;

            // ∂cos/∂b2 = (fh1 × b1) + (b3 × fh2)
            const double db2_x = (fh1_y * b1z - fh1_z * b1y) + (b3y * fh2_z - b3z * fh2_y);
            const double db2_y = (fh1_z * b1x - fh1_x * b1z) + (b3z * fh2_x - b3x * fh2_z);
            const double db2_z = (fh1_x * b1y - fh1_y * b1x) + (b3x * fh2_y - b3y * fh2_x);

            const double f1x =  dV_dcos * db1_x;
            const double f1y =  dV_dcos * db1_y;
            const double f1z =  dV_dcos * db1_z;

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

            fxt[i1] += f1x; fyt[i1] += f1y; fzt[i1] += f1z;
            fxt[i2] += f2x; fyt[i2] += f2y; fzt[i2] += f2z;
            fxt[i3] += f3x; fyt[i3] += f3y; fzt[i3] += f3z;
            fxt[i4] += f4x; fyt[i4] += f4y; fzt[i4] += f4z;

            // Virial per dihedral (see comment above loop).
            v_xx += -b1x*f1x + b2x*(f3x+f4x) + b3x*f4x;
            v_yy += -b1y*f1y + b2y*(f3y+f4y) + b3y*f4y;
            v_zz += -b1z*f1z + b2z*(f3z+f4z) + b3z*f4z;
            v_xy += -b1x*f1y + b2x*(f3y+f4y) + b3x*f4y;
            v_xz += -b1x*f1z + b2x*(f3z+f4z) + b3x*f4z;
            v_yz += -b1y*f1z + b2y*(f3z+f4z) + b3y*f4z;
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
    f.stress[0] += v_xx;
    f.stress[1] += v_yy;
    f.stress[2] += v_zz;
    f.stress[3] += v_xy;
    f.stress[4] += v_xz;
    f.stress[5] += v_yz;
}

void DihedralMultiHarmonic::destroy(const DestroyContext& /*ctx*/) {
    fx_thr_.clear();
    fy_thr_.clear();
    fz_thr_.clear();
    // Owned-tuple cache is owned by the cpu::Dihedral base.
}

} // namespace cpu
} // namespace ATV
