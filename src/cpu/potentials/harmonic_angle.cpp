/// @file cpu/potentials/harmonic_angle.cpp
/// @brief CPU harmonic angle bending implementation.

#include "atomiverse/cpu/potentials/harmonic_angle.h"
#include "atomiverse/core/mat3.h"
#include "atomiverse/cpu/md/mpi/domain.h"

#include <algorithm>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace ATV {
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
    const int n_angles = static_cast<int>(f.angles.atomi.size());
    if (n_angles == 0) return;

    const auto* dd = ctx.runtime.dd;
    const bool multi_rank = (dd != nullptr) && (ctx.runtime.mpi_size > 1);
    int n_atoms = f.atoms.size();
    if (multi_rank) {
        n_atoms = dd->n_local + dd->n_ghost;
        if (local_angles_generation_ != dd->neighbor_generation) {
            local_angles_.clear();
            local_angles_.reserve(static_cast<std::size_t>(n_angles));
            for (int a = 0; a < n_angles; ++a) {
                const std::int64_t owner_tag =
                    static_cast<std::int64_t>(f.angles.atomj[a]);
                auto it_owner = dd->tag_to_local.find(owner_tag);
                if (it_owner != dd->tag_to_local.end() &&
                    it_owner->second < dd->n_local) {
                    local_angles_.push_back(a);
                }
            }
            local_angles_generation_ = dd->neighbor_generation;
        }
    }
    const int n_work = multi_rank ? static_cast<int>(local_angles_.size())
                                  : n_angles;
    if (n_work == 0) return;

    double energy_total = 0;
    // Angle virial: Σ_atoms r_α · F_β = d_ji_α · F_i_β + d_jk_α · F_k_β
    // (centre-atom force F_j = -(F_i + F_k) drops out by Newton's 3rd).
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
        for (int ia = 0; ia < n_work; ++ia) {
            const int a = multi_rank ? local_angles_[static_cast<std::size_t>(ia)]
                                     : ia;
            int i = f.angles.atomi[a];
            int j = f.angles.atomj[a]; // center atom
            int k = f.angles.atomk[a];
            int type = f.angles.type[a];

            // Multi-rank: angles.atom{i,j,k} store global tags. Owner =
            // rank holding centre atom (j). Skip if not owned, or if any
            // of the three atoms is beyond this rank's halo skin.
            if (multi_rank) {
                auto it_j = dd->tag_to_local.find(static_cast<std::int64_t>(j));
                if (it_j == dd->tag_to_local.end()) continue;
                auto it_i = dd->tag_to_local.find(static_cast<std::int64_t>(i));
                auto it_k = dd->tag_to_local.find(static_cast<std::int64_t>(k));
                if (it_i == dd->tag_to_local.end() ||
                    it_k == dd->tag_to_local.end()) continue;
                i = it_i->second;
                j = it_j->second;
                k = it_k->second;
            }

            // Pack the angle tuple's positions into stack-resident triples
            // (B.3): three contiguous loads per atom feed the math below
            // from registers, removing the per-component xtensor strided
            // indirection. Bit-identical: the FP operations and order are
            // unchanged, only the address computation moves outward.
            const double ri[3] = { x_ptr[i], y_ptr[i], z_ptr[i] };
            const double rj[3] = { x_ptr[j], y_ptr[j], z_ptr[j] };
            const double rk[3] = { x_ptr[k], y_ptr[k], z_ptr[k] };

            // r_ji = r_i - r_j
            double dx_ji = ri[0] - rj[0];
            double dy_ji = ri[1] - rj[1];
            double dz_ji = ri[2] - rj[2];
            f.box.min_image(dx_ji, dy_ji, dz_ji);

            // r_jk = r_k - r_j
            double dx_jk = rk[0] - rj[0];
            double dy_jk = rk[1] - rj[1];
            double dz_jk = rk[2] - rj[2];
            f.box.min_image(dx_jk, dy_jk, dz_jk);

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

            // Virial contribution (d_ji · F_i + d_jk · F_k). Equivalent to
            // Σ_atoms r · F once the centre-atom contribution cancels, so
            // this is translation-invariant and safe under PBC.
            v_xx += dx_ji * fi_x + dx_jk * fk_x;
            v_yy += dy_ji * fi_y + dy_jk * fk_y;
            v_zz += dz_ji * fi_z + dz_jk * fk_z;
            v_xy += dx_ji * fi_y + dx_jk * fk_y;
            v_xz += dx_ji * fi_z + dx_jk * fk_z;
            v_yz += dy_ji * fi_z + dy_jk * fk_z;
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

void HarmonicAngle::destroy(const cpu::DestroyContext& ctx) {
    fx_thr_.clear();
    fy_thr_.clear();
    fz_thr_.clear();
    local_angles_.clear();
    local_angles_generation_ = ~std::uint64_t{0};
}

} // namespace ATV
