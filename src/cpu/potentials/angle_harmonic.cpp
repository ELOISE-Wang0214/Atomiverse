/// @file cpu/potentials/angle_harmonic.cpp
/// @brief CPU AngleHarmonic — V = K(θ − θ₀)² implementation.

#include "atomiverse/cpu/potentials/angle_harmonic.h"
#include "atomiverse/core/mat3.h"
#include "atomiverse/cpu/md/mpi/domain.h"

#include <algorithm>
#include <cmath>
#include <cstdint>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace ATV {
namespace cpu {

AngleHarmonic::AngleHarmonic(std::vector<double> K, std::vector<double> theta0)
    : K_(std::move(K)), theta0_(std::move(theta0)) {}

void AngleHarmonic::build(const BuildContext& ctx) {
    n_threads_ = 1;
#ifdef _OPENMP
    n_threads_ = omp_get_max_threads();
#endif
    int n = ctx.frame.atoms.nmax;
    fx_thr_.resize(n_threads_);
    fy_thr_.resize(n_threads_);
    fz_thr_.resize(n_threads_);
    for (int t = 0; t < n_threads_; ++t) {
        fx_thr_[t].resize(n, double(0));
        fy_thr_[t].resize(n, double(0));
        fz_thr_[t].resize(n, double(0));
    }
}

void AngleHarmonic::compute(const LaunchContext& ctx) {
    auto& f = ctx.frame;
    auto& a = f.atoms;
    const auto* dd = ctx.runtime.dd;
    const bool multi_rank = (dd != nullptr) && (ctx.runtime.mpi_size > 1);
    const int nlocal = multi_rank ? dd->n_local : a.nlocal;

    // Quick exit if no angles anywhere.
    int max_na = 0;
    for (int i = 0; i < nlocal; ++i)
        max_na = std::max(max_na, a.num_angle[i]);
    if (max_na == 0)
        return;

    auto resolve = [&](std::int64_t tag) -> int {
        if (multi_rank) {
            auto it = dd->tag_to_local.find(tag);
            return (it == dd->tag_to_local.end()) ? -1 : it->second;
        }
        auto it = a.map.find(tag);
        return (it == a.map.end()) ? -1 : it->second;
    };

    const int n_atoms = a.n_total();
    double energy_total = 0;
    // Angle virial: Σ_atoms r_α · F_β = d_ji_α · F_i_β + d_jk_α · F_k_β
    // (centre-atom force F_j = -(F_i + F_k) drops out by Newton's 3rd).
    double v_xx = 0, v_yy = 0, v_zz = 0, v_xy = 0, v_xz = 0, v_yz = 0;

#pragma omp parallel reduction(+ : energy_total) reduction(+ : v_xx)           \
    reduction(+ : v_yy) reduction(+ : v_zz) reduction(+ : v_xy)                \
    reduction(+ : v_xz) reduction(+ : v_yz)
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

        const double* __restrict x_ptr = a.x.data();
        const double* __restrict y_ptr = a.y.data();
        const double* __restrict z_ptr = a.z.data();

        // LAMMPS-style per-atom iteration: each owned atom `j` (center)
        // owns the angles where it sits in the middle. The two non-center
        // partners (i, k) are stored at angle_atom(j, 2*m + 0/1).
#pragma omp for schedule(static)
        for (int j_owner = 0; j_owner < nlocal; ++j_owner) {
            const int na = a.num_angle[j_owner];
            const int j = j_owner;
            for (int m = 0; m < na; ++m) {
                const std::int64_t ti = a.angle_atom(j, 2 * m + 0);
                const std::int64_t tk = a.angle_atom(j, 2 * m + 1);
                const int i = resolve(ti);
                const int k = resolve(tk);
                if (i < 0 || k < 0)
                    continue;
                const int type = a.angle_type(j, m);

            // Pack the angle tuple's positions into stack-resident triples
            // (B.3): three contiguous loads per atom feed the math below
            // from registers, removing the per-component xtensor strided
            // indirection. Bit-identical: the FP operations and order are
            // unchanged, only the address computation moves outward.
            const double ri[3] = {x_ptr[i], y_ptr[i], z_ptr[i]};
            const double rj[3] = {x_ptr[j], y_ptr[j], z_ptr[j]};
            const double rk[3] = {x_ptr[k], y_ptr[k], z_ptr[k]};

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

            double r_ji =
                std::sqrt(dx_ji * dx_ji + dy_ji * dy_ji + dz_ji * dz_ji);
            double r_jk =
                std::sqrt(dx_jk * dx_jk + dy_jk * dy_jk + dz_jk * dz_jk);

            double cos_theta =
                (dx_ji * dx_jk + dy_ji * dy_jk + dz_ji * dz_jk) / (r_ji * r_jk);
            cos_theta = std::max(double(-1), std::min(double(1), cos_theta));
            double theta = std::acos(cos_theta);

            double dtheta = theta - theta0_[type];
            energy_total += K_[type] * dtheta * dtheta;

            // dV/dtheta = 2 * K * dtheta
            double dV = double(2) * K_[type] * dtheta;

            // sin(theta) for force projection
            double sin_theta = std::sqrt(double(1) - cos_theta * cos_theta);
            if (sin_theta < double(1e-10))
                continue;

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
            fxt[i] += fi_x;
            fyt[i] += fi_y;
            fzt[i] += fi_z;
            fxt[k] += fk_x;
            fyt[k] += fk_y;
            fzt[k] += fk_z;
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
            }  // for m (per-atom angle entry)
        }      // for j_owner (per-atom)

#pragma omp for schedule(static)
        for (int i = 0; i < n_atoms; ++i) {
            for (int t = 0; t < n_threads_; ++t) {
                f.atoms.fx[i] += fx_thr_[t][i];
                f.atoms.fy[i] += fy_thr_[t][i];
                f.atoms.fz[i] += fz_thr_[t][i];
            }
        }
    }

    f.scalars.energy += energy_total;
    f.scalars.energy_angle += energy_total;
    f.scalars.stress[0] += v_xx;
    f.scalars.stress[1] += v_yy;
    f.scalars.stress[2] += v_zz;
    f.scalars.stress[3] += v_xy;
    f.scalars.stress[4] += v_xz;
    f.scalars.stress[5] += v_yz;
}

void AngleHarmonic::destroy(const DestroyContext& /*ctx*/) {
    fx_thr_.clear();
    fy_thr_.clear();
    fz_thr_.clear();
    // Owned-tuple cache is owned by the cpu::Angle base.
}

} // namespace cpu
} // namespace ATV
