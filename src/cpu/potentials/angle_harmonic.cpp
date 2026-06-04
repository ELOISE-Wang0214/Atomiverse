/// @file cpu/potentials/angle_harmonic.cpp
/// @brief CPU AngleHarmonic — V = K(θ − θ₀)² implementation.

#include "atomiverse/cpu/potentials/angle_harmonic.h"
#include "atomiverse/cpu/md/resolve_view.h"

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
    int n_threads = 1;
#ifdef _OPENMP
    n_threads = omp_get_max_threads();
#endif
    scratch_.build(n_threads, ctx.frame.atoms.nmax);
}

void AngleHarmonic::compute(const LaunchContext& ctx) {
    auto& f = ctx.frame;
    auto& a = f.atoms;

    const md::ResolveView rv(ctx);
    const int nlocal = rv.n_local();
    const int n_force = rv.n_force();

    // Quick exit if no angles anywhere.
    int max_na = 0;
    for (int i = 0; i < nlocal; ++i)
        max_na = std::max(max_na, a.num_angle[i]);
    if (max_na == 0)
        return;

    double energy_total = 0;
    // Angle virial: Σ d_ji ⊗ F_i + d_jk ⊗ F_k (centre force cancels by
    // Newton's 3rd, so this is translation-invariant and PBC-safe).
    md::Voigt6 vir;

    const double* __restrict x_ptr = a.x.data();
    const double* __restrict y_ptr = a.y.data();
    const double* __restrict z_ptr = a.z.data();

#pragma omp declare reduction(vadd : md::Voigt6 : omp_out += omp_in)            \
    initializer(omp_priv = md::Voigt6{})

#pragma omp parallel reduction(+ : energy_total) reduction(vadd : vir)
    {
        int tid = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num();
#endif
        const md::ThreadForce tf = scratch_.begin_thread(tid, n_force);

        // LAMMPS-style per-atom iteration: each owned atom `j` (centre) owns
        // the angles where it sits in the middle. The two non-centre partners
        // (i, k) are stored at angle_atom(j, 2*m + 0/1).
#pragma omp for schedule(static)
        for (int j = 0; j < nlocal; ++j) {
            const int na = a.num_angle[j];
            for (int m = 0; m < na; ++m) {
                const int i = rv.resolve(a.angle_atom(j, 2 * m + 0));
                const int k = rv.resolve(a.angle_atom(j, 2 * m + 1));
                if (i < 0 || k < 0)
                    continue;
                const int type = a.angle_type(j, m);

                const double ri[3] = {x_ptr[i], y_ptr[i], z_ptr[i]};
                const double rj[3] = {x_ptr[j], y_ptr[j], z_ptr[j]};
                const double rk[3] = {x_ptr[k], y_ptr[k], z_ptr[k]};

                double dx_ji = ri[0] - rj[0];
                double dy_ji = ri[1] - rj[1];
                double dz_ji = ri[2] - rj[2];
                f.box.min_image(dx_ji, dy_ji, dz_ji);

                double dx_jk = rk[0] - rj[0];
                double dy_jk = rk[1] - rj[1];
                double dz_jk = rk[2] - rj[2];
                f.box.min_image(dx_jk, dy_jk, dz_jk);

                const double r_ji =
                    std::sqrt(dx_ji * dx_ji + dy_ji * dy_ji + dz_ji * dz_ji);
                const double r_jk =
                    std::sqrt(dx_jk * dx_jk + dy_jk * dy_jk + dz_jk * dz_jk);

                double cos_theta =
                    (dx_ji * dx_jk + dy_ji * dy_jk + dz_ji * dz_jk)
                    / (r_ji * r_jk);
                cos_theta = std::max(double(-1), std::min(double(1), cos_theta));
                const double theta = std::acos(cos_theta);

                const double dtheta = theta - theta0_[type];
                energy_total += K_[type] * dtheta * dtheta;

                const double dV = double(2) * K_[type] * dtheta;

                const double sin_theta =
                    std::sqrt(double(1) - cos_theta * cos_theta);
                if (sin_theta < double(1e-10))
                    continue;

                const double a_coeff = dV / (r_ji * sin_theta);
                const double c_coeff = dV / (r_jk * sin_theta);

                const double fi_x =
                    a_coeff * (dx_jk / r_jk - cos_theta * dx_ji / r_ji);
                const double fi_y =
                    a_coeff * (dy_jk / r_jk - cos_theta * dy_ji / r_ji);
                const double fi_z =
                    a_coeff * (dz_jk / r_jk - cos_theta * dz_ji / r_ji);

                const double fk_x =
                    c_coeff * (dx_ji / r_ji - cos_theta * dx_jk / r_jk);
                const double fk_y =
                    c_coeff * (dy_ji / r_ji - cos_theta * dy_jk / r_jk);
                const double fk_z =
                    c_coeff * (dz_ji / r_ji - cos_theta * dz_jk / r_jk);

                // Force on j: -(fi + fk) (Newton's 3rd law).
                tf.fx[i] += fi_x;
                tf.fy[i] += fi_y;
                tf.fz[i] += fi_z;
                tf.fx[k] += fk_x;
                tf.fy[k] += fk_y;
                tf.fz[k] += fk_z;
                tf.fx[j] -= fi_x + fk_x;
                tf.fy[j] -= fi_y + fk_y;
                tf.fz[j] -= fi_z + fk_z;

                vir.add_pair(dx_ji, dy_ji, dz_ji, fi_x, fi_y, fi_z);
                vir.add_pair(dx_jk, dy_jk, dz_jk, fk_x, fk_y, fk_z);
            }
        }
    }

    scratch_.reduce_into(f, n_force);

    f.scalars.energy += energy_total;
    f.scalars.energy_angle += energy_total;
    vir.flush_to(f.scalars);
}

void AngleHarmonic::destroy(const DestroyContext& /*ctx*/) {
    scratch_.clear();
    // Owned-tuple cache is owned by the cpu::Angle base.
}

} // namespace cpu
} // namespace ATV
