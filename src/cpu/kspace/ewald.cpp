/// @file src/cpu/kspace/ewald.cpp
/// @brief Reciprocal-space Ewald summation implementation (CPU).
#include "atomiverse/cpu/kspace/ewald.h"

#include <cmath>
#include <stdexcept>

#include "atomiverse/core/units.h"
#include "atomiverse/cpu/frame.h"
#include "atomiverse/cpu/md/owned_range.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef ATV_HAS_MPI
#include <mpi.h>

#include "atomiverse/cpu/md/mpi/domain.h"
#endif

namespace ATV {
namespace cpu {

KSpaceEwald::KSpaceEwald(double alpha, double accuracy)
    : alpha_(alpha), accuracy_(accuracy) {
    if (alpha <= 0.0)
        throw std::invalid_argument("KSpaceEwald: alpha must be > 0");
}

double KSpaceEwald::recommend_alpha(double cutoff, double accuracy) {
    // LAMMPS heuristic: erfc(g·rc)/rc ≈ accuracy.
    return (1.35 - 0.15 * std::log(accuracy)) / cutoff;
}

void KSpaceEwald::ensure_kvectors(const cpu::Frame& frame) {
    const auto& box = frame.box;
    std::array<double, 9> cell{box(0, 0), box(0, 1), box(0, 2),
                               box(1, 0), box(1, 1), box(1, 2),
                               box(2, 0), box(2, 1), box(2, 2)};
    if (have_kvectors_ && cell == cell_cache_)
        return;
    cell_cache_ = cell;

    // Reciprocal basis b_i = 2π · (column i of cell⁻¹); g = Σ m_i b_i.
    const auto& ic = box.inv_cell();
    const double tp = 2.0 * ATV::kPI;
    double b[3][3];
    for (int i = 0; i < 3; ++i) {
        b[i][0] = tp * ic(0, i);
        b[i][1] = tp * ic(1, i);
        b[i][2] = tp * ic(2, i);
    }
    auto blen = [&](int i) {
        return std::sqrt(b[i][0] * b[i][0] + b[i][1] * b[i][1] +
                         b[i][2] * b[i][2]);
    };
    const double gcut = 2.0 * alpha_ * std::sqrt(-std::log(accuracy_));
    const double gcut2 = gcut * gcut;
    int kmax[3];
    for (int i = 0; i < 3; ++i)
        kmax[i] = std::max(1, static_cast<int>(std::ceil(gcut / blen(i))));

    const double inv4a2 = 1.0 / (4.0 * alpha_ * alpha_);
    kx_.clear();
    ky_.clear();
    kz_.clear();
    kw_.clear();
    kvterm_.clear();
    for (int mx = -kmax[0]; mx <= kmax[0]; ++mx)
        for (int my = -kmax[1]; my <= kmax[1]; ++my)
            for (int mz = -kmax[2]; mz <= kmax[2]; ++mz) {
                if (mx == 0 && my == 0 && mz == 0)
                    continue;
                const double gx = mx * b[0][0] + my * b[1][0] + mz * b[2][0];
                const double gy = mx * b[0][1] + my * b[1][1] + mz * b[2][1];
                const double gz = mx * b[0][2] + my * b[1][2] + mz * b[2][2];
                const double g2 = gx * gx + gy * gy + gz * gz;
                if (g2 > gcut2)
                    continue;
                kx_.push_back(gx);
                ky_.push_back(gy);
                kz_.push_back(gz);
                kw_.push_back(std::exp(-g2 * inv4a2) / g2);
                kvterm_.push_back(-2.0 * (1.0 / g2 + inv4a2));
            }
    sr_.assign(kx_.size(), 0.0);
    si_.assign(kx_.size(), 0.0);
    have_kvectors_ = true;
}

void KSpaceEwald::build(const BuildContext& ctx) { ensure_kvectors(ctx.frame); }

void KSpaceEwald::launch(const LaunchContext& ctx) {
    auto& f = ctx.frame;
    ensure_kvectors(f);
    const int ng = static_cast<int>(kx_.size());

    const auto* dd = ctx.runtime.dd;
    const bool mpi = dd && ctx.runtime.mpi_size > 1;
    const int n_local = mpi ? dd->n_local : f.atoms.size();
    const bool first_rank = !mpi || ctx.runtime.mpi_rank == 0;

    const double* __restrict x = f.atoms.x.data();
    const double* __restrict y = f.atoms.y.data();
    const double* __restrict z = f.atoms.z.data();
    const double* __restrict q = f.atoms.charge.data();

    const double* __restrict kx = kx_.data();
    const double* __restrict ky = ky_.data();
    const double* __restrict kz = kz_.data();
    const double* __restrict kw = kw_.data();

    // ── Partial structure factor S(k) over owned atoms ─────────────────
#pragma omp parallel for schedule(static)
    for (int k = 0; k < ng; ++k) {
        double sr = 0, si = 0;
        const double kxk = kx[k], kyk = ky[k], kzk = kz[k];
        for (int i = 0; i < n_local; ++i) {
            const double phase = kxk * x[i] + kyk * y[i] + kzk * z[i];
            sr += q[i] * std::cos(phase);
            si += q[i] * std::sin(phase);
        }
        sr_[k] = sr;
        si_[k] = si;
    }

#ifdef ATV_HAS_MPI
    if (mpi) {
        MPI_Allreduce(MPI_IN_PLACE, sr_.data(), ng, MPI_DOUBLE, MPI_SUM,
                      dd->cart_comm);
        MPI_Allreduce(MPI_IN_PLACE, si_.data(), ng, MPI_DOUBLE, MPI_SUM,
                      dd->cart_comm);
    }
#endif

    const double vol = f.box.volume();
    const double K = ATV::kCoulomb * 2.0 * ATV::kPI / vol;

    // ── Reciprocal forces on owned atoms from global S(k) ──────────────
#pragma omp parallel for schedule(static)
    for (int i = 0; i < n_local; ++i) {
        const double xi = x[i], yi = y[i], zi = z[i], qi = q[i];
        if (qi == 0.0)
            continue;
        double fx = 0, fy = 0, fz = 0;
        for (int k = 0; k < ng; ++k) {
            const double phase = kx[k] * xi + ky[k] * yi + kz[k] * zi;
            const double c = std::cos(phase), s = std::sin(phase);
            const double fac = 2.0 * K * kw[k] * qi * (s * sr_[k] - c * si_[k]);
            fx += fac * kx[k];
            fy += fac * ky[k];
            fz += fac * kz[k];
        }
        f.atoms.fx[i] += fx;
        f.atoms.fy[i] += fy;
        f.atoms.fz[i] += fz;
    }

    // ── Self energy (partitions across ranks: owned atoms only) ────────
    double qsq_local = 0.0, q_local = 0.0;
    for (int i = 0; i < n_local; ++i) {
        qsq_local += q[i] * q[i];
        q_local += q[i];
    }
    const double e_self = -ATV::kCoulomb * alpha_ / std::sqrt(ATV::kPI) * qsq_local;
    f.scalars.energy += e_self;
    f.scalars.energy_pair += e_self;

    // ── Reciprocal + net-charge energy/stress: rank 0 only ─────────────
    if (first_rank) {
        double e_recip = 0;
        double s_xx = 0, s_yy = 0, s_zz = 0, s_xy = 0, s_xz = 0, s_yz = 0;
        for (int k = 0; k < ng; ++k) {
            const double sk2 = sr_[k] * sr_[k] + si_[k] * si_[k];
            const double ug = K * kw[k] * sk2;
            e_recip += ug;
            const double vt = kvterm_[k];
            s_xx += ug * (1.0 + vt * kx[k] * kx[k]);
            s_yy += ug * (1.0 + vt * ky[k] * ky[k]);
            s_zz += ug * (1.0 + vt * kz[k] * kz[k]);
            s_xy += ug * (vt * kx[k] * ky[k]);
            s_xz += ug * (vt * kx[k] * kz[k]);
            s_yz += ug * (vt * ky[k] * kz[k]);
        }

        double q_total = q_local;
#ifdef ATV_HAS_MPI
        if (mpi)
            MPI_Allreduce(MPI_IN_PLACE, &q_total, 1, MPI_DOUBLE, MPI_SUM,
                          dd->cart_comm);
#endif
        const double e_net =
            -ATV::kCoulomb * ATV::kPI / (2.0 * vol * alpha_ * alpha_) *
            q_total * q_total;

        f.scalars.energy += e_recip + e_net;
        f.scalars.energy_pair += e_recip + e_net;
        f.scalars.stress[0] += s_xx;
        f.scalars.stress[1] += s_yy;
        f.scalars.stress[2] += s_zz;
        f.scalars.stress[3] += s_xy;
        f.scalars.stress[4] += s_xz;
        f.scalars.stress[5] += s_yz;
    }
#ifdef ATV_HAS_MPI
    else if (mpi) {
        // Non-root ranks must still participate in the net-charge Allreduce
        // to match the collective issued on rank 0.
        double q_total = q_local;
        MPI_Allreduce(MPI_IN_PLACE, &q_total, 1, MPI_DOUBLE, MPI_SUM,
                      dd->cart_comm);
    }
#endif
}

} // namespace cpu
} // namespace ATV
