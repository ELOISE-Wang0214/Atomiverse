/// @file src/cpu/kspace/ewald.cpp
/// @brief Reciprocal-space Ewald summation implementation (CPU).
#include "atomiverse/cpu/kspace/ewald.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/core/units.h"
#include "atomiverse/cpu/component/owned_range.h"
#include "atomiverse/cpu/frame/frame.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef ATV_HAS_MPI
#include <mpi.h>

#include "atomiverse/cpu/comm/domain.h"
#endif

namespace ATV {
namespace cpu {

KSpaceEwald::KSpaceEwald(double alpha, double accuracy)
    : alpha_(alpha), accuracy_(accuracy) {
    if (alpha <= 0.0)
        throw std::invalid_argument("KSpaceEwald: alpha must be > 0");
    if (!(accuracy > 0.0 && accuracy < 1.0))
        throw std::invalid_argument("KSpaceEwald: accuracy must be in (0, 1)");
}

double KSpaceEwald::recommend_alpha(double cutoff, double accuracy) {
    // LAMMPS heuristic: erfc(g·rc)/rc ≈ accuracy.
    return (1.35 - 0.15 * std::log(accuracy)) / cutoff;
}

void KSpaceEwald::ensure_kvectors(const cpu::Frame& frame, bool initial_build) {
    const auto& box = frame.simbox;
    if (have_kvectors_ && box_identity_ == &box
        && box.revision() == cell_revision_)
        return;

    // Reciprocal basis b_i = 2π · (column i of cell⁻¹); g = Σ m_i b_i.
    const auto& ic = box.inv_cell();
    const double tp = 2.0 * ATV::units::pi;
    double b[3][3];
    for (int i = 0; i < 3; ++i) {
        b[i][0] = tp * ic(0, i);
        b[i][1] = tp * ic(1, i);
        b[i][2] = tp * ic(2, i);
    }
    auto blen = [&](int i) {
        return std::sqrt(b[i][0] * b[i][0] + b[i][1] * b[i][1]
                         + b[i][2] * b[i][2]);
    };
    const double gcut = 2.0 * alpha_ * std::sqrt(-std::log(accuracy_));
    const double gcut2 = gcut * gcut;
    int kmax[3];
    for (int i = 0; i < 3; ++i)
        kmax[i] = std::max(1, static_cast<int>(std::ceil(gcut / blen(i))));

    const std::size_t max_vectors =
        static_cast<std::size_t>(2 * kmax[0] + 1)
            * static_cast<std::size_t>(2 * kmax[1] + 1)
            * static_cast<std::size_t>(2 * kmax[2] + 1)
        - 1;
    const bool capacity_grows =
        kx_.capacity() < max_vectors || ky_.capacity() < max_vectors
        || kz_.capacity() < max_vectors || kw_.capacity() < max_vectors
        || kvterm_.capacity() < max_vectors || sr_.capacity() < max_vectors
        || si_.capacity() < max_vectors;
    if (capacity_grows) {
        kx_.reserve(max_vectors);
        ky_.reserve(max_vectors);
        kz_.reserve(max_vectors);
        kw_.reserve(max_vectors);
        kvterm_.reserve(max_vectors);
        sr_.reserve(max_vectors);
        si_.reserve(max_vectors);
    }
    if (have_kvectors_ && !initial_build) {
        if (capacity_grows)
            note_capacity_rebuild();
        else
            note_cell_refresh();
    }

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
    sr_.resize(kx_.size());
    si_.resize(kx_.size());
    std::fill(sr_.begin(), sr_.end(), 0.0);
    std::fill(si_.begin(), si_.end(), 0.0);
    box_identity_ = &box;
    cell_revision_ = box.revision();
    have_kvectors_ = true;
}

void KSpaceEwald::build(const BuildContext& ctx) {
    note_build();
    ensure_kvectors(ctx.frame, true);
}

void KSpaceEwald::compute(const LaunchContext& ctx) {
    auto& f = ctx.frame;
    ensure_kvectors(f, false);
    const int ng = static_cast<int>(kx_.size());

    const auto* dd = ctx.runtime.md.dd;
    const bool mpi = dd && ctx.runtime.mpi_size > 1;
    const int n_local = mpi ? dd->n_local : f.nlocal();
    const bool first_rank = !mpi || ctx.runtime.mpi_rank == 0;

    const StorageReal* __restrict x = f.span(kAtom_x);
    const StorageReal* __restrict y = f.span(kAtom_y);
    const StorageReal* __restrict z = f.span(kAtom_z);
    const StorageReal* __restrict q = f.span(kAtom_charge);
    AccReal* __restrict sp_fx = f.span(kAtom_fx);
    AccReal* __restrict sp_fy = f.span(kAtom_fy);
    AccReal* __restrict sp_fz = f.span(kAtom_fz);

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

    const double vol = f.simbox.volume();
    const double K = ATV::units::coulomb * 2.0 * ATV::units::pi / vol;

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
        sp_fx[i] += fx;
        sp_fy[i] += fy;
        sp_fz[i] += fz;
    }

    // ── Self energy (partitions across ranks: owned atoms only) ────────
    double qsq_local = 0.0, q_local = 0.0;
    for (int i = 0; i < n_local; ++i) {
        qsq_local += q[i] * q[i];
        q_local += q[i];
    }
    const double e_self =
        -ATV::units::coulomb * alpha_ / std::sqrt(ATV::units::pi) * qsq_local;
    f.meta.get(kEnergy) += e_self;
    f.meta.get(kEnergyPair) += e_self;

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
        const double e_net = -ATV::units::coulomb * ATV::units::pi
                             / (2.0 * vol * alpha_ * alpha_) * q_total
                             * q_total;

        f.meta.get(kEnergy) += e_recip + e_net;
        f.meta.get(kEnergyPair) += e_recip + e_net;
        f.meta.get(kStress)[0] += s_xx;
        f.meta.get(kStress)[1] += s_yy;
        f.meta.get(kStress)[2] += s_zz;
        f.meta.get(kStress)[3] += s_xy;
        f.meta.get(kStress)[4] += s_xz;
        f.meta.get(kStress)[5] += s_yz;
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

void KSpaceEwald::destroy(const DestroyContext&) {
    std::vector<double>().swap(kx_);
    std::vector<double>().swap(ky_);
    std::vector<double>().swap(kz_);
    std::vector<double>().swap(kw_);
    std::vector<double>().swap(kvterm_);
    std::vector<double>().swap(sr_);
    std::vector<double>().swap(si_);
    box_identity_ = nullptr;
    cell_revision_ = 0;
    have_kvectors_ = false;
}

std::size_t KSpaceEwald::memory_bytes() const {
    return (kx_.capacity() + ky_.capacity() + kz_.capacity() + kw_.capacity()
            + kvterm_.capacity() + sr_.capacity() + si_.capacity())
           * sizeof(double);
}

} // namespace cpu
} // namespace ATV
