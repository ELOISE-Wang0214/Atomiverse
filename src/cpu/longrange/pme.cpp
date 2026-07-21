/// @file src/cpu/longrange/pme.cpp
/// @brief Smooth Particle-Mesh Ewald — Frame/MPI/FFTW wrapper over longrange/.
///
/// Built only when ATV_BUILD_KSPACE is defined (links FFTW3).
#include "atomiverse/cpu/longrange/pme.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

#include <fftw3.h>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/core/units.h"
#include "atomiverse/cpu/component/owned_range.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/longrange/pme_influence.h"
#include "atomiverse/longrange/pme_spline.h"

#ifdef ATV_HAS_MPI
#include <mpi.h>

#include "atomiverse/cpu/comm/domain.h"
#endif

namespace ATV {
namespace cpu {

KSpacePME::KSpacePME(double alpha, double grid_spacing, int order)
    : alpha_(alpha), spacing_(grid_spacing), order_(order) {
    if (alpha <= 0.0)
        throw std::invalid_argument("KSpacePME: alpha must be > 0");
    if (grid_spacing <= 0.0)
        throw std::invalid_argument("KSpacePME: grid_spacing must be > 0");
    if (order < 3)
        throw std::invalid_argument("KSpacePME: order must be ≥ 3");
}

KSpacePME::~KSpacePME() { release_storage(); }

KSpacePME::KSpacePME(KSpacePME&& other) noexcept
    : KSpace(std::move(other)), alpha_(other.alpha_), spacing_(other.spacing_),
      order_(other.order_), box_identity_(other.box_identity_),
      cell_revision_(other.cell_revision_), have_grid_(other.have_grid_),
      grid_(std::move(other.grid_)), plan_fwd_(other.plan_fwd_),
      plan_bwd_(other.plan_bwd_), atom_capacity_(other.atom_capacity_),
      base_(std::move(other.base_)), theta_(std::move(other.theta_)),
      dtheta_(std::move(other.dtheta_)) {
    for (int d = 0; d < 3; ++d) {
        K_[d] = other.K_[d];
        bsp_mod_[d] = std::move(other.bsp_mod_[d]);
        other.K_[d] = 0;
    }
    other.box_identity_ = nullptr;
    other.cell_revision_ = 0;
    other.have_grid_ = false;
    other.plan_fwd_ = nullptr;
    other.plan_bwd_ = nullptr;
    other.atom_capacity_ = 0;
}

void KSpacePME::ensure_grid(const cpu::Frame& frame, bool initial_build) {
    const auto& box = frame.simbox;
    if (have_grid_ && box_identity_ == &box && box.revision() == cell_revision_)
        return;

    int next_K[3];
    for (int d = 0; d < 3; ++d) {
        const double len = box.length(d);
        next_K[d] = longrange::good_fft_size(
            static_cast<int>(std::ceil(len / spacing_)));
        if (next_K[d] < order_ + 1)
            next_K[d] = longrange::good_fft_size(order_ + 1);
    }

    const bool shape_changed = !have_grid_ || next_K[0] != K_[0]
                               || next_K[1] != K_[1] || next_K[2] != K_[2];
    if (have_grid_ && !initial_build) {
        if (shape_changed)
            note_capacity_rebuild();
        else
            note_cell_refresh();
    }

    if (shape_changed) {
        if (plan_fwd_)
            fftw_destroy_plan(static_cast<fftw_plan>(plan_fwd_));
        if (plan_bwd_)
            fftw_destroy_plan(static_cast<fftw_plan>(plan_bwd_));
        plan_fwd_ = nullptr;
        plan_bwd_ = nullptr;

        for (int d = 0; d < 3; ++d) {
            K_[d] = next_K[d];
            bsp_mod_[d] = longrange::spline_moduli(K_[d], order_);
        }
        const int ntot = K_[0] * K_[1] * K_[2];
        grid_.resize(static_cast<std::size_t>(ntot));
        auto* data = reinterpret_cast<fftw_complex*>(grid_.data());
        // fftw_plan_dft_3d(n0,n1,n2): n0 slowest, n2 fastest. Lay out grid
        // with axis 2 fastest: idx = k2 + K2*(k1 + K1*k0).
        plan_fwd_ = fftw_plan_dft_3d(K_[0], K_[1], K_[2], data, data,
                                     FFTW_FORWARD, FFTW_ESTIMATE);
        plan_bwd_ = fftw_plan_dft_3d(K_[0], K_[1], K_[2], data, data,
                                     FFTW_BACKWARD, FFTW_ESTIMATE);
        if (!plan_fwd_ || !plan_bwd_) {
            if (plan_fwd_)
                fftw_destroy_plan(static_cast<fftw_plan>(plan_fwd_));
            if (plan_bwd_)
                fftw_destroy_plan(static_cast<fftw_plan>(plan_bwd_));
            plan_fwd_ = nullptr;
            plan_bwd_ = nullptr;
            throw std::runtime_error("KSpacePME: FFTW plan creation failed");
        }
    }

    box_identity_ = &box;
    cell_revision_ = box.revision();
    have_grid_ = true;
}

void KSpacePME::ensure_atom_capacity(int required, bool initial_build) {
    if (required <= atom_capacity_)
        return;
    const int grown = atom_capacity_ > 0 ? atom_capacity_ * 2 : 1;
    const int next = std::max(required, grown);
    base_.resize(static_cast<std::size_t>(3 * next));
    theta_.resize(static_cast<std::size_t>(3 * order_ * next));
    dtheta_.resize(static_cast<std::size_t>(3 * order_ * next));
    atom_capacity_ = next;
    if (!initial_build)
        note_capacity_rebuild();
}

void KSpacePME::build(const BuildContext& ctx) {
    note_build();
    ensure_grid(ctx.frame, true);
    ensure_atom_capacity(ctx.frame.capacity(), true);
}

void KSpacePME::compute(const LaunchContext& ctx) {
    auto& f = ctx.frame;
    ensure_grid(f, false);

    const auto* dd = ctx.runtime.md.dd;
    const bool mpi = dd && ctx.runtime.mpi_size > 1;
    const int n_local = mpi ? dd->n_local : f.nlocal();
    const bool first_rank = !mpi || ctx.runtime.mpi_rank == 0;
    ensure_atom_capacity(n_local, false);

    const int ntot = K_[0] * K_[1] * K_[2];

    const StorageReal* __restrict x = f.span(kAtom_x);
    const StorageReal* __restrict y = f.span(kAtom_y);
    const StorageReal* __restrict z = f.span(kAtom_z);
    const StorageReal* __restrict q = f.span(kAtom_charge);
    AccReal* __restrict sp_fx = f.span(kAtom_fx);
    AccReal* __restrict sp_fy = f.span(kAtom_fy);
    AccReal* __restrict sp_fz = f.span(kAtom_fz);

    const auto& ic = f.simbox.inv_cell();
    const double vol = f.simbox.volume();

    std::fill(grid_.begin(), grid_.end(), std::complex<double>(0, 0));

    // ── Spread owned charges onto the grid ─────────────────────────────
    longrange::pme_fill_atom_splines(x, y, z, n_local, ic, K_, order_,
                                     base_.data(), theta_.data(),
                                     dtheta_.data());
    longrange::pme_spread_charges(q, n_local, order_, K_, base_.data(),
                                  theta_.data(), grid_.data());

#ifdef ATV_HAS_MPI
    if (mpi)
        MPI_Allreduce(MPI_IN_PLACE, reinterpret_cast<double*>(grid_.data()),
                      2 * ntot, MPI_DOUBLE, MPI_SUM, dd->cart_comm);
#endif

    // ── Forward FFT ────────────────────────────────────────────────────
    fftw_execute_dft(static_cast<fftw_plan>(plan_fwd_),
                     reinterpret_cast<fftw_complex*>(grid_.data()),
                     reinterpret_cast<fftw_complex*>(grid_.data()));

    // ── Influence function: energy + virial + convolution ──────────────
    double e_recip = 0;
    double stress[6];
    longrange::pme_apply_influence(grid_.data(), K_, bsp_mod_[0].data(),
                                   bsp_mod_[1].data(), bsp_mod_[2].data(), ic,
                                   alpha_, vol, e_recip, stress);

    // ── Inverse FFT → real-space convolution grid ──────────────────────
    fftw_execute_dft(static_cast<fftw_plan>(plan_bwd_),
                     reinterpret_cast<fftw_complex*>(grid_.data()),
                     reinterpret_cast<fftw_complex*>(grid_.data()));

    // ── Interpolate reciprocal forces onto owned atoms ─────────────────
    longrange::pme_interpolate_forces(q, n_local, order_, K_, ic, vol,
                                      base_.data(), theta_.data(),
                                      dtheta_.data(), grid_.data(), sp_fx,
                                      sp_fy, sp_fz);

    // ── Self energy (per-rank owned) + net-charge (rank 0) ─────────────
    double qsq_local = 0, q_local = 0;
    for (int i = 0; i < n_local; ++i) {
        qsq_local += q[i] * q[i];
        q_local += q[i];
    }
    // Same formulas as longrange::EwaldSum::{self_energy,net_charge_energy}.
    const double e_self =
        -ATV::units::coulomb * alpha_ / std::sqrt(ATV::units::pi) * qsq_local;
    f.meta.get(kEnergy) += e_self;
    f.meta.get(kEnergyPair) += e_self;

    double q_total = q_local;
#ifdef ATV_HAS_MPI
    if (mpi)
        MPI_Allreduce(MPI_IN_PLACE, &q_total, 1, MPI_DOUBLE, MPI_SUM,
                      dd->cart_comm);
#endif
    if (first_rank) {
        const double e_net = -ATV::units::coulomb * ATV::units::pi
                             / (2.0 * vol * alpha_ * alpha_) * q_total
                             * q_total;
        f.meta.get(kEnergy) += e_recip + e_net;
        f.meta.get(kEnergyPair) += e_recip + e_net;
        // Reciprocal virial only (same as KSpaceEwald): the net-charge term
        // is pure 1/V energy with zero contribution on a neutral cell, and
        // is intentionally omitted from the stress tally for parity.
        f.meta.get(kStress)[0] += stress[0];
        f.meta.get(kStress)[1] += stress[1];
        f.meta.get(kStress)[2] += stress[2];
        f.meta.get(kStress)[3] += stress[3];
        f.meta.get(kStress)[4] += stress[4];
        f.meta.get(kStress)[5] += stress[5];
    }
}

void KSpacePME::release_storage() noexcept {
    if (plan_fwd_)
        fftw_destroy_plan(static_cast<fftw_plan>(plan_fwd_));
    if (plan_bwd_)
        fftw_destroy_plan(static_cast<fftw_plan>(plan_bwd_));
    plan_fwd_ = nullptr;
    plan_bwd_ = nullptr;

    std::vector<std::complex<double>>().swap(grid_);
    for (auto& moduli : bsp_mod_)
        std::vector<double>().swap(moduli);
    std::vector<int>().swap(base_);
    std::vector<double>().swap(theta_);
    std::vector<double>().swap(dtheta_);

    K_[0] = K_[1] = K_[2] = 0;
    atom_capacity_ = 0;
    box_identity_ = nullptr;
    cell_revision_ = 0;
    have_grid_ = false;
}

void KSpacePME::destroy(const DestroyContext&) { release_storage(); }

std::size_t KSpacePME::memory_bytes() const {
    std::size_t bytes = grid_.capacity() * sizeof(std::complex<double>);
    for (const auto& moduli : bsp_mod_)
        bytes += moduli.capacity() * sizeof(double);
    bytes += base_.capacity() * sizeof(int);
    bytes += (theta_.capacity() + dtheta_.capacity()) * sizeof(double);
    return bytes;
}

} // namespace cpu
} // namespace ATV
