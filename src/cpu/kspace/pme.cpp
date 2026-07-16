/// @file src/cpu/kspace/pme.cpp
/// @brief Smooth Particle-Mesh Ewald reciprocal sum implementation (CPU).
///
/// Built only when ATV_BUILD_KSPACE is defined (links FFTW3).
#include "atomiverse/cpu/kspace/pme.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

#include <fftw3.h>

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

namespace {

// Cardinal B-spline weights and derivatives on the n-point stencil for a
// fractional offset `frac` ∈ [0,1). Fills theta[k]=M_n(frac+k) and
// dtheta[k]=dM_n/d(frac) at frac+k, k = 0..n-1. Canonical Essmann/Amber
// `fill_bspline` (J. Chem. Phys. 103, 8577 (1995), appendix).
void bspline(double frac, int n, double* theta, double* dtheta) {
    theta[n - 1] = 0.0;
    theta[1] = frac;
    theta[0] = 1.0 - frac;
    // Standard B-spline up-recursion to order n-1.
    for (int k = 3; k < n; ++k) {
        const double div = 1.0 / (k - 1);
        theta[k - 1] = div * frac * theta[k - 2];
        for (int j = 1; j < k - 1; ++j)
            theta[k - 1 - j] = div
                               * ((frac + j) * theta[k - 2 - j]
                                  + (k - j - frac) * theta[k - 1 - j]);
        theta[0] = div * (1.0 - frac) * theta[0];
    }
    // Derivative taken from the order n-1 values.
    dtheta[0] = -theta[0];
    for (int j = 1; j < n; ++j)
        dtheta[j] = theta[j - 1] - theta[j];
    // Final up-recursion: order n-1 → n.
    const double div = 1.0 / (n - 1);
    theta[n - 1] = div * frac * theta[n - 2];
    for (int j = 1; j < n - 1; ++j)
        theta[n - 1 - j] = div
                           * ((frac + j) * theta[n - 2 - j]
                              + (n - j - frac) * theta[n - 1 - j]);
    theta[0] = div * (1.0 - frac) * theta[0];
}

// |b_d(m)|² Euler spline modulus for grid size K, order n.
std::vector<double> spline_moduli(int K, int n) {
    // Spline values at integers: data[k] = M_n(k), k=0..n-1 (frac=0).
    std::vector<double> data(n, 0.0), ddummy(n, 0.0);
    bspline(0.0, n, data.data(), ddummy.data());
    // fill_bspline(0) returns [M_n(1), M_n(2), …, M_n(n-1), 0], so the spline
    // value M_n(k+1) used in the Euler modulus is data[k].
    std::vector<double> mod(K, 0.0);
    for (int i = 0; i < K; ++i) {
        double sr = 0, si = 0;
        for (int k = 0; k < n - 1; ++k) {
            const double arg = 2.0 * ATV::units::pi * i * k / K;
            sr += data[k] * std::cos(arg);
            si += data[k] * std::sin(arg);
        }
        const double denom = sr * sr + si * si;
        mod[i] = denom > 1e-12 ? 1.0 / denom : 0.0;
    }
    // Fix exact zeros (even orders at i=K/2) by neighbour averaging.
    for (int i = 0; i < K; ++i) {
        if (mod[i] == 0.0) {
            const double lo = mod[(i - 1 + K) % K];
            const double hi = mod[(i + 1) % K];
            mod[i] = 0.5 * (lo + hi);
        }
    }
    return mod;
}

int good_fft_size(int n) {
    // Smallest m ≥ n that factorises into {2,3,5,7}.
    auto ok = [](int m) {
        for (int p : {2, 3, 5, 7})
            while (m % p == 0)
                m /= p;
        return m == 1;
    };
    if (n < 4)
        n = 4;
    while (!ok(n))
        ++n;
    return n;
}

} // namespace

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
        next_K[d] = good_fft_size(static_cast<int>(std::ceil(len / spacing_)));
        if (next_K[d] < order_ + 1)
            next_K[d] = good_fft_size(order_ + 1);
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
            bsp_mod_[d] = spline_moduli(K_[d], order_);
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

    const int K0 = K_[0], K1 = K_[1], K2 = K_[2];
    const int ntot = K0 * K1 * K2;
    const int n = order_;

    const StorageReal* __restrict x = f.span(kAtom_x);
    const StorageReal* __restrict y = f.span(kAtom_y);
    const StorageReal* __restrict z = f.span(kAtom_z);
    const StorageReal* __restrict q = f.span(kAtom_charge);
    AccReal* __restrict sp_fx = f.span(kAtom_fx);
    AccReal* __restrict sp_fy = f.span(kAtom_fy);
    AccReal* __restrict sp_fz = f.span(kAtom_fz);

    const auto& ic =
        f.simbox.inv_cell(); // fractional coord: s_d = Σ_b ic(b,d) r_b
    const double vol = f.simbox.volume();

    // Per-atom spline tables live in persistent component-owned workspace.
    auto& base = base_;
    auto& th = theta_;
    auto& dth = dtheta_;

    std::fill(grid_.begin(), grid_.end(), std::complex<double>(0, 0));

    // ── Spread owned charges onto the grid ─────────────────────────────
    for (int i = 0; i < n_local; ++i) {
        const double qi = q[i];
        double s[3];
        for (int d = 0; d < 3; ++d)
            s[d] = ic(0, d) * x[i] + ic(1, d) * y[i] + ic(2, d) * z[i];
        for (int d = 0; d < 3; ++d) {
            double u = s[d] * K_[d];
            u -= std::floor(u / K_[d]) * K_[d]; // wrap into [0,K)
            const int u0 = static_cast<int>(std::floor(u));
            const double frac = u - u0;
            bspline(frac, n, &th[(3 * i + d) * n], &dth[(3 * i + d) * n]);
            base[3 * i + d] = u0;
        }
        if (qi == 0.0)
            continue;
        const double* thx = &th[(3 * i + 0) * n];
        const double* thy = &th[(3 * i + 1) * n];
        const double* thz = &th[(3 * i + 2) * n];
        // theta[k] = M_n(frac + n-1-k) pairs with grid point u0 - (n-1) + k.
        const int nm1 = n - 1;
        for (int a = 0; a < n; ++a) {
            const int g0 = ((base[3 * i + 0] - nm1 + a) % K0 + K0) % K0;
            const double wa = qi * thx[a];
            for (int b = 0; b < n; ++b) {
                const int g1 = ((base[3 * i + 1] - nm1 + b) % K1 + K1) % K1;
                const double wab = wa * thy[b];
                for (int c = 0; c < n; ++c) {
                    const int g2 = ((base[3 * i + 2] - nm1 + c) % K2 + K2) % K2;
                    const int idx = g2 + K2 * (g1 + K1 * g0);
                    grid_[idx] += wab * thz[c];
                }
            }
        }
    }

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
    // Reciprocal energy (Essmann et al.):
    //   E = (k_e/(2πV)) Σ_{m≠0} exp(−π²|m|²/α²)/|m|² · B(m) · |Q(m)|²
    // with m built from inv_cell (no explicit 2π; the 2π is absorbed into
    // the prefactor). The per-mode virial matches KSpaceEwald after the
    // k=2πm substitution:
    //   V_ab = ug · (δ_ab + vt · m_a m_b),  vt = −2(1/|m|² + π²/α²).
    // Rank 0 alone accumulates energy/stress; every rank still applies the
    // influence function so the inverse FFT yields forces on owned atoms.
    const double inv_pi2_a2 =
        ATV::units::pi * ATV::units::pi / (alpha_ * alpha_);
    const double vt_const = inv_pi2_a2; // π²/α² term inside vt
    double e_recip = 0;
    double s_xx = 0, s_yy = 0, s_zz = 0, s_xy = 0, s_xz = 0, s_yz = 0;
    for (int m0 = 0; m0 < K0; ++m0) {
        const int h0 = (m0 <= K0 / 2) ? m0 : m0 - K0;
        const double bm0 = bsp_mod_[0][m0];
        for (int m1 = 0; m1 < K1; ++m1) {
            const int h1 = (m1 <= K1 / 2) ? m1 : m1 - K1;
            const double bm01 = bm0 * bsp_mod_[1][m1];
            for (int m2 = 0; m2 < K2; ++m2) {
                const int idx = m2 + K2 * (m1 + K1 * m0);
                if (m0 == 0 && m1 == 0 && m2 == 0) {
                    grid_[idx] = 0;
                    continue;
                }
                const int h2 = (m2 <= K2 / 2) ? m2 : m2 - K2;
                // Reciprocal vector m = Σ_d h_d · (column d of cell⁻¹).
                double mx = 0, my = 0, mz = 0;
                const int hh[3] = {h0, h1, h2};
                for (int d = 0; d < 3; ++d) {
                    mx += hh[d] * ic(0, d);
                    my += hh[d] * ic(1, d);
                    mz += hh[d] * ic(2, d);
                }
                const double msq = mx * mx + my * my + mz * mz;
                const double B = bm01 * bsp_mod_[2][m2];
                const double H = std::exp(-inv_pi2_a2 * msq) / msq * B;
                const std::complex<double> Qm = grid_[idx];
                const double q2 = Qm.real() * Qm.real() + Qm.imag() * Qm.imag();
                e_recip += H * q2;
                // Virial geometric factor in m-space (see header comment).
                const double vt = -2.0 * (1.0 / msq + vt_const);
                const double ug = H * q2; // scaled by escale below
                s_xx += ug * (1.0 + vt * mx * mx);
                s_yy += ug * (1.0 + vt * my * my);
                s_zz += ug * (1.0 + vt * mz * mz);
                s_xy += ug * (vt * mx * my);
                s_xz += ug * (vt * mx * mz);
                s_yz += ug * (vt * my * mz);
                grid_[idx] = H * Qm;
            }
        }
    }
    const double escale = ATV::units::coulomb / (2.0 * ATV::units::pi * vol);
    e_recip *= escale;
    s_xx *= escale;
    s_yy *= escale;
    s_zz *= escale;
    s_xy *= escale;
    s_xz *= escale;
    s_yz *= escale;

    // ── Inverse FFT → real-space convolution grid ──────────────────────
    fftw_execute_dft(static_cast<fftw_plan>(plan_bwd_),
                     reinterpret_cast<fftw_complex*>(grid_.data()),
                     reinterpret_cast<fftw_complex*>(grid_.data()));

    // ── Interpolate reciprocal forces onto owned atoms ─────────────────
    // F_{i,a} = -(units::coulomb/πV) q_i Σ_stencil C(grid) Σ_d dθ_d Π_{d'≠d} θ_{d'}
    //           · K_d · ic(a,d)
    const double fscale = ATV::units::coulomb / (ATV::units::pi * vol);
    for (int i = 0; i < n_local; ++i) {
        const double qi = q[i];
        if (qi == 0.0)
            continue;
        const double* thx = &th[(3 * i + 0) * n];
        const double* thy = &th[(3 * i + 1) * n];
        const double* thz = &th[(3 * i + 2) * n];
        const double* dthx = &dth[(3 * i + 0) * n];
        const double* dthy = &dth[(3 * i + 1) * n];
        const double* dthz = &dth[(3 * i + 2) * n];
        double g_d[3] = {0, 0, 0}; // gradient in (scaled fractional) dims
        const int nm1 = n - 1;
        for (int a = 0; a < n; ++a) {
            const int g0 = ((base[3 * i + 0] - nm1 + a) % K0 + K0) % K0;
            for (int b = 0; b < n; ++b) {
                const int g1 = ((base[3 * i + 1] - nm1 + b) % K1 + K1) % K1;
                for (int c = 0; c < n; ++c) {
                    const int g2 = ((base[3 * i + 2] - nm1 + c) % K2 + K2) % K2;
                    const int idx = g2 + K2 * (g1 + K1 * g0);
                    const double cval = grid_[idx].real();
                    g_d[0] += cval * dthx[a] * thy[b] * thz[c];
                    g_d[1] += cval * thx[a] * dthy[b] * thz[c];
                    g_d[2] += cval * thx[a] * thy[b] * dthz[c];
                }
            }
        }
        // Chain rule: ∂E/∂r_a = fscale·q_i·Σ_d g_d·K_d·ic(a,d), with
        // dθ = dM/d(frac) and ∂(frac_d)/∂r_a = K_d·ic(a,d). Force = −∂E/∂r.
        double fx = 0, fy = 0, fz = 0;
        for (int d = 0; d < 3; ++d) {
            const double gk = g_d[d] * K_[d];
            fx += gk * ic(0, d);
            fy += gk * ic(1, d);
            fz += gk * ic(2, d);
        }
        sp_fx[i] += -fscale * qi * fx;
        sp_fy[i] += -fscale * qi * fy;
        sp_fz[i] += -fscale * qi * fz;
    }

    // ── Self energy (per-rank owned) + net-charge (rank 0) ─────────────
    double qsq_local = 0, q_local = 0;
    for (int i = 0; i < n_local; ++i) {
        qsq_local += q[i] * q[i];
        q_local += q[i];
    }
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
        f.meta.get(kStress)[0] += s_xx;
        f.meta.get(kStress)[1] += s_yy;
        f.meta.get(kStress)[2] += s_zz;
        f.meta.get(kStress)[3] += s_xy;
        f.meta.get(kStress)[4] += s_xz;
        f.meta.get(kStress)[5] += s_yz;
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
