/// @file src/cpu/gpw/g0w0.cpp
/// @brief GpwDft::g0w0 — Γ-point spectral G0W0@KS on the converged RKS
///        orbitals. dRPA screening + MT exchange + local v_xc.

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

#include "atomiverse/core/units.h"
#include "atomiverse/cpu/gpw/gpw_solver.h"
#include "atomiverse/cpu/gpw/internal/mp2_kernels.h"
#include "atomiverse/cpu/gpw/internal/xc_grid.h"
#include "atomiverse/cpu/gpw/internal/xc_grid_gga.h"

namespace ATV {
namespace cpu {
namespace gpw {
namespace {

constexpr double kEta = 1e-8;    ///< Spectral broadening [Hartree].
constexpr double kQpeTol = 1e-8; ///< Newton residual on ω [Hartree].
constexpr int kQpeMaxIter = 100;

/// Diagonal Re Σ_c,p(ω) [Hartree]. tdm row-major [s, m, p] with
/// strides (n_mo * n_mo, n_mo, 1).
double sigma_c_diag(double omega, int p, int n_occ, int n_mo,
                    const std::vector<double>& eps,
                    const std::vector<double>& Om,
                    const std::vector<double>& tdm) {
    const int n_exc = static_cast<int>(Om.size());
    const std::size_t nmo = static_cast<std::size_t>(n_mo);
    double re = 0.0;
    for (int s = 0; s < n_exc; ++s) {
        const double Os = Om[static_cast<std::size_t>(s)];
        const std::size_t s0 =
            static_cast<std::size_t>(s) * nmo * nmo
            + static_cast<std::size_t>(p);
        for (int m = 0; m < n_occ; ++m) {
            const double t = tdm[s0 + static_cast<std::size_t>(m) * nmo];
            const double den_r = omega + Os - eps[static_cast<std::size_t>(m)];
            const double den2 = den_r * den_r + kEta * kEta;
            re += (t * t) * den_r / den2;
        }
        for (int m = n_occ; m < n_mo; ++m) {
            const double t = tdm[s0 + static_cast<std::size_t>(m) * nmo];
            const double den_r = omega - Os - eps[static_cast<std::size_t>(m)];
            const double den2 = den_r * den_r + kEta * kEta;
            re += (t * t) * den_r / den2;
        }
    }
    return re;
}

bool solve_qp(double eps_p, double vx_m_vxc, int p, int n_occ, int n_mo,
              const std::vector<double>& eps, const std::vector<double>& Om,
              const std::vector<double>& tdm, double& omega_out) {
    double w = eps_p;
    constexpr double kH = 1e-6;
    for (int it = 0; it < kQpeMaxIter; ++it) {
        const double sig = sigma_c_diag(w, p, n_occ, n_mo, eps, Om, tdm);
        const double f = w - eps_p - (sig + vx_m_vxc);
        if (std::abs(f) < kQpeTol) {
            omega_out = w;
            return true;
        }
        const double sig_p =
            sigma_c_diag(w + kH, p, n_occ, n_mo, eps, Om, tdm);
        const double dsig = (sig_p - sig) / kH;
        const double df = 1.0 - dsig;
        if (!(std::abs(df) > 1e-14))
            return false;
        w -= f / df;
    }
    omega_out = w;
    return std::abs(w - eps_p
                    - (sigma_c_diag(w, p, n_occ, n_mo, eps, Om, tdm)
                       + vx_m_vxc))
           < 10.0 * kQpeTol;
}

} // namespace

const G0W0Result& GpwDft::g0w0(const std::vector<GpwAtom>& atoms) {
    if (!built_)
        throw std::logic_error(
            "GpwDft::g0w0 called before build() — validate the "
            "configuration first (no silent auto-build)");
    if (is_uks())
        throw std::invalid_argument(
            "GpwDft::g0w0: closed-shell RKS only (UGW is a follow-on)");
    if (use_ot)
        throw std::invalid_argument(
            "GpwDft::g0w0: OT stores occupied orbitals only — canonical "
            "virtuals/eigenvalues require the diagonalization path "
            "(use_ot = false)");
    if (!result_.converged)
        throw std::logic_error(
            "GpwDft::g0w0 needs a converged run() first (no correlation "
            "on an unconverged reference)");
    if (result_.n_kpts > 1)
        throw std::invalid_argument(
            "GpwDft::g0w0: multi-k SCF reference not supported "
            "(k-point GW is a residual; re-run with k_mesh = {1,1,1} for "
            "Γ G0W0)");

    const int nb_atoms = [&] {
        int nb = 0;
        for (const auto& a : atoms)
            nb += static_cast<int>(a.basis.size());
        return nb;
    }();
    const int nb = static_cast<int>(result_.C.shape(0));
    const int n_mo = static_cast<int>(result_.C.shape(1));
    if (nb_atoms != nb)
        throw std::invalid_argument(
            "GpwDft::g0w0: atom basis count (" + std::to_string(nb_atoms)
            + ") does not match the converged orbitals (" + std::to_string(nb)
            + ") — same geometry/basis as the last run() required");
    if (static_cast<int>(result_.eps.size()) != n_mo)
        throw std::logic_error(
            "GpwDft::g0w0: eps/C shape mismatch in the converged result");
    if (n_mo <= n_occ)
        throw std::invalid_argument(
            "GpwDft::g0w0: no virtual orbitals in the converged result "
            "(n_mo = "
            + std::to_string(n_mo) + ", n_occ = " + std::to_string(n_occ)
            + ")");

    int n_virt = n_mo - n_occ;
    if (mp2_n_virt > 0)
        n_virt = std::min(n_virt, mp2_n_virt);

    constexpr double kMinGapHa = 1e-6;
    const double gap = result_.eps[static_cast<std::size_t>(n_occ)]
                       - result_.eps[static_cast<std::size_t>(n_occ - 1)];
    if (gap < kMinGapHa)
        throw std::invalid_argument(
            "GpwDft::g0w0: vanishing HOMO-LUMO gap (" + std::to_string(gap)
            + " Ha < 1e-6) — dRPA (A−B) is singular");

    const bool periodic = !nonperiodic;
    const UniformGrid grid = make_uniform_grid(cell, cutoff_ry);
    std::vector<GridBasisFunction> basis;
    basis.reserve(static_cast<std::size_t>(nb));
    for (const auto& a : atoms)
        for (const auto& bf : a.basis)
            basis.push_back(bf);
    const CompactBasisGrid cbg =
        build_compact_basis_grid(grid, basis, periodic, eps_default);

    // All MOs on the density grid (screening virtuals are a prefix of the
    // virtual block; QP levels cover the full spectrum).
    const std::vector<double> psi_all =
        internal::mo_on_grid(cbg, result_.C, 0, n_mo);
    const long npts = grid.npts();
    const double gw = grid.weight();

    std::vector<double> psi_occ(
        static_cast<std::size_t>(n_occ) * static_cast<std::size_t>(npts));
    std::vector<double> psi_virt(
        static_cast<std::size_t>(n_virt) * static_cast<std::size_t>(npts));
    for (int i = 0; i < n_occ; ++i)
        for (long g = 0; g < npts; ++g)
            psi_occ[static_cast<std::size_t>(i)
                        * static_cast<std::size_t>(npts)
                    + static_cast<std::size_t>(g)] =
                psi_all[static_cast<std::size_t>(i)
                            * static_cast<std::size_t>(npts)
                        + static_cast<std::size_t>(g)];
    for (int a = 0; a < n_virt; ++a)
        for (long g = 0; g < npts; ++g)
            psi_virt[static_cast<std::size_t>(a)
                         * static_cast<std::size_t>(npts)
                     + static_cast<std::size_t>(g)] =
                psi_all[static_cast<std::size_t>(n_occ + a)
                            * static_cast<std::size_t>(npts)
                        + static_cast<std::size_t>(g)];

    const bool use_ri = mp2_ri && periodic;
    const std::vector<double> v_ovov =
        use_ri ? internal::mo_eri_ovov_ri(grid, psi_occ, n_occ, psi_virt, n_virt)
               : internal::mo_eri_ovov(grid, psi_occ, n_occ, psi_virt, n_virt,
                                       nonperiodic);

    const std::vector<double> eps_occ(result_.eps.begin(),
                                      result_.eps.begin() + n_occ);
    const std::vector<double> eps_virt(result_.eps.begin() + n_occ,
                                       result_.eps.begin() + n_occ + n_virt);

    internal::DirectRpaAmplitudes rpa;
    internal::direct_rpa_rhf(v_ovov, eps_occ, eps_virt, rpa);
    const int n_exc = static_cast<int>(rpa.omega.size());
    const int nov = n_occ * n_virt;

    // V_ia(r) intermediates (periodic Coulomb; ρ_ia neutral).
    const std::vector<double> V_ov = internal::mo_pair_potentials_ov(
        grid, psi_occ, n_occ, psi_virt, n_virt);

    // Transition densities tdm[s, m, p] = Σ_ia √2 (X+Y)_ia^s (ia|mp)
    // Layout: [s][m][p] row-major, strides (n_mo*n_mo, n_mo, 1).
    // For each s form W_s(g) = Σ_ia √2 (X+Y)_ia^s V_ia(g), then
    // tdm[s,m,p] = w Σ_g W_s(g) ψ_m(g) ψ_p(g).
    const double sqrt2 = std::sqrt(2.0);
    std::vector<double> tdm(static_cast<std::size_t>(n_exc)
                                * static_cast<std::size_t>(n_mo)
                                * static_cast<std::size_t>(n_mo),
                            0.0);
    std::vector<double> Ws(static_cast<std::size_t>(npts));
    for (int s = 0; s < n_exc; ++s) {
        std::fill(Ws.begin(), Ws.end(), 0.0);
        for (int ia = 0; ia < nov; ++ia) {
            const double z =
                sqrt2
                * (rpa.X[static_cast<std::size_t>(s)
                             * static_cast<std::size_t>(nov)
                         + static_cast<std::size_t>(ia)]
                   + rpa.Y[static_cast<std::size_t>(s)
                               * static_cast<std::size_t>(nov)
                           + static_cast<std::size_t>(ia)]);
            if (z == 0.0)
                continue;
            const double* Via =
                V_ov.data()
                + static_cast<std::size_t>(ia) * static_cast<std::size_t>(npts);
            for (long g = 0; g < npts; ++g)
                Ws[static_cast<std::size_t>(g)] += z * Via[g];
        }
        const std::size_t s0 = static_cast<std::size_t>(s)
                               * static_cast<std::size_t>(n_mo)
                               * static_cast<std::size_t>(n_mo);
        for (int m = 0; m < n_mo; ++m) {
            const double* pm =
                psi_all.data()
                + static_cast<std::size_t>(m) * static_cast<std::size_t>(npts);
            for (int p = 0; p < n_mo; ++p) {
                const double* pp =
                    psi_all.data()
                    + static_cast<std::size_t>(p)
                          * static_cast<std::size_t>(npts);
                double acc = 0.0;
                for (long g = 0; g < npts; ++g)
                    acc += Ws[static_cast<std::size_t>(g)] * pm[g] * pp[g];
                tdm[s0 + static_cast<std::size_t>(m)
                             * static_cast<std::size_t>(n_mo)
                         + static_cast<std::size_t>(p)] = gw * acc;
            }
        }
    }

    // Exact-exchange diagonals (MT Coulomb).
    const std::vector<double> sigma_x =
        internal::exchange_diag_mt(grid, psi_all, n_mo, n_occ);

    // Local v_xc diagonals from the same functional / closed-shell density.
    std::vector<double> rho(static_cast<std::size_t>(npts), 0.0);
    for (int i = 0; i < n_occ; ++i) {
        const double* pi =
            psi_all.data()
            + static_cast<std::size_t>(i) * static_cast<std::size_t>(npts);
        for (long g = 0; g < npts; ++g)
            rho[static_cast<std::size_t>(g)] += 2.0 * pi[g] * pi[g];
    }
    const XcOnGrid xc_grid =
        (this->xc == XcKind::PBE)
            ? xc_grid_pbe(grid, rho)
            : (this->xc == XcKind::PBEsol
                   ? xc_grid_pbesol(grid, rho)
                   : (this->xc == XcKind::RevPBE
                          ? xc_grid_revpbe(grid, rho)
                          : (this->xc == XcKind::BLYP
                                 ? xc_grid_blyp(grid, rho)
                                 : xc_grid_pade(grid, rho))));

    std::vector<double> vxc_diag(static_cast<std::size_t>(n_mo), 0.0);
    for (int p = 0; p < n_mo; ++p) {
        const double* pp =
            psi_all.data()
            + static_cast<std::size_t>(p) * static_cast<std::size_t>(npts);
        double acc = 0.0;
        for (long g = 0; g < npts; ++g)
            acc +=
                pp[g] * pp[g] * xc_grid.vxc_r[static_cast<std::size_t>(g)];
        vxc_diag[static_cast<std::size_t>(p)] = gw * acc;
    }

    std::vector<double> eps = result_.eps; // full spectrum

    g0w0_result_ = G0W0Result{};
    g0w0_result_.n_occ = n_occ;
    g0w0_result_.n_virt = n_virt;
    g0w0_result_.ri = use_ri;
    g0w0_result_.mf_energy_ha = eps;
    g0w0_result_.qp_energy_ha.resize(static_cast<std::size_t>(n_mo));
    g0w0_result_.sigma_x_ha = sigma_x;
    g0w0_result_.vxc_ha = vxc_diag;
    g0w0_result_.converged = true;

    for (int p = 0; p < n_mo; ++p) {
        const double vx_m_vxc = sigma_x[static_cast<std::size_t>(p)]
                                - vxc_diag[static_cast<std::size_t>(p)];
        double qp = eps[static_cast<std::size_t>(p)];
        const bool ok =
            solve_qp(eps[static_cast<std::size_t>(p)], vx_m_vxc, p, n_occ,
                     n_mo, eps, rpa.omega, tdm, qp);
        if (!ok)
            g0w0_result_.converged = false;
        g0w0_result_.qp_energy_ha[static_cast<std::size_t>(p)] = qp;
    }

    const int homo = n_occ - 1;
    const int lumo = n_occ;
    const double e_h =
        g0w0_result_.qp_energy_ha[static_cast<std::size_t>(homo)];
    const double e_l =
        g0w0_result_.qp_energy_ha[static_cast<std::size_t>(lumo)];
    const double mf_h = eps[static_cast<std::size_t>(homo)];
    const double mf_l = eps[static_cast<std::size_t>(lumo)];
    g0w0_result_.ip_ev = -e_h * units::kHartreeInEv;
    g0w0_result_.ea_ev = e_l * units::kHartreeInEv;
    g0w0_result_.gap_ev = (e_l - e_h) * units::kHartreeInEv;
    g0w0_result_.mf_gap_ev = (mf_l - mf_h) * units::kHartreeInEv;
    return g0w0_result_;
}

} // namespace gpw
} // namespace cpu
} // namespace ATV
