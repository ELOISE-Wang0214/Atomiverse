/// @file cpu/scf/gradient.cpp
/// @brief Unified analytic nuclear gradient for all mean-field methods.
///
/// Single function `scf_gradient()` handles RHF, UHF, RKS, and UKS.
/// Method dispatch:
///   HF: one-electron + Coulomb + exchange + Pulay + nuclear
///   KS: one-electron + Coulomb + XC_grid  + Pulay + nuclear

#include "atomiverse/cpu/scf/gradient.h"
#include "atomiverse/cpu/scf/integral_derivatives.h"
#include "atomiverse/cpu/scf/basis_eval.h"
#include "atomiverse/cpu/scf/grid.h"
#include "atomiverse/cpu/scf/xc_functional.h"
#include "atomiverse/core/log.h"
#include "atomiverse/core/units.h"
#include "atomiverse/core/constants.h"

#include <xtensor/views/xview.hpp>
#include <xtensor/io/xio.hpp>
#include <xtensor-blas/xlinalg.hpp>
#include <cmath>
#include <functional>
#include <utility>

namespace ATV {
namespace cpu {
namespace scf {

// ── Shared helpers ───────────────────────────────────────────────────

static std::pair<double, double> gradient_metrics(const xt::xtensor<double, 2>& grad)
{
    const int natm = static_cast<int>(grad.shape()[0]);
    double max_abs = 0.0;
    double max_atom_norm = 0.0;
    for (int A = 0; A < natm; ++A) {
        double atom_norm2 = 0.0;
        for (int d = 0; d < 3; ++d) {
            const double v = grad(A, d);
            max_abs = std::max(max_abs, std::abs(v));
            atom_norm2 += v * v;
        }
        max_atom_norm = std::max(max_atom_norm, std::sqrt(atom_norm2));
    }
    return {max_abs, max_atom_norm};
}

static void accumulate_W(
    xt::xtensor<double, 2>& W,
    const xt::xtensor<double, 2>& C,
    const xt::xtensor<double, 1>& eigenvalues,
    int n_occ, int /*n_basis*/)
{
    // W += Σ_i ε_i · c_i · c_i^T   (energy-weighted density matrix)
    for (int i = 0; i < n_occ; ++i)
        for (int mu = 0; mu < static_cast<int>(C.shape()[0]); ++mu)
            for (int nu = 0; nu < static_cast<int>(C.shape()[0]); ++nu)
                W(mu, nu) += eigenvalues(i) * C(mu, i) * C(nu, i);
}

static void accumulate_common_gradient(
    xt::xtensor<double, 2>& grad,
    const xt::xtensor<double, 2>& P,
    const xt::xtensor<double, 2>& W,
    const IntegralDerivPack& dpack,
    int /*n_basis*/, int natm)
{
    // grad(A,d) = Tr[P · (dT + dV)] − Tr[W · dS] + dV_nuc
    // Each (A, d) slice of dpack.dT/dV/dS is an [n×n] matrix.
    for (int A = 0; A < natm; ++A) {
        for (int d = 0; d < 3; ++d) {
            auto dH = xt::view(dpack.dT, A, d, xt::all(), xt::all())
                     + xt::view(dpack.dV, A, d, xt::all(), xt::all());
            auto dS = xt::view(dpack.dS, A, d, xt::all(), xt::all());
            double val = xt::sum(P * dH)() - xt::sum(W * dS)();
            val += dpack.dV_nuc(A, d);
            grad(A, d) += val;
        }
    }
}

static GradientResult make_result(
    const xt::xtensor<double, 2>& grad, double energy_hartree, int natm)
{
    // SCF gradient lives in Hartree/Bohr; Frame / MD drivers expect eV/Å.
    using units::convert;
    GradientResult out;
    out.forces = xt::zeros<double>({(int)natm, (int)3});
    for (int A = 0; A < natm; ++A)
        for (int d = 0; d < 3; ++d)
            out.forces(A, d) =
                -convert(grad(A, d), units::Force::HartreePerBohr, units::Force::eVPerAngstrom);
    out.energy = convert(energy_hartree, units::Energy::Hartree, units::Energy::eV);
    return out;
}

/// @brief Accumulate the exchange-correlation contribution dE_xc/dR_A.
///
/// For the LDA-style (vrho) part:
///   dE_xc/dR_{A,d} = sum_g w * [vrho_a * drho_a/dR_{A,d}
///                             + vrho_b * drho_b/dR_{A,d}]
///   drho_a/dR_{A,d} = -2 * sum_{mu in A} sum_nu P_a(mu,nu) dphi_mu_d phi_nu
/// In matrix form: drho_a/dR_{A,d}[g] = -2 * sum_{mu in A} T_a(g,mu) dPhi_d(g,mu)
/// with T_a = Phi @ P_a (precomputed via GEMM). For a GGA the per-point block
/// additionally contracts vsigma with d(grad_rho)/dR_A, whose nuclear
/// derivative needs the basis-function Hessian (Term2 below).
///
/// @par Diagnosed defect (spec phase1-productization-02-ks-dshell-gradient)
/// The KS gradient was wrong by ~4 orders of magnitude. The defect was NOT in
/// the d-shell polynomial-lowering term — the `poly_grad`/`poly_hess` tables in
/// `basis_eval.h` are FD-verified correct (`test_basis_eval_dshell`), and it
/// was not l-dependent at all. Two independent root causes were found, both
/// restricted-KS-wide; they only surfaced on d-bearing bases because the sole
/// restricted-KS gradient exercise (`default_config`) happened to use one:
///
///  1. **Spin-density convention.** A converged restricted SCF stores
///     `result.P_alpha == result.P_beta == P_total` (the full double-occupied
///     density, factor 2 folded in). The XC functional is spin-resolved and
///     must be fed `rho_alpha = rho_beta = rho_total/2`. The energy path
///     (`hamiltonian_builder.cpp`) already halves it; this gradient assembly
///     fed the unhalved `P_total` straight into the density build, so
///     `vrho`/`vsigma` were evaluated at twice the density. Fixed by the
///     `xc_spin_scale` factor below — the dominant ~24 eV/Å error.
///
///  2. **GGA XC-potential factor of 2.** The GGA contribution to the XC
///     potential matrix in `xc_grid_engine.h` carried a spurious `* 2.0`,
///     so the SCF converged at a density where `Vxc != dE_xc/dP`. The
///     analytic nuclear gradient relies on the variational identity
///     `dE/dP = 0` (occ-virt) to drop the density-relaxation terms; with an
///     inconsistent `Vxc` that cancellation failed, leaving a residual
///     ~0.15 eV/Å on the KS gradient even after fix (1). Fixed in
///     `xc_grid_engine.h` (the `ga_d`/`gb_d` coefficients). LSDA was never
///     affected (no `vsigma`), which is how the GGA-only nature was isolated.
///
/// @param grad_out      [natm,3] gradient accumulator, units Hartree/Bohr.
/// @param integrals     Integral pack at the current geometry.
/// @param result        Converged SCF result (density, MO coefficients).
/// @param config        SCF configuration (spin, xc_type, grid sizes).
/// @param shell_to_atom Map from basis-function index to owning atom.
/// @param atom_shells   Per-atom list of owned basis-function indices.
/// @param natm          Number of atoms.
/// @param n             Number of basis functions.

static void accumulate_xc_gradient(
    xt::xtensor<double, 2>& grad_out,
    const IntegralPack& integrals,
    const SCFResult& result,
    const SCFConfig& config,
    const std::vector<int>& shell_to_atom,
    const std::vector<std::vector<int>>& atom_shells,
    int natm, int n)
{
    bool is_gga = (config.xc_type == XCType::PBE);
    const int N = static_cast<int>(n);
    const size_t Ns = static_cast<size_t>(N);

    auto grid = build_molecular_grid(integrals.atoms.data(), natm,
                                     config.grid_level,
                                     config.n_radial, config.n_angular);
    const int G = static_cast<int>(grid.points.size());
    const size_t Gs = static_cast<size_t>(G);

    // ── 1. Build Phi[G×N], dPhi[G×N×3], (Hess[G×N×9] for GGA) ─────
    xt::xtensor<double, 2> Phi = xt::zeros<double>({Gs, Ns});
    xt::xtensor<double, 3> dPhi = xt::zeros<double>({Gs, Ns, size_t(3)});
    xt::xtensor<double, 3> Hess;
    if (is_gga) Hess = xt::zeros<double>({Gs, Ns, size_t(9)});

    for (int mu = 0; mu < N; ++mu) {
        const auto& sh = integrals.shells[mu];
        double min_alpha = sh.alpha[0];
        for (int k = 1; k < sh.n_prim; ++k)
            if (sh.alpha[k] < min_alpha) min_alpha = sh.alpha[k];

        for (int g = 0; g < G; ++g) {
            double rx = static_cast<double>(grid.points[g].x);
            double ry = static_cast<double>(grid.points[g].y);
            double rz = static_cast<double>(grid.points[g].z);
            double dx = rx - sh.center[0], dy = ry - sh.center[1], dz = rz - sh.center[2];
            double r2 = dx*dx + dy*dy + dz*dz;
            if (min_alpha * r2 > 36.0) continue;

            if (is_gga) {
                double val, gx, gy, gz, hxx, hxy, hxz, hyy, hyz, hzz;
                eval_basis_all(sh, rx, ry, rz, val, gx, gy, gz,
                               hxx, hxy, hxz, hyy, hyz, hzz);
                Phi(g, mu) = val;
                dPhi(g, mu, 0) = gx; dPhi(g, mu, 1) = gy; dPhi(g, mu, 2) = gz;
                // Hess stored as flat [xx,xy,xz, yx,yy,yz, zx,zy,zz]
                Hess(g, mu, 0) = hxx; Hess(g, mu, 1) = hxy; Hess(g, mu, 2) = hxz;
                Hess(g, mu, 3) = hxy; Hess(g, mu, 4) = hyy; Hess(g, mu, 5) = hyz;
                Hess(g, mu, 6) = hxz; Hess(g, mu, 7) = hyz; Hess(g, mu, 8) = hzz;
            } else {
                double val, gx, gy, gz;
                eval_basis_val_grad(sh, rx, ry, rz, val, gx, gy, gz);
                Phi(g, mu) = val;
                dPhi(g, mu, 0) = gx; dPhi(g, mu, 1) = gy; dPhi(g, mu, 2) = gz;
            }
        }
    }

    // ── 2. Density via GEMM ────────────────────────────────────────────
    //
    // Spin-density convention (must match hamiltonian_builder.cpp's XC path):
    //   - restricted:  result.P_alpha == result.P_beta == P_total (factor 2
    //                  folded in). The XC functional is spin-resolved, so the
    //                  per-spin density blocks fed to it are P_total/2.
    //   - unrestricted: result.P_alpha / result.P_beta already carry one spin
    //                   each — no scaling.
    // Feeding the unhalved P_total into a restricted-KS XC gradient evaluates
    // vrho/vsigma at twice the correct density and doubles drho/dR — the
    // ~4-orders-of-magnitude defect this fix corrects.
    const bool restricted_ks = (config.spin == SpinType::Restricted);
    const double xc_spin_scale = restricted_ks ? 0.5 : 1.0;
    xt::xtensor<double, 2> Pa({Ns, Ns}), Pb({Ns, Ns});
    for (int mu = 0; mu < N; ++mu)
        for (int nu = 0; nu < N; ++nu) {
            Pa(mu, nu) = xc_spin_scale
                       * static_cast<double>(result.P_alpha(mu, nu));
            Pb(mu, nu) = xc_spin_scale
                       * static_cast<double>(result.P_beta(mu, nu));
        }

    auto T_a = xt::linalg::dot(Phi, Pa);  // [G × N]
    auto T_b = xt::linalg::dot(Phi, Pb);

    xt::xtensor<double, 1> rho_a_arr = xt::maximum(xt::sum(T_a * Phi, {1}), 0.0);
    xt::xtensor<double, 1> rho_b_arr = xt::maximum(xt::sum(T_b * Phi, {1}), 0.0);

    // Gradient density for GGA.
    xt::xtensor<double, 2> grad_rho_a, grad_rho_b;  // [G × 3]
    xt::xtensor<double, 1> sigma_aa_arr, sigma_ab_arr, sigma_bb_arr;
    if (is_gga) {
        grad_rho_a = xt::zeros<double>({Gs, size_t(3)});
        grad_rho_b = xt::zeros<double>({Gs, size_t(3)});
        for (int d = 0; d < 3; ++d) {
            auto dPhi_d = xt::view(dPhi, xt::all(), xt::all(), d);
            xt::view(grad_rho_a, xt::all(), d) = 2.0 * xt::sum(T_a * dPhi_d, {1});
            xt::view(grad_rho_b, xt::all(), d) = 2.0 * xt::sum(T_b * dPhi_d, {1});
        }
        sigma_aa_arr = xt::sum(grad_rho_a * grad_rho_a, {1});
        sigma_ab_arr = xt::sum(grad_rho_a * grad_rho_b, {1});
        sigma_bb_arr = xt::sum(grad_rho_b * grad_rho_b, {1});
    }

    // ── 3. Batch XC evaluation ─────────────────────────────────────────
    xt::xtensor<double, 1> exc_arr = xt::zeros<double>({Gs});
    xt::xtensor<double, 1> vrho_a_arr = xt::zeros<double>({Gs});
    xt::xtensor<double, 1> vrho_b_arr = xt::zeros<double>({Gs});
    xt::xtensor<double, 1> vsigma_aa_arr, vsigma_ab_arr, vsigma_bb_arr;
    if (is_gga) {
        vsigma_aa_arr = xt::zeros<double>({Gs});
        vsigma_ab_arr = xt::zeros<double>({Gs});
        vsigma_bb_arr = xt::zeros<double>({Gs});
    }

    evaluate_xc_batch(config.xc_type, G,
                      rho_a_arr.data(), rho_b_arr.data(),
                      is_gga ? sigma_aa_arr.data() : nullptr,
                      is_gga ? sigma_ab_arr.data() : nullptr,
                      is_gga ? sigma_bb_arr.data() : nullptr,
                      exc_arr.data(), vrho_a_arr.data(), vrho_b_arr.data(),
                      is_gga ? vsigma_aa_arr.data() : nullptr,
                      is_gga ? vsigma_ab_arr.data() : nullptr,
                      is_gga ? vsigma_bb_arr.data() : nullptr);

    // ── 4. Weights ─────────────────────────────────────────────────────
    xt::xtensor<double, 1> w_arr = xt::zeros<double>({Gs});
    for (int g = 0; g < G; ++g) w_arr(g) = static_cast<double>(grid.points[g].weight);

    // ── 5. XC gradient accumulation (per atom, per direction) ──────────
    //
    // This computes the explicit basis-function-derivative contribution to
    // dE_xc/dR_A on a grid held FIXED in space (Johnson-Gill-Pople
    // convention, JCP 98, 5612 (1993)): only the basis-function centers
    // move with their atoms. The Becke partition-weight / grid-coordinate
    // response is the small omitted piece (sub-mHa/Bohr on a converged
    // grid); the matching FD oracle (scf_finite_difference_gradient) pins
    // the grid likewise. The density-relaxation contribution cancels
    // because the SCF is variational and Vxc == dE_xc/dP (see the GGA Vxc
    // fix in xc_grid_engine.h).
    //
    // LDA term:
    //   dE_xc/dR_{A,d} = -2 * sum_g w * (vrho_a * sum_{mu in A} T_a(g,mu)*dPhi_d
    //                                    + vrho_b * sum_{mu in A} T_b(g,mu)*dPhi_d)
    // Pre-weight T by w*vrho: WT_a(g,mu) = w(g)*vrho_a(g)*T_a(g,mu).
    auto w_vrho_a_col = xt::reshape_view(w_arr * vrho_a_arr, {Gs, size_t(1)});
    auto w_vrho_b_col = xt::reshape_view(w_arr * vrho_b_arr, {Gs, size_t(1)});
    xt::xtensor<double, 2> WT_a = T_a * w_vrho_a_col;
    xt::xtensor<double, 2> WT_b = T_b * w_vrho_b_col;

    for (int A = 0; A < natm; ++A) {
        const auto& shells_A = atom_shells[A];
        if (shells_A.empty()) continue;

        for (int dd = 0; dd < 3; ++dd) {
            double xc_grad = 0;

            // LDA contribution via pre-weighted T
            for (int mu : shells_A) {
                for (int g = 0; g < G; ++g) {
                    double dphi_d = dPhi(g, mu, dd);
                    xc_grad -= 2.0 * (WT_a(g, mu) + WT_b(g, mu)) * dphi_d;
                }
            }

            // GGA contribution (per-point for correctness)
            if (is_gga) {
                for (int g = 0; g < G; ++g) {
                    double rho = rho_a_arr(g) + rho_b_arr(g);
                    if (rho < 1e-20) continue;
                    double w = w_arr(g);
                    double vsa = vsigma_aa_arr(g), vsab = vsigma_ab_arr(g), vsb = vsigma_bb_arr(g);

                    // grad_rho_a(g,i) = 2*sum_{mu,nu} Phi(g,mu)*Pa(mu,nu)*dPhi(g,nu,i)
                    // d(grad_rho_a(g,i))/dR_{A,d} has exactly 2 terms:
                    //   Term1 (dPhi_mu/dR_A): 2*sum_{mu in A, nu} Pa(mu,nu)*(-dPhi(g,mu,d))*dPhi(g,nu,i)
                    //   Term2 (dHess_nu/dR_A): 2*sum_{nu in A, mu} Pa(mu,nu)*Phi(g,mu)*(-Hess(g,nu,d*i))
                    double dgrad_rho_a[3] = {0,0,0}, dgrad_rho_b[3] = {0,0,0};

                    // Term1: derivative of the Phi(mu) factor when mu is on atom A
                    for (int mu : shells_A) {
                        double neg_dphi_mu_d = -dPhi(g, mu, dd);
                        for (int nu = 0; nu < N; ++nu) {
                            double pa = Pa(mu, nu), pb = Pb(mu, nu);
                            for (int i = 0; i < 3; ++i) {
                                dgrad_rho_a[i] += 2.0 * pa * neg_dphi_mu_d * dPhi(g, nu, i);
                                dgrad_rho_b[i] += 2.0 * pb * neg_dphi_mu_d * dPhi(g, nu, i);
                            }
                        }
                    }
                    // Term2: derivative of the dPhi(nu)/dr_i factor when nu is on atom A
                    for (int nu : shells_A) {
                        for (int mu = 0; mu < N; ++mu) {
                            double pa = Pa(mu, nu), pb = Pb(mu, nu);
                            double phi_mu = Phi(g, mu);
                            for (int i = 0; i < 3; ++i) {
                                double neg_h_nu_di = -Hess(g, nu, dd*3 + i);
                                dgrad_rho_a[i] += 2.0 * pa * phi_mu * neg_h_nu_di;
                                dgrad_rho_b[i] += 2.0 * pb * phi_mu * neg_h_nu_di;
                            }
                        }
                    }

                    double dsigma_aa = 0, dsigma_ab = 0, dsigma_bb = 0;
                    for (int i = 0; i < 3; ++i) {
                        dsigma_aa += 2.0 * grad_rho_a(g, i) * dgrad_rho_a[i];
                        dsigma_bb += 2.0 * grad_rho_b(g, i) * dgrad_rho_b[i];
                        dsigma_ab += grad_rho_a(g, i) * dgrad_rho_b[i]
                                   + grad_rho_b(g, i) * dgrad_rho_a[i];
                    }
                    xc_grad += w * (vsa * dsigma_aa + vsab * dsigma_ab + vsb * dsigma_bb);
                }
            }
            grad_out(A, dd) += static_cast<double>(xc_grad);
        }
    }
}

// ── Unified gradient ─────────────────────────────────────────────────

GradientResult scf_gradient(
    const Frame& frame,
    const IntegralPack& integrals,
    const SCFResult& result,
    const SCFConfig& config)
{
    const auto t_total = LogClock::now();
    int natm = static_cast<int>(frame.atoms.size());
    int n = integrals.n_basis;
    bool is_ks = (config.method == MeanFieldType::KS);
    bool restricted = (config.spin == SpinType::Restricted);
    const char* method_label = is_ks ? "KS" : "HF";
    const char* spin_label = restricted ? "restricted" : "unrestricted";

    ATV_INFO(
        "[gradient] start method={} spin={} natm={} basis={}",
        method_label, spin_label, natm, n);

    const auto t_dpack = LogClock::now();
    auto dpack = build_integral_derivatives(integrals);
    ATV_DEBUG("[gradient] derivative_pack_ms={:.1f}", elapsed_ms(t_dpack));

    // Restricted: P_alpha IS P_total (factor 2 included in scf_solve)
    // Unrestricted: P_total = P_alpha + P_beta (per-spin)
    xt::xtensor<double, 2> P = restricted
        ? xt::xtensor<double, 2>(result.P_alpha)
        : xt::xtensor<double, 2>(result.P_alpha + result.P_beta);

    xt::xtensor<double, 2> W = xt::zeros<double>({n, n});
    accumulate_W(W, result.C_alpha, result.eigenvalues_alpha, result.n_alpha, n);
    accumulate_W(W, result.C_beta, result.eigenvalues_beta, result.n_beta, n);

    xt::xtensor<double, 2> grad = xt::zeros<double>({(int)natm, (int)3});

    const auto t_common = LogClock::now();
    xt::xtensor<double, 2> grad_common = xt::zeros<double>({(int)natm, (int)3});
    accumulate_common_gradient(grad_common, P, W, dpack, n, natm);
    grad += grad_common;
    const auto [common_max_abs, common_max_atom_norm] = gradient_metrics(grad_common);
    ATV_DEBUG(
        "[gradient] common_terms_ms={:.1f}",
        elapsed_ms(t_common));
    ATV_INFO(
        "[gradient] common_terms max_abs_hb={:.6f} max_atom_norm_hb={:.6f}",
        common_max_abs, common_max_atom_norm);

    // Two-electron term: contracted on the fly to avoid materializing dERI.
    const auto t_two_e = LogClock::now();
    xt::xtensor<double, 2> grad_two_e = xt::zeros<double>({(int)natm, (int)3});
    accumulate_two_electron_gradient(grad_two_e, frame, integrals, P, result, config);
    grad += grad_two_e;
    const auto [two_e_max_abs, two_e_max_atom_norm] = gradient_metrics(grad_two_e);
    ATV_DEBUG(
        "[gradient] two_e_ms={:.1f}",
        elapsed_ms(t_two_e));
    ATV_INFO(
        "[gradient] two_e_terms max_abs_hb={:.6f} max_atom_norm_hb={:.6f}",
        two_e_max_abs, two_e_max_atom_norm);

    // XC grid contribution (KS only)
    if (is_ks) {
        const auto t_xc = LogClock::now();
        std::vector<int> shell_to_atom(n);
        std::vector<std::vector<int>> atom_shells(natm);
        {
            int idx = 0;
            for (int a = 0; a < natm; ++a) {
                for (int mu = idx; mu < n; ++mu) {
                    double dx = integrals.shells[mu].center[0] - integrals.atoms[a].pos[0];
                    double dy = integrals.shells[mu].center[1] - integrals.atoms[a].pos[1];
                    double dz = integrals.shells[mu].center[2] - integrals.atoms[a].pos[2];
                    if (dx*dx + dy*dy + dz*dz < 1e-20) {
                        shell_to_atom[mu] = a;
                        atom_shells[a].push_back(mu);
                    } else {
                        break;
                    }
                }
                idx += static_cast<int>(atom_shells[a].size());
            }
        }
        xt::xtensor<double, 2> grad_xc = xt::zeros<double>({(int)natm, (int)3});
        accumulate_xc_gradient(grad_xc, integrals, result, config,
                               shell_to_atom, atom_shells, natm, n);

        // Enforce translational invariance of the XC gradient.
        //
        // The exact dE_xc/dR is translation-invariant (a rigid shift of the
        // molecule leaves E_xc unchanged), so sum_A dE_xc/dR_A = 0. The
        // fixed-grid term computed above approximates dE_xc/dR with a
        // numerical-quadrature error whose net (sum-over-atoms) component
        // is a spurious uniform force on the whole molecule — the residual
        // of integrating grad_r f_xc, which the Becke partition-weight
        // response would otherwise cancel. Projecting out that net
        // component (subtracting the per-atom mean) removes the spurious
        // translation while leaving every *relative* (internal) force —
        // the physically meaningful part for geometry optimization and MD
        // — unchanged. This is the standard translational-invariance
        // enforcement for incomplete integration grids; on the SG-0-ish
        // (50,110) grid the removed residual is ~4e-4 Ha/Bohr.
        if (natm > 0) {
            for (int d = 0; d < 3; ++d) {
                double net = 0.0;
                for (int A = 0; A < natm; ++A) net += grad_xc(A, d);
                const double shift = net / static_cast<double>(natm);
                for (int A = 0; A < natm; ++A) grad_xc(A, d) -= shift;
            }
        }
        grad += grad_xc;
        const auto [xc_max_abs, xc_max_atom_norm] = gradient_metrics(grad_xc);
        ATV_DEBUG(
            "[gradient] xc_ms={:.1f}",
            elapsed_ms(t_xc));
        ATV_INFO(
            "[gradient] xc_terms max_abs_hb={:.6f} max_atom_norm_hb={:.6f}",
            xc_max_abs, xc_max_atom_norm);
    }

    auto out = make_result(grad, result.energy, natm);
    double max_force = 0.0;
    for (int A = 0; A < natm; ++A)
        for (int d = 0; d < 3; ++d)
            max_force = std::max(max_force, std::abs(out.forces(A, d)));

    ATV_INFO(
        "[gradient] done total_ms={:.1f} energy_eV={:.6f} max_force_component_eVA={:.6f}",
        elapsed_ms(t_total), out.energy, max_force);

    return out;
}

// ── Finite-difference gradient ───────────────────────────────────────

GradientResult scf_finite_difference_gradient(
    const Frame& frame,
    const SCFConfig& config,
    int charge,
    double h,
    BasisSet basis)
{
    int natm = static_cast<int>(frame.atoms.size());

    // The analytic KS XC gradient (accumulate_xc_gradient) is the
    // fixed-grid gradient (Johnson-Gill-Pople convention, JCP 98, 5612
    // (1993)): the integration grid is held fixed in space and only the
    // basis-function centers move with their atoms. For the FD oracle to
    // be a valid check of that analytic gradient, the +/- displaced
    // energies must SHARE one grid pinned to the reference geometry —
    // otherwise the moving-grid (Becke partition-weight) response, which
    // the analytic gradient deliberately omits, contaminates the FD
    // result. Build the reference grid once and pin it via
    // SCFConfig::fixed_grid for every displacement.
    SCFConfig fd_config = config;
    MolecularGrid fixed_grid;
    if (config.method == MeanFieldType::KS) {
        const auto ref_integrals = build_integrals(frame, basis, charge);
        if (!ref_integrals.atoms.empty()) {
            fixed_grid = build_molecular_grid(
                ref_integrals.atoms.data(),
                static_cast<int>(ref_integrals.atoms.size()),
                config.grid_level,
                config.n_radial, config.n_angular);
            fd_config.fixed_grid = &fixed_grid;
        }
    }

    auto energy_fn = [&](const Frame& f) {
        return scf_solve(build_integrals(f, basis, charge), fd_config).energy;
    };
    double ref_energy = energy_fn(frame);

    xt::xtensor<double, 2> forces = xt::zeros<double>({(int)natm, (int)3});

    for (int A = 0; A < natm; ++A) {
        for (int d = 0; d < 3; ++d) {
            Frame fwd = frame;
            if (d == 0) fwd.atoms.x(A) += h;
            else if (d == 1) fwd.atoms.y(A) += h;
            else fwd.atoms.z(A) += h;
            double e_fwd = energy_fn(fwd);

            Frame bwd = frame;
            if (d == 0) bwd.atoms.x(A) -= h;
            else if (d == 1) bwd.atoms.y(A) -= h;
            else bwd.atoms.z(A) -= h;
            double e_bwd = energy_fn(bwd);

            // FD uses Å step, so dE/dr is already in eV/Å once the Hartree
            // energy is converted — no extra 1/Bohr factor.
            forces(A, d) =
                -units::convert((e_fwd - e_bwd) / (2.0 * h),
                                units::Energy::Hartree, units::Energy::eV);
        }
    }

    GradientResult out;
    out.forces = forces;
    out.energy = units::convert(ref_energy, units::Energy::Hartree, units::Energy::eV);
    return out;
}

} // namespace scf
} // namespace cpu
} // namespace ATV
