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
#include "atomiverse/core/units.h"
#include "atomiverse/core/constants.h"

#include <xtensor/xview.hpp>
#include <xtensor/xio.hpp>
#include <cmath>
#include <functional>

namespace ATV {
namespace cpu {
namespace scf {

// ── Shared helpers ───────────────────────────────────────────────────

static void accumulate_W(
    xt::xtensor<Real, 2>& W,
    const xt::xtensor<Real, 2>& C,
    const xt::xtensor<Real, 1>& eigenvalues,
    Int n_occ, Int n_basis)
{
    for (Int i = 0; i < n_occ; ++i)
        for (Int mu = 0; mu < n_basis; ++mu)
            for (Int nu = 0; nu < n_basis; ++nu)
                W(mu, nu) += eigenvalues(i) * C(mu, i) * C(nu, i);
}

static void accumulate_common_gradient(
    xt::xtensor<Real, 2>& grad,
    const xt::xtensor<Real, 2>& P,
    const xt::xtensor<Real, 2>& W,
    const IntegralDerivPack& dpack,
    Int n_basis, int natm)
{
    for (int A = 0; A < natm; ++A) {
        for (int d = 0; d < 3; ++d) {
            Real val = 0;
            for (Int mu = 0; mu < n_basis; ++mu)
                for (Int nu = 0; nu < n_basis; ++nu) {
                    val += P(mu, nu) * (dpack.dT(A, d, mu, nu) + dpack.dV(A, d, mu, nu));
                    val -= W(mu, nu) * dpack.dS(A, d, mu, nu);
                }
            val += dpack.dV_nuc(A, d);
            grad(A, d) += val;
        }
    }
}

static GradientResult make_result(
    const xt::xtensor<Real, 2>& grad, Real energy_hartree, int natm)
{
    GradientResult out;
    out.forces = xt::zeros<Real>({(Int)natm, (Int)3});
    for (int A = 0; A < natm; ++A)
        for (int d = 0; d < 3; ++d)
            out.forces(A, d) = -grad(A, d) * kForceAU2Metal;
    out.energy = energy_hartree * kHartree2eV;
    return out;
}

// ── XC grid contribution (KS only) ──────────────────────────────────

static void accumulate_xc_gradient(
    xt::xtensor<Real, 2>& grad,
    const IntegralPack& integrals,
    const SCFResult& result,
    const SCFConfig& config,
    const std::vector<int>& shell_to_atom,
    const std::vector<std::vector<Int>>& atom_shells,
    int natm, Int n)
{
    bool is_gga = (config.xc_type == XCType::PBE);

    auto grid = build_molecular_grid(integrals.atoms.data(), natm,
                                     config.n_radial, config.n_angular);

    std::vector<double> Pa(n*n), Pb(n*n);
    for (Int mu = 0; mu < n; ++mu)
        for (Int nu = 0; nu < n; ++nu) {
            Pa[mu*n + nu] = static_cast<double>(result.P_alpha(mu, nu));
            Pb[mu*n + nu] = static_cast<double>(result.P_beta(mu, nu));
        }

    std::vector<double> phi(n);
    std::vector<double> dphi(n * 3);
    std::vector<double> hess(is_gga ? n * 9 : 0);

    for (const auto& gp : grid.points) {
        double rx = static_cast<double>(gp.x);
        double ry = static_cast<double>(gp.y);
        double rz = static_cast<double>(gp.z);
        double w  = static_cast<double>(gp.weight);

        bool any_nonzero = false;
        for (Int mu = 0; mu < n; ++mu) {
            if (is_gga) {
                eval_basis_all(integrals.shells[mu], rx, ry, rz,
                    phi[mu], dphi[mu*3], dphi[mu*3+1], dphi[mu*3+2],
                    hess[mu*9+0], hess[mu*9+1], hess[mu*9+2],
                    hess[mu*9+4], hess[mu*9+5], hess[mu*9+8]);
                hess[mu*9+3] = hess[mu*9+1];
                hess[mu*9+6] = hess[mu*9+2];
                hess[mu*9+7] = hess[mu*9+5];
            } else {
                eval_basis_val_grad(integrals.shells[mu], rx, ry, rz,
                    phi[mu], dphi[mu*3], dphi[mu*3+1], dphi[mu*3+2]);
            }
            if (std::abs(phi[mu]) > 1e-20) any_nonzero = true;
        }
        if (!any_nonzero) continue;

        double rho_a = 0, rho_b = 0;
        for (Int mu = 0; mu < n; ++mu)
            for (Int nu = 0; nu < n; ++nu) {
                double pp = phi[mu] * phi[nu];
                rho_a += Pa[mu*n + nu] * pp;
                rho_b += Pb[mu*n + nu] * pp;
            }
        rho_a = std::max(0.0, rho_a);
        rho_b = std::max(0.0, rho_b);
        if (rho_a + rho_b < 1e-20) continue;

        double grad_rho_a[3] = {0,0,0}, grad_rho_b[3] = {0,0,0};
        double sigma_aa = 0, sigma_ab = 0, sigma_bb = 0;
        if (is_gga) {
            for (Int mu = 0; mu < n; ++mu)
                for (Int nu = 0; nu < n; ++nu) {
                    double pa = Pa[mu*n + nu];
                    double pb = Pb[mu*n + nu];
                    for (int i = 0; i < 3; ++i) {
                        grad_rho_a[i] += pa * dphi[mu*3+i] * phi[nu];
                        grad_rho_b[i] += pb * dphi[mu*3+i] * phi[nu];
                    }
                }
            for (int i = 0; i < 3; ++i) { grad_rho_a[i] *= 2.0; grad_rho_b[i] *= 2.0; }
            sigma_aa = grad_rho_a[0]*grad_rho_a[0] + grad_rho_a[1]*grad_rho_a[1] + grad_rho_a[2]*grad_rho_a[2];
            sigma_ab = grad_rho_a[0]*grad_rho_b[0] + grad_rho_a[1]*grad_rho_b[1] + grad_rho_a[2]*grad_rho_b[2];
            sigma_bb = grad_rho_b[0]*grad_rho_b[0] + grad_rho_b[1]*grad_rho_b[1] + grad_rho_b[2]*grad_rho_b[2];
        }

        XCOutput xc = evaluate_xc(config.xc_type,
            static_cast<Real>(rho_a), static_cast<Real>(rho_b),
            static_cast<Real>(sigma_aa), static_cast<Real>(sigma_ab), static_cast<Real>(sigma_bb));
        double vrho_a = static_cast<double>(xc.vrho_alpha);
        double vrho_b = static_cast<double>(xc.vrho_beta);

        for (int A = 0; A < natm; ++A) {
            const auto& shells_A = atom_shells[A];
            if (shells_A.empty()) continue;

            for (int dd = 0; dd < 3; ++dd) {
                double xc_contrib = 0;

                double drho_a = 0, drho_b = 0;
                for (Int mu : shells_A) {
                    double dphi_d = dphi[mu*3 + dd];
                    for (Int nu = 0; nu < n; ++nu) {
                        drho_a += Pa[mu*n + nu] * dphi_d * phi[nu];
                        drho_b += Pb[mu*n + nu] * dphi_d * phi[nu];
                    }
                }
                drho_a *= -2.0;
                drho_b *= -2.0;
                xc_contrib += vrho_a * drho_a + vrho_b * drho_b;

                if (is_gga) {
                    double vsa = static_cast<double>(xc.vsigma_aa);
                    double vsab = static_cast<double>(xc.vsigma_ab);
                    double vsb = static_cast<double>(xc.vsigma_bb);

                    double dgrad_rho_a[3] = {0,0,0}, dgrad_rho_b[3] = {0,0,0};

                    for (Int mu : shells_A) {
                        double dphi_mu_d = dphi[mu*3 + dd];
                        for (Int nu = 0; nu < n; ++nu) {
                            double pa = Pa[mu*n + nu];
                            double pb = Pb[mu*n + nu];
                            for (int i = 0; i < 3; ++i) {
                                double h_mu_di = hess[mu*9 + dd*3 + i];
                                double contrib = (-h_mu_di) * phi[nu]
                                               + (-dphi_mu_d) * dphi[nu*3+i];
                                dgrad_rho_a[i] += pa * contrib;
                                dgrad_rho_b[i] += pb * contrib;
                            }
                        }
                    }
                    for (Int nu : shells_A) {
                        double dphi_nu_d = dphi[nu*3 + dd];
                        for (Int mu = 0; mu < n; ++mu) {
                            if (shell_to_atom[mu] == A) continue;
                            double pa = Pa[mu*n + nu];
                            double pb = Pb[mu*n + nu];
                            for (int i = 0; i < 3; ++i) {
                                double h_nu_di = hess[nu*9 + dd*3 + i];
                                double contrib = phi[mu] * (-h_nu_di)
                                               + dphi[mu*3+i] * (-dphi_nu_d);
                                dgrad_rho_a[i] += pa * contrib;
                                dgrad_rho_b[i] += pb * contrib;
                            }
                        }
                    }

                    double dsigma_aa = 0, dsigma_ab = 0, dsigma_bb = 0;
                    for (int i = 0; i < 3; ++i) {
                        dsigma_aa += 2.0 * grad_rho_a[i] * dgrad_rho_a[i];
                        dsigma_bb += 2.0 * grad_rho_b[i] * dgrad_rho_b[i];
                        dsigma_ab += grad_rho_a[i] * dgrad_rho_b[i] + grad_rho_b[i] * dgrad_rho_a[i];
                    }

                    xc_contrib += vsa * dsigma_aa + vsab * dsigma_ab + vsb * dsigma_bb;
                }

                grad(A, dd) += static_cast<Real>(w * xc_contrib);
            }
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
    int natm = static_cast<int>(frame.atoms.size());
    Int n = integrals.n_basis;
    bool is_ks = (config.method == MeanFieldType::KS);
    bool restricted = (config.spin == SpinType::Restricted);

    auto dpack = build_integral_derivatives_sto3g(frame);

    // Restricted: P_alpha IS P_total (factor 2 included in scf_solve)
    // Unrestricted: P_total = P_alpha + P_beta (per-spin)
    xt::xtensor<Real, 2> P = restricted
        ? xt::xtensor<Real, 2>(result.P_alpha)
        : xt::xtensor<Real, 2>(result.P_alpha + result.P_beta);

    xt::xtensor<Real, 2> W = xt::zeros<Real>({n, n});
    accumulate_W(W, result.C_alpha, result.eigenvalues_alpha, result.n_alpha, n);
    accumulate_W(W, result.C_beta, result.eigenvalues_beta, result.n_beta, n);

    xt::xtensor<Real, 2> grad = xt::zeros<Real>({(Int)natm, (Int)3});

    accumulate_common_gradient(grad, P, W, dpack, n, natm);

    // Two-electron term: dispatch based on method + spin
    // RHF uses a fused Coulomb-exchange formula with matching index patterns;
    // UHF separates Coulomb and per-spin exchange with different index patterns.
    // These are NOT equivalent for restricted P_alpha = P_beta = P/2.
    bool is_restricted_hf = (!is_ks && config.spin == SpinType::Restricted);

    for (int A = 0; A < natm; ++A)
        for (int d = 0; d < 3; ++d)
            for (Int mu = 0; mu < n; ++mu)
                for (Int nu = 0; nu < n; ++nu)
                    for (Int k = 0; k < n; ++k)
                        for (Int l = 0; l < n; ++l) {
                            if (is_restricted_hf) {
                                // RHF fused: 0.5*P*P*[dERI_coulomb - 0.5*dERI_exchange]
                                grad(A, d) += 0.5 * P(mu, nu) * P(k, l)
                                     * (dpack.dERI(A, d, mu, nu, k, l)
                                        - 0.5 * dpack.dERI(A, d, mu, k, nu, l));
                            } else if (is_ks) {
                                // KS: Coulomb only
                                grad(A, d) += 0.5 * P(mu, nu) * P(k, l)
                                     * dpack.dERI(A, d, mu, nu, k, l);
                            } else {
                                // UHF: separate Coulomb + per-spin exchange
                                grad(A, d) += 0.5 * P(mu, nu) * P(k, l)
                                     * dpack.dERI(A, d, mu, nu, k, l);
                                grad(A, d) -= 0.5 * (result.P_alpha(mu, k) * result.P_alpha(nu, l)
                                            + result.P_beta(mu, k) * result.P_beta(nu, l))
                                     * dpack.dERI(A, d, mu, k, nu, l);
                            }
                        }

    // XC grid contribution (KS only)
    if (is_ks) {
        std::vector<int> shell_to_atom(n);
        std::vector<std::vector<Int>> atom_shells(natm);
        {
            int idx = 0;
            for (int a = 0; a < natm; ++a) {
                for (Int mu = idx; mu < n; ++mu) {
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
        accumulate_xc_gradient(grad, integrals, result, config,
                               shell_to_atom, atom_shells, natm, n);
    }

    return make_result(grad, result.energy, natm);
}

// ── Finite-difference gradient ───────────────────────────────────────

GradientResult scf_finite_difference_gradient(
    const Frame& frame,
    const SCFConfig& config,
    int charge,
    Real h)
{
    int natm = static_cast<int>(frame.atoms.size());

    auto energy_fn = [&](const Frame& f) {
        return scf_solve(build_integrals_sto3g(f, charge), config).energy;
    };
    Real ref_energy = energy_fn(frame);

    xt::xtensor<Real, 2> forces = xt::zeros<Real>({(Int)natm, (Int)3});

    for (int A = 0; A < natm; ++A) {
        for (int d = 0; d < 3; ++d) {
            Frame fwd = frame;
            if (d == 0) fwd.atoms.x(A) += h;
            else if (d == 1) fwd.atoms.y(A) += h;
            else fwd.atoms.z(A) += h;
            Real e_fwd = energy_fn(fwd);

            Frame bwd = frame;
            if (d == 0) bwd.atoms.x(A) -= h;
            else if (d == 1) bwd.atoms.y(A) -= h;
            else bwd.atoms.z(A) -= h;
            Real e_bwd = energy_fn(bwd);

            forces(A, d) = -(e_fwd - e_bwd) / (2.0 * h) * kHartree2eV;
        }
    }

    GradientResult out;
    out.forces = forces;
    out.energy = ref_energy * kHartree2eV;
    return out;
}

} // namespace scf
} // namespace cpu
} // namespace ATV
