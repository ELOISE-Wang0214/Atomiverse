/// @file cpu/scf/gradient.cpp
/// @brief Analytic nuclear gradients for RHF and UHF.
///
/// RHF gradient:
///   dE/dR_A = sum_{mu,nu} P_{mu,nu} * dHcore_{mu,nu}/dR_A
///           + (1/2) sum_{mu,nu,kl} P_{mu,nu} * P_{kl}
///             * [d(mu nu|kl)/dR_A - (1/2) d(mu k|nu l)/dR_A]
///           - sum_{mu,nu} W_{mu,nu} * dS_{mu,nu}/dR_A
///           + dV_nuc/dR_A
///
/// where W_{mu,nu} = sum_i^{occ} epsilon_i * C_{mu,i} * C_{nu,i}
/// is the energy-weighted density matrix.

#include "atomiverse/cpu/scf/gradient.h"
#include "atomiverse/cpu/scf/integral_derivatives.h"
#include "atomiverse/core/units.h"
#include "atomiverse/core/constants.h"

#include <xtensor/xview.hpp>
#include <xtensor/xio.hpp>
#include <cmath>

namespace ATV {
namespace cpu {
namespace scf {

// ── RHF gradient ─────────────────────────────────────────────────────

GradientResult rhf_gradient(
    const Frame& frame,
    const IntegralPack& integrals,
    const RHFResult& result)
{
    int natm = static_cast<int>(frame.n_atoms());
    Int n = integrals.n_basis;
    Int n_occ = integrals.n_occ();

    // Build integral derivatives
    auto dpack = build_integral_derivatives_sto3g(frame);

    // Energy-weighted density matrix: W = 2 * C_occ * diag(epsilon_occ) * C_occ^T
    // (factor 2 for RHF closed-shell double occupancy)
    xt::xtensor<Real, 2> W = xt::zeros<Real>({n, n});
    for (Int i = 0; i < n_occ; ++i)
        for (Int mu = 0; mu < n; ++mu)
            for (Int nu = 0; nu < n; ++nu)
                W(mu, nu) += 2.0 * result.eigenvalues(i) * result.C(mu, i) * result.C(nu, i);

    const auto& P = result.P;

    // Gradient in atomic units (Hartree/Bohr)
    xt::xtensor<Real, 2> grad = xt::zeros<Real>({(Int)natm, (Int)3});

    for (int A = 0; A < natm; ++A) {
        for (int d = 0; d < 3; ++d) {
            Real val = 0;

            // (1) One-electron: sum_{mu,nu} P_{mu,nu} * dHcore_{mu,nu}/dR_A
            //     dHcore = dT + dV
            for (Int mu = 0; mu < n; ++mu)
                for (Int nu = 0; nu < n; ++nu)
                    val += P(mu, nu) * (dpack.dT(A, d, mu, nu) + dpack.dV(A, d, mu, nu));

            // (2) Two-electron: (1/2) sum_{mu,nu,kl} P_{mu,nu} * P_{kl}
            //     * [d(mu nu|kl)/dR_A - (1/2) d(mu k|nu l)/dR_A]
            for (Int mu = 0; mu < n; ++mu)
                for (Int nu = 0; nu < n; ++nu)
                    for (Int k = 0; k < n; ++k)
                        for (Int l = 0; l < n; ++l)
                            val += 0.5 * P(mu, nu) * P(k, l)
                                 * (dpack.dERI(A, d, mu, nu, k, l)
                                    - 0.5 * dpack.dERI(A, d, mu, k, nu, l));

            // (3) Pulay force: -sum_{mu,nu} W_{mu,nu} * dS_{mu,nu}/dR_A
            for (Int mu = 0; mu < n; ++mu)
                for (Int nu = 0; nu < n; ++nu)
                    val -= W(mu, nu) * dpack.dS(A, d, mu, nu);

            // (4) Nuclear repulsion
            val += dpack.dV_nuc(A, d);

            grad(A, d) = val;
        }
    }

    // Convert to forces in eV/Angstrom: F = -grad, with unit conversion
    GradientResult out;
    out.forces = xt::zeros<Real>({(Int)natm, (Int)3});
    for (int A = 0; A < natm; ++A)
        for (int d = 0; d < 3; ++d)
            out.forces(A, d) = -grad(A, d) * kForceAU2Metal;

    out.energy = result.energy * kHartree2eV;
    return out;
}

// ── UHF gradient ─────────────────────────────────────────────────────

GradientResult uhf_gradient(
    const Frame& frame,
    const IntegralPack& integrals,
    const UHFResult& result)
{
    int natm = static_cast<int>(frame.n_atoms());
    Int n = integrals.n_basis;

    auto dpack = build_integral_derivatives_sto3g(frame);

    // Total density matrix P = P_alpha + P_beta
    xt::xtensor<Real, 2> P = result.P_alpha + result.P_beta;

    // Energy-weighted density matrices (no factor 2 -- UHF single occupancy)
    // W_alpha = C_alpha_occ * diag(epsilon_alpha_occ) * C_alpha_occ^T
    // W = W_alpha + W_beta
    xt::xtensor<Real, 2> W = xt::zeros<Real>({n, n});
    for (Int i = 0; i < result.n_alpha; ++i)
        for (Int mu = 0; mu < n; ++mu)
            for (Int nu = 0; nu < n; ++nu)
                W(mu, nu) += result.eigenvalues_alpha(i)
                           * result.C_alpha(mu, i) * result.C_alpha(nu, i);
    for (Int i = 0; i < result.n_beta; ++i)
        for (Int mu = 0; mu < n; ++mu)
            for (Int nu = 0; nu < n; ++nu)
                W(mu, nu) += result.eigenvalues_beta(i)
                           * result.C_beta(mu, i) * result.C_beta(nu, i);

    xt::xtensor<Real, 2> grad = xt::zeros<Real>({(Int)natm, (Int)3});

    for (int A = 0; A < natm; ++A) {
        for (int d = 0; d < 3; ++d) {
            Real val = 0;

            // (1) One-electron: P_total * dHcore
            for (Int mu = 0; mu < n; ++mu)
                for (Int nu = 0; nu < n; ++nu)
                    val += P(mu, nu) * (dpack.dT(A, d, mu, nu) + dpack.dV(A, d, mu, nu));

            // (2) Two-electron: UHF form
            // Coulomb: (1/2) P_total * P_total * d(mu nu|kl)/dR_A
            // Exchange: -(1/2) [P_alpha * P_alpha + P_beta * P_beta] * d(mu k|nu l)/dR_A
            for (Int mu = 0; mu < n; ++mu)
                for (Int nu = 0; nu < n; ++nu)
                    for (Int k = 0; k < n; ++k)
                        for (Int l = 0; l < n; ++l) {
                            val += 0.5 * P(mu, nu) * P(k, l)
                                 * dpack.dERI(A, d, mu, nu, k, l);
                            val -= 0.5 * (result.P_alpha(mu, k) * result.P_alpha(nu, l)
                                        + result.P_beta(mu, k) * result.P_beta(nu, l))
                                 * dpack.dERI(A, d, mu, k, nu, l);
                        }

            // (3) Pulay force
            for (Int mu = 0; mu < n; ++mu)
                for (Int nu = 0; nu < n; ++nu)
                    val -= W(mu, nu) * dpack.dS(A, d, mu, nu);

            // (4) Nuclear repulsion
            val += dpack.dV_nuc(A, d);

            grad(A, d) = val;
        }
    }

    GradientResult out;
    out.forces = xt::zeros<Real>({(Int)natm, (Int)3});
    for (int A = 0; A < natm; ++A)
        for (int d = 0; d < 3; ++d)
            out.forces(A, d) = -grad(A, d) * kForceAU2Metal;

    out.energy = result.energy * kHartree2eV;
    return out;
}

// ── Finite-difference gradients ──────────────────────────────────────

GradientResult rhf_finite_difference_gradient(
    const Frame& frame,
    int charge,
    Real scf_tol,
    int max_iter,
    Real h)
{
    int natm = static_cast<int>(frame.n_atoms());

    // Reference energy at current geometry
    auto integrals = build_integrals_sto3g(frame, charge);
    auto ref = rhf_solve(integrals, scf_tol, max_iter);

    xt::xtensor<Real, 2> forces = xt::zeros<Real>({(Int)natm, (Int)3});

    // Central differences: F_A,d = -[E(R+h) - E(R-h)] / (2h)
    for (int A = 0; A < natm; ++A) {
        for (int d = 0; d < 3; ++d) {
            // Displace +h
            Frame fwd = frame;
            if (d == 0) fwd.x(A) += h;
            else if (d == 1) fwd.y(A) += h;
            else fwd.z(A) += h;

            auto int_fwd = build_integrals_sto3g(fwd, charge);
            auto res_fwd = rhf_solve(int_fwd, scf_tol, max_iter);

            // Displace -h
            Frame bwd = frame;
            if (d == 0) bwd.x(A) -= h;
            else if (d == 1) bwd.y(A) -= h;
            else bwd.z(A) -= h;

            auto int_bwd = build_integrals_sto3g(bwd, charge);
            auto res_bwd = rhf_solve(int_bwd, scf_tol, max_iter);

            // F = -dE/dR, convert Hartree -> eV, Angstrom already
            Real dE_dR = (res_fwd.energy - res_bwd.energy) / (2.0 * h);
            forces(A, d) = -dE_dR * kHartree2eV;
        }
    }

    GradientResult out;
    out.forces = forces;
    out.energy = ref.energy * kHartree2eV;
    return out;
}

GradientResult uhf_finite_difference_gradient(
    const Frame& frame,
    int n_alpha,
    int n_beta,
    int charge,
    Real scf_tol,
    int max_iter,
    Real h)
{
    int natm = static_cast<int>(frame.n_atoms());

    auto integrals = build_integrals_sto3g(frame, charge);
    auto ref = uhf_solve(integrals, n_alpha, n_beta, scf_tol, max_iter);

    xt::xtensor<Real, 2> forces = xt::zeros<Real>({(Int)natm, (Int)3});

    for (int A = 0; A < natm; ++A) {
        for (int d = 0; d < 3; ++d) {
            Frame fwd = frame;
            if (d == 0) fwd.x(A) += h;
            else if (d == 1) fwd.y(A) += h;
            else fwd.z(A) += h;
            auto int_fwd = build_integrals_sto3g(fwd, charge);
            auto res_fwd = uhf_solve(int_fwd, n_alpha, n_beta, scf_tol, max_iter);

            Frame bwd = frame;
            if (d == 0) bwd.x(A) -= h;
            else if (d == 1) bwd.y(A) -= h;
            else bwd.z(A) -= h;
            auto int_bwd = build_integrals_sto3g(bwd, charge);
            auto res_bwd = uhf_solve(int_bwd, n_alpha, n_beta, scf_tol, max_iter);

            Real dE_dR = (res_fwd.energy - res_bwd.energy) / (2.0 * h);
            forces(A, d) = -dE_dR * kHartree2eV;
        }
    }

    GradientResult out;
    out.forces = forces;
    out.energy = ref.energy * kHartree2eV;
    return out;
}

} // namespace scf
} // namespace cpu
} // namespace ATV
