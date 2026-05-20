/// @file cpu/scf/xc_grid_engine.h
/// @brief Batched XC grid integration engine — GEMM-based density/Vxc assembly.
///
/// Replaces per-point XC integration with matrix operations:
///   1. Build Phi[n_grid × n_basis] — batched basis function evaluation
///   2. rho = rowsum((Phi @ P) * Phi) — density via GEMM
///   3. XC batch evaluation — vectorized over grid points
///   4. Vxc = Phi_w^T @ Phi_w — potential matrix via GEMM
///
/// Shared by scf_solve.cpp (Fock build) and gradient.cpp (XC gradient).
///
/// @see PySCF numint.py for the reference algorithm.
#pragma once

#include <vector>

#include <xtensor/containers/xtensor.hpp>
#include <xtensor/views/xview.hpp>
#include <xtensor-blas/xlinalg.hpp>

#include "atomiverse/cpu/scf/grid.h"
#include "atomiverse/cpu/scf/native_integrals.h"
#include "atomiverse/cpu/scf/xc_functional.h"
#include "atomiverse/cpu/scf/basis_eval.h"

namespace ATV {
namespace cpu {
namespace scf {

/// @brief Result of batched XC grid integration.
struct XCGridResult {
    xt::xtensor<double, 2> Vxc_alpha;   ///< XC potential matrix (alpha) [n_basis × n_basis].
    xt::xtensor<double, 2> Vxc_beta;    ///< XC potential matrix (beta) [n_basis × n_basis].
    double E_xc;                         ///< Total XC energy [Hartree].
    double n_alpha_integrated;           ///< Grid-integrated alpha electron count.
    double n_beta_integrated;            ///< Grid-integrated beta electron count.
};

/// @brief Batched XC grid integration engine.
///
/// Builds AO value matrix Phi and gradient matrices dPhi once per call,
/// then uses BLAS GEMM for density contraction and Vxc assembly.
///
/// @param grid    Molecular quadrature grid.
/// @param shells  Basis function shell descriptors [n_basis].
/// @param n_basis Number of basis functions.
/// @param P_alpha Alpha density matrix [n_basis × n_basis].
/// @param P_beta  Beta density matrix [n_basis × n_basis].
/// @param xc_type XC functional type (LSDA or PBE).
/// @return XCGridResult with Vxc matrices, E_xc, and integration diagnostics.
inline XCGridResult evaluate_xc_batched(
    const MolecularGrid& grid,
    const std::vector<native::Shell>& shells,
    int n_basis,
    const xt::xtensor<double, 2>& P_alpha,
    const xt::xtensor<double, 2>& P_beta,
    XCType xc_type)
{
    const int G = static_cast<int>(grid.points.size());
    const int N = static_cast<int>(n_basis);
    const bool is_gga = (xc_type == XCType::PBE);
    const size_t Gs = static_cast<size_t>(G);
    const size_t Ns = static_cast<size_t>(N);

    // ── 1. Build Phi[G × N] (and dPhi_d[G × N] for GGA) ──────────────
    // Row-major for element-wise ops; will convert for BLAS as needed.
    xt::xtensor<double, 2> Phi = xt::zeros<double>({Gs, Ns});
    xt::xtensor<double, 2> dPhi_x, dPhi_y, dPhi_z;
    if (is_gga) {
        dPhi_x = xt::zeros<double>({Gs, Ns});
        dPhi_y = xt::zeros<double>({Gs, Ns});
        dPhi_z = xt::zeros<double>({Gs, Ns});
    }

    // Screen: track which grid points have nonzero basis functions.
    std::vector<bool> active(G, false);

    // Loop: outer over shells (cache-friendly for exp), inner over grid points.
    // Delegate the per-point polynomial logic to basis_eval.h so d-shells
    // (and higher l) inherit the same code path as the cube writers.
    for (int mu = 0; mu < N; ++mu) {
        const auto& sh = shells[mu];
        for (int g = 0; g < G; ++g) {
            const double rx = static_cast<double>(grid.points[g].x);
            const double ry = static_cast<double>(grid.points[g].y);
            const double rz = static_cast<double>(grid.points[g].z);
            if (is_gga) {
                double v, gx, gy, gz;
                eval_basis_val_grad(sh, rx, ry, rz, v, gx, gy, gz);
                Phi(g, mu) = v;
                dPhi_x(g, mu) = gx;
                dPhi_y(g, mu) = gy;
                dPhi_z(g, mu) = gz;
                if (std::abs(v) > 1e-20) active[g] = true;
            } else {
                const double v = eval_basis(sh, rx, ry, rz);
                Phi(g, mu) = v;
                if (std::abs(v) > 1e-20) active[g] = true;
            }
        }
    }

    // ── 2. Density via GEMM: T = Phi @ P, rho = rowsum(T * Phi) ───────
    // Convert P to double for computation
    xt::xtensor<double, 2> Pa({Ns, Ns}), Pb({Ns, Ns});
    for (int mu = 0; mu < N; ++mu)
        for (int nu = 0; nu < N; ++nu) {
            Pa(mu, nu) = static_cast<double>(P_alpha(mu, nu));
            Pb(mu, nu) = static_cast<double>(P_beta(mu, nu));
        }

    // T_a = Phi @ Pa  [G × N]  — GEMM
    auto T_a = xt::linalg::dot(Phi, Pa);
    auto T_b = xt::linalg::dot(Phi, Pb);

    // rho_a[g] = sum_nu T_a(g,nu) * Phi(g,nu)  — element-wise multiply + row sum
    xt::xtensor<double, 1> rho_a_arr = xt::sum(T_a * Phi, {1});
    xt::xtensor<double, 1> rho_b_arr = xt::sum(T_b * Phi, {1});

    // Clamp densities >= 0
    rho_a_arr = xt::maximum(rho_a_arr, 0.0);
    rho_b_arr = xt::maximum(rho_b_arr, 0.0);

    // ── 2b. Gradient density for GGA ───────────────────────────────────
    // grad_rho_a[d] = 2 * sum_nu T_a(g,nu) * dPhi_d(g,nu)
    xt::xtensor<double, 1> sigma_aa_arr, sigma_ab_arr, sigma_bb_arr;
    xt::xtensor<double, 2> grad_rho_ax, grad_rho_ay, grad_rho_az;
    xt::xtensor<double, 2> grad_rho_bx, grad_rho_by, grad_rho_bz;
    if (is_gga) {
        // grad_rho_a_d[g] = 2 * sum_nu T_a(g,nu) * dPhi_d(g,nu)
        auto gax = 2.0 * xt::sum(T_a * dPhi_x, {1});
        auto gay = 2.0 * xt::sum(T_a * dPhi_y, {1});
        auto gaz = 2.0 * xt::sum(T_a * dPhi_z, {1});
        auto gbx = 2.0 * xt::sum(T_b * dPhi_x, {1});
        auto gby = 2.0 * xt::sum(T_b * dPhi_y, {1});
        auto gbz = 2.0 * xt::sum(T_b * dPhi_z, {1});

        sigma_aa_arr = gax*gax + gay*gay + gaz*gaz;
        sigma_ab_arr = gax*gbx + gay*gby + gaz*gbz;
        sigma_bb_arr = gbx*gbx + gby*gby + gbz*gbz;

        // Store for GGA Vxc assembly
        grad_rho_ax = xt::xtensor<double, 2>(xt::reshape_view(gax, {Gs, size_t(1)}));
        grad_rho_ay = xt::xtensor<double, 2>(xt::reshape_view(gay, {Gs, size_t(1)}));
        grad_rho_az = xt::xtensor<double, 2>(xt::reshape_view(gaz, {Gs, size_t(1)}));
        grad_rho_bx = xt::xtensor<double, 2>(xt::reshape_view(gbx, {Gs, size_t(1)}));
        grad_rho_by = xt::xtensor<double, 2>(xt::reshape_view(gby, {Gs, size_t(1)}));
        grad_rho_bz = xt::xtensor<double, 2>(xt::reshape_view(gbz, {Gs, size_t(1)}));
    }

    // ── 3. Batch XC evaluation ─────────────────────────────────────────
    xt::xtensor<double, 1> exc_arr = xt::zeros<double>({Gs});
    xt::xtensor<double, 1> vrho_a_arr = xt::zeros<double>({Gs});
    xt::xtensor<double, 1> vrho_b_arr = xt::zeros<double>({Gs});
    xt::xtensor<double, 1> vsigma_aa_out, vsigma_ab_out, vsigma_bb_out;
    if (is_gga) {
        vsigma_aa_out = xt::zeros<double>({Gs});
        vsigma_ab_out = xt::zeros<double>({Gs});
        vsigma_bb_out = xt::zeros<double>({Gs});
    }

    evaluate_xc_batch(xc_type, G,
                      rho_a_arr.data(), rho_b_arr.data(),
                      is_gga ? sigma_aa_arr.data() : nullptr,
                      is_gga ? sigma_ab_arr.data() : nullptr,
                      is_gga ? sigma_bb_arr.data() : nullptr,
                      exc_arr.data(), vrho_a_arr.data(), vrho_b_arr.data(),
                      is_gga ? vsigma_aa_out.data() : nullptr,
                      is_gga ? vsigma_ab_out.data() : nullptr,
                      is_gga ? vsigma_bb_out.data() : nullptr);

    // ── 4. Integrate E_xc and electron counts ──────────────────────────
    xt::xtensor<double, 1> w_arr = xt::zeros<double>({Gs});
    for (int g = 0; g < G; ++g)
        w_arr(g) = static_cast<double>(grid.points[g].weight);

    double E_xc = xt::sum(w_arr * exc_arr * (rho_a_arr + rho_b_arr))();
    double n_alpha_int = xt::sum(w_arr * rho_a_arr)();
    double n_beta_int  = xt::sum(w_arr * rho_b_arr)();

    // ── 5. Vxc via GEMM ───────────────────────────────────────────────
    // LDA: Vxc = Phi_w^T @ Phi  where Phi_w(g,mu) = w(g)*vrho(g)*Phi(g,mu)
    // This gives: Vxc(mu,nu) = sum_g w(g) * vrho(g) * Phi(g,mu) * Phi(g,nu)
    //
    // For GGA, add gradient terms:
    //   Vxc += sum_d (Phi_gd_a)^T @ dPhi_d + dPhi_d^T @ (Phi_gd_a)
    // where Phi_gd_a(g,mu) = w(g) * ga_d(g) * Phi(g,mu)
    //   and ga_d = 2*(2*vsigma_aa*grad_rho_a_d + vsigma_ab*grad_rho_b_d)

    // LDA contribution (common to both LDA and GGA)
    // Phi_wa(g, mu) = w(g) * vrho_a(g) * Phi(g, mu)
    auto w_vrho_a = xt::reshape_view(w_arr * vrho_a_arr, {Gs, size_t(1)});
    auto w_vrho_b = xt::reshape_view(w_arr * vrho_b_arr, {Gs, size_t(1)});
    xt::xtensor<double, 2> Phi_wa = Phi * w_vrho_a;
    xt::xtensor<double, 2> Phi_wb = Phi * w_vrho_b;

    // Vxc_alpha = Phi_wa^T @ Phi  [N × N]
    auto Vxc_a_d = xt::linalg::dot(xt::transpose(Phi_wa), Phi);
    auto Vxc_b_d = xt::linalg::dot(xt::transpose(Phi_wb), Phi);

    if (is_gga) {
        // GGA contribution to the XC potential matrix.
        //
        //   Vxc_a^sigma(mu,nu) = sum_g sum_d ga_d(g)
        //                            * [Phi(mu)dPhi(nu,d) + dPhi(mu,d)Phi(nu)]
        //   ga_d = w * (2*vsigma_aa*grad_rho_a_d + vsigma_ab*grad_rho_b_d)
        //   gb_d = w * (2*vsigma_bb*grad_rho_b_d + vsigma_ab*grad_rho_a_d)
        //
        // E_xc = sum_g w f(sigma), sigma_aa = grad_rho_a . grad_rho_a, with the
        // P-symmetric form grad_rho_a_d = sum_kl P_a(kl)[dPhi(k,d)Phi(l) +
        // Phi(k)dPhi(l,d)]. Then dE/dP_a(mu,nu) = sum_g w [2*vsigma_aa*gra_d
        // + vsigma_ab*grb_d] * [dPhi(mu,d)Phi(nu)+Phi(mu)dPhi(nu,d)]. The
        // symmetric [Phi*dPhi + dPhi*Phi] pair is already supplied by the two
        // GEMMs below, so ga_d carries NO extra factor of 2. An earlier
        // spurious `* 2.0` here doubled the GGA potential, leaving the SCF
        // converged at a density where Vxc != dE_xc/dP; that broke the
        // density-relaxation cancellation in the analytic KS nuclear gradient
        // (spec phase1-productization-02-ks-dshell-gradient).
        auto vsa = xt::reshape_view(vsigma_aa_out, {Gs, size_t(1)});
        auto vsab = xt::reshape_view(vsigma_ab_out, {Gs, size_t(1)});
        auto vsb = xt::reshape_view(vsigma_bb_out, {Gs, size_t(1)});
        auto w_col = xt::reshape_view(w_arr, {Gs, size_t(1)});

        // For each Cartesian direction d:
        //   Vxc_a += w * (ga_d * Phi)^T @ dPhi_d + w * (ga_d * dPhi_d)^T @ Phi
        // which is: Vxc_a += (Phi * w*ga_d)^T @ dPhi_d + (dPhi_d * w*ga_d)^T @ Phi
        // Symmetrized: Vxc_a += (Phi * w*ga_d)^T @ dPhi_d + dPhi_d^T @ (Phi * w*ga_d)
        // But Phi * w*ga_d transposed-dotted with dPhi_d gives a non-symmetric matrix.
        // The correct formula for LDA+GGA combined Vxc(mu,nu) = sum_g w *
        //   [vrho * phi_mu * phi_nu + sum_d ga_d * (phi_mu * dphi_d_nu + dphi_d_mu * phi_nu)]
        // = LDA part + sum_d [ga_d * Phi]^T @ dPhi_d + sum_d [ga_d * dPhi_d]^T @ Phi

        const xt::xtensor<double, 2>* grad_rho_a[3] = {&grad_rho_ax, &grad_rho_ay, &grad_rho_az};
        const xt::xtensor<double, 2>* grad_rho_b[3] = {&grad_rho_bx, &grad_rho_by, &grad_rho_bz};
        const xt::xtensor<double, 2>* dPhi_arr[3] = {&dPhi_x, &dPhi_y, &dPhi_z};

        for (int d = 0; d < 3; ++d) {
            auto ga_d = w_col * (2.0 * vsa * (*grad_rho_a[d]) + vsab * (*grad_rho_b[d]));
            auto gb_d = w_col * (2.0 * vsb * (*grad_rho_b[d]) + vsab * (*grad_rho_a[d]));
            const auto& dPhi_d = *dPhi_arr[d];

            // Vxc_a += (Phi * ga_d)^T @ dPhi_d + dPhi_d^T @ (Phi * ga_d)
            auto Phi_ga = Phi * ga_d;
            Vxc_a_d += xt::linalg::dot(xt::transpose(Phi_ga), dPhi_d);
            Vxc_a_d += xt::linalg::dot(xt::transpose(dPhi_d), Phi_ga);

            auto Phi_gb = Phi * gb_d;
            Vxc_b_d += xt::linalg::dot(xt::transpose(Phi_gb), dPhi_d);
            Vxc_b_d += xt::linalg::dot(xt::transpose(dPhi_d), Phi_gb);
        }
    }

    // Symmetrize (GEMM may produce slight asymmetry from rounding)
    Vxc_a_d = 0.5 * (Vxc_a_d + xt::transpose(Vxc_a_d));
    Vxc_b_d = 0.5 * (Vxc_b_d + xt::transpose(Vxc_b_d));

    // Convert to double
    XCGridResult result;
    result.Vxc_alpha = xt::xtensor<double, 2>::from_shape({Ns, Ns});
    result.Vxc_beta  = xt::xtensor<double, 2>::from_shape({Ns, Ns});
    for (int mu = 0; mu < N; ++mu)
        for (int nu = 0; nu < N; ++nu) {
            result.Vxc_alpha(mu, nu) = static_cast<double>(Vxc_a_d(mu, nu));
            result.Vxc_beta(mu, nu)  = static_cast<double>(Vxc_b_d(mu, nu));
        }
    result.E_xc = static_cast<double>(E_xc);
    result.n_alpha_integrated = static_cast<double>(n_alpha_int);
    result.n_beta_integrated  = static_cast<double>(n_beta_int);
    return result;
}

} // namespace scf
} // namespace cpu
} // namespace ATV
