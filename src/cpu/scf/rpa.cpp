/// @file cpu/scf/rpa.cpp
/// @brief Closed-shell direct RPA (dRPA) correlation energy — Furche plasmon
///        formula on the Coulomb-only Casida spectrum (α_hf = 0).

#include "atomiverse/cpu/scf/rpa.h"

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/containers/xtensor.hpp>

#include "atomiverse/cpu/scf/dense_matrix.h"
#include "atomiverse/cpu/scf/internal/integrals.h"
#include "atomiverse/cpu/scf/internal/orbital_hessian.h"
#include "atomiverse/cpu/scf/rhf.h"
#include "atomiverse/cpu/scf/slater_determinant.h"

namespace ATV {
namespace cpu {
namespace scf {
namespace {

constexpr double kAmBFloor = 1e-10;  ///< (A−B) positivity floor [Hartree].
constexpr double kOmega2Tol = 1e-10; ///< |ω²| zero / instability guard.

/// Symmetric M^{p} via eigh (p = ±½ for the RPA reduction).
RealMatrix sym_power(const RealMatrix& M, double p, const char* what) {
    auto [s, U] = xt::linalg::eigh(M);
    const std::size_t m = M.shape()[0];
    if (s(0) <= kAmBFloor)
        throw std::runtime_error(
            std::string(what) + ": (A - B) is not positive definite "
            "(lowest eigenvalue = " + std::to_string(s(0))
            + " Ha) — dRPA reduction undefined (no fallback)");
    RealMatrix Us = U;
    for (std::size_t c = 0; c < m; ++c) {
        const double w = std::pow(s(c), p);
        for (std::size_t r = 0; r < m; ++r)
            Us(r, c) *= w;
    }
    RealMatrix out = xt::linalg::dot(Us, xt::transpose(U));
    for (std::size_t r = 0; r < m; ++r)
        for (std::size_t c = r + 1; c < m; ++c) {
            const double v = 0.5 * (out(r, c) + out(c, r));
            out(r, c) = v;
            out(c, r) = v;
        }
    return out;
}

/// Dense dRPA eigenvalues only (plasmon needs Ω_s and Tr A, not X/Y).
void solve_drpa_omega(const internal::OrbitalHessian& hess, RealVector& omega,
                      double& tr_A) {
    const RealMatrix A =
        hess.dense_A(internal::Channel::singlet, /*alpha_hf=*/0.0);
    const RealMatrix B =
        hess.dense_B(internal::Channel::singlet, /*alpha_hf=*/0.0);
    const std::size_t n = A.shape()[0];
    tr_A = 0.0;
    for (std::size_t i = 0; i < n; ++i)
        tr_A += A(i, i);

    const RealMatrix Mp = A + B;
    const RealMatrix Mm = A - B;
    const RealMatrix T = sym_power(Mm, 0.5, "RHF::direct_rpa/dRPA");
    RealMatrix M = xt::linalg::dot(T, RealMatrix(xt::linalg::dot(Mp, T)));
    for (std::size_t r = 0; r < n; ++r)
        for (std::size_t c = r + 1; c < n; ++c) {
            const double v = 0.5 * (M(r, c) + M(c, r));
            M(r, c) = v;
            M(c, r) = v;
        }
    auto [w2, Z] = xt::linalg::eigh(M);
    (void)Z;

    omega = xt::zeros<double>({n});
    for (std::size_t s = 0; s < n; ++s) {
        if (w2(s) < -kOmega2Tol)
            throw std::runtime_error(
                "RHF::direct_rpa: dRPA omega^2 = " + std::to_string(w2(s))
                + " Ha^2 (root " + std::to_string(s)
                + ") is negative — imaginary excitation (no fallback)");
        if (w2(s) <= kOmega2Tol)
            throw std::runtime_error(
                "RHF::direct_rpa: dRPA omega^2 = " + std::to_string(w2(s))
                + " Ha^2 (root " + std::to_string(s)
                + ") is numerically zero (no fallback)");
        omega(s) = std::sqrt(w2(s));
    }
}

} // namespace

DirectRPAResult RHF::direct_rpa() const {
    if (!has_result())
        throw std::runtime_error(
            "RHF::direct_rpa: no SCF result — call run(frame) first");

    const SlaterDeterminant& det = determinant();
    if (!det.restricted() || det.generalized())
        throw std::runtime_error(
            "RHF::direct_rpa: closed-shell RHF reference required "
            "(no fallback)");

    const IntegralPack& pack = integral_workspace_();
    const xt::xtensor<double, 4> eri_ao =
        (pack.ERI.size() > 0) ? pack.ERI : pack.eri_working();

    using Spin = SlaterDeterminant::Spin;
    const RealMatrix& C = det.coefficients(Spin::alpha);
    const RealVector& eps_raw = det.orbital_energies(Spin::alpha);
    const int n_mo = static_cast<int>(C.shape()[1]);
    if (n_mo != static_cast<int>(eps_raw.size()))
        throw std::runtime_error("RHF::direct_rpa: MO dimension mismatch");

    int n_occ = 0;
    for (int p = 0; p < n_mo; ++p) {
        const double f =
            det.occupations(Spin::alpha)(static_cast<std::size_t>(p));
        if (f == 1.0)
            ++n_occ;
        else if (f != 0.0)
            throw std::invalid_argument(
                "RHF::direct_rpa: fractional occupation — integer "
                "occupations only");
    }
    const int n_vir = n_mo - n_occ;
    if (n_occ < 1 || n_vir < 1)
        throw std::invalid_argument(
            "RHF::direct_rpa: need at least one occupied and one virtual "
            "orbital (no fallback)");

    const internal::OrbitalHessian hess{det, eri_ao};
    if (hess.n_occ() != n_occ || hess.n_virt() != n_vir)
        throw std::runtime_error(
            "RHF::direct_rpa: Hessian occupation mismatch");

    RealVector Om;
    double tr_A = 0.0;
    solve_drpa_omega(hess, Om, tr_A);

    double sum_om = 0.0;
    DirectRPAResult out;
    out.tr_A = tr_A;
    out.omega_ha.resize(Om.size());
    for (std::size_t s = 0; s < Om.size(); ++s) {
        out.omega_ha[s] = Om(s);
        sum_om += Om(s);
    }
    // Furche plasmon formula (closed-shell spin-adapted singlet dRPA).
    out.e_corr = 0.5 * (sum_om - tr_A);
    out.e_tot = energy() + out.e_corr;
    return out;
}

} // namespace scf
} // namespace cpu
} // namespace ATV
