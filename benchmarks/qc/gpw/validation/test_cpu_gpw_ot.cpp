/// @file tests/cpu/gpw/test_cpu_gpw_ot.cpp
/// @brief OT inner-minimizer unit tests, decoupled from the SCF/grid:
///        (1) parameterization keeps C^T S C = I exactly,
///        (2) the analytic gradient matches finite differences (guards the
///            factor-4 / S-projector / sign conventions),
///        (3) fixed-H OT minimization reproduces the generalized
///            diagonalization band energy and density to ~1e-9.

#include <cmath>
#include <cstdio>
#include <random>
#include <tuple>

#include <gtest/gtest.h>

#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/containers/xtensor.hpp>

#include "atomiverse/cpu/gpw/gpw_solver.h"
#include "atomiverse/cpu/gpw/internal/cp2k_basis.h"
#include "atomiverse/cpu/gpw/internal/cp2k_io.h"
#include "atomiverse/cpu/gpw/internal/ot.h"
#include "atomiverse/cpu/gpw/basis_data/gpw_data.h"
#include "atomiverse/cpu/linalg/density.h"
#include "atomiverse/cpu/linalg/eigensolve.h"

using ATV::cpu::linalg::RealMatrix;
using ATV::cpu::gpw::internal::ot_gradient;
using ATV::cpu::gpw::internal::ot_minimize_fixed;
using ATV::cpu::gpw::internal::ot_orthonormalize;
using ATV::cpu::gpw::internal::ot_parameterize;
using ATV::cpu::gpw::internal::ot_tangent;

namespace {

// Deterministic small SPD overlap, symmetric H, and an orthonormal guess.
struct Problem {
    RealMatrix S, H, T, C0;
    int nb, no;
};

Problem make_problem(int nb, int no, unsigned seed) {
    std::mt19937 rng(seed);
    std::normal_distribution<double> nd(0.0, 1.0);
    RealMatrix A = xt::zeros<double>({nb, nb});
    RealMatrix H = xt::zeros<double>({nb, nb});
    for (int i = 0; i < nb; ++i)
        for (int j = 0; j < nb; ++j) {
            A(i, j) = nd(rng);
            H(i, j) = nd(rng);
        }
    RealMatrix S = xt::linalg::dot(A, xt::transpose(A));
    for (int i = 0; i < nb; ++i) S(i, i) += static_cast<double>(nb); // well-conditioned SPD
    RealMatrix Hs = 0.5 * (H + xt::transpose(H));                    // symmetric
    RealMatrix Cr = xt::zeros<double>({nb, no});
    for (int i = 0; i < nb; ++i)
        for (int j = 0; j < no; ++j) Cr(i, j) = nd(rng);
    RealMatrix C0 = ot_orthonormalize(Cr, S);
    RealMatrix T = xt::zeros<double>({nb, nb}); // identity preconditioner basis
    for (int i = 0; i < nb; ++i) T(i, i) = 1.0;
    return {S, Hs, T, C0, nb, no};
}

double max_abs(const RealMatrix& m) {
    double r = 0.0;
    for (double v : m) r = std::max(r, std::abs(v));
    return r;
}

} // namespace

TEST(GpwOt, ParameterizationKeepsOrthonormal) {
    const Problem p = make_problem(10, 4, 1);
    std::mt19937 rng(7);
    std::normal_distribution<double> nd(0.0, 1.0);
    RealMatrix d = xt::zeros<double>({p.nb, p.no});
    for (int i = 0; i < p.nb; ++i)
        for (int j = 0; j < p.no; ++j) d(i, j) = nd(rng);
    const RealMatrix x = ot_tangent(p.C0, p.S, d); // S-tangent

    for (double alpha : {0.05, 0.5, 1.7, 4.0}) {
        const RealMatrix C = ot_parameterize(p.C0, alpha * x, p.S);
        RealMatrix G = xt::linalg::dot(xt::transpose(C), xt::linalg::dot(p.S, C));
        for (int i = 0; i < p.no; ++i) G(i, i) -= 1.0;
        EXPECT_LT(max_abs(G), 1e-11) << "C^T S C != I at alpha=" << alpha;
    }
}

TEST(GpwOt, GradientMatchesFiniteDifference) {
    const Problem p = make_problem(9, 3, 2);
    auto band = [&](const RealMatrix& C) {
        return 2.0 * xt::sum(C * xt::linalg::dot(p.H, C))(); // 2 Tr(C^T H C)
    };
    const RealMatrix g = ot_gradient(p.H, p.S, p.C0);

    std::mt19937 rng(11);
    std::normal_distribution<double> nd(0.0, 1.0);
    RealMatrix dd = xt::zeros<double>({p.nb, p.no});
    for (int i = 0; i < p.nb; ++i)
        for (int j = 0; j < p.no; ++j) dd(i, j) = nd(rng);
    const RealMatrix d = ot_tangent(p.C0, p.S, dd);

    const double eps = 1e-6;
    const double e_plus = band(ot_parameterize(p.C0, eps * d, p.S));
    const double e_minus = band(ot_parameterize(p.C0, -eps * d, p.S));
    const double fd = (e_plus - e_minus) / (2.0 * eps);
    const double analytic = xt::sum(g * d)(); // directional derivative g . d
    std::printf("[GPW OT] g.d analytic=%.10f  fd=%.10f  d=%.2e\n", analytic, fd,
                std::abs(analytic - fd));
    EXPECT_NEAR(analytic, fd, 1e-5 * std::max(1.0, std::abs(fd)));
}

TEST(GpwOt, FixedHReproducesDiagonalization) {
    const Problem p = make_problem(12, 4, 3);
    // Reference: generalized diagonalization, occupy the no lowest states.
    auto [eig, Call] = ATV::cpu::linalg::diagonalize_fock_generalized(p.H, p.S);
    double band_ref = 0.0;
    for (int i = 0; i < p.no; ++i) band_ref += 2.0 * eig(i);
    const RealMatrix P_ref =
        ATV::cpu::linalg::build_density_matrix(Call, p.no, p.nb, 2.0);

    const auto r = ot_minimize_fixed(p.H, p.S, p.T, p.C0, 1e-9, 500);
    const RealMatrix P_ot =
        ATV::cpu::linalg::build_density_matrix(r.C, p.no, p.nb, 2.0);

    std::printf("[GPW OT] fixed-H  E_ot=%.10f  E_diag=%.10f  dE=%.2e  iters=%d "
                "conv=%d  ||dP||=%.2e\n",
                r.band_energy, band_ref, std::abs(r.band_energy - band_ref),
                r.iters, (int)r.converged, max_abs(P_ot - P_ref));
    EXPECT_TRUE(r.converged);
    EXPECT_NEAR(r.band_energy, band_ref, 1e-8);
    EXPECT_LT(max_abs(P_ot - P_ref), 1e-6); // same occupied subspace (density)
}

// ── Self-consistent equivalence: full gpw_scf, OT vs diagonalization ────────
namespace {
using ATV::cpu::gpw::GpwConfig;
using ATV::cpu::gpw::GpwResult;
using ATV::cpu::gpw::gpw_scf;
using ATV::cpu::gpw::make_gpw_atom;
using ATV::cpu::gpw::builtin_gth_basis;
using ATV::cpu::gpw::builtin_gth_potential;
} // namespace

TEST(GpwOtScf, H2MatchesDiagonalization) {
    const auto basis = builtin_gth_basis("H", "DZV-GTH-PADE");
    const auto pot = builtin_gth_potential("H", "GTH-PADE-q1");
    const double rA[3] = {3.5, 2.0, 2.0};
    const double rB[3] = {4.5, 2.0, 2.0};
    std::vector<ATV::cpu::gpw::GpwAtom> atoms = {make_gpw_atom(rA, basis, pot),
                                                 make_gpw_atom(rB, basis, pot)};
    GpwConfig cfg;
    cfg.cell[0] = 8.0; cfg.cell[1] = 4.0; cfg.cell[2] = 4.0;
    cfg.cutoff_ry = 140.0; cfg.n_occ = 1; cfg.max_iter = 300; cfg.tol = 1e-8;

    cfg.use_ot = false; const GpwResult diag = gpw_scf(atoms, cfg);
    cfg.use_ot = true;  const GpwResult ot = gpw_scf(atoms, cfg);
    std::printf("[GPW OT-SCF H2] diag E=%.10f (%d it)  ot E=%.10f (%d it)  dE=%.2e\n",
                diag.e_total, diag.iters, ot.e_total, ot.iters,
                std::abs(ot.e_total - diag.e_total));
    EXPECT_TRUE(diag.converged);
    EXPECT_TRUE(ot.converged);
    EXPECT_NEAR(ot.e_total, diag.e_total, 1e-7);
}

TEST(GpwOtScf, ZnOMatchesDiagonalization) {
    constexpr double kAng = 1.8897259885789;
    const auto bZn = builtin_gth_basis("Zn", "DZVP-MOLOPT-SR-GTH");
    const auto bO = builtin_gth_basis("O", "DZVP-MOLOPT-SR-GTH");
    const auto pZn = builtin_gth_potential("Zn", "GTH-PBE-q12");
    const auto pO = builtin_gth_potential("O", "GTH-PBE-q6");
    const double rZn[3] = {4.0 * kAng, 4.0 * kAng, 4.0 * kAng};
    const double rO[3] = {4.0 * kAng, 4.0 * kAng, 5.72 * kAng};
    std::vector<ATV::cpu::gpw::GpwAtom> atoms = {make_gpw_atom(rZn, bZn, pZn),
                                                 make_gpw_atom(rO, bO, pO)};
    GpwConfig cfg;
    cfg.cell[0] = cfg.cell[1] = cfg.cell[2] = 8.0 * kAng;
    cfg.cutoff_ry = 120.0; // cheap grid (consistency, not absolute accuracy)
    cfg.n_occ = 9; cfg.max_iter = 400; cfg.tol = 1e-8; cfg.xc = ATV::cpu::gpw::XcKind::PBE;

    cfg.use_ot = false; const GpwResult diag = gpw_scf(atoms, cfg);
    cfg.use_ot = true;  const GpwResult ot = gpw_scf(atoms, cfg);
    std::printf("[GPW OT-SCF ZnO] diag E=%.10f (%d it)  ot E=%.10f (%d it)  dE=%.2e\n",
                diag.e_total, diag.iters, ot.e_total, ot.iters,
                std::abs(ot.e_total - diag.e_total));
    EXPECT_TRUE(diag.converged);
    EXPECT_TRUE(ot.converged);
    EXPECT_NEAR(ot.e_total, diag.e_total, 2e-7);
}
