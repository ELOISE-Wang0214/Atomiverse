/// @file tests/cpu/scf/test_ks_grid_response_gradient.cpp
/// @brief Tests for the Kohn-Sham XC grid-response nuclear gradient.
///
/// Spec: ks-xc-grid-response-gradient.
///
/// The CPU KS analytic nuclear gradient previously computed the XC
/// contribution `dE_xc/dR_A` on a grid held fixed in space and projected
/// out the residual net-translational quadrature error to restore
/// translational invariance. This file is the RED-before-GREEN harness for
/// the genuine grid-response term (Becke weight-derivative + grid-
/// coordinate), a verbatim port of PySCF's `grids_response_becke`.
///
/// Acceptance criteria covered:
///   - ac-001: becke_weight_derivative matches a central difference of
///     becke_partition_weights_for_test (h=1e-4 Bohr) to <= 1e-6.
///   - ac-002: per-point sum_A dw_g/dR_A == 0 (rigid-shift identity).
///   - ac-003: with xc_grad_project_translation=false, the net
///     translational XC force on water PBE/def2-SVP is <= 1e-8 Ha/Bohr.
///   - ac-004: the analytic KS gradient matches a moving-grid central-
///     difference FD (grid regenerated per displaced geometry, h=1e-3
///     Bohr, SCF tol 1e-9) to <= 1e-5 Ha/Bohr for water PBE/def2-SVP and
///     PBE/6-31G*, RKS and UKS.

#ifdef ATV_HAS_RHF

#include <gtest/gtest.h>

#include <cmath>
#include <cstdio>
#include <vector>

#include "atomiverse/core/mat3.h"
#include "atomiverse/core/units.h"
#include "atomiverse/cpu/frame.h"
#include "atomiverse/cpu/scf/basis_set.h"
#include "atomiverse/cpu/scf/gradient.h"
#include "atomiverse/cpu/scf/grid.h"
#include "atomiverse/cpu/scf/integrals.h"
#include "atomiverse/cpu/scf/native_integrals.h"
#include "atomiverse/cpu/scf/scf_solve.h"

using namespace ATV;
using namespace ATV::cpu::scf;

namespace {

/// @brief Convert a force component from eV/Angstrom to Hartree/Bohr.
inline double ev_a_to_ha_bohr(double f_ev_a) {
    return units::convert(f_ev_a, units::Force::eVPerAngstrom,
                          units::Force::HartreePerBohr);
}

/// @brief Water near its def2-SVP / 6-31G* RHF equilibrium (C2v, z-axis).
cpu::Frame make_h2o() {
    cpu::Frame f(3);
    f.atoms.type[0] = 8;
    f.atoms.mass[0] = 15.999;
    f.atoms.type[1] = 1;
    f.atoms.mass[1] = 1.008;
    f.atoms.type[2] = 1;
    f.atoms.mass[2] = 1.008;
    f.atoms.x[0] = 0.000000;
    f.atoms.y[0] = 0.000000;
    f.atoms.z[0] = 0.117790;
    f.atoms.x[1] = 0.000000;
    f.atoms.y[1] = 0.755453;
    f.atoms.z[1] = -0.471161;
    f.atoms.x[2] = 0.000000;
    f.atoms.y[2] = -0.755453;
    f.atoms.z[2] = -0.471161;
    f.box = Mat3::cubic(100.0);
    return f;
}

/// @brief Asymmetric three-atom geometry in Bohr for weight-derivative
///        unit checks (O + 2 H, no symmetry to mask sign errors).
std::vector<native::Atom> make_asymmetric_h2o_bohr() {
    return {
        native::Atom{8, {0.10, 0.05, 0.20}},
        native::Atom{1, {1.43, 0.31, -0.89}},
        native::Atom{1, {-1.21, -1.10, -0.77}},
    };
}

} // namespace

// ── ac-001 — weight derivative vs finite difference ──────────────────────

TEST(KsGridResponseGradient, BeckeWeightDerivativeMatchesFD) {
    const auto atoms = make_asymmetric_h2o_bohr();
    const int natm = static_cast<int>(atoms.size());
    const double h = 1e-4; // Bohr

    // A scattered set of probe points spanning the molecular region.
    const double probes[][3] = {
        {0.30, 0.40, 0.10},  {-0.50, 0.20, 0.60},  {0.90, -0.40, -0.20},
        {0.05, 1.10, -0.30}, {-0.80, -0.60, 0.40}, {1.20, 0.10, 0.55},
        {0.00, 0.00, 1.30},  {-1.40, 0.50, -0.10},
    };

    double max_err = 0.0;
    std::vector<double> dw(static_cast<size_t>(natm) * 3);
    std::vector<double> w_plus(natm), w_minus(natm);

    for (const auto& p : probes) {
        const double rx = p[0], ry = p[1], rz = p[2];
        for (int parent = 0; parent < natm; ++parent) {
            // Analytic: total d(normalized weight of `parent`)/dR_A. The
            // grid point rides rigidly with its parent atom, so the
            // returned derivative is the TOTAL derivative -- partial w.r.t.
            // R_A plus, when A == parent, the grid-point displacement.
            becke_weight_derivative(atoms.data(), natm, rx, ry, rz, parent,
                                    GridLevel::SG0, dw.data());

            for (int A = 0; A < natm; ++A) {
                for (int d = 0; d < 3; ++d) {
                    // Central difference of the SAME total derivative:
                    // displace atom A by +/- h and, when A is the parent,
                    // displace the probe point identically (the grid point
                    // rides with its parent).
                    auto a_plus = atoms;
                    auto a_minus = atoms;
                    a_plus[A].pos[d] += h;
                    a_minus[A].pos[d] -= h;
                    double rp[3] = {rx, ry, rz};
                    double rm[3] = {rx, ry, rz};
                    if (A == parent) {
                        rp[d] += h;
                        rm[d] -= h;
                    }
                    becke_partition_weights_for_test(a_plus.data(), natm, rp[0],
                                                     rp[1], rp[2],
                                                     w_plus.data());
                    becke_partition_weights_for_test(a_minus.data(), natm,
                                                     rm[0], rm[1], rm[2],
                                                     w_minus.data());
                    const double fd =
                        (w_plus[parent] - w_minus[parent]) / (2.0 * h);
                    const double an = dw[A * 3 + d];
                    max_err = std::max(max_err, std::abs(an - fd));
                    EXPECT_NEAR(an, fd, 1e-6)
                        << "parent=" << parent << " A=" << A << " d=" << d
                        << " at (" << rx << "," << ry << "," << rz << ")";
                }
            }
        }
    }
    std::printf("[grid-response] ac-001 max|dw_analytic - dw_FD| = %.3e\n",
                max_err);
}

// ── ac-002 — rigid-shift identity sum_A dw_g/dR_A == 0 ───────────────────

TEST(KsGridResponseGradient, BeckeWeightDerivativeSumsToZero) {
    const auto atoms = make_asymmetric_h2o_bohr();
    const int natm = static_cast<int>(atoms.size());

    const double probes[][3] = {
        {0.30, 0.40, 0.10},  {-0.50, 0.20, 0.60},  {0.90, -0.40, -0.20},
        {0.05, 1.10, -0.30}, {-0.80, -0.60, 0.40}, {1.20, 0.10, 0.55},
    };

    double max_residual = 0.0;
    std::vector<double> dw(static_cast<size_t>(natm) * 3);
    for (const auto& p : probes) {
        for (int parent = 0; parent < natm; ++parent) {
            becke_weight_derivative(atoms.data(), natm, p[0], p[1], p[2],
                                    parent, GridLevel::SG0, dw.data());
            for (int d = 0; d < 3; ++d) {
                double s = 0.0;
                for (int A = 0; A < natm; ++A)
                    s += dw[A * 3 + d];
                max_residual = std::max(max_residual, std::abs(s));
                EXPECT_NEAR(s, 0.0, 1e-10)
                    << "rigid-shift identity broken: parent=" << parent
                    << " axis=" << d;
            }
        }
    }
    std::printf("[grid-response] ac-002 max|sum_A dw_g/dR_A| = %.3e\n",
                max_residual);
}

// ── ac-003 — net translational XC force vanishes without projection ──────

TEST(KsGridResponseGradient, NetTranslationalForceVanishesWithoutProjection) {
    const auto frame = make_h2o();
    const auto basis = cpu::scf::BasisSet::Def2SVP;

    auto ints = cpu::scf::build_integrals(frame, basis);
    auto cfg = cpu::scf::default_config(ints); // KS / PBE / MatchPySCF
    cfg.tol = 1e-9;
    cfg.max_iter = 200;
    cfg.xc_grad_project_translation = false; // genuine first-principles path

    auto res = cpu::scf::scf_solve(ints, cfg);
    ASSERT_TRUE(res.converged) << "PBE RKS should converge";

    auto g = cpu::scf::scf_gradient(frame, ints, res, cfg);

    double max_net = 0.0;
    for (int d = 0; d < 3; ++d) {
        double s = 0.0;
        for (int A = 0; A < 3; ++A)
            s += ev_a_to_ha_bohr(g.forces(A, d));
        max_net = std::max(max_net, std::abs(s));
        EXPECT_LE(std::abs(s), 1e-8) << "net translational force on axis " << d
                                     << " must vanish without projection";
    }
    std::printf("[grid-response] ac-003 ||sum_A dE_xc/dR_A|| = %.3e Ha/Bohr\n",
                max_net);
}

// ── ac-004 — analytic KS gradient vs moving-grid FD ──────────────────────

namespace {

/// @brief Run the moving-grid FD parity check for one config.
///
/// The FD oracle regenerates the integration grid at each displaced
/// geometry (fixed_grid left null), so it probes the genuine moving-grid
/// gradient — exactly what the projection shortcut cannot reproduce.
void check_moving_grid_fd(const cpu::Frame& frame, cpu::scf::BasisSet basis,
                          cpu::scf::SCFConfig cfg, const char* label) {
    cfg.tol = 1e-9;
    cfg.max_iter = 200;
    cfg.xc_grad_project_translation = false;

    auto ints = cpu::scf::build_integrals(frame, basis);
    auto res = cpu::scf::scf_solve(ints, cfg);
    ASSERT_TRUE(res.converged) << label << ": SCF should converge";
    auto analytic = cpu::scf::scf_gradient(frame, ints, res, cfg);

    // Moving-grid FD: fixed_grid stays null, so scf_finite_difference_gradient
    // regenerates the grid at each displacement.
    const double h_ang = 1e-3 / units::kAngstromInBohr;
    auto fd = cpu::scf::scf_finite_difference_gradient(frame, cfg, /*charge=*/0,
                                                       h_ang, basis);

    double max_err = 0.0;
    const int natm = static_cast<int>(frame.atoms.size());
    for (int A = 0; A < natm; ++A) {
        for (int d = 0; d < 3; ++d) {
            const double e = ev_a_to_ha_bohr(
                std::abs(analytic.forces(A, d) - fd.forces(A, d)));
            max_err = std::max(max_err, e);
            EXPECT_NEAR(ev_a_to_ha_bohr(analytic.forces(A, d)),
                        ev_a_to_ha_bohr(fd.forces(A, d)), 1e-5)
                << label << " atom=" << A << " axis=" << d;
        }
    }
    std::printf("[grid-response] ac-004 %s max|analytic-FD| = %.3e Ha/Bohr\n",
                label, max_err);
}

} // namespace

TEST(KsGridResponseGradient, MovingGridFD_RKS_PBE_def2SVP) {
    const auto frame = make_h2o();
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::Def2SVP);
    auto cfg = cpu::scf::default_config(ints); // RKS / PBE / MatchPySCF
    check_moving_grid_fd(frame, cpu::scf::BasisSet::Def2SVP, cfg,
                         "RKS PBE/def2-SVP");
}

TEST(KsGridResponseGradient, MovingGridFD_RKS_PBE_631Gs) {
    const auto frame = make_h2o();
    auto ints =
        cpu::scf::build_integrals(frame, cpu::scf::BasisSet::Pople631Gs);
    auto cfg = cpu::scf::default_config(ints); // RKS / PBE / MatchPySCF
    check_moving_grid_fd(frame, cpu::scf::BasisSet::Pople631Gs, cfg,
                         "RKS PBE/6-31G*");
}

TEST(KsGridResponseGradient, MovingGridFD_UKS_PBE_def2SVP) {
    const auto frame = make_h2o();
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::Def2SVP);
    auto cfg = cpu::scf::default_config(ints);
    cfg.spin = cpu::scf::SpinType::Unrestricted;
    cfg.n_alpha = ints.n_occ();
    cfg.n_beta = ints.n_occ();
    check_moving_grid_fd(frame, cpu::scf::BasisSet::Def2SVP, cfg,
                         "UKS PBE/def2-SVP");
}

TEST(KsGridResponseGradient, MovingGridFD_UKS_PBE_631Gs) {
    const auto frame = make_h2o();
    auto ints =
        cpu::scf::build_integrals(frame, cpu::scf::BasisSet::Pople631Gs);
    auto cfg = cpu::scf::default_config(ints);
    cfg.spin = cpu::scf::SpinType::Unrestricted;
    cfg.n_alpha = ints.n_occ();
    cfg.n_beta = ints.n_occ();
    check_moving_grid_fd(frame, cpu::scf::BasisSet::Pople631Gs, cfg,
                         "UKS PBE/6-31G*");
}

#endif // ATV_HAS_RHF
