/// @file tests/cpu/scf/test_hybrid_gradient.cpp
/// @brief B3LYP / PBE0 analytic nuclear gradient + hybrid geometry-opt tests.
///
/// Acceptance criteria for spec `phase3-hybrid-dft-03-hybrid-gradient`:
///   - ac-001: B3LYP/6-31G* analytic force vs PySCF grad.RKS, <= 1e-5 Ha/Bohr
///             on H2O, CH4, NH3.
///   - ac-002: PBE0/6-31G* analytic force vs PySCF grad.RKS, <= 1e-5 Ha/Bohr
///             on H2O, CH4, NH3.
///   - ac-003: B3LYP/PBE0 analytic force vs moving-grid finite difference,
///             <= 1e-5 Ha/Bohr (independent oracle; does not trust PySCF).
///   - ac-004: translational invariance |sum_A F_A| <= 1e-8 Ha/Bohr.
///   - ac-005: B3LYP/6-31G* geometry optimization converges and the
///             optimized water RMS <= 1e-3 Angstrom + final energy
///             <= 1e-6 Ha vs the PySCF optimized reference.
///
/// Reference data are GENUINE PySCF runs (not invented), generated with:
///   PySCF 2.13.0
///   gto.M(basis="6-31g*", unit="Angstrom", cart=True)
///   dft.RKS(mol), mf.xc = "b3lyp" / "pbe0", conv_tol = 1e-12,
///   grids.level = 3
///   forces  = mf.nuc_grad_method().kernel()   # dE/dR, Hartree/Bohr
///   geom-opt: scipy.optimize.minimize(BFGS, gtol=1e-6) on PySCF E + grad,
///             converged to max|grad| = 2.7e-6 Ha/Bohr.
///
/// Geometries are verbatim from `benchmarks/qc/molecules/*.xyz`, identical
/// to the make_* fixtures in `test_hybrid_singlepoint.cpp` — the exact
/// geometries the PySCF reference runs used.

#ifdef ATV_HAS_RHF

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "atomiverse/core/mat3.h"
#include "atomiverse/core/units.h"
#include "atomiverse/cpu/frame.h"
#include "atomiverse/cpu/minimizer/lbfgs.h"
#include "atomiverse/cpu/minimizer/minimize.h"
#include "atomiverse/cpu/scf/basis_set.h"
#include "atomiverse/cpu/scf/gradient.h"
#include "atomiverse/cpu/scf/integrals.h"
#include "atomiverse/cpu/scf/scf_evaluator.h"
#include "atomiverse/cpu/scf/scf_solve.h"
#include "atomiverse/cpu/scf/xc_functional.h"

using namespace ATV;
using ATV::cpu::scf::BasisSet;
using ATV::cpu::scf::XCType;

namespace {

/// @brief Convert a force component from eV/Angstrom to Hartree/Bohr.
inline double ev_a_to_ha_bohr(double f_ev_a) {
    return units::convert(f_ev_a, units::Force::eVPerAngstrom,
                          units::Force::HartreePerBohr);
}

// ── Molecule geometries (verbatim from benchmarks/qc/molecules/*.xyz) ────

cpu::Frame make_h2o() {
    cpu::Frame f(3);
    f.atoms.type = {8, 1, 1};
    f.atoms.x = {0.000000, 0.000000, 0.000000};
    f.atoms.y = {0.000000, 0.755453, -0.755453};
    f.atoms.z = {0.117790, -0.471161, -0.471161};
    f.atoms.mass = {15.999, 1.008, 1.008};
    f.box = Mat3::cubic(100.0);
    return f;
}

cpu::Frame make_ch4() {
    cpu::Frame f(5);
    f.atoms.type = {6, 1, 1, 1, 1};
    f.atoms.x = {0.000000, 0.627640, -0.627640, 0.627640, -0.627640};
    f.atoms.y = {0.000000, 0.627640, -0.627640, -0.627640, 0.627640};
    f.atoms.z = {0.000000, 0.627640, 0.627640, -0.627640, -0.627640};
    f.atoms.mass = {12.011, 1.008, 1.008, 1.008, 1.008};
    f.box = Mat3::cubic(100.0);
    return f;
}

cpu::Frame make_nh3() {
    cpu::Frame f(4);
    f.atoms.type = {7, 1, 1, 1};
    f.atoms.x = {0.000000, 0.000000, 0.812328, -0.812328};
    f.atoms.y = {0.000000, 0.938035, -0.469017, -0.469017};
    f.atoms.z = {0.113880, -0.265825, -0.265825, -0.265825};
    f.atoms.mass = {14.007, 1.008, 1.008, 1.008};
    f.box = Mat3::cubic(100.0);
    return f;
}

cpu::Frame make_molecule(const std::string& stem) {
    if (stem == "h2o")
        return make_h2o();
    if (stem == "ch4")
        return make_ch4();
    if (stem == "nh3")
        return make_nh3();
    throw std::runtime_error("unknown molecule: " + stem);
}

// ── PySCF 2.13.0 reference forces, grad.RKS(xc=...).kernel() ─────────────
//
// Stored as the PySCF gradient dE/dR in Hartree/Bohr ([natm][3]). The
// Atomiverse analytic gradient is compared after converting its eV/Angstrom
// forces back to Hartree/Bohr with the FORCE = -dE/dR sign flip.
//
// B3LYP/6-31G*/cart, grids.level = 3, conv_tol = 1e-12.
struct RefForces {
    const char* stem;
    std::vector<std::array<double, 3>> grad; // dE/dR, Ha/Bohr
};

const std::vector<RefForces>& b3lyp_reference() {
    static const std::vector<RefForces> ref = {
        {"h2o",
         {
             {-1.612833381966e-15, -1.705742921477e-15, -1.310086217001e-02},
             {-2.213629307622e-16, -8.249004083599e-03, 6.554047210348e-03},
             {4.765752619961e-17, 8.249004083600e-03, 6.554047210347e-03},
         }},
        {"ch4",
         {
             {-3.195908022222e-15, 1.773291293473e-15, 1.313315824143e-15},
             {-2.431483947497e-03, -2.431483947495e-03, -2.431483947496e-03},
             {2.431483947497e-03, 2.431483947495e-03, -2.431483947496e-03},
             {-2.431483947497e-03, 2.431483947495e-03, 2.431483947495e-03},
             {2.431483947498e-03, -2.431483947495e-03, 2.431483947496e-03},
         }},
        {"nh3",
         {
             {-3.981440475394e-14, -1.742668852297e-05, -8.209616933049e-03},
             {9.831397196930e-16, -4.527993911985e-03, 2.730620284688e-03},
             {-3.944416598231e-03, 2.273199951145e-03, 2.737743167405e-03},
             {3.944416598234e-03, 2.273199951144e-03, 2.737743167405e-03},
         }},
    };
    return ref;
}

// PBE0/6-31G*/cart, grids.level = 3, conv_tol = 1e-12.
const std::vector<RefForces>& pbe0_reference() {
    static const std::vector<RefForces> ref = {
        {"h2o",
         {
             {-1.338316571841e-15, -8.456672528774e-16, -8.415435826656e-03},
             {2.356922883545e-16, -5.231182081665e-03, 4.209926730038e-03},
             {-2.289012244022e-16, 5.231182081664e-03, 4.209926730037e-03},
         }},
        {"ch4",
         {
             {7.002495363078e-16, -2.078354458644e-15, -1.989664157781e-15},
             {-2.176037678373e-03, -2.176037678373e-03, -2.176037678372e-03},
             {2.176037678373e-03, 2.176037678375e-03, -2.176037678372e-03},
             {-2.176037678374e-03, 2.176037678373e-03, 2.176037678373e-03},
             {2.176037678372e-03, -2.176037678373e-03, 2.176037678373e-03},
         }},
        {"nh3",
         {
             {2.013460145031e-14, -1.763242409458e-05, -5.684077768786e-03},
             {-1.068314767823e-16, -2.605231774153e-03, 1.888906795949e-03},
             {-2.279577859184e-03, 1.311792297317e-03, 1.896415878683e-03},
             {2.279577859185e-03, 1.311792297318e-03, 1.896415878681e-03},
         }},
    };
    return ref;
}

// ── B3LYP/6-31G* optimized water reference ──────────────────────────────
//
// scipy BFGS on PySCF 2.13.0 B3LYP/6-31G* E + analytic gradient, converged
// to max|grad| = 2.679e-6 Ha/Bohr. Optimized geometry (Angstrom):
//   O   0.0000000000   0.0000000000   0.1242918872
//   H   0.0000000000   0.7615569969  -0.4744358138
//   H  -0.0000000000  -0.7615569969  -0.4744358138
// OH bond = 0.968733 A, HOH angle = 103.6519 deg.
constexpr double kB3LYPOptEnergyHa = -76.4089538451;
constexpr double kB3LYPOptGeom[3][3] = {
    {0.0000000000, 0.0000000000, 0.1242918872},
    {0.0000000000, 0.7615569969, -0.4744358138},
    {0.0000000000, -0.7615569969, -0.4744358138},
};

// Run the analytic hybrid gradient for one molecule at 6-31G*.
cpu::scf::GradientResult run_hybrid_gradient(const std::string& stem,
                                             XCType xc) {
    auto frame = make_molecule(stem);
    auto integrals = cpu::scf::build_integrals(frame, BasisSet::Pople631Gs);
    cpu::scf::SCFConfig cfg = (xc == XCType::B3LYP)
                                  ? cpu::scf::b3lyp_config(integrals)
                                  : cpu::scf::pbe0_config(integrals);
    cfg.tol = 1e-9;
    cfg.max_iter = 300;
    cfg.xc_grad_project_translation = false;
    auto result = cpu::scf::scf_solve(integrals, cfg);
    EXPECT_TRUE(result.converged) << "hybrid SCF did not converge: " << stem;
    return cpu::scf::scf_gradient(frame, integrals, result, cfg);
}

// Compare an Atomiverse GradientResult to a PySCF dE/dR reference.
void check_against_pyscf(const std::string& stem, XCType xc,
                         const std::vector<std::array<double, 3>>& ref,
                         const char* label) {
    auto g = run_hybrid_gradient(stem, xc);
    const int natm = static_cast<int>(ref.size());
    double max_err = 0.0;
    for (int A = 0; A < natm; ++A) {
        for (int d = 0; d < 3; ++d) {
            // Atomiverse stores FORCE (= -dE/dR) in eV/Angstrom; convert
            // back to the PySCF gradient convention (dE/dR, Ha/Bohr).
            const double grad_atv = -ev_a_to_ha_bohr(g.forces(A, d));
            const double err = std::abs(grad_atv - ref[A][d]);
            max_err = std::max(max_err, err);
            EXPECT_NEAR(grad_atv, ref[A][d], 1e-5)
                << label << " " << stem << " atom=" << A << " axis=" << d;
        }
    }
    std::printf("[hybrid-grad] %s %s max|analytic-PySCF| = %.3e Ha/Bohr\n",
                label, stem.c_str(), max_err);
}

} // namespace

// ── ac-001: B3LYP analytic force vs PySCF ───────────────────────────────

TEST(HybridGradient, B3LYPForceMatchesPySCF) {
    for (const auto& rf : b3lyp_reference())
        check_against_pyscf(rf.stem, XCType::B3LYP, rf.grad, "B3LYP/6-31G*");
}

// ── ac-002: PBE0 analytic force vs PySCF ────────────────────────────────

TEST(HybridGradient, PBE0ForceMatchesPySCF) {
    for (const auto& rf : pbe0_reference())
        check_against_pyscf(rf.stem, XCType::PBE0, rf.grad, "PBE0/6-31G*");
}

// ── ac-003: B3LYP/PBE0 analytic force vs moving-grid finite difference ───
//
// Independent oracle: scf_finite_difference_gradient regenerates the grid at
// each displaced geometry (fixed_grid left null because
// xc_grad_project_translation stays false), so it probes the genuine
// moving-grid gradient without trusting any PySCF convention.

namespace {

void check_moving_grid_fd(const std::string& stem, XCType xc,
                          const char* label) {
    auto frame = make_molecule(stem);
    auto integrals = cpu::scf::build_integrals(frame, BasisSet::Pople631Gs);
    cpu::scf::SCFConfig cfg = (xc == XCType::B3LYP)
                                  ? cpu::scf::b3lyp_config(integrals)
                                  : cpu::scf::pbe0_config(integrals);
    cfg.tol = 1e-9;
    cfg.max_iter = 300;
    cfg.xc_grad_project_translation = false;

    auto result = cpu::scf::scf_solve(integrals, cfg);
    ASSERT_TRUE(result.converged) << label << " " << stem;
    auto analytic = cpu::scf::scf_gradient(frame, integrals, result, cfg);

    // Moving-grid FD: h = 1e-3 Bohr expressed as an Angstrom step.
    const double h_ang = 1e-3 / units::kAngstromInBohr;
    auto fd = cpu::scf::scf_finite_difference_gradient(
        frame, cfg, /*charge=*/0, h_ang, BasisSet::Pople631Gs);

    const int natm = static_cast<int>(frame.atoms.size());
    double max_err = 0.0;
    for (int A = 0; A < natm; ++A) {
        for (int d = 0; d < 3; ++d) {
            const double a = ev_a_to_ha_bohr(analytic.forces(A, d));
            const double f = ev_a_to_ha_bohr(fd.forces(A, d));
            max_err = std::max(max_err, std::abs(a - f));
            EXPECT_NEAR(a, f, 1e-5)
                << label << " " << stem << " atom=" << A << " axis=" << d;
        }
    }
    std::printf("[hybrid-grad] %s %s max|analytic-FD| = %.3e Ha/Bohr\n", label,
                stem.c_str(), max_err);
}

} // namespace

TEST(HybridGradient, B3LYPForceMatchesFiniteDifference) {
    check_moving_grid_fd("h2o", XCType::B3LYP, "B3LYP/6-31G*");
}

TEST(HybridGradient, PBE0ForceMatchesFiniteDifference) {
    check_moving_grid_fd("h2o", XCType::PBE0, "PBE0/6-31G*");
}

// ── ac-004: translational invariance sum_A F_A = 0 ──────────────────────

TEST(HybridGradient, HybridGradientIsTranslationallyInvariant) {
    for (XCType xc : {XCType::B3LYP, XCType::PBE0}) {
        const char* label = (xc == XCType::B3LYP) ? "B3LYP" : "PBE0";
        for (const char* stem : {"h2o", "ch4", "nh3"}) {
            auto g = run_hybrid_gradient(stem, xc);
            const int natm = static_cast<int>(g.forces.shape()[0]);
            double max_net = 0.0;
            for (int d = 0; d < 3; ++d) {
                double s = 0.0;
                for (int A = 0; A < natm; ++A)
                    s += ev_a_to_ha_bohr(g.forces(A, d));
                max_net = std::max(max_net, std::abs(s));
                EXPECT_LE(std::abs(s), 1e-8)
                    << label << " " << stem << " axis " << d
                    << ": net translational force must vanish";
            }
            std::printf("[hybrid-grad] %s %s ||sum_A F_A|| = %.3e Ha/Bohr\n",
                        label, stem, max_net);
        }
    }
}

// ── ac-005: B3LYP/6-31G* geometry optimization converges + matches PySCF ─
//
// Composition-only proof: b3lyp_config flows into the existing
// cpu::minimizer::Minimize + SCFEvaluator + LBFGS stack with no new
// production code. Water is started slightly off the B3LYP minimum; the
// optimizer must converge to the checked-in PySCF optimized geometry.

TEST(HybridGradient, B3LYPGeometryOptimizationMatchesPySCF) {
    // Start water slightly compressed from the B3LYP minimum so the
    // optimizer has real work to do.
    cpu::Frame frame(3);
    frame.atoms.type = {8, 1, 1};
    frame.atoms.x = {0.000000, 0.000000, 0.000000};
    frame.atoms.y = {0.000000, 0.745000, -0.745000};
    frame.atoms.z = {0.130000, -0.480000, -0.480000};
    frame.atoms.mass = {15.999, 1.008, 1.008};
    frame.box = Mat3::cubic(100.0);

    auto integrals = cpu::scf::build_integrals(frame, BasisSet::Pople631Gs);
    auto config = cpu::scf::b3lyp_config(integrals);
    config.tol = 1e-9;
    config.max_iter = 300;

    cpu::minimizer::Minimize opt;
    // alpha0 = 1.0 takes the full L-BFGS quasi-Newton step (the textbook
    // choice); the Armijo backtracking line search handles any overshoot.
    // A small alpha0 would cap every step at a tiny fraction of the Newton
    // step and stall progress with a ~(1-alpha0) geometric decay. Force
    // tolerances are tight so the optimizer drives onto the true minimum.
    cpu::minimizer::LBFGS step{/*memory=*/8,       /*alpha0=*/1.0,
                               /*ftol=*/0.01, // eV/Angstrom (~2e-4 Ha/Bohr)
                               /*max_ls=*/20,
                               /*rms_ftol=*/0.005,
                               /*max_dtol=*/0.005,
                               /*rms_dtol=*/0.002,
                               /*etol=*/1e-5};
    step.line_search();
    opt.add(std::move(step));
    opt.on(ATV::minimizer::Phase::force)
        .add(cpu::scf::SCFEvaluator{config, BasisSet::Pople631Gs});
    opt.build();
    const int max_iter = 60;
    const int iters = opt.run(frame, max_iter);

    // Converged means the optimizer stopped early on its convergence check.
    EXPECT_LT(iters, max_iter)
        << "B3LYP geometry optimization did not converge in " << max_iter
        << " iterations";

    // Optimized-geometry RMS vs PySCF reference. Water has a translational
    // gauge freedom (the optimizer does not pin the centroid); remove the
    // centroid offset on both geometries before the RMS comparison, then
    // also account for a possible mirror flip along x (all-zero by C2v).
    double com_atv[3] = {0, 0, 0}, com_ref[3] = {0, 0, 0};
    for (int A = 0; A < 3; ++A) {
        com_atv[0] += frame.atoms.x(A);
        com_atv[1] += frame.atoms.y(A);
        com_atv[2] += frame.atoms.z(A);
        com_ref[0] += kB3LYPOptGeom[A][0];
        com_ref[1] += kB3LYPOptGeom[A][1];
        com_ref[2] += kB3LYPOptGeom[A][2];
    }
    for (int d = 0; d < 3; ++d) {
        com_atv[d] /= 3.0;
        com_ref[d] /= 3.0;
    }

    double sum_sq = 0.0;
    for (int A = 0; A < 3; ++A) {
        const double ax = frame.atoms.x(A) - com_atv[0];
        const double ay = frame.atoms.y(A) - com_atv[1];
        const double az = frame.atoms.z(A) - com_atv[2];
        const double rx = kB3LYPOptGeom[A][0] - com_ref[0];
        const double ry = kB3LYPOptGeom[A][1] - com_ref[1];
        const double rz = kB3LYPOptGeom[A][2] - com_ref[2];
        sum_sq += (ax - rx) * (ax - rx) + (ay - ry) * (ay - ry)
                  + (az - rz) * (az - rz);
    }
    const double rms = std::sqrt(sum_sq / 3.0);
    std::printf("[hybrid-grad] ac-005 geometry RMS vs PySCF = %.3e Angstrom\n",
                rms);
    EXPECT_LE(rms, 1e-3)
        << "optimized B3LYP/6-31G* water geometry RMS exceeds 1e-3 Angstrom";

    // Final energy vs PySCF optimized reference. SCFEvaluator writes the
    // KS total energy (eV) into frame.scalars.energy; convert back to Hartree.
    const double e_final_ha = units::convert(
        frame.scalars.energy, units::Energy::eV, units::Energy::Hartree);
    std::printf("[hybrid-grad] ac-005 final E = %.10f Ha (PySCF %.10f Ha)\n",
                e_final_ha, kB3LYPOptEnergyHa);
    EXPECT_NEAR(e_final_ha, kB3LYPOptEnergyHa, 1e-6)
        << "optimized B3LYP/6-31G* water energy mismatch vs PySCF";
}

#endif // ATV_HAS_RHF
