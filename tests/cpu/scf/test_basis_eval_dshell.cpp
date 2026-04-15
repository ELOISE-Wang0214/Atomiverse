/// @file tests/cpu/scf/test_basis_eval_dshell.cpp
/// @brief Consistency tests for cpu::scf::basis_eval with d-shells (l = 2).
///
/// Covers:
///   1) Finite-difference checks of the gradient and Hessian entries
///      for every Cartesian d AO (xx, xy, xz, yy, yz, zz).
///   2) Density-normalization smoke test through the full pipeline:
///      build 6-31G* integrals on water, run RHF, integrate ρ on a fine
///      Cartesian grid, expect ≈ n_electrons.
///
/// Failure of (1) flags algebra regressions in poly_{val,grad,hess};
/// failure of (2) flags an end-to-end break in the KS-grid path once
/// d-shells enter the AO set.

#ifdef ATV_HAS_RHF

#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <vector>

#include "atomiverse/core/mat3.h"
#include "atomiverse/cpu/frame.h"
#include "atomiverse/cpu/scf/basis_eval.h"
#include "atomiverse/cpu/scf/basis_set.h"
#include "atomiverse/cpu/scf/integrals.h"
#include "atomiverse/cpu/scf/scf_solve.h"

using namespace ATV;

namespace {

cpu::Frame make_water() {
    cpu::Frame f(3);
    f.atoms.x[0]=0.0; f.atoms.y[0]=0.0;    f.atoms.z[0]= 0.117; f.atoms.type[0]=8;
    f.atoms.x[1]=0.0; f.atoms.y[1]= 0.757; f.atoms.z[1]=-0.468; f.atoms.type[1]=1;
    f.atoms.x[2]=0.0; f.atoms.y[2]=-0.757; f.atoms.z[2]=-0.468; f.atoms.type[2]=1;
    f.box = Mat3::cubic(100.0);
    return f;
}

} // namespace

// ── 1) Finite-difference checks on a single d-shell. ──────────────────

TEST(BasisEvalDShell, GradientAndHessianMatchFiniteDifference) {
    // Use the first d-shell of oxygen in 6-31G* so we inherit a real
    // primitive spec rather than hand-crafting one.
    auto frame = make_water();
    auto integrals = cpu::scf::build_integrals(
        frame, cpu::scf::BasisSet::Pople631Gs);

    // Find a d-shell in the integrals.
    const cpu::scf::native::Shell* d_shell = nullptr;
    for (const auto& sh : integrals.shells) {
        if (sh.l == 2) { d_shell = &sh; break; }
    }
    ASSERT_NE(d_shell, nullptr) << "6-31G* on O should include at least one d shell";

    // Evaluate at a grid point offset from the shell center so none of the
    // Cartesian polynomial components vanish trivially.
    const double rx = d_shell->center[0] + 0.43;
    const double ry = d_shell->center[1] - 0.17;
    const double rz = d_shell->center[2] + 0.29;

    // Iterate over the 6 Cartesian AOs sharing this primitive spec.
    for (int cart = 0; cart < 6; ++cart) {
        cpu::scf::native::Shell sh = *d_shell;
        sh.cart = cart;

        double val, gx, gy, gz, hxx, hxy, hxz, hyy, hyz, hzz;
        cpu::scf::eval_basis_all(sh, rx, ry, rz, val,
                                 gx, gy, gz, hxx, hxy, hxz, hyy, hyz, hzz);

        // FD gradient
        const double h = 1e-5;
        double vpx = cpu::scf::eval_basis(sh, rx+h, ry, rz);
        double vmx = cpu::scf::eval_basis(sh, rx-h, ry, rz);
        double vpy = cpu::scf::eval_basis(sh, rx, ry+h, rz);
        double vmy = cpu::scf::eval_basis(sh, rx, ry-h, rz);
        double vpz = cpu::scf::eval_basis(sh, rx, ry, rz+h);
        double vmz = cpu::scf::eval_basis(sh, rx, ry, rz-h);
        double fd_gx = (vpx - vmx) / (2 * h);
        double fd_gy = (vpy - vmy) / (2 * h);
        double fd_gz = (vpz - vmz) / (2 * h);

        EXPECT_NEAR(gx, fd_gx, 1e-6) << "d-cart=" << cart << " gx";
        EXPECT_NEAR(gy, fd_gy, 1e-6) << "d-cart=" << cart << " gy";
        EXPECT_NEAR(gz, fd_gz, 1e-6) << "d-cart=" << cart << " gz";

        // FD diagonal Hessian from second-differences of eval_basis.
        double v0 = val;
        double fd_hxx = (vpx - 2 * v0 + vmx) / (h * h);
        double fd_hyy = (vpy - 2 * v0 + vmy) / (h * h);
        double fd_hzz = (vpz - 2 * v0 + vmz) / (h * h);
        EXPECT_NEAR(hxx, fd_hxx, 1e-4) << "d-cart=" << cart << " hxx";
        EXPECT_NEAR(hyy, fd_hyy, 1e-4) << "d-cart=" << cart << " hyy";
        EXPECT_NEAR(hzz, fd_hzz, 1e-4) << "d-cart=" << cart << " hzz";

        // FD off-diagonal Hessian: central five-point mix.
        double vpxpy = cpu::scf::eval_basis(sh, rx+h, ry+h, rz);
        double vmxmy = cpu::scf::eval_basis(sh, rx-h, ry-h, rz);
        double vpxmy = cpu::scf::eval_basis(sh, rx+h, ry-h, rz);
        double vmxpy = cpu::scf::eval_basis(sh, rx-h, ry+h, rz);
        double fd_hxy = (vpxpy - vpxmy - vmxpy + vmxmy) / (4 * h * h);
        EXPECT_NEAR(hxy, fd_hxy, 1e-4) << "d-cart=" << cart << " hxy";
    }
}

// ── 2) End-to-end density integration on a d-shell basis. ────────────

TEST(BasisEvalDShell, DensityIntegralRecoversElectronCountOn631GsWater) {
    auto frame = make_water();
    auto integrals = cpu::scf::build_integrals(
        frame, cpu::scf::BasisSet::Pople631Gs);
    auto config = cpu::scf::baseline_config(integrals);
    config.max_iter = 200;
    auto result = cpu::scf::scf_solve(integrals, config);
    ASSERT_TRUE(result.converged) << "RHF/6-31G* on water should converge";

    // Bounding-box Cartesian grid in Bohr.  Fine enough spacing to recover
    // ~10 electrons within a few 0.1 e absolute tolerance.
    const double pad = 6.0;
    const double step = 0.25;
    double xmin = integrals.atoms[0].pos[0], xmax = xmin;
    double ymin = integrals.atoms[0].pos[1], ymax = ymin;
    double zmin = integrals.atoms[0].pos[2], zmax = zmin;
    for (const auto& a : integrals.atoms) {
        xmin = std::min(xmin, a.pos[0]); xmax = std::max(xmax, a.pos[0]);
        ymin = std::min(ymin, a.pos[1]); ymax = std::max(ymax, a.pos[1]);
        zmin = std::min(zmin, a.pos[2]); zmax = std::max(zmax, a.pos[2]);
    }
    xmin -= pad; xmax += pad;
    ymin -= pad; ymax += pad;
    zmin -= pad; zmax += pad;
    const int nx = static_cast<int>(std::ceil((xmax - xmin) / step)) + 1;
    const int ny = static_cast<int>(std::ceil((ymax - ymin) / step)) + 1;
    const int nz = static_cast<int>(std::ceil((zmax - zmin) / step)) + 1;

    const auto& shells = integrals.shells;
    const int n_basis = static_cast<int>(shells.size());
    std::vector<double> phi(n_basis);

    // RHF: P_alpha already holds the full double-occupied density.
    const auto& P = result.P_alpha;

    double integral = 0.0;
    for (int ix = 0; ix < nx; ++ix) {
        double rx = xmin + ix * step;
        for (int iy = 0; iy < ny; ++iy) {
            double ry = ymin + iy * step;
            for (int iz = 0; iz < nz; ++iz) {
                double rz = zmin + iz * step;
                for (int mu = 0; mu < n_basis; ++mu)
                    phi[mu] = cpu::scf::eval_basis(shells[mu], rx, ry, rz);
                double rho = 0;
                for (int mu = 0; mu < n_basis; ++mu) {
                    double row = 0;
                    for (int nu = 0; nu < n_basis; ++nu)
                        row += P(mu, nu) * phi[nu];
                    rho += phi[mu] * row;
                }
                integral += rho;
            }
        }
    }
    integral *= step * step * step;

    // Water = 10 electrons.  0.25 Bohr grid / 6-Bohr pad → ~1 % accuracy.
    EXPECT_NEAR(integral, 10.0, 0.3)
        << "grid-integrated RHF/6-31G* density should recover 10 electrons";
}

// ── 3) KS smoke test on a d-shell basis. ─────────────────────────────

TEST(BasisEvalDShell, PBE_631Gs_WaterConvergesAndIntegratesNElectrons) {
    auto frame = make_water();
    auto integrals = cpu::scf::build_integrals(
        frame, cpu::scf::BasisSet::Pople631Gs);

    // Unrestricted KS / PBE on water — exercises both the d-shell basis_eval
    // path (through xc_grid_engine + GGA gradient) and the native::eri
    // general-l dispatch inside build_j_direct.
    auto config = cpu::scf::uks_config(
        integrals, integrals.n_occ(), integrals.n_occ(), cpu::scf::XCType::PBE);
    config.max_iter = 100;
    config.tol = 1e-8;

    auto result = cpu::scf::scf_solve(integrals, config);
    ASSERT_TRUE(result.converged) << "PBE/6-31G* on water should converge";

    // n_radial=50, n_angular=110 is SG-0-ish; alpha/beta electron count
    // recovery is limited to a few 1e-3.
    EXPECT_NEAR(result.integrated_n_alpha, static_cast<double>(config.n_alpha), 5e-3);
    EXPECT_NEAR(result.integrated_n_beta,  static_cast<double>(config.n_beta),  5e-3);
    EXPECT_LT(result.xc_energy, 0.0);
}

#endif // ATV_HAS_RHF
