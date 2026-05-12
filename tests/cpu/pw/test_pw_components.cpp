/// @file test_pw_components.cpp
/// @brief Unit tests for PW SCF components.
#include <gtest/gtest.h>
#include "atomiverse/core/box.h"
#include "atomiverse/core/vec3.h"
#include "atomiverse/cpu/context.h"
#include "atomiverse/cpu/component.h"
#include "atomiverse/cpu/pw/internal/g_vector.h"
#include "atomiverse/cpu/pw/internal/crystal.h"
#include "atomiverse/cpu/pw/pw_workspace.h"
#include "atomiverse/cpu/pw/pw_poisson.h"
#include "atomiverse/cpu/pw/pw_lda.h"
#include "atomiverse/cpu/pw/pw_hamiltonian.h"
#include "atomiverse/cpu/pw/pw_veff.h"
#include "atomiverse/cpu/pw/pw_davidson.h"
#include "atomiverse/cpu/pw/pw_density.h"
#include "atomiverse/cpu/pw/pw_convergence.h"

using namespace ATV;
using namespace ATV::cpu;
using namespace ATV::cpu::pw;
using namespace ATV::cpu::pw::internal;

// Helper: create a minimal PW test fixture
struct PWFixture {
    cpu::Context ctx;
    cpu::Frame frame{2};
    GVectorList gv;
    PWWorkspace ws;
    cpu::BuildContext bctx{ctx, frame};
    cpu::LaunchContext lctx{ctx, frame, 0, 0};

    PWFixture() {
        frame = build_si_diamond(5.43);
        Box box;
        box.cell = frame.box;
        gv = enumerate_g_vectors(box, Real(15), Vec3::zeros());

        ws.n_g = gv.count();
        ws.gvec = &gv;
        auto ng = static_cast<std::size_t>(ws.n_g);
        ws.rho_G.resize({ng});
        ws.rho_G_old.resize({ng});
        ws.V_H_G.resize({ng});
        ws.V_xc_G.resize({ng});
        ws.V_eff_G.resize({ng});
        ws.V_loc_G_real.resize({ng});
        ws.rho_r.resize({ng});
        ws.V_xc_r.resize({ng});
        ws.T_G.resize({ng});

        for (int i = 0; i < ws.n_g; ++i) {
            ws.rho_G[i] = (i == 0) ? std::complex<Real>(Real(4)/box.volume(), 0) : std::complex<Real>(0, 0);
            ws.V_loc_G_real[i] = (i == 0) ? Real(-0.5) : Real(0);
        }

        bctx.pw_workspace = &ws;
    }
};

TEST(PWComponents, Poisson) {
    PWFixture f;
    PW_Poisson p;
    p.build(f.bctx);
    p.launch(f.lctx);
    EXPECT_NEAR(std::abs(f.ws.V_H_G[0]), 0.0, 1e-12); // G=0 zeroed
}

TEST(PWComponents, VeffAssembly) {
    PWFixture f;
    for (int i = 0; i < f.ws.n_g; ++i) {
        f.ws.V_H_G[i] = std::complex<Real>(Real(0.1), 0);
        f.ws.V_xc_G[i] = std::complex<Real>(Real(-0.05), 0);
    }
    PW_Veff v;
    v.build(f.bctx);
    v.launch(f.lctx);
    for (int i = 0; i < f.ws.n_g; ++i) {
        Real expected = f.ws.V_loc_G_real[i] + Real(0.1) + Real(-0.05);
        EXPECT_NEAR(f.ws.V_eff_G[i].real(), expected, Real(1e-12));
    }
}

TEST(PWComponents, HamiltonianKinetic) {
    PWFixture f;
    PW_Hamiltonian h;
    h.build(f.bctx);
    h.launch(f.lctx);
    for (int i = 0; i < f.ws.n_g; ++i)
        EXPECT_GE(static_cast<double>(f.ws.T_G[i]), 0.0);
}

TEST(PWComponents, DavidsonEigenvalues) {
    PWFixture f;
    PW_Hamiltonian h;
    h.build(f.bctx);
    PW_Veff v;
    v.build(f.bctx);
    v.launch(f.lctx);
    h.launch(f.lctx);

    PW_Davidson d(4, Real(1e-8), 20);
    d.build(f.bctx);
    d.launch(f.lctx);
    EXPECT_EQ(f.ws.eigenvalues.size(), 4u);
    for (size_t i = 1; i < f.ws.eigenvalues.size(); ++i)
        EXPECT_LE(static_cast<double>(f.ws.eigenvalues[i-1]),
                  static_cast<double>(f.ws.eigenvalues[i]) + 1e-10);
}

TEST(PWComponents, ConvergenceCheck) {
    PWFixture f;
    PW_Convergence c(Real(1e-8), 0);
    c.build(f.bctx);

    for (int i = 0; i < f.ws.n_g; ++i) {
        f.ws.rho_G[i] = std::complex<Real>(Real(1e-7) * Real(i+1), 0);
        f.ws.rho_G_old[i] = std::complex<Real>(0, 0);
    }
    c.launch(f.lctx);
    EXPECT_FALSE(f.ws.converged);

    for (int i = 0; i < f.ws.n_g; ++i) {
        f.ws.rho_G[i] = std::complex<Real>(Real(1e-12), 0);
        f.ws.rho_G_old[i] = std::complex<Real>(0, 0);
    }
    c.launch(f.lctx);
    EXPECT_TRUE(f.ws.converged);
}

TEST(PWComponents, WarmupSuppression) {
    PWFixture f;
    PW_Convergence c(Real(1e-20), 3);
    c.build(f.bctx);
    for (int step = 0; step < 3; ++step) {
        f.lctx.step = step;
        c.launch(f.lctx);
        EXPECT_FALSE(f.ws.converged) << "warmup step " << step;
    }
}
