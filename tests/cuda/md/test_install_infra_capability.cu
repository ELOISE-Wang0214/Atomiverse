/// @file tests/cuda/md/test_install_infra_capability.cu
/// @brief Capability-typed auto-injection tests for `cuda::md::MD::install_infra_`.
///
/// Covers `cuda-md-review-cleanup-01-component-capability` ac-006/ac-008:
/// - `cuda::md::MD::install_infra_` consults `reads()` / `needs_neighbor_list()`
///   / `capabilities()` instead of `name() == "..."` string proxies and
///   `dynamic_cast<const cuda::Pair*>`.
/// - A pending component carrying `(reads() & Capability::Stress) != None`
///   auto-injects `cuda::Virial`. This is the key behavior change
///   vs the old `role() == Barostat` proxy: `cuda::md::NoseHooverNPT`
///   reports `role() == Integrator` but consumes Stress, so the new
///   typed check picks it up where the old string/role proxy could not.
/// - A pending component reporting `needs_neighbor_list() == true`
///   without a registered `cuda::NeighborList` throws at `build()`.
///
/// These tests *cannot* introspect `pending_` directly (it's private to
/// the Backend). Instead they observe injection through externally
/// visible effects: a successful `build()` + `run()` for the NPT path
/// (would-throw if Virial wasn't auto-injected because NPT's
/// Stress-reading branch would fault) and a `runtime_error` for the
/// missing-NeighborList path.

#include <cmath>
#include <cstdlib>
#include <gtest/gtest.h>
#include <string>

#include "atomiverse/core/mat3.h"
#include "atomiverse/cpu/frame.h"
#include "atomiverse/cuda/md/md.cuh"
#include "atomiverse/cuda/md/nose_hoover_npt.cuh"
#include "atomiverse/cuda/md/velocity_verlet.cuh"
#include "atomiverse/cuda/neighbor_list.cuh"
#include "atomiverse/cuda/observer/stress.cuh"
#include "atomiverse/cuda/potentials/lj.cuh"
#include "atomiverse/facade/component_role.h"

using namespace ATV;

namespace {

/// Same builder as `test_cuda_md_components.cu::make_lj_system`.
cpu::Frame make_lj_system(int n_side, double density) {
    int n = n_side * n_side * n_side;
    double box_len = std::cbrt(double(n) / density);
    double spacing = box_len / double(n_side);
    cpu::Frame frame(n);
    int idx = 0;
    for (int ix = 0; ix < n_side; ++ix)
        for (int iy = 0; iy < n_side; ++iy)
            for (int iz = 0; iz < n_side; ++iz) {
                frame.atoms.x[idx] = (double(ix) + 0.5) * spacing;
                frame.atoms.y[idx] = (double(iy) + 0.5) * spacing;
                frame.atoms.z[idx] = (double(iz) + 0.5) * spacing;
                frame.atoms.mass[idx] = 1.0;
                frame.atoms.type[idx] = 0;
                ++idx;
            }
    frame.box = Mat3::cubic(box_len);
    std::srand(42);
    for (int i = 0; i < n; ++i) {
        frame.atoms.vx[i] = 1e-3 * (double(std::rand()) / RAND_MAX - 0.5);
        frame.atoms.vy[i] = 1e-3 * (double(std::rand()) / RAND_MAX - 0.5);
        frame.atoms.vz[i] = 1e-3 * (double(std::rand()) / RAND_MAX - 0.5);
    }
    double mvx = 0, mvy = 0, mvz = 0;
    for (int i = 0; i < n; ++i) {
        mvx += frame.atoms.vx[i];
        mvy += frame.atoms.vy[i];
        mvz += frame.atoms.vz[i];
    }
    mvx /= n; mvy /= n; mvz /= n;
    for (int i = 0; i < n; ++i) {
        frame.atoms.vx[i] -= mvx;
        frame.atoms.vy[i] -= mvy;
        frame.atoms.vz[i] -= mvz;
    }
    return frame;
}

cuda::md::NoseHooverNPT::Config default_npt_cfg() {
    cuda::md::NoseHooverNPT::Config cfg{};
    cfg.T_target = 1.0;
    cfg.tau_T = 100.0;
    cfg.P_target = 0.0;
    cfg.tau_P = 1000.0;
    cfg.dt = 0.001;
    return cfg;
}

} // namespace

// ── ac-001 + ac-006: capability-typed virtuals exist and route ──────

TEST(InstallInfraCapability, ComponentBaseDefaultsAreNone) {
    // Smoke: a default-constructed bare ComponentBase subclass (here we
    // use cuda::md::VelocityVerlet, which doesn't override reads /
    // needs_neighbor_list) returns the None / false defaults.
    cuda::md::VelocityVerlet vv{0.001};
    EXPECT_EQ(vv.reads(), Capability::None);
    EXPECT_FALSE(vv.needs_neighbor_list());
}

TEST(InstallInfraCapability, PairDeclaresNeedsNeighborList) {
    // PairLJCut (cuda::LJ) inherits from cuda::Pair, which overrides
    // needs_neighbor_list()==true. The OLD behavior was to detect pair
    // styles via dynamic_cast<const cuda::Pair*>; the new behavior is
    // a typed virtual on ComponentBase.
    cuda::LJ lj{1.0, 1.0, 3.0};
    EXPECT_TRUE(lj.needs_neighbor_list());
}

TEST(InstallInfraCapability, NPTDeclaresReadsStress) {
    // cuda::md::NoseHooverNPT is the canonical Integrator that *also*
    // acts as a barostat — it consumes Stress to drive ε̇ updates. With
    // the old install_infra_ proxy (role==Barostat) this was invisible
    // and required users to add cuda::Virial by hand. The new
    // reads() virtual exposes the dependency directly.
    cuda::md::NoseHooverNPT npt{default_npt_cfg()};
    EXPECT_EQ(npt.role(), ComponentRole::Integrator);
    EXPECT_NE((npt.reads() & Capability::Stress), Capability::None);
}

TEST(InstallInfraCapability, VirialDeclaresStressCapability) {
    cuda::Virial vs;
    EXPECT_EQ(vs.role(), ComponentRole::Observer);
    EXPECT_NE((vs.capabilities() & Capability::Stress), Capability::None);
}

// ── ac-008: end-to-end auto-injection for the NPT/Stress path ───────

TEST(InstallInfraCapability, NPTAutoInjectsVirial) {
    // Before the rewrite, `cuda::md::NoseHooverNPT` had role==Integrator
    // and install_infra_ checked role==Barostat — so this configuration
    // would NOT auto-inject Virial and the NPT stress-feedback
    // branch would fault. After the rewrite, install_infra_ checks
    // `reads() & Stress != None`, so NPT pulls in Virial
    // automatically.
    auto frame = make_lj_system(4, 0.8);
    const double cutoff = 3.0;
    cuda::md::MD sim;
    sim.set_neighbor(cuda::NeighborList{cutoff});
    sim.set_pair(cuda::LJ{1.0, 1.0, cutoff});
    sim.set_integrator(cuda::md::NoseHooverNPT{default_npt_cfg()});
    EXPECT_NO_THROW(sim.build());
    EXPECT_NO_THROW(sim.run(frame, 1));
}

// ── ac-006: missing-NeighborList detection via needs_neighbor_list ──

TEST(InstallInfraCapability, MissingNeighborThrowsViaTypedQuery) {
    // The old install_infra_ identified pair evaluators via
    // dynamic_cast<const cuda::Pair*>. The new install_infra_ uses
    // ComponentBase::needs_neighbor_list(). Behavior is identical;
    // this test guards the typed contract.
    cuda::md::MD sim;
    sim.set_pair(cuda::LJ{1.0, 1.0, 3.0});
    sim.set_integrator(cuda::md::VelocityVerlet{0.001});
    EXPECT_THROW(
        {
            try {
                sim.build();
            } catch (const std::runtime_error& e) {
                EXPECT_NE(std::string(e.what()).find("set_neighbor"),
                          std::string::npos);
                throw;
            }
        },
        std::runtime_error);
}
