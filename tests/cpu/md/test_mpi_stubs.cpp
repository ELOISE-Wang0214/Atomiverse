/// @file tests/cpu/md/test_mpi_stubs.cpp
/// @brief RED tests for the four Phase-1 no-op MPI components.
///
/// Spec: .claude/specs/parallel-md-engine.md § A.4 (components table).
/// Expected bindings (Phase enum as int):
///   HaloExchangeForward → pre_force @ -20
///   HaloExchangeReverse → post_force @ -1000
///   AtomMigration       → pre_force @ -30
///   ReduceEnergy        → end_of_step @ -50
///
/// In single-rank mode every launch() is a no-op: frame positions must be
/// byte-identical before and after.

#include <gtest/gtest.h>

#include <cstring>
#include <string>
#include <vector>

#include "atomiverse/cpu/context.h"
#include "atomiverse/cpu/frame.h"
#include "atomiverse/cpu/md/md.h"
#include "atomiverse/cpu/md/velocity_verlet.h"
#include "atomiverse/cpu/md/mpi/atom_migration.h"
#include "atomiverse/cpu/md/mpi/halo_exchange.h"
#include "atomiverse/cpu/md/mpi/reduce.h"
#include "atomiverse/core/binding.h"
#include "atomiverse/core/mat3.h"
#include "atomiverse/facade/stage.h"

using namespace ATV;

namespace {

cpu::Frame make_fcc_argon(int n_cells, double a) {
    int n = n_cells * n_cells * n_cells * 4;
    cpu::Frame frame(n);
    double basis[4][3] = {
        {0.0, 0.0, 0.0}, {0.5, 0.5, 0.0},
        {0.5, 0.0, 0.5}, {0.0, 0.5, 0.5},
    };
    int idx = 0;
    for (int ix = 0; ix < n_cells; ++ix)
        for (int iy = 0; iy < n_cells; ++iy)
            for (int iz = 0; iz < n_cells; ++iz)
                for (int b = 0; b < 4; ++b) {
                    frame.atoms.x[idx] = (ix + basis[b][0]) * a;
                    frame.atoms.y[idx] = (iy + basis[b][1]) * a;
                    frame.atoms.z[idx] = (iz + basis[b][2]) * a;
                    frame.atoms.mass[idx] = 39.948;
                    frame.atoms.type[idx] = 18;
                    ++idx;
                }
    frame.box = Mat3::cubic(n_cells * a);
    return frame;
}

bool has_binding(const std::vector<DefaultBinding>& bindings, int phase, int order) {
    if (bindings.size() != 1) return false;
    return bindings[0].phase == phase && bindings[0].order == order;
}

} // namespace

// ── bindings ─────────────────────────────────────────────────────────

TEST(MpiStubBindings, HaloExchangeForwardBindings) {
    cpu::md::mpi::HaloExchangeForward hx{};
    auto b = hx.default_bindings();
    EXPECT_TRUE(has_binding(b, static_cast<int>(ATV::md::Phase::pre_force), -20))
        << "HaloExchangeForward must bind pre_force @ -20";
}

TEST(MpiStubBindings, HaloExchangeReverseBindings) {
    cpu::md::mpi::HaloExchangeReverse hx{};
    auto b = hx.default_bindings();
    EXPECT_TRUE(has_binding(b, static_cast<int>(ATV::md::Phase::post_force), -1000))
        << "HaloExchangeReverse must bind post_force @ -1000";
}

TEST(MpiStubBindings, AtomMigrationBinding) {
    cpu::md::mpi::AtomMigration am{};
    auto b = am.default_bindings();
    EXPECT_TRUE(has_binding(b, static_cast<int>(ATV::md::Phase::pre_force), -30))
        << "AtomMigration must bind pre_force @ -30";
}

TEST(MpiStubBindings, ReduceEnergyBinding) {
    cpu::md::mpi::ReduceEnergy re{};
    auto b = re.default_bindings();
    EXPECT_TRUE(has_binding(b, static_cast<int>(ATV::md::Phase::end_of_step), -50))
        << "ReduceEnergy must bind end_of_step @ -50";
}

// ── single-rank no-op: positions byte-identical ──────────────────────

TEST(MpiStubBindings, SingleRankLaunchNoOp) {
    auto frame = make_fcc_argon(3, 5.26);
    const int n = frame.atoms.size();

    // Snapshot positions before.
    std::vector<double> x_before(n), y_before(n), z_before(n);
    std::memcpy(x_before.data(), frame.atoms.x.data(), n * sizeof(double));
    std::memcpy(y_before.data(), frame.atoms.y.data(), n * sizeof(double));
    std::memcpy(z_before.data(), frame.atoms.z.data(), n * sizeof(double));

    // Wire a bare MD driver with the four stubs only (plus auto infra + a cheap
    // integrator to have a full step). Drive 10 steps, verify positions unchanged
    // *by the stubs* — since VelocityVerlet only moves if velocities are nonzero.
    {
        cpu::md::MD sim;
        sim.add(cpu::md::VelocityVerlet{1e-6});  // negligible dt
        sim.add(cpu::md::mpi::HaloExchangeForward{});
        sim.add(cpu::md::mpi::HaloExchangeReverse{});
        sim.add(cpu::md::mpi::AtomMigration{});
        sim.add(cpu::md::mpi::ReduceEnergy{});
        sim.build();
        sim.run(frame, 10);
    }

    // With no forces or velocities, the stubs must be bit-identical no-ops.
    EXPECT_EQ(0, std::memcmp(x_before.data(), frame.atoms.x.data(),
                             n * sizeof(double)))
        << "MPI stubs perturbed frame.atoms.x in single-rank mode";
    EXPECT_EQ(0, std::memcmp(y_before.data(), frame.atoms.y.data(),
                             n * sizeof(double)))
        << "MPI stubs perturbed frame.atoms.y in single-rank mode";
    EXPECT_EQ(0, std::memcmp(z_before.data(), frame.atoms.z.data(),
                             n * sizeof(double)))
        << "MPI stubs perturbed frame.atoms.z in single-rank mode";
}

// ── class-name hygiene ───────────────────────────────────────────────

TEST(MpiStubBindings, ClassNamesAreDistinct) {
    cpu::md::mpi::HaloExchangeForward fwd{};
    cpu::md::mpi::HaloExchangeReverse rev{};
    const std::string n_fwd = fwd.name();
    const std::string n_rev = rev.name();
    EXPECT_NE(n_fwd, n_rev);
    EXPECT_NE(n_fwd.find("halo"), std::string::npos)
        << "HaloExchangeForward::name() should contain 'halo': got '" << n_fwd << "'";
    EXPECT_NE(n_rev.find("halo"), std::string::npos)
        << "HaloExchangeReverse::name() should contain 'halo': got '" << n_rev << "'";
}
