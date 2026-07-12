/// @file tests/cpu/md/test_write_sets.cpp
/// @brief Unit: every component that mutates x or v declares it in capabilities().
///
/// `Capability` is documented as the write-set — "what frame fields may launch() write".
/// Before constraints landed, that documentation was a lie for the entire integrator
/// family: VelocityVerlet rewrote every position and velocity in the system and reported
/// `Capability::None`. Nothing broke, because nothing read the write-set. That is exactly
/// what makes it dangerous — the first engine feature to trust the declaration (a
/// constraint asking "who moved my atoms?") would have been silently misinformed.
///
/// This test is the thing that keeps the declaration honest. It is deliberately dumb: it
/// does not check the physics, only that a component which is KNOWN to write x or v says
/// so. A new integrator that forgets fails here, at unit-test speed, and not three
/// months later inside whatever finally starts reading the bits.

#include <gtest/gtest.h>

#include "atomiverse/cpu/md/deform.h"
#include "atomiverse/cpu/md/fix_csvr.h"
#include "atomiverse/cpu/md/fix_langevin.h"
#include "atomiverse/cpu/md/fix_rattle.h"
#include "atomiverse/cpu/md/init.h"
#include "atomiverse/cpu/md/nose_hoover.h"
#include "atomiverse/cpu/md/nose_hoover_npt.h"
#include "atomiverse/cpu/md/nve_limit.h"
#include "atomiverse/cpu/md/velocity_verlet.h"

using namespace ATV;

namespace {
/// NPT's Config is an aggregate with no defaults -- named once, used twice.
cpu::md::NoseHooverNPT::Config npt_cfg() {
    return {/*T_target=*/300.0, /*tau_T=*/100.0, /*P_target=*/0.0,
            /*tau_P=*/1000.0, /*dt=*/1.0};
}
} // namespace

/// Integrators move atoms AND velocities. Both bits, every one of them.
TEST(WriteSets, IntegratorsDeclarePositionAndVelocity) {
    const auto both = Capability::Position | Capability::Velocity;

    EXPECT_TRUE(has_capability(cpu::md::VelocityVerlet{1.0}.capabilities(), both));
    EXPECT_TRUE(has_capability(cpu::md::NVELimit{1.0, 0.1}.capabilities(), both));
    EXPECT_TRUE(has_capability(cpu::md::NoseHooverNPT{npt_cfg()}.capabilities(),
                               both));
}

/// Velocity-only rescalers: they touch v and must not claim Position they never write.
TEST(WriteSets, VelocityRescalersDeclareVelocityOnly) {
    const auto nh = cpu::md::NoseHoover{300.0, 100.0, 1.0}.capabilities();
    EXPECT_TRUE(has_capability(nh, Capability::Velocity));
    EXPECT_FALSE(has_capability(nh, Capability::Position));

    const auto csvr = cpu::md::FixCSVR{300.0, 100.0, 1.0, 7u}.capabilities();
    EXPECT_TRUE(has_capability(csvr, Capability::Velocity));
    EXPECT_FALSE(has_capability(csvr, Capability::Position));

    const auto init = cpu::md::InitVelocity{300.0, 12345u}.capabilities();
    EXPECT_TRUE(has_capability(init, Capability::Velocity));
}

/// Deform moves atoms with the box; it does not touch velocities.
TEST(WriteSets, DeformDeclaresPositionOnly) {
    const auto d = cpu::md::Deform{cpu::md::Deform::Mode::erate, 1e-4, 1.0}.capabilities();
    EXPECT_TRUE(has_capability(d, Capability::Position));
    EXPECT_FALSE(has_capability(d, Capability::Velocity));
}

/// A Langevin thermostat is NOT a velocity writer — it adds friction + noise to the
/// FORCE and lets the integrator do the moving. Declaring Velocity here would be as
/// wrong as declaring None on VelocityVerlet, in the other direction.
TEST(WriteSets, LangevinDeclaresForceNotVelocity) {
    const auto l = cpu::md::FixLangevin{300.0, 300.0, 100.0, 1.0, 42u}.capabilities();
    EXPECT_TRUE(has_capability(l, Capability::Force));
    EXPECT_FALSE(has_capability(l, Capability::Velocity));
}

/// RATTLE corrects both stages: positions after the drift, velocities after the kick.
TEST(WriteSets, RattleDeclaresPositionAndVelocity) {
    cpu::md::FixRattle r{1.0};
    r.constrain_bond(0, 0.9572);
    EXPECT_TRUE(has_capability(r.capabilities(),
                               Capability::Position | Capability::Velocity));
}

/// The whole point, stated once: none of these may report an empty write-set.
TEST(WriteSets, NoStateMutatorReportsAnEmptyWriteSet) {
    const Capability caps[] = {
        cpu::md::VelocityVerlet{1.0}.capabilities(),
        cpu::md::NVELimit{1.0, 0.1}.capabilities(),
        cpu::md::NoseHoover{300.0, 100.0, 1.0}.capabilities(),
        cpu::md::NoseHooverNPT{npt_cfg()}.capabilities(),
        cpu::md::FixCSVR{300.0, 100.0, 1.0, 7u}.capabilities(),
        cpu::md::FixLangevin{300.0, 300.0, 100.0, 1.0, 42u}.capabilities(),
        cpu::md::InitVelocity{300.0, 12345u}.capabilities(),
        cpu::md::Deform{cpu::md::Deform::Mode::erate, 1e-4, 1.0}.capabilities(),
    };
    for (const auto c : caps)
        EXPECT_NE(c, Capability::None) << "a component mutates frame state but declares "
                                          "Capability::None — the write-set is a lie";
}
