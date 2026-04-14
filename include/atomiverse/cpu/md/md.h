/// @file cpu/md/md.h
/// @brief cpu::md::MD driver alias with auto-added infrastructure components.
#pragma once

#include "atomiverse/core/device.h"
#include "atomiverse/cpu/backend_traits.h"
#include "atomiverse/facade/driver.h"
#include "atomiverse/facade/stage.h"
#include "atomiverse/cpu/md/init.h"

namespace ATV {
namespace cpu {
namespace md {

/// @brief CPU MD simulation driver.
///
/// Convenience class that auto-adds infrastructure components:
/// ZeroForces (pre_force, -1000), ZeroEnergy (pre_force, -999),
/// ZeroKinetic (end_of_step, -1000), KineticEnergy (end_of_step, -100).
///
/// Usage:
/// @code{.cpp}
/// cpu::md::MD sim;
/// sim.add(cpu::md::VelocityVerlet{0.001});
/// sim.on(ATV::md::Phase::force).add(cpu::md::NeighborList{6.8}, -10);
/// sim.build();
/// sim.run(frame, 1000);
/// @endcode
class MD : public Driver<ATV::md::Phase, cpu::device_tag> {
public:
    MD() {
        add(ZeroForces{});
        add(ZeroEnergy{});
        add(ZeroKinetic{});
        add(KineticEnergy{});
    }
};

} // namespace md
} // namespace cpu
} // namespace ATV
