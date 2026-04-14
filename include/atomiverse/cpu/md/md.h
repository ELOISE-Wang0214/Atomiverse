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

/// @brief CPU MD specification with auto-added infrastructure components.
///
/// Auto-adds infrastructure components:
/// ZeroForces (pre_force, -1000), ZeroEnergy (pre_force, -999),
/// ZeroKinetic (end_of_step, -1000), KineticEnergy (end_of_step, -100).
///
/// Usage:
/// @code{.cpp}
/// cpu::md::Spec spec;
/// spec.add(cpu::md::VelocityVerlet{0.001});
/// spec.on(ATV::md::Phase::force).add(cpu::md::NeighborList{6.8}, -10);
/// auto plan = spec.compile();
/// cpu::md::Driver driver(std::move(plan));
/// driver.run(frame, 1000);
/// @endcode
class Spec : public ATV::DriverSpec<ATV::md::Phase, cpu::device_tag> {
    using Base = ATV::DriverSpec<ATV::md::Phase, cpu::device_tag>;
public:
    Spec() {
        Base::add(ZeroForces{});
        Base::add(ZeroEnergy{});
        Base::add(ZeroKinetic{});
        Base::add(KineticEnergy{});
    }
};

/// @brief Pure runtime CPU MD driver executing a compiled plan.
class Driver : public ATV::Driver<ATV::md::Phase, cpu::device_tag> {
    using Base = ATV::Driver<ATV::md::Phase, cpu::device_tag>;
public:
    using Base::Base;
};

/// @brief Transitional legacy CPU MD facade.
class MD : public ATV::detail::LegacyDriverFacade<Spec, Driver> {};

} // namespace md
} // namespace cpu
} // namespace ATV
