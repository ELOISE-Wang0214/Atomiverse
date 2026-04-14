/// @file cuda/minimizer/driver.cuh
/// @brief Shared minimizer Spec and Driver types for CUDA minimizers.
#pragma once

#include "atomiverse/facade/driver.h"
#include "atomiverse/cuda/md/init.cuh"
#include "atomiverse/facade/stage.h"

namespace ATV {
namespace minimizer {

/// @brief CUDA minimizer specification with auto-added zeroing infrastructure.
class Spec : public ATV::DriverSpec<Phase> {
    using Base = ATV::DriverSpec<Phase>;
public:
    Spec() {
        Base::on(Phase::force).add(md::ZeroForces{}, -1000);
        Base::on(Phase::force).add(md::ZeroEnergy{}, -999);
    }
};

/// @brief Pure runtime CUDA minimizer driver executing a compiled plan.
class Driver : public ATV::Driver<Phase> {
    using Base = ATV::Driver<Phase>;
public:
    using Base::Base;
};

} // namespace minimizer
} // namespace ATV
