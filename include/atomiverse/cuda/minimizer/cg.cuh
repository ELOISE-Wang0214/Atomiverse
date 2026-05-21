/// @file cuda/minimizer/cg.cuh
/// @brief minimizer::CG -- conjugate-gradient energy minimizer driver.
///
/// Automatically adds infrastructure components:
/// - ZeroForces before the force phase (order -1000)
/// - ZeroEnergy before the force phase (order -999)
///
/// Force evaluators (e.g., LJ) can be reused from the md namespace
/// by binding to minimizer::Phase::force.
///
/// @see minimizer::Phase for the minimizer phase sequence
/// @see Driver for the inherited fluent API
#pragma once

#include "atomiverse/cuda/md/init.cuh"
#include "atomiverse/facade/driver.h"
#include "atomiverse/facade/stage.h"

namespace ATV {
namespace minimizer {

/// @brief CG minimizer driver with auto-added zeroing infrastructure.
///
/// Automatically adds:
/// - ZeroForces at force phase (order -1000).
/// - ZeroEnergy at force phase (order -999).
///
/// @see md::MD for the molecular dynamics driver
class CG : public Driver<Phase> {
  public:
    /// @brief Construct CG driver with auto-added zeroing infrastructure.
    CG() {
        on(Phase::force).add(md::ZeroForces{}, -1000);
        on(Phase::force).add(md::ZeroEnergy{}, -999);
    }
};

} // namespace minimizer
} // namespace ATV
