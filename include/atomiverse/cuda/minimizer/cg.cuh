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

#include "atomiverse/cuda/minimizer/driver.cuh"

namespace ATV {
namespace minimizer {

/// @brief Transitional conjugate-gradient minimizer facade.
///
/// @see md::MD for the molecular dynamics driver
class CG : public ATV::detail::LegacyDriverFacade<Spec, Driver> {};

} // namespace minimizer
} // namespace ATV
