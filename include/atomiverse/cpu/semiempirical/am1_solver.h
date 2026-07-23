/// @file cpu/semiempirical/am1_solver.h
/// @brief AM1 concrete method — public members + Frame build/run.
#pragma once

#include <memory>
#include <utility>

#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/semiempirical/internal/nddo_solver.h"

namespace ATV {
namespace cpu {
namespace semiempirical {

using AM1Result = NDDOResult;

/// @brief AM1 concrete semiempirical method (RHF pattern: members + run).
///
/// Owns an internal NDDO engine. Parameters are public members — there is no
/// public `NDDOConfig` / `NDDOSolver` bag. `AM1Evaluator` wraps this type.
class AM1 {
  public:
    /// @brief Total molecular charge [e].
    int charge = 0;
    /// @brief Spin multiplicity 2S+1 (dimensionless).
    int multiplicity = 1;
    /// @brief SCF iteration cap (dimensionless).
    int max_iter = 200;
    /// @brief Density residual tolerance (dimensionless).
    double density_tol = 1.0e-8;
    /// @brief Energy residual tolerance [eV].
    double energy_tol = 1.0e-10;

    AM1();
    ~AM1();
    AM1(AM1&&) noexcept;
    AM1& operator=(AM1&&) noexcept;
    /// @brief Copy parameters only (does not copy SCF engine state).
    AM1(const AM1& other);
    AM1& operator=(const AM1& other);

    void build(const cpu::Frame& reference);
    AM1Result& run(const cpu::Frame& frame);
    AM1Result& compute_forces(cpu::Frame& frame, double step_angstrom = 1.0e-4);

    const AM1Result& last_result() const;
    bool built() const noexcept;
    bool has_density() const noexcept;
    const NDDOSystem& system() const;
    const NDDOScfState& state() const;

  private:
    void sync_engine_();
    std::unique_ptr<NDDOSolver> engine_;
};

/// @brief Historical alias of the concrete AM1 method.
using AM1Solver = AM1;

/// @brief Transitional parameter bag (evaluator / task / python).
struct AM1Config {
    int charge = 0;
    int multiplicity = 1;
    int max_iter = 200;
    double density_tol = 1.0e-8;
    double energy_tol = 1.0e-10;
};

} // namespace semiempirical
} // namespace cpu
} // namespace ATV
