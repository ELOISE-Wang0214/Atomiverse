/// @file cpu/response_driver.h
/// @brief cpu::ResponseDriver — finite-difference response (Hessian) driver.
///
/// `ResponseDriver` is the fourth CPU outer-driver category, sitting
/// alongside `EvaluateDriver` / `DynamicsDriver` / `OptimizeDriver`. It
/// is the category CLAUDE.md's driver-taxonomy section explicitly
/// sanctions as a genuinely distinct outer-loop semantic ("Response for
/// finite-difference Hessian").
///
/// Termination semantic: a fixed `6N` displacement budget — every
/// displacement must complete for the Hessian to be defined — followed
/// by one finalize pass. There is no convergence success path and no
/// early-stop request, so `ResponseDriver` performs **no**
/// `outer_status` polling; `facade/iteration_status.h` plays no role
/// in this category.
///
/// Loop (`run(frame)`):
///   - `setup_pass(frame)`
///   - `step_once(d)` for `d` in `0 … 6N-1`  (the displacement passes)
///   - `step_once(6N)`                       (the finalize pass)
///   - `teardown_pass()`
/// `N = frame.atoms.size()` and `6N` are derived inside `run()` — never a
/// `run()` argument and never a Config field.
///
/// Validation rule (build-time):
///   - >= 1 component carrying `Capability::Force` (the force evaluator)
///   - >= 1 component carrying `Capability::Hessian` (the accumulator)
///
/// Concrete user-facing class: `cpu::vib::Hessian` (in cpu/vib/hessian.h).
///
/// Inheritance: ResponseDriver → CPUDriver<vib::Phase> →
/// Driver<vib::Phase, cpu::device_tag>.
#pragma once

#include "atomiverse/cpu/backend_traits.h"
#include "atomiverse/cpu/component.h"
#include "atomiverse/cpu/cpu_driver.h"
#include "atomiverse/cpu/frame.h"
#include "atomiverse/facade/component_role.h"
#include "atomiverse/facade/stage.h"

namespace ATV {
namespace cpu {

/// @brief Finite-difference response (numerical-Hessian) driver base.
class ResponseDriver : public CPUDriver<ATV::vib::Phase> {
public:
    /// @brief Run the full `6N`-displacement finite-difference response loop.
    ///
    /// Derives `N = frame.atoms.size()` internally. Executes one Plan
    /// pass per displacement (`6N` passes total) plus a final
    /// finalize pass where the Symmetrize component runs. The frame is
    /// restored to its equilibrium geometry by the final pass.
    ///
    /// @param frame Host-resident system state. Coordinates are
    ///              displaced and restored in place; on return the frame
    ///              carries the equilibrium geometry again.
    void run(cpu::Frame& frame) {
        const int n_atoms = frame.atoms.size();
        const int n_disp  = 6 * n_atoms; ///< Two displacements per Cartesian coord.

        setup_pass(frame);
        for (int d = 0; d < n_disp; ++d) {
            // No outer_status reset / poll: the Response category has no
            // convergence success path and no early-stop request — every
            // displacement must complete for the Hessian to be defined.
            step_once(d);
        }
        // Final finalize pass: Symmetrize and any future diagonalizer run
        // here; DisplaceStep restores the equilibrium geometry.
        step_once(n_disp);
        teardown_pass();
    }

protected:
    void validate_components() const override {
        bool has_force_writer  = false;
        bool has_hessian_writer = false;
        for_each_pending([&](const cpu::ComponentBase& c) {
            if (has_capability(c.capabilities(), Capability::Force))
                has_force_writer = true;
            if (has_capability(c.capabilities(), Capability::Hessian))
                has_hessian_writer = true;
        });
        if (!has_force_writer)
            throw_missing("Component with Capability::Force");
        if (!has_hessian_writer)
            throw_missing("Component with Capability::Hessian");
    }
};

} // namespace cpu
} // namespace ATV
