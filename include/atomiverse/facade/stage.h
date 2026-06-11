/// @file facade/stage.h
/// @brief Phase enums that define the staged execution skeleton for each
///        driver type. CUDA-free header.
#pragma once

#include <string>
#include <vector>

namespace ATV {

/// @brief MD simulation phases -- LAMMPS-style naming.
namespace md {

/// @brief Phase order for one molecular-dynamics step — LAMMPS-Verlet-granular
///        so the whole timestep is legible from the enum (see
///        `.claude/notes/cpu-md-execution-model.md`). Every phase is always
///        walked even when empty; an empty phase is an explicit "nothing
///        happens here" for this run, not an absence to infer. Within a phase
///        components order by a small integer (negative = infrastructure first,
///        0 = the op, positive = consumers) so a new op slots in without adding
///        a phase.
enum class Phase {
    initial_integrate, ///< First half-step: v += 0.5*dt*f/m, x += dt*v
    post_integrate,    ///< Hook right after kick/drift (usually empty)
    pre_exchange,      ///< Re-home atoms that left this rank (atom migration)
    pre_neighbor,      ///< Build ghost atoms (halo send) before the neighbor list
    neighbor,          ///< (Re)build the cell / pair neighbor list
    pre_force,         ///< Zero force/energy/stress accumulators before eval
    force,             ///< Force/energy/virial accumulation (pair + bonded)
    post_force,        ///< Reverse ghost-force comm, virial finalize, constraints
    final_integrate,   ///< Second half-step: v += 0.5*dt*f/m
    end_of_step,       ///< Global reductions, thermostat/barostat, observers, I/O
    COUNT              ///< Sentinel -- total number of phases
};

/// @brief Human-readable names for the MD phases.
/// @param phase_tag Tag value used only for overload selection.
/// @return Phase names in the same order as the enum values, excluding `COUNT`.
inline std::vector<std::string> phase_names(Phase phase_tag) {
    (void)phase_tag;
    return {"initial_integrate", "post_integrate", "pre_exchange",
            "pre_neighbor",      "neighbor",       "pre_force",
            "force",             "post_force",     "final_integrate",
            "end_of_step"};
}

} // namespace md

/// @brief Energy minimizer phases.
namespace minimizer {

/// @brief Phase order for one minimization iteration.
///
/// @note `line_search` is a reserved slot with no current binder: the
///       built-in optimizers (LBFGS Armijo backtracking) implement their
///       line search inline across successive driver iterations, because
///       each trial energy needs a full `force` pass. The slot stays for
///       future exact / dedicated line-search components.
enum class Phase {
    search_direction,  ///< Compute descent direction
    line_search,       ///< Line search along descent direction (reserved)
    force,             ///< Force evaluation (shared with MD evaluators)
    convergence_check, ///< Check stopping criteria
    end_of_step,       ///< Per-iteration observers, recorders
    COUNT              ///< Sentinel -- total number of phases
};

/// @brief Human-readable names for the minimizer phases.
/// @param phase_tag Tag value used only for overload selection.
/// @return Phase names in the same order as the enum values, excluding `COUNT`.
inline std::vector<std::string> phase_names(Phase phase_tag) {
    (void)phase_tag;
    return {"search_direction", "line_search", "force", "convergence_check",
            "end_of_step"};
}

} // namespace minimizer

/// @brief One-shot evaluation phases.
///
/// Used by `cpu::EvaluateDriver`. There is no `end_of_step` — Evaluate
/// does not have a step concept. The three phases give callers slots to
/// place pre-evaluation setup (zeroing forces / preparing buffers),
/// the evaluator components themselves, and post-evaluation
/// finalization (recorders, observers).
namespace evaluate {

/// @brief Phase order for one Evaluate pass.
enum class Phase {
    prepare,  ///< Pre-evaluation setup (Zero*, buffer prep)
    evaluate, ///< Evaluator components run here
    finalize, ///< Post-evaluation finalization (recorders, observers)
    COUNT     ///< Sentinel -- total number of phases
};

/// @brief Human-readable names for the Evaluate phases.
/// @param phase_tag Tag value used only for overload selection.
/// @return Phase names in the same order as the enum values, excluding `COUNT`.
inline std::vector<std::string> phase_names(Phase phase_tag) {
    (void)phase_tag;
    return {"prepare", "evaluate", "finalize"};
}

} // namespace evaluate

// NOTE: there is deliberately no `scf::Phase`. SCF — gaussian-basis
// (`cpu::scf::SCFSolver`, `cuda::scf::SCFSolver`) and plane-wave
// (`cpu::pw::PWSolver`) — is an *inner* solver behind an Evaluator
// component, never an outer driver category, so it carries no outer
// phase enum. The historical `scf::Phase` was deleted with the
// `cpu::pw` migration onto the PWSolver/PWEvaluator shape.

/// @brief Vibrational-analysis (Response category) phases.
///
/// Used by `cpu::ResponseDriver` / `cpu::vib::Hessian`. The Response
/// category runs a fixed `6N` displacement budget (no convergence /
/// stop-request polling) followed by one finalize pass. The phase
/// ordering reflects what runs within a single displacement pass
/// (`setup` → `displace` → `force` → `accumulate`) plus the
/// finalize-only phases (`diagonalize`, `thermochem`, `end_of_step`).
///
/// @note `force` is the slot where the user-registered force evaluator
///       runs — exactly as AIMD registers an evaluator into `md::Phase::force`.
///       `thermochem` is reserved for the thermochemistry link and
///       carries no component in the numerical-Hessian link.
namespace vib {

/// @brief Phase order for one Response-category pass.
enum class Phase {
    setup, ///< Cache the equilibrium geometry (ReferenceGeometry).
    displace, ///< Apply +/-h displacement for the current coordinate (DisplaceStep).
    force,      ///< Force evaluation (user-registered evaluator runs here).
    accumulate, ///< Central-difference Hessian-row assembly (HessianAccumulator).
    diagonalize, ///< Symmetrization / future diagonalization (Symmetrize).
    thermochem,  ///< Reserved — thermochemistry link; no component here.
    end_of_step, ///< Per-pass observers / recorders.
    COUNT        ///< Sentinel -- total number of phases.
};

/// @brief Human-readable names for the vibrational-analysis phases.
/// @param phase_tag Tag value used only for overload selection.
/// @return Phase names in the same order as the enum values, excluding `COUNT`.
inline std::vector<std::string> phase_names(Phase phase_tag) {
    (void)phase_tag;
    return {"setup",       "displace",   "force",      "accumulate",
            "diagonalize", "thermochem", "end_of_step"};
}

} // namespace vib
} // namespace ATV
