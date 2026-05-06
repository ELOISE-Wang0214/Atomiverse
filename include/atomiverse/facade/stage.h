/// @file facade/stage.h
/// @brief Phase enums that define the staged execution skeleton for each
///        driver type. CUDA-free header.
#pragma once

#include <string>
#include <vector>

namespace ATV {

/// @brief MD simulation phases -- LAMMPS-style naming.
namespace md {

/// @brief Phase order for one molecular-dynamics step.
enum class Phase {
    initial_integrate, ///< First half-step: v += 0.5*dt*f/m, x += dt*v
    pre_force,         ///< Pre-force hooks (pair list rebuild, constraint setup)
    force,             ///< Force evaluation (evaluators accumulate here)
    post_force,        ///< Post-force corrections (constraints, long-range)
    final_integrate,   ///< Second half-step: v += 0.5*dt*f/m
    end_of_step,       ///< Per-step observers, recorders, controllers
    COUNT              ///< Sentinel -- total number of phases
};

/// @brief Human-readable names for the MD phases.
/// @param phase_tag Tag value used only for overload selection.
/// @return Phase names in the same order as the enum values, excluding `COUNT`.
inline std::vector<std::string> phase_names(Phase phase_tag) {
    (void)phase_tag;
    return {"initial_integrate", "pre_force", "force",
            "post_force", "final_integrate", "end_of_step"};
}

} // namespace md

/// @brief Energy minimizer phases.
namespace minimizer {

/// @brief Phase order for one minimization iteration.
enum class Phase {
    search_direction,  ///< Compute descent direction
    line_search,       ///< Line search along descent direction
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
    return {"search_direction", "line_search", "force",
            "convergence_check", "end_of_step"};
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

/// @brief Self-consistent field phases.
///
/// @note Legacy CUDA-only enum. CPU SCF is encapsulated by
///       `cpu::scf::SCFSolver` (an inner solver, not a Driver subclass)
///       and does not use this enum. This namespace will retire once
///       the CUDA SCF migration ships.
namespace scf {

/// @brief Phase order for one SCF iteration.
enum class Phase {
    build_hamiltonian,    ///< Construct effective Hamiltonian H(P) = H0 + alpha*P
    diagonalize,          ///< Eigensolve H -> eigenvalues + eigenvectors
    build_density,        ///< Rebuild density P = C_occ * C_occ^T
    convergence_check,    ///< Check ||P_new - P_old||_F < tol
    end_of_step,          ///< Per-iteration observers, recorders
    COUNT                 ///< Sentinel -- total number of phases
};

/// @brief Human-readable names for the SCF phases.
/// @param phase_tag Tag value used only for overload selection.
/// @return Phase names in the same order as the enum values, excluding `COUNT`.
inline std::vector<std::string> phase_names(Phase phase_tag) {
    (void)phase_tag;
    return {"build_hamiltonian", "diagonalize", "build_density",
            "convergence_check", "end_of_step"};
}

} // namespace scf
} // namespace ATV
