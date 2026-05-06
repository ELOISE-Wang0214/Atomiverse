/// @file cpu/scf/scf_solve.cpp
/// @brief Thin shim — constructs a transient `cpu::scf::SCFSolver` and
///        runs it. The SCF cycle body now lives in
///        `cpu/scf/scf_solver.cpp` (see Phase 3 of the driver-taxonomy
///        spec).
///
/// Retained as a transitional public entry point for `rhf_solve` /
/// `uhf_solve` / `uks_solve` and for legacy callers. New code should
/// construct an `SCFSolver` directly (or use `cpu::Evaluate +
/// cpu::scf::SCFEvaluator` for the outer-facade single-point pattern).

#include "atomiverse/cpu/scf/scf_solve.h"
#include "atomiverse/cpu/scf/scf_solver.h"

namespace ATV {
namespace cpu {
namespace scf {

SCFResult scf_solve(const IntegralPack& integrals, const SCFConfig& config,
                    const SCFInitialGuess* guess) {
    SCFSolver solver{config};
    solver.build(integrals);
    return solver.run(integrals, guess);
}

} // namespace scf
} // namespace cpu
} // namespace ATV
