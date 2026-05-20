/// @file cpu/scf/scf_solve.cpp
/// @brief Function-style SCF entry point — constructs a transient
///        `cpu::scf::SCFSolver` and runs it. The SCF cycle body lives in
///        `cpu/scf/scf_solver.cpp` (see Phase 3 of the driver-taxonomy
///        spec).
///
/// `scf_solve()` is the low-level functional API used by `rhf_solve` /
/// `uhf_solve` / `uks_solve` and by the finite-difference gradient path.
/// Outer-facade single-point energy uses `cpu::Evaluate +
/// cpu::scf::SCFEvaluator` instead.

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
