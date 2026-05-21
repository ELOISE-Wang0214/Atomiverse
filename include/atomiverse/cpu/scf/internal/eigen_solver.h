/// @file cpu/scf/internal/eigen_solver.h
/// @brief internal::EigenSolver — generalized symmetric eigensolve (F C = S C ε).
///
/// **Internal SCF helper.** Does NOT inherit `cpu::ComponentBase`.
///
/// Thin stateful wrapper around `diagonalize_fock_generalized` (LAPACK
/// dsygvd via xtensor-blas). Held as a member of `cpu::scf::SCFSolver`
/// so future LAPACK workspace pre-allocation can land in `prepare()`
/// without touching call sites.
#pragma once

#include <utility>

#include "atomiverse/cpu/scf/scf_utils.h"

namespace ATV {
namespace cpu {
namespace scf {
namespace internal {

class EigenSolver {
  public:
    EigenSolver() = default;

    /// @brief Pre-allocate workspaces sized for `n_basis`. Currently a
    /// no-op (xtensor-blas allocates LAPACK workspaces internally per
    /// call); reserved for Phase 4 in-place dsygvd integration.
    void prepare(int /*n_basis*/) {}

    /// @brief Solve F C = S C ε.
    /// @return (eigenvalues [n], eigenvectors [n×n]).
    std::pair<RealVector, RealMatrix> step(const RealMatrix& F,
                                           const RealMatrix& S) const {
        return diagonalize_fock_generalized(F, S);
    }
};

} // namespace internal
} // namespace scf
} // namespace cpu
} // namespace ATV
