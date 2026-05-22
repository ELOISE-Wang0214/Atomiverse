/// @file cuda/scf/internal/rhf_density.cuh
/// @brief internal::RHFDensity -- back-transform C, rebuild P, compute RHF energy.
///
/// **Internal SCF helper.** Does NOT inherit `cuda::ComponentBase` and is
/// never registered into a `Driver` Plan. Owned by `cuda::scf::SCFSolver`,
/// which drives its `build` / `step` / `destroy` lifecycle directly.
///
/// In `step()`:
///   1. C = X * C' (back-transform from orthogonal basis)
///   2. P_old = P; P = 2 * C_occ * C_occ^T (RHF closed-shell)
///   3. E = 0.5 * Tr[(Hcore + F) * P] + V_nuc
///   4. Write E to frame.scalars->energy
///
/// @see cuda::scf::SCFSolver for the inner solver that owns this helper
/// @see internal::RHFFock for the Fock matrix builder
#pragma once

#include "atomiverse/cuda/component.cuh"
#include "atomiverse/cuda/scf/scf_workspace.cuh"

namespace ATV {
namespace cuda {
namespace scf {
namespace internal {

#ifdef ATV_HAS_RHF

/// @brief Back-transform MO coefficients, rebuild density, compute RHF energy.
class RHFDensity {
  public:
    /// @brief Construct with occupation count.
    /// @param n_occ Number of doubly-occupied orbitals [dimensionless].
    explicit RHFDensity(int n_occ) : n_occ_(n_occ) {}

    /// @brief Allocate cuBLAS handle and the Hcore+F scratch buffer.
    /// @param ctx Build context (CUDA stream + device frame).
    /// @param ws SCF electronic-structure workspace.
    void build(const BuildContext& ctx, SCFWorkspace& ws);

    /// @brief Back-transform C, rebuild P, compute and write RHF energy.
    /// @param ctx Launch context (CUDA stream + device frame + step index).
    /// @param ws SCF electronic-structure workspace.
    void step(const LaunchContext& ctx, SCFWorkspace& ws);

    /// @brief Free cuBLAS handle and scratch buffer.
    /// @param ctx Destroy context (CUDA stream + device frame).
    /// @param ws SCF electronic-structure workspace.
    void destroy(const DestroyContext& ctx, SCFWorkspace& ws);

  private:
    int n_occ_;                     ///< Doubly-occupied orbital count.
    void* cublas_handle_ = nullptr; ///< cuBLAS handle.
    double* d_HplusF_ =
        nullptr; ///< Temporary Hcore+F for energy [device, n*n].
};

#endif // ATV_HAS_RHF

} // namespace internal
} // namespace scf
} // namespace cuda
} // namespace ATV
