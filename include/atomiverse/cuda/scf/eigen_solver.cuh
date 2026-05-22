/// @file cuda/scf/eigen_solver.cuh
/// @brief EigenSolver -- dense symmetric eigendecomposition via cuSOLVER.
///
/// Diagonalizes the effective Hamiltonian H each SCF iteration:
///   H * C = C * diag(eig)
///
/// Uses cusolverDnSsyevd (float) or cusolverDnDsyevd (double) depending
/// on compile-time double precision.
///
/// @see scf::Phase::diagonalize
/// @see Hamiltonian for the component that builds H
#pragma once

#include <string>
#include <vector>

#include "atomiverse/cuda/component.cuh"
#include "atomiverse/facade/stage.h"

namespace ATV {
namespace cuda {
namespace scf {

/// @brief Dense symmetric eigendecomposition via cuSOLVER syevd.
///
/// In build(), allocates the cuSOLVER workspace.
/// In launch(), copies H -> C then solves in-place.
/// Eigenvalues written to frame.scf.eig, eigenvectors to frame.scf.C.
class EigenSolver : public ComponentBase {
  public:
    /// @brief Component name.
    /// @return "eigen_solver".
    std::string name() const override { return "eigen_solver"; }

    /// @brief Default binding: diagonalize phase at order 0.
    std::vector<DefaultBinding> default_bindings() const override {
        return {{static_cast<int>(ATV::scf::Phase::diagonalize), 0}};
    }

    /// @brief Allocate cuSOLVER handle and workspace.
    /// @param ctx Build context.
    void build(const BuildContext& ctx) override;

    /// @brief Diagonalize H: copy H -> C, solve C in-place, eigenvalues to eig.
    /// @param ctx Launch context.
    void launch(const LaunchContext& ctx) override;

    /// @brief Free cuSOLVER handle and workspace.
    /// @param ctx Destroy context.
    void destroy(const DestroyContext& ctx) override;

  private:
    void* cusolver_handle_ = nullptr; ///< cuSOLVER handle.
    double* d_work_ = nullptr;        ///< cuSOLVER workspace [device].
    int work_size_ = 0;               ///< Workspace size in elements.
    int* d_info_ = nullptr; ///< cuSOLVER convergence info [device, 1 element].
};

} // namespace scf
} // namespace cuda
} // namespace ATV
