/// @file cuda/scf/density_builder.cuh
/// @brief DensityBuilder -- rebuild density matrix P = C_occ * C_occ^T.
///
/// Each SCF iteration, saves the current density to P_old, then rebuilds P
/// from the first n_occ eigenvectors using cuBLAS syrk.
///
/// @see scf::Phase::build_density
/// @see EigenSolver for the component that produces eigenvectors C
#pragma once

#include <string>
#include <vector>

#include "atomiverse/cuda/component.cuh"
#include "atomiverse/facade/stage.h"

namespace ATV {
namespace cuda {
namespace scf {

/// @brief Rebuild density matrix P from occupied eigenvectors.
///
/// P = C[:, 0:n_occ] * C[:, 0:n_occ]^T  (cuBLAS ssyrk/dsyrk)
///
/// Before rebuilding, copies P -> P_old for convergence checking.
class DensityBuilder : public ComponentBase {
  public:
    /// @brief Component name.
    /// @return "density_builder".
    std::string name() const override { return "density_builder"; }

    /// @brief Default binding: build_density phase at order 0.
    std::vector<DefaultBinding> default_bindings() const override {
        return {{static_cast<int>(ATV::scf::Phase::build_density), 0}};
    }

    /// @brief Allocate cuBLAS handle.
    /// @param ctx Build context.
    void build(const BuildContext& ctx) override;

    /// @brief Copy P -> P_old, then rebuild P = C_occ * C_occ^T.
    /// @param ctx Launch context.
    void launch(const LaunchContext& ctx) override;

    /// @brief Free cuBLAS handle.
    /// @param ctx Destroy context.
    void destroy(const DestroyContext& ctx) override;

  private:
    void* cublas_handle_ = nullptr; ///< cuBLAS handle.
};

} // namespace scf
} // namespace cuda
} // namespace ATV
