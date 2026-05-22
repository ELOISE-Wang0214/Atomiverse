/// @file cuda/scf/convergence_check.cuh
/// @brief ConvergenceCheck -- SCF convergence via Frobenius norm of P - P_old.
///
/// Computes ||P - P_old||_F using cuBLAS (axpy + nrm2) and compares against
/// the specified tolerance. Sets frame.scf.converged when satisfied.
///
/// When compute_energy is true (default), also computes E = dot(H, P) and
/// writes it to frame.energy. Set to false when another component (e.g.
/// RHFDensity) computes the proper energy.
///
/// @see scf::Phase::convergence_check
/// @see DensityBuilder for the component that produces P and P_old
#pragma once

#include <string>
#include <vector>

#include "atomiverse/cuda/component.cuh"
#include "atomiverse/facade/stage.h"

namespace ATV {
namespace cuda {
namespace scf {

/// @brief SCF convergence check: ||P_new - P_old||_F < tol.
///
/// Optionally computes SCF energy E = trace(H * P) via cuBLAS dot.
class ConvergenceCheck : public ComponentBase {
  public:
    /// @brief Construct convergence checker.
    /// @param tol Frobenius norm tolerance for convergence.
    /// @param compute_energy If true, compute E=dot(H,P) and write to frame.energy.
    explicit ConvergenceCheck(double tol, bool compute_energy = true)
        : tol_(tol), compute_energy_(compute_energy) {}

    /// @brief Component name.
    /// @return "convergence_check".
    std::string name() const override { return "convergence_check"; }

    /// @brief Default binding: convergence_check phase at order 0.
    std::vector<DefaultBinding> default_bindings() const override {
        return {{static_cast<int>(ATV::scf::Phase::convergence_check), 0}};
    }

    /// @brief Allocate cuBLAS handle.
    /// @param ctx Build context.
    void build(const BuildContext& ctx) override;

    /// @brief Compute ||P - P_old||_F, check convergence, optionally compute energy.
    /// @param ctx Launch context.
    void launch(const LaunchContext& ctx) override;

    /// @brief Free cuBLAS handle.
    /// @param ctx Destroy context.
    void destroy(const DestroyContext& ctx) override;

  private:
    double tol_;                    ///< Convergence tolerance (Frobenius norm).
    bool compute_energy_;           ///< Whether to compute and write energy.
    void* cublas_handle_ = nullptr; ///< cuBLAS handle.
};

} // namespace scf
} // namespace cuda
} // namespace ATV
