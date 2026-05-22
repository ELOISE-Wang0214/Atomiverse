/// @file cuda/scf/hamiltonian.cuh
/// @brief Hamiltonian -- construct effective Hamiltonian H = H0 + alpha * P.
///
/// Allocates and initializes the SCF workspace matrices (H0, H, P, P_old, C, eig)
/// in build(). Constructs the effective Hamiltonian each SCF iteration in launch().
///
/// H0 is a deterministic diagonally-dominant symmetric matrix seeded from
/// atom indices. The coupling term alpha*P creates a self-consistent cycle.
///
/// @see scf::Phase::build_hamiltonian
/// @see DensityBuilder for the component that rebuilds P
#pragma once

#include <string>
#include <vector>

#include "atomiverse/cuda/component.cuh"
#include "atomiverse/facade/stage.h"

namespace ATV {
namespace cuda {
namespace scf {

/// @brief Construct effective Hamiltonian H(P) = H0 + alpha * P.
///
/// Also responsible for SCF workspace allocation in build():
/// allocates H0, H, P, P_old, C, eig, work on device and initializes H0
/// as a deterministic symmetric matrix.
///
/// @par Parameters
/// - n: matrix dimension (= frame.n_atoms())
/// - n_occ: number of occupied orbitals
/// - alpha: coupling strength for density feedback
class Hamiltonian : public ComponentBase {
  public:
    /// @brief Construct Hamiltonian builder.
    /// @param n_occ Number of occupied orbitals.
    /// @param alpha Coupling parameter for H = H0 + alpha*P.
    Hamiltonian(int n_occ, double alpha) : n_occ_(n_occ), alpha_(alpha) {}

    /// @brief Component name.
    /// @return "hamiltonian".
    std::string name() const override { return "hamiltonian"; }

    /// @brief Default binding: build_hamiltonian phase at order 0.
    std::vector<DefaultBinding> default_bindings() const override {
        return {{static_cast<int>(ATV::scf::Phase::build_hamiltonian), 0}};
    }

    /// @brief Allocate SCF workspace and initialize H0 on device.
    /// @param ctx Build context.
    void build(const BuildContext& ctx) override;

    /// @brief Compute H = H0 + alpha * P via cuBLAS.
    /// @param ctx Launch context.
    void launch(const LaunchContext& ctx) override;

    /// @brief Free cuBLAS handle.
    /// @param ctx Destroy context.
    void destroy(const DestroyContext& ctx) override;

  private:
    int n_occ_;    ///< Number of occupied orbitals.
    double alpha_; ///< Density coupling parameter.
    void* cublas_handle_ =
        nullptr; ///< cuBLAS handle (cast to cublasHandle_t in .cu).
};

} // namespace scf
} // namespace cuda
} // namespace ATV
