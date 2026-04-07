/// @file cuda/scf/rhf_fock.cuh
/// @brief RHFFock -- build Fock matrix F = Hcore + G(P), transform F' = X^T F X.
///
/// In build(), uploads Hcore, ERI, and X = S^{-1/2} to device.
/// In launch(), builds two-electron matrix G(P) from ERI and current density,
/// forms Fock matrix F, and transforms to orthogonal basis F'.
///
/// @see scf::Phase::build_hamiltonian
/// @see RHFDensity for back-transform and density rebuild
#pragma once

#include <string>
#include <vector>

#include "atomiverse/cuda/component.cuh"
#include "atomiverse/facade/stage.h"

#ifdef ATV_HAS_RHF
#include "atomiverse/cpu/scf/integrals.h"
#endif

namespace ATV {
namespace scf {

#ifdef ATV_HAS_RHF

/// @brief Build Fock matrix and transform to orthogonal basis (CPU integral path).
///
/// build(): compute X = S^{-1/2} on host, upload Hcore/ERI/X to device.
/// launch(): G(P) from ERI, F = Hcore + G, F' = X^T F X.
class RHFFock : public ComponentBase {
public:
    /// @brief Construct with precomputed CPU integrals.
    explicit RHFFock(const cpu::scf::IntegralPack& integrals)
        : integrals_(integrals) {}

    std::string name() const override { return "rhf_fock"; }

    std::vector<DefaultBinding> default_bindings() const override {
        return {{static_cast<int>(Phase::build_hamiltonian), 0}};
    }

    void build(const BuildContext& ctx) override;
    void launch(const LaunchContext& ctx) override;
    void destroy(const DestroyContext& ctx) override;

private:
    cpu::scf::IntegralPack integrals_;
    void* cublas_handle_ = nullptr;
    double* d_temp_ = nullptr;
};

/// @brief Build Fock matrix — CUDA-only path (GPU integral kernels).
///
/// build(): upload basis params to device, launch GPU integral kernels to
///          compute S, T, V, ERI directly on device, then Hcore = T + V
///          and X = S^{-1/2} via cuSOLVER — zero CPU integral computation.
/// launch(): same as RHFFock (G(P) from ERI, F = Hcore + G, F' = X^T F X).
///
/// Usage:
/// @code{.cpp}
/// solver.add(scf::RHFFockDevice{n_prim, exponents, coefficients, n_occ});
/// @endcode
class RHFFockDevice : public ComponentBase {
public:
    /// @brief Construct for CUDA-only integral computation.
    /// @param n_prim Primitives per shell (3 for STO-3G).
    /// @param exponents Gaussian exponents [n_prim].
    /// @param coefficients Raw contraction coefficients [n_prim] (NOT normalized).
    /// @param n_occ Number of doubly-occupied orbitals.
    RHFFockDevice(int n_prim, const double* exponents,
                  const double* coefficients, int n_occ);

    std::string name() const override { return "rhf_fock_device"; }

    std::vector<DefaultBinding> default_bindings() const override {
        return {{static_cast<int>(Phase::build_hamiltonian), 0}};
    }

    void build(const BuildContext& ctx) override;
    void launch(const LaunchContext& ctx) override;
    void destroy(const DestroyContext& ctx) override;

private:
    int n_prim_;
    int n_occ_;
    std::vector<double> exponents_;
    std::vector<double> coefficients_;
    void* cublas_handle_ = nullptr;
    void* cusolver_handle_ = nullptr;
    double* d_temp_ = nullptr;
    double* d_S_ = nullptr;  ///< Overlap (device, for X computation).
    double* d_T_ = nullptr;  ///< Kinetic (device, temp).
    double* d_V_ = nullptr;  ///< Nuclear (device, temp).
};

#endif // ATV_HAS_RHF

} // namespace scf
} // namespace ATV
