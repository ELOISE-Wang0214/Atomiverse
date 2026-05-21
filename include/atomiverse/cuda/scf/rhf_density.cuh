/// @file cuda/scf/rhf_density.cuh
/// @brief RHFDensity -- back-transform C, rebuild P, compute RHF energy.
///
/// In launch():
///   1. C = X * C' (back-transform from orthogonal basis)
///   2. P_old = P; P = 2 * C_occ * C_occ^T (RHF closed-shell)
///   3. E = 0.5 * Tr[(Hcore + F) * P] + V_nuc
///   4. Write E to frame.energy
///
/// @see scf::Phase::build_density
/// @see RHFFock for the Fock matrix builder
#pragma once

#include <string>
#include <vector>

#include "atomiverse/cuda/component.cuh"
#include "atomiverse/facade/stage.h"

namespace ATV {
namespace scf {

#ifdef ATV_HAS_RHF

/// @brief Back-transform MO coefficients, rebuild density, compute RHF energy.
class RHFDensity : public ComponentBase {
  public:
    /// @brief Construct with occupation count.
    /// @param n_occ Number of doubly-occupied orbitals.
    explicit RHFDensity(int n_occ) : n_occ_(n_occ) {}

    std::string name() const override { return "rhf_density"; }

    std::vector<DefaultBinding> default_bindings() const override {
        return {{static_cast<int>(Phase::build_density), 0}};
    }

    void build(const BuildContext& ctx) override;
    void launch(const LaunchContext& ctx) override;
    void destroy(const DestroyContext& ctx) override;

  private:
    int n_occ_;
    void* cublas_handle_ = nullptr;
    double* d_HplusF_ =
        nullptr; ///< Temporary Hcore+F for energy [device, n*n].
};

#endif // ATV_HAS_RHF

} // namespace scf
} // namespace ATV
