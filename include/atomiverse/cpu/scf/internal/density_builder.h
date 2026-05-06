/// @file cpu/scf/internal/density_builder.h
/// @brief internal::DensityBuilder — P = factor · C_occ · C_occᵀ.
///
/// **Internal SCF helper.** Does NOT inherit `cpu::ComponentBase`.
///
/// Stateless thin wrapper around `build_density_matrix`. Phase 4
/// introduces an in-place variant that writes to a pre-allocated
/// workspace; the public surface stays the same.
#pragma once

#include "atomiverse/cpu/scf/scf_utils.h"

namespace ATV {
namespace cpu {
namespace scf {
namespace internal {

class DensityBuilder {
public:
    DensityBuilder() = default;

    /// @brief Pre-allocate workspaces sized for `n_basis`. No-op today.
    void prepare(int /*n_basis*/) {}

    /// @brief Build the density matrix from MO coefficients.
    /// @param C        MO coefficient matrix [n × n].
    /// @param n_occ    Number of occupied orbitals.
    /// @param n_basis  Basis dimension.
    /// @param factor   Occupancy factor (2 for closed-shell RHF, 1 per spin
    ///                 for unrestricted).
    /// @return Density matrix [n × n].
    RealMatrix step(const RealMatrix& C, int n_occ, int n_basis,
                    double factor = 1.0) const {
        return build_density_matrix(C, n_occ, n_basis, factor);
    }
};

} // namespace internal
} // namespace scf
} // namespace cpu
} // namespace ATV
