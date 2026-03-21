/// @file cpu/scf/operators.h
/// @brief Factored two-electron operator builders: J (Coulomb), K (exchange), G (fused).
///
/// These replace the fused build_G in RHF, build_uhf_G in UHF, and build_coulomb_J
/// in UKS with composable, independent operators that any method can use.
///
/// @see cpu::scf::scf_solve for the unified SCF driver that uses these
#pragma once

#include <xtensor/xtensor.hpp>

#include "atomiverse/core/real.h"
#include "atomiverse/core/int.h"
#include "atomiverse/cpu/scf/scf_utils.h"

namespace ATV {
namespace cpu {
namespace scf {

/// @brief Build Coulomb matrix J from ERI and total density.
/// @details J_ij = Σ_kl P_kl (ij|kl). Exploits (i,j) symmetry.
/// @param P_total Total density matrix P_alpha + P_beta [n×n].
/// @param ERI Two-electron integrals in chemist's notation [n^4].
/// @param n Basis dimension.
/// @param[out] J Coulomb matrix [n×n].
void build_coulomb(const RealMatrix& P_total,
                   const xt::xtensor<Real, 4>& ERI, Int n,
                   RealMatrix& J);

/// @brief Build exchange matrix K from ERI and per-spin density.
/// @details K_ij = Σ_kl P_kl (ik|jl). Exploits (i,j) symmetry.
/// @param P_spin Per-spin density matrix [n×n].
/// @param ERI Two-electron integrals in chemist's notation [n^4].
/// @param n Basis dimension.
/// @param[out] K Exchange matrix [n×n].
void build_exchange(const RealMatrix& P_spin,
                    const xt::xtensor<Real, 4>& ERI, Int n,
                    RealMatrix& K);

/// @brief Build fused Coulomb + exchange G = J - 0.5*K (RHF shortcut).
/// @details G_ij = Σ_kl P_kl [(ij|kl) - 0.5*(ik|jl)].
///          P is the total RHF density (includes factor of 2).
/// @param P_total Total RHF density matrix [n×n].
/// @param ERI Two-electron integrals in chemist's notation [n^4].
/// @param n Basis dimension.
/// @param[out] G Fused two-electron matrix [n×n].
void build_coulomb_exchange_fused(const RealMatrix& P_total,
                                  const xt::xtensor<Real, 4>& ERI, Int n,
                                  RealMatrix& G);

} // namespace scf
} // namespace cpu
} // namespace ATV
