/// @file cpu/gpw/internal/gpw_forces.h
/// @brief Analytic ionic forces for the closed-shell GPW SCF (Γ and k-points).
///
/// Replaces the 6N+1-SCF finite-difference oracle (`gpw_forces_fd`) with a
/// single SCF solve followed by one gradient pass. At the SCF stationary point
/// the force on atom A is
///
///   F_A = -Σ_μν P_μν d(T+V_pp+V_grid)_μν/dR_A
///         + Σ_μν W_μν dS_μν/dR_A
///         + ∫ V_es(r) ∂n_c^A(r)/∂R_A dr,
///
/// with the energy-weighted density matrix W = Σ_i 2 ε_i c_{μi} c_{νi} (closed
/// shell). All AO-derivative terms use the Gaussian center-derivative identity
/// ∂φ_a/∂A_x = 2α φ_{a+1_x} − a_x φ_{a−1_x}, so they reuse the existing overlap
/// / kinetic / GTH / collocation kernels with raised/lowered Cartesian powers.
/// The GTH operator-center derivative is obtained from the basis derivatives via
/// translational invariance (∂_C = −(∂_μ + ∂_ν)). The grid term holds the
/// converged effective potential V_eff fixed (Hellmann–Feynman for the density).
/// V_es held fixed likewise for the moving Gaussian ionic core n_c^A.
///
/// On a Monkhorst–Pack / explicit k mesh the same assembly is k-summed with
/// Bloch-phased complex matrices and complex χ^k (see `internal/kpoints.h`);
/// `GpwResult::C_occ_k` supplies the occupied MOs used to rebuild P_k / W_k.
///
/// Reference: VandeVondele, Krack, Mohamed, Parrinello, Chassaing, Hutter,
/// Comput. Phys. Commun. 167, 103 (2005) (Quickstep GPW forces).
#pragma once

#include <array>
#include <vector>

#include "atomiverse/cpu/gpw/gpw_solver.h" // GpwAtom, GpwDft

namespace ATV {
namespace cpu {
namespace gpw {

/// @brief Analytic GPW forces F_a = -dE/dR_a (one SCF + one gradient pass).
///
/// Runs `GpwDft`, then assembles the Pulay (W·dS), core-Hamiltonian (P·dT,
/// P·dV_pp), grid (P·dV_grid against the converged V_eff) and ionic-core
/// Hellmann–Feynman (∫V_es ∂n_c) terms analytically.
///
/// @param atoms GPW atoms (basis + GTH PP), as for `GpwDft::run`.
/// @param cfg   SCF configuration (grid, XC, electrostatics).
/// @param guess_C Optional warm-start MO coefficients for the internal SCF
///        (e.g. the energy pass just converged at this geometry — XL-BOMD
///        path). Default nullptr keeps the cold solve bit-for-bit.
/// @return Forces [Hartree/Bohr], one xyz triple per atom, in input order.
/// @throws std::runtime_error if the SCF does not converge.
std::vector<std::array<double, 3>>
gpw_forces_analytic(const std::vector<GpwAtom>& atoms, const GpwDft& cfg,
                    const linalg::RealMatrix* guess_C = nullptr);

/// @brief Analytic ionic forces at a **fixed** electronic state (no SCF).
///
/// Car–Parrinello path: nuclei feel \(-\partial E/\partial R\) at the current
/// density \(P = 2 C C^\top\) with energy-weighted
/// \(W = 2\,C\,(C^\top H C)\,C^\top\) (constraint contribution off the BO
/// surface). Same force assembly as `gpw_forces_analytic` after its SCF.
///
/// @param atoms Geometry + basis/PP [Bohr].
/// @param cfg   Grid / XC / electrostatics (same functional as the energy).
/// @param P     Closed-shell AO density [e], nb × nb.
/// @param C_occ Occupied MO coefficients, nb × n_occ (S-orthonormal).
/// @param H     KS matrix at this density [Hartree], nb × nb.
/// @return Forces [Hartree/Bohr], one xyz triple per atom.
std::vector<std::array<double, 3>>
gpw_forces_at_density(const std::vector<GpwAtom>& atoms, const GpwDft& cfg,
                      const linalg::RealMatrix& P,
                      const linalg::RealMatrix& C_occ,
                      const linalg::RealMatrix& H);

} // namespace gpw
} // namespace cpu
} // namespace ATV
