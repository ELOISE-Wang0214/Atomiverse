/// @file cpu/gpw/internal/mp2_kernels.h
/// @brief Private GPW MP2 kernels: MO values on the density grid, the
///        occupied–virtual MO repulsion tensor via per-pair Poisson solves,
///        and the canonical closed-shell doubles sum.
///
/// Private workspace behind `GpwDft::mp2` — never on the public include
/// path (unit tests may include it). The `(ia|jb)` route is CP2K's
/// Quickstep wavefunction-correlation integral scheme: collocate
/// \(\rho_{ia}=\psi_i\psi_a\), solve Poisson, contract with
/// \(\rho_{jb}\) by grid quadrature. The same tensor + denominators feed
/// the periodic RPA follow-on (`direct_rpa_rhf_energy`).
#pragma once

#include <vector>

#include "atomiverse/cpu/gpw/internal/collocate.h" // CompactBasisGrid, UniformGrid
#include "atomiverse/cpu/linalg/matrix.h"          // RealMatrix

namespace ATV {
namespace cpu {
namespace gpw {
namespace internal {

/// @brief Evaluate MO values ψ_p(r_g) on the grid from the compact
///        (lattice-summed, CSR-by-point) AO representation.
/// @param cbg Compact AO grid from `build_compact_basis_grid` (periodic
///        Γ-image-summed or isolated, matching the SCF).
/// @param C  MO coefficients (nbasis × n_mo, columns = orbitals).
/// @param p0 First MO column to evaluate.
/// @param np Number of MO columns.
/// @return Row-major (np × npts): value of MO `p0+p` at grid point g is
///         `out[p*npts + g]`.
std::vector<double> mo_on_grid(const CompactBasisGrid& cbg,
                               const linalg::RealMatrix& C, int p0, int np);

/// @brief Occupied–virtual MO repulsion tensor (ia|jb) via per-pair Poisson.
///
/// For each pair density ρ_ia(r)=ψ_i(r)ψ_a(r) one Poisson solve gives
/// V_ia(r); then (ia|jb) = Σ_g w V_ia(g) ρ_jb(g). Uses the (ia|jb)=(jb|ia)
/// symmetry — n_occ·n_virt Poisson solves total. G=0 is dropped by the
/// periodic solver (neutral pair densities, S_ia = 0). Canonical foothold
/// path (`mp2_ri = false`); product path uses `mo_eri_ovov_ri`.
/// @param grid     Uniform density grid (as in the SCF).
/// @param psi_occ  Occupied MO values, row-major (n_occ × npts) [mo_on_grid].
/// @param n_occ    Occupied count.
/// @param psi_virt Virtual MO values, row-major (n_virt × npts).
/// @param n_virt   Virtual count.
/// @param nonperiodic true → truncated-Coulomb (MT-class) isolated solver.
/// @return Row-major ((n_occ·n_virt) × (n_occ·n_virt)) tensor,
///         `(ia|jb) = out[(i*n_virt+a)*n_occ*n_virt + j*n_virt+b]` [Hartree].
std::vector<double> mo_eri_ovov(const UniformGrid& grid,
                                const std::vector<double>& psi_occ, int n_occ,
                                const std::vector<double>& psi_virt, int n_virt,
                                bool nonperiodic);

/// @brief Periodic RI-GPW (ia|jb) via plane-wave Coulomb intermediates.
///
/// Del Ben / Quickstep product path: for each pair density ρ_ia form the
/// G-space Coulomb potential (4π/|G|², G=0 dropped) via the shared FFT
/// Poisson kernel, store the three-index real-space intermediate V_ia
/// (memory \(\mathcal O(n_{\mathrm{ov}} N_g)\)), then contract
/// \((ia|jb)=\sum_g w\, V_{ia}(g)\,\rho_{jb}(g)\). Numerically identical to
/// `mo_eri_ovov(..., nonperiodic=false)` on the same grid — the PW
/// resolution of the identity for the pair density. Periodic only (isolated
/// cells keep the MT Poisson foothold via `mo_eri_ovov`).
/// @param grid     Uniform density grid (as in the SCF).
/// @param psi_occ  Occupied MO values, row-major (n_occ × npts).
/// @param n_occ    Occupied count.
/// @param psi_virt Virtual MO values, row-major (n_virt × npts).
/// @param n_virt   Virtual count.
/// @return Same layout as `mo_eri_ovov` [Hartree].
std::vector<double> mo_eri_ovov_ri(const UniformGrid& grid,
                                   const std::vector<double>& psi_occ,
                                   int n_occ,
                                   const std::vector<double>& psi_virt,
                                   int n_virt);

/// @brief Canonical closed-shell MP2 doubles sum (real orbitals).
///
/// E_OS = Σ (ia|jb)² / (ε_i+ε_j−ε_a−ε_b),
/// E_SS = Σ (ia|jb)[(ia|jb)−(ib|ja)] / (ε_i+ε_j−ε_a−ε_b); both ≤ 0 for a
/// gapped spectrum. E_corr(unscaled) = E_OS + E_SS.
/// @param v_ovov   (ia|jb) tensor from `mo_eri_ovov` [Hartree].
/// @param eps_occ  Occupied orbital energies (length n_occ) [Hartree].
/// @param eps_virt Virtual orbital energies (length n_virt) [Hartree].
/// @param e_os Out: unscaled opposite-spin part [Hartree].
/// @param e_ss Out: unscaled same-spin part [Hartree].
void mp2_rhf_energy(const std::vector<double>& v_ovov,
                    const std::vector<double>& eps_occ,
                    const std::vector<double>& eps_virt, double& e_os,
                    double& e_ss);

/// @brief Closed-shell direct RPA (dRPA) correlation via the Furche
///        plasmon formula on the Coulomb-only (\(\alpha_{\mathrm{hf}}=0\))
///        singlet Casida spectrum.
///
/// Builds
/// \(A_{ia,jb}=\delta_{ij}\delta_{ab}(\varepsilon_a-\varepsilon_i)+2(ia|jb)\),
/// \(B_{ia,jb}=2(ia|jb)\) (real orbitals; composite index
/// \(p = i\cdot n_{\mathrm{virt}}+a\) matching `mo_eri_ovov` layout),
/// reduces \((A-B)^{1/2}(A+B)(A-B)^{1/2}\), and returns
/// \(E_c=\tfrac12(\sum_s\Omega_s-\mathrm{Tr}\,A)\). Same science as
/// molecular `OrbitalHessian::direct_rpa(singlet, /*alpha_hf=*/0)`.
///
/// @param v_ovov   (ia|jb) tensor from `mo_eri_ovov` / `mo_eri_ovov_ri` [Ha].
/// @param eps_occ  Occupied orbital energies (length n_occ) [Hartree].
/// @param eps_virt Virtual orbital energies (length n_virt) [Hartree].
/// @param e_corr   Out: dRPA correlation energy [Hartree].
/// @param tr_A     Out: \(\mathrm{Tr}\,A\) [Hartree].
/// @param omega    Out: positive \(\Omega_s\) ascending [Hartree], length
///                 \(n_{\mathrm{occ}}\cdot n_{\mathrm{virt}}\).
/// @throws std::runtime_error when (A−B) is not PD or any \(\Omega^2\le 0\).
void direct_rpa_rhf_energy(const std::vector<double>& v_ovov,
                           const std::vector<double>& eps_occ,
                           const std::vector<double>& eps_virt, double& e_corr,
                           double& tr_A, std::vector<double>& omega);

} // namespace internal
} // namespace gpw
} // namespace cpu
} // namespace ATV
