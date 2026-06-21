/// @file cpu/gpw/gpw_solver.h
/// @brief Self-contained GPW (Gaussian and Plane Waves) closed-shell SCF solver.
///
/// Structurally parallel to `cpu::scf::SCFSolver` and `cpu::pw::PWSolver`: a
/// self-contained inner solver for one electronic-structure method. The Kohn-Sham
/// problem is solved in the contracted Gaussian AO basis (small dense generalized
/// eigenproblem H C = S C ε), while the Hartree + XC + long-range local
/// pseudopotential are evaluated on a uniform real-space grid (collocate / FFT
/// Poisson / integrate). It reuses only shared *physics/algebra* from public
/// headers — `cpu::linalg` (generalized eigensolve, density build) and
/// `cpu::xc::pade` (functional) — never another domain's `internal`.
///
/// Scope: closed-shell, Gamma point, LDA (PADE). Periodic neutral electrostatics
/// via the GTH Gaussian core charge. Forces / GGA / nonlocal-for-all-l and full
/// CP2K bit-parity are follow-on milestones.
#pragma once

#include <array>
#include <vector>

#include "atomiverse/cpu/gpw/internal/collocate.h"
#include "atomiverse/cpu/gpw/internal/cp2k_basis.h" // GaussianBasis
#include "atomiverse/cpu/gpw/internal/cp2k_io.h"    // GthPotential, GthNonlocalChannel
#include "atomiverse/cpu/linalg/matrix.h"           // RealMatrix

namespace ATV {
namespace cpu {
namespace gpw {

/// @brief One atom's GPW data: valence charge, basis, and GTH pseudopotential.
struct GpwAtom {
    double z_val = 0.0;                       ///< valence charge (electrons it contributes).
    double pos[3] = {0, 0, 0};                ///< position [Bohr].
    std::vector<GridBasisFunction> basis;     ///< contracted GTOs, already centered at pos.
    double rloc = 0.0;                        ///< GTH local radius.
    std::vector<double> cloc;                 ///< GTH local short-range C coefficients.
    std::vector<GthNonlocalChannel> nonlocal; ///< GTH nonlocal projector channels.
};

/// @brief Exchange-correlation functional选择.
enum class XcKind { PADE, PBE };

/// @brief GPW SCF configuration.
struct GpwConfig {
    double cell[3] = {0, 0, 0}; ///< orthorhombic box lengths [Bohr].
    double cutoff_ry = 300.0;   ///< density grid cutoff [Ry] (finest multigrid level).
    int ngrids = 1;             ///< CP2K multigrid levels; 1 = single grid (default, bit-identical).
    double rel_cutoff = 0.0;    ///< CP2K REL_CUTOFF [Ry]; exponent->level map when ngrids>1.
    int n_occ = 0;              ///< number of doubly-occupied orbitals (N_elec / 2).
    int max_iter = 100;         ///< max SCF iterations.
    double tol = 1e-8;          ///< SCF energy convergence [Hartree].
    double mix = 0.5;           ///< linear density-mixing factor (used only when diis_depth == 0).
    int diis_depth = 6;         ///< Pulay DIIS history depth; 0 disables DIIS (falls back to linear `mix`).
    XcKind xc = XcKind::PADE;   ///< LDA (PADE) or GGA (PBE).
    bool nonperiodic = false;   ///< use the isolated (MT) Poisson solver.
    /// @brief Use the OT (orbital transformation) minimizer instead of dense
    /// diagonalization. Diagonalization-free SCF (gradient O(nb²·n_occ)); the
    /// path to large systems where the O(nb³) eigensolve is prohibitive.
    /// Converges to the same gapped ground state as the diagonalization path.
    bool use_ot = false;
};

/// @brief GPW SCF result with an energy decomposition.
struct GpwResult {
    double e_total = 0.0;   ///< total energy [Hartree].
    double e_kinetic = 0.0; ///< Tr(P T).
    double e_loc_sr = 0.0;  ///< Tr(P V_loc^sr).
    double e_nonlocal = 0.0;///< Tr(P V_nl).
    double e_es = 0.0;      ///< total electrostatic (Hartree + e-ion-lr + ion-ion).
    double e_xc = 0.0;      ///< XC energy.
    bool converged = false;
    int iters = 0;

    /// @name Converged wavefunction (populated on convergence; empty otherwise).
    /// Consistent set from a final diagonalization of the converged KS matrix:
    /// `P == build_density_matrix(C, n_occ)`. Used to assemble analytic forces
    /// (energy-weighted density matrix W = Σ_i 2 ε_i c_i c_iᵀ).
    /// @{
    linalg::RealMatrix C;      ///< MO coefficients, columns = orbitals (nb×nb).
    linalg::RealMatrix P;      ///< Converged density matrix (nb×nb).
    std::vector<double> eps;   ///< Orbital energies ε_i [Hartree], length nb.
    /// @}
};

/// @brief Build a GpwAtom from parsed CP2K data (basis + GTH potential).
///
/// Converts each contracted shell of `basis` into Cartesian GridBasisFunctions
/// centered at `pos`, folding the per-primitive Gaussian normalization
/// N(alpha,L) = (2 alpha/pi)^{3/4} (4 alpha)^{L/2} into the coefficient (the L
/// alpha-power is what fixes the spanned space to match CP2K; the constant
/// factor is irrelevant to the energy). Angular functions are the (2L+1)
/// pure-spherical real solid harmonics from real_solid_harmonics(L), for any L.
/// @param pos Atom position [Bohr]; z_val is taken from `pot.zion`.
GpwAtom make_gpw_atom(const double pos[3], const GaussianBasis& basis,
                      const GthPotential& pot);

/// @brief Run the closed-shell GPW SCF for a set of atoms.
GpwResult gpw_scf(const std::vector<GpwAtom>& atoms, const GpwConfig& cfg);

/// @brief Finite-difference forces F_a = -dE/dR_a via central differences.
///
/// Runs `gpw_scf` at +/- h displacements of each atom along each axis
/// (6*N_atoms SCF solves). Correct but expensive; analytic GPW forces are a
/// follow-on. Returns forces [Hartree/Bohr], one xyz triple per atom.
/// @param h Displacement step [Bohr].
std::vector<std::array<double, 3>>
gpw_forces_fd(const std::vector<GpwAtom>& atoms, const GpwConfig& cfg, double h = 1e-3);

} // namespace gpw
} // namespace cpu
} // namespace ATV
