/// @file cpu/gpw/gpw_solver.h
/// @brief Self-contained GPW (Gaussian and Plane Waves) SCF solver.
///
/// Structurally parallel to `cpu::scf::SCFSolver` and `cpu::pw::PW`: a
/// self-contained inner solver for one electronic-structure method. The Kohn-Sham
/// problem is solved in the contracted Gaussian AO basis (small dense generalized
/// eigenproblem H C = S C ε), while the Hartree + XC + long-range local
/// pseudopotential are evaluated on a uniform real-space grid (collocate / FFT
/// Poisson / integrate). It reuses only shared *physics/algebra* from public
/// headers — `cpu::linalg` (generalized eigensolve, density build) and
/// `cpu::xc::pade` (functional) — never another domain's `internal`.
///
/// Scope: Γ-point closed-shell (RKS, `n_occ`) **and** unrestricted (UKS,
/// `n_alpha`/`n_beta` with PADE LSD). Periodic neutral electrostatics via the
/// GTH Gaussian core charge. UKS PBE-GGA, k-points, and UKS analytic forces
/// are follow-on.
#pragma once

#include <array>
#include <string>
#include <vector>

#include "atomiverse/cpu/gpw/gaussian_basis.h" // GaussianBasis
#include "atomiverse/cpu/gpw/gth_potential.h" // GthPotential, GthNonlocalChannel
#include "atomiverse/cpu/gpw/internal/collocate.h"
#include "atomiverse/cpu/linalg/matrix.h" // RealMatrix

namespace ATV {
namespace cpu {
namespace gpw {

/// @brief One atom's GPW data: valence charge, basis, and GTH pseudopotential.
struct GpwAtom {
    double z_val = 0.0;        ///< valence charge (electrons it contributes).
    int z = 0;                 ///< atomic number Z (for DFT-D3; 0 = unset).
    double pos[3] = {0, 0, 0}; ///< position [Bohr].
    std::vector<GridBasisFunction>
        basis;                ///< contracted GTOs, already centered at pos.
    double rloc = 0.0;        ///< GTH local radius.
    std::vector<double> cloc; ///< GTH local short-range C coefficients.
    std::vector<GthNonlocalChannel>
        nonlocal; ///< GTH nonlocal projector channels.

    GpwAtom() = default;

    /// @brief Construct from a contracted basis + GTH potential.
    ///
    /// Converts each contracted shell of `basis` into Cartesian
    /// GridBasisFunctions centered at `pos`, folding the per-primitive Gaussian
    /// normalization N(alpha,L) = (2 alpha/pi)^{3/4} (4 alpha)^{L/2} into the
    /// coefficient (the L alpha-power is what fixes the spanned space to match
    /// CP2K; the constant factor is irrelevant to the energy). Angular
    /// functions are the (2L+1) pure-spherical real solid harmonics from
    /// real_solid_harmonics(L), for any L.
    /// @param pos   Atom position [Bohr]; z_val is taken from `pot.zion`.
    /// @param basis Contracted Gaussian basis to center at `pos`.
    /// @param pot   GTH pseudopotential (local + nonlocal channels).
    /// @param z     Atomic number Z (needed only for DFT-D3, GpwDft.d3).
    GpwAtom(const double pos[3], const GaussianBasis& basis,
            const GthPotential& pot, int z = 0);
};

/// @brief Exchange-correlation functional choice (RKS).
///
/// `PADE` — LDA Padé; `PBE` / `PBEsol` — pure GGA (shared `evaluate_xc`).
/// UKS currently supports `PADE` only.
enum class XcKind { PADE, PBE, PBEsol };

/// @brief GPW SCF result with an energy decomposition.
struct GpwResult {
    double e_total = 0.0;    ///< total energy [Hartree].
    double e_kinetic = 0.0;  ///< Tr(P T)  (RKS: closed P; UKS: Tr((Pα+Pβ)T)).
    double e_loc_sr = 0.0;   ///< Tr(P V_loc^sr).
    double e_nonlocal = 0.0; ///< Tr(P V_nl).
    double e_es = 0.0; ///< total electrostatic (Hartree + e-ion-lr + ion-ion).
    double e_xc = 0.0; ///< XC energy.
    double e_disp = 0.0; ///< DFT-D3 dispersion [Hartree] (0 when d3/d3bj off).
    double s2 = 0.0;   ///< ⟨S²⟩ (Löwdin; RKS → 0; UKS after convergence).
    bool converged = false;
    int iters = 0;

    /// @name Converged wavefunction (populated on convergence; empty otherwise).
    /// RKS: `P == build_density_matrix(C, n_occ, ·, 2)` from the closed H.
    /// UKS: `C`/`P`/`eps` are α; `C_beta`/`P_beta`/`eps_beta` are β.
    /// @{
    linalg::RealMatrix C;    ///< MO coefficients α (or RKS), columns = orbitals.
    linalg::RealMatrix P;    ///< Density matrix α (or RKS closed-shell P).
    std::vector<double> eps; ///< Orbital energies α (or RKS) [Hartree].
    linalg::RealMatrix C_beta;    ///< UKS β MO coefficients (empty for RKS).
    linalg::RealMatrix P_beta;    ///< UKS β density (empty for RKS).
    std::vector<double> eps_beta; ///< UKS β orbital energies (empty for RKS).
    /// @}
};

/// @brief Concrete GPW DFT method — members + build/run
///        (mirrors `cpu::pw::PW` / `cpu::scf::SCFSolver`).
///
/// **NOT a `Driver<...>` subclass** — this is the inner solver behind
/// `GPWEvaluator`. Construction fixes the configuration; `build()` runs
/// the one-time fail-fast validation against a reference atom set (the
/// no-fallback rule: invalid configuration throws instead of degrading);
/// `run()` is re-entrant and rebuilds every geometry-dependent structure
/// per call (grids, lattice images, S/T/Vpp — same policy as `PW`).
/// Move-only, like the sibling solvers.
class GpwDft {
  public:
    // ── Method parameters (were GpwConfig fields) ─────────────────────
    double cell[3] = {0, 0, 0}; ///< orthorhombic box lengths [Bohr].
    double cutoff_ry = 300.0; ///< density grid cutoff [Ry].
    int ngrids = 1; ///< CP2K multigrid levels.
    double rel_cutoff = 0.0; ///< CP2K REL_CUTOFF [Ry].
    int n_occ = 0; ///< RKS: doubly-occupied orbitals (N_elec / 2).
    /// UKS occupations (valence electrons per spin). When
    /// `n_alpha + n_beta > 0`, the unrestricted path is used and `n_occ` is
    /// ignored. High-spin O₂ (GTH-q6, 12 valence e⁻, multiplicity 3) → (7, 5).
    int n_alpha = 0;
    int n_beta = 0;
    int max_iter = 100; ///< max SCF iterations.
    double tol = 1e-8; ///< SCF energy convergence [Hartree].
    double eps_default = 1e-8; ///< Gaussian-tail cutoff.
    double mix = 0.5; ///< linear density-mixing (diis_depth == 0).
    int diis_depth = 6; ///< Pulay DIIS history depth.
    double lindep_threshold = 1e-8; ///< S-eigenvalue keep threshold.
    double level_shift = 0.0; ///< virtual level shift [Hartree].
    XcKind xc = XcKind::PADE; ///< LDA (PADE) or GGA (PBE / PBEsol).
    bool nonperiodic = false; ///< isolated (MT) Poisson solver.
    bool local_pp_on_grid = true; ///< collocate GTH local SR PP on grid.
    bool d3 = false; ///< Grimme DFT-D3 zero damping (PBE). Mutually exclusive with d3bj.
    /// @brief Grimme DFT-D3 Becke–Johnson damping (PBE params, CP2K CN form).
    /// Mutually exclusive with `d3` — both true throws at run.
    bool d3bj = false; ///< Grimme DFT-D3 dispersion.
    bool use_ot = false; ///< OT minimizer instead of diagonalization.

    GpwDft() = default;
    GpwDft(const GpwDft&) = delete;
    GpwDft& operator=(const GpwDft&) = delete;
    GpwDft(GpwDft&&) = default;
    GpwDft& operator=(GpwDft&&) = default;

    /// @brief One-time fail-fast validation against a reference system.
    ///
    /// @param reference Atoms the solver will be run on (used for
    ///        structural validation only; `run()` takes the live set).
    /// @throws std::invalid_argument when the atom set is empty, when
    ///         neither RKS (`n_occ >= 1`) nor UKS (`n_alpha+n_beta >= 1`) is
    ///         set, UKS+PBE is requested, or a periodic cell edge is
    ///         non-positive.
    void build(const std::vector<GpwAtom>& reference);

    /// @return true when this method is configured for UKS (`n_alpha+n_beta>0`).
    bool is_uks() const { return n_alpha + n_beta > 0; }

    /// @brief Run the GPW SCF (RKS or UKS) at the given geometry.
    ///
    /// Re-entrant: geometry-dependent structures are rebuilt per call.
    /// @param atoms Atom set at the current geometry [Bohr].
    /// @param guess_C Optional warm-start MO coefficients from an earlier run.
    /// @return The owned result of this run (valid until the next `run()`).
    /// @throws std::logic_error when called before `build()`.
    GpwResult& run(const std::vector<GpwAtom>& atoms,
                   const linalg::RealMatrix* guess_C = nullptr);

    /// @brief Result of the most recent `run()` (default-constructed —
    ///        `converged == false` — before the first run).
    const GpwResult& last_result() const { return result_; }

  private:
    GpwResult result_;   ///< Owned result buffer (last run).
    bool built_ = false; ///< Whether build() validated the configuration.
};

/// @brief Finite-difference forces F_a = -dE/dR_a via central differences.
///
/// Runs a `GpwDft` at +/- h displacements of each atom along each axis
/// (6*N_atoms SCF solves). Correct but expensive; analytic GPW forces are a
/// follow-on. Returns forces [Hartree/Bohr], one xyz triple per atom.
/// @param h Displacement step [Bohr].
std::vector<std::array<double, 3>>
gpw_forces_fd(const std::vector<GpwAtom>& atoms, const GpwDft& method,
              double h = 1e-3);

} // namespace gpw
} // namespace cpu
} // namespace ATV
