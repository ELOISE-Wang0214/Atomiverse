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
/// Scope: closed-shell (RKS, `n_occ`) at Γ **or** on a Monkhorst-Pack k mesh
/// (`k_mesh` / explicit `k_points`; complex Bloch matrices per k, k-summed
/// density), and Γ-point unrestricted (UKS, `n_alpha`/`n_beta` with PADE
/// LSD). Periodic neutral electrostatics via the GTH Gaussian core charge.
/// UKS PBE-GGA, k-point UKS, and analytic k-point forces are follow-on.
#pragma once

#include <array>
#include <string>
#include <vector>

#include "atomiverse/cpu/gpw/gaussian_basis.h" // GaussianBasis
#include "atomiverse/cpu/gpw/gth_potential.h" // GthPotential, GthNonlocalChannel
#include "atomiverse/cpu/gpw/internal/collocate.h"
#include "atomiverse/cpu/gpw/mp2.h" // Mp2Result
#include "atomiverse/cpu/linalg/matrix.h" // RealMatrix
#include "atomiverse/molrs.hpp"

namespace ATV {
namespace cpu {
namespace gpw {

/// @brief Dudarev DFT+U (CP2K PLUS_U) on one element / angular-momentum shell.
///
/// Parameters are constructor args / members — no `PlusUConfig` bag.
/// @param element Correlated element (e.g. Ni).
/// @param l       Angular momentum of the correlated shell (d → 2).
/// @param u_minus_j_ev Effective U−J [eV] (public surface; converted to Ha inside).
struct PlusU {
    molrs::Element element = molrs::Element::H;
    int l = 2;
    double u_minus_j_ev = 0.0;

    PlusU() = default;
    PlusU(molrs::Element element, int l, double u_minus_j_ev)
        : element(element), l(l), u_minus_j_ev(u_minus_j_ev) {}
};

/// @brief One atom's GPW data: valence charge, basis, and GTH pseudopotential.
struct GpwAtom {
    double z_val = 0.0;        ///< valence charge (electrons it contributes).
    int z = 0;                 ///< atomic number Z (for DFT-D3; 0 = unset).
    double pos[3] = {0, 0, 0}; ///< position [Bohr].
    std::vector<GridBasisFunction>
        basis;                ///< contracted GTOs, already centered at pos.
    /// @brief Angular momentum ℓ of each entry in `basis` (same length; for DFT+U).
    std::vector<int> basis_l;
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
    double e_disp = 0.0; ///< DFT-D3/D4 dispersion [Hartree] (0 when off).
    double e_dftu = 0.0; ///< Dudarev DFT+U energy [Hartree] (0 when plus_u empty).
    /// @name k-point sampling payload (populated by the k-point RKS path only;
    /// empty for Γ / UKS runs). One entry per sampled k.
    /// @{
    std::vector<std::array<double, 3>> k_frac; ///< fractional k coordinates.
    std::vector<double> k_weights;             ///< quadrature weights (Σ = 1).
    /// @brief Band energies per k [Hartree], ascending; one vector per k
    ///        (length = linearly independent basis directions at that k).
    std::vector<std::vector<double>> eps_k;
    /// @brief Total sampled k-points after mesh / TR fold / explicit list.
    int n_kpts = 0;
    /// @brief k-points owned by this MPI rank under round-robin
    ///        (`ik = rank; ik < n_kpts; ik += nprocs`). Equals `n_kpts` when
    ///        MPI is inactive or `nprocs == 1`.
    int n_kpts_local = 0;
    /// @}
    double s2 = 0.0;   ///< ⟨S²⟩ (Löwdin; RKS → 0; UKS after convergence).
    /// @brief Mulliken spin populations m_A = Nα−Nβ [e] per atom (empty if RKS).
    std::vector<double> mulliken_spin;
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
    bool d3 = false; ///< Grimme DFT-D3 zero damping (PBE). Exclusive vs d3bj/d4.
    /// @brief Grimme DFT-D3 Becke–Johnson damping (PBE params, CP2K CN form).
    /// Mutually exclusive with `d3` / `d4` — any two true throws at run.
    bool d3bj = false; ///< Grimme DFT-D3(BJ) dispersion.
    /// @brief Grimme DFT-D4 BJ damping (PBE params, two-body; periodic EEQ).
    /// Mutually exclusive with `d3` / `d3bj`.
    bool d4 = false; ///< Grimme DFT-D4 dispersion.
    bool use_ot = false; ///< OT minimizer instead of diagonalization.
    /// @brief Monkhorst–Pack k mesh (CP2K `&KPOINTS SCHEME MONKHORST-PACK`
    ///        convention: κ_d = (2r − N_d − 1)/(2N_d), so even N is a
    ///        Γ-avoiding shifted mesh). `{1,1,1}` keeps the Γ-point path.
    ///        The k-point path is closed-shell RKS, periodic, single-grid,
    ///        `local_pp_on_grid` only — other combinations throw at `build()`.
    int k_mesh[3] = {1, 1, 1};
    /// @brief Fold ±k pairs by time reversal (exact for RKS; weights doubled).
    bool k_time_reversal = true;
    /// @brief Explicit fractional k list (CP2K `SCHEME GENERAL`); non-empty
    ///        overrides `k_mesh`. Same length as `k_point_weights`
    ///        (weights must sum to 1; no implicit normalization).
    std::vector<std::array<double, 3>> k_points;
    std::vector<double> k_point_weights; ///< weights for `k_points` [Σ = 1].
    /// @brief Dudarev PLUS_U entries (empty = off). U−J in eV on each entry.
    std::vector<PlusU> plus_u;
    /// @brief MP2 virtual-space cap for `mp2()` (0 = all kept virtuals).
    int mp2_n_virt = 0;
    /// @brief Grimme opposite-spin scale for `mp2()` (1.0 = plain MP2;
    ///        SCS-MP2 1.2 — same member name as the molecular methods).
    double mp2_os_scale = 1.0;
    /// @brief Grimme same-spin scale for `mp2()` (1.0 = plain MP2; SCS 1/3).
    double mp2_ss_scale = 1.0;
    /// @brief Optional per-atom AFM/FM spin seed for the initial UKS density
    ///        (units-free sign/magnitude; size 0 or n_atoms). Applied only to
    ///        the cold-start densities, not as an external potential.
    std::vector<double> atom_spin_guess;

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
    /// @param guess_P Optional RKS initial-density override (n_ao × n_ao) —
    ///        the XL-BOMD auxiliary density. Diagonalization/DIIS RKS path
    ///        only; non-null with UKS or `use_ot` throws (no-fallback), and a
    ///        shape mismatch throws inside the solve.
    /// @return The owned result of this run (valid until the next `run()`).
    /// @throws std::logic_error when called before `build()`.
    /// @throws std::invalid_argument for `guess_P` with UKS or `use_ot`.
    GpwResult& run(const std::vector<GpwAtom>& atoms,
                   const linalg::RealMatrix* guess_C = nullptr,
                   const linalg::RealMatrix* guess_P = nullptr);

    /// @brief Result of the most recent `run()` (default-constructed —
    ///        `converged == false` — before the first run).
    const GpwResult& last_result() const { return result_; }

    /// @brief Γ-point canonical GPW MP2 correlation on the converged
    ///        closed-shell orbitals of the last `run()`.
    ///
    /// The occupied–virtual (ia|jb) integrals are computed on the SCF's
    /// density grid: collocate ψ_i ψ_a, one Poisson solve per pair
    /// (periodic FFT, G=0 dropped — divergence-free since S_ia = 0 —
    /// or the truncated-Coulomb isolated solver when `nonperiodic`),
    /// grid quadrature against ψ_j ψ_b (CP2K `&WF_CORRELATION` MP2
    /// `DIRECT_CANONICAL`; Del Ben, Hutter, VandeVondele, JCTC 8, 4177
    /// (2012)). Virtual space capped by `mp2_n_virt`; `e_corr` applies
    /// `mp2_os_scale` / `mp2_ss_scale`.
    /// @param atoms Same geometry [Bohr] the last `run()` converged on.
    /// @return The owned MP2 result (valid until the next `mp2()`).
    /// @throws std::logic_error before `build()` or a converged `run()`.
    /// @throws std::invalid_argument for UKS, for `use_ot` (no canonical
    ///         virtuals), when no virtuals remain, or for a vanishing
    ///         HOMO–LUMO gap (< 1e-6 Ha; degenerate denominators).
    const Mp2Result& mp2(const std::vector<GpwAtom>& atoms);

  private:
    GpwResult result_;    ///< Owned result buffer (last run).
    Mp2Result mp2_result_; ///< Owned MP2 result buffer (last mp2()).
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

/// @brief Cell stress σ [Hartree/Bohr³], Voigt order xx,yy,zz,xy,xz,yz.
///
/// v1: orthorhombic strain-FD of the SCF energy at fixed fractional
/// coordinates (diagonal only; off-diagonal components are 0). Same strategy
/// as `cpu::pw::pw_stress` strain-FD foothold. Analytic stress is a follow-on.
///
/// @param atoms Geometry [Bohr] + basis/PP (as for `GpwDft::run`).
/// @param method SCF configuration; `cell` must match the atom box [Bohr].
/// @param h Dimensionless strain step (default 2e-3).
/// @param axes Per-axis mask (x,y,z); a masked axis skips its two FD SCF
///        solves and reports exactly 0 — for variable-cell relaxation of a
///        subset of the box edges (CP2K `CELL_OPT` `CONSTRAINT` twin).
/// @param guess_C Optional converged unstrained MO coefficients used to warm
///        start every strained SCF (pure cost cut; the strain is O(h)).
/// @return Voigt stress [Hartree/Bohr³].
std::array<double, 6>
gpw_stress(const std::vector<GpwAtom>& atoms, const GpwDft& method,
           double h = 2e-3, const std::array<bool, 3>& axes = {true, true, true},
           const linalg::RealMatrix* guess_C = nullptr);

} // namespace gpw
} // namespace cpu
} // namespace ATV
