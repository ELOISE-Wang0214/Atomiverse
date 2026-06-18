/// @file cpu/scf/scf_solve.h
/// @brief Unified mean-field SCF solver -- RHF, UHF, UKS through one engine.
///
/// Factorized design:
///   - SpinType: Restricted (1 channel) or Unrestricted (2 channels)
///   - MeanFieldType: HF (J - K) or KS (J + V_xc)
///   - SCFConfig: all parameters in one struct
///   - SCFResult: full energy decomposition, spin-resolved orbitals, diagnostics
///
/// The method-specific convenience solvers (rhf_solve, uhf_solve, uks_solve)
/// are thin wrappers around scf_solve() that preset the matching SCFConfig
/// and repack into the method-specific result struct.
///
/// @see cpu/scf/scf_linalg.h and cpu/linalg/dense.h for shared LA
/// @see cpu::scf::operators.h for J/K operator builders
/// @see cpu::scf::diis.h for DIIS accelerator
#pragma once

#include <vector>
#include <xtensor/containers/xtensor.hpp>

#include "atomiverse/cpu/scf/grid.h"
#include "atomiverse/cpu/scf/integrals.h"
#include "atomiverse/cpu/scf/rhf_solve.h"
#include "atomiverse/cpu/scf/xc_functional.h"

namespace ATV {
namespace cpu {
namespace scf {

/// @brief Mean-field method type.
enum class MeanFieldType { HF, KS };

/// @brief Spin treatment.
enum class SpinType { Restricted, Unrestricted };

/// @brief Initial-density guess type.
///   Hcore: diagonalize H = T + V (core Hamiltonian). Simple but poor for
///          molecules -- typically needs ~3x more SCF iterations.
///   MINAO: project atomic fragment densities in an ANO-trimmed minimal basis
///          onto the target AO basis. Matches PySCF's default and yields
///          fast, reliable convergence. Requires all atoms be in the MINAO
///          table (H, C, N, O, F at present); unsupported atoms fall back
///          to Hcore automatically.
enum class InitGuess { Hcore, MINAO };

/// @brief SCF configuration -- all method/convergence/grid parameters.
struct SCFConfig {
    SpinType spin =
        SpinType::Unrestricted; ///< Restricted or unrestricted spin treatment.
    MeanFieldType method =
        MeanFieldType::HF;         ///< Hartree-Fock or Kohn-Sham mean field.
    XCType xc_type = XCType::LSDA; ///< Only used when method == KS.
    int n_alpha = 0;               ///< Alpha electrons.
    int n_beta = 0;                ///< Beta electrons.

    // Convergence (PySCF-aligned defaults for parameter values;
    // damping retained at 0.5 until MINAO lands — Hcore guess + damping=0
    // + DIIS-from-iter-1 diverges on aspirin-scale systems. Drop to 0.0
    // once MINAO is available.)
    double tol = 1e-9;         ///< Energy convergence tolerance |dE| (Hartree).
    double conv_tol_grad = -1; ///< DIIS-error tolerance. If <0, uses sqrt(tol).
    int max_iter = 100;        ///< Maximum SCF iterations.
    bool use_diis = true;      ///< Enable DIIS acceleration.
    int diis_max_vectors = 8;  ///< DIIS history size (PySCF default = 8).
    int diis_start_iter =
        1; ///< Start DIIS at this iteration (PySCF default = 1).
    double damping_factor =
        0.5; ///< Linear mixing for pre-DIIS iterations. Auto-dropped
             ///< to 0 when MINAO guess is available (matches PySCF).
    InitGuess init_guess = InitGuess::MINAO; ///< Initial-density guess type.

    // Grid (KS only)
    /// @brief Named DFT integration-grid quality level.
    /// @details Selects the radial scheme + angular pruning. `SG0` (the
    ///          struct default) is the legacy uniform-order path that
    ///          honours `n_radial` / `n_angular` verbatim. `SG1` is the
    ///          Gill-Johnson-Pople standard grid; `MatchPySCF` reproduces
    ///          PySCF's `level=3` grid and is the production default set
    ///          by `default_config`. When `grid_level` is not `SG0`,
    ///          `n_radial` / `n_angular` are ignored.
    GridLevel grid_level = GridLevel::SG0;
    int n_radial = 50;   ///< Radial grid points per atom (SG0 only).
    int n_angular = 110; ///< Angular grid points per atom (SG0 only).

    /// @brief Force the native 2e gradient engine even when libcint is available.
    /// @details libcint (when compiled in) is the STO-3G path for
    ///          `accumulate_two_electron_gradient`. Setting this to `true`
    ///          routes through the native Rys-quadrature-based derivative
    ///          instead. Reserved for tests that compare native-vs-libcint
    ///          outputs (see `test_gradient_vs_libcint.cpp`) and benchmarks
    ///          that time the native path. Has no effect on the Fock build,
    ///          which is always native.
    bool force_native_two_electron_deriv = false;

    /// @brief Debug flag: project out the net translational XC gradient.
    /// @details Historically the KS XC nuclear gradient was computed on a
    ///          grid held fixed in space, then the residual net-translational
    ///          quadrature error (~4e-4 Ha/Bohr) was projected out to restore
    ///          `sum_A dE_xc/dR_A = 0`. With the genuine Becke grid-response
    ///          term (weight-derivative + grid-coordinate) now assembled in
    ///          `accumulate_xc_gradient`, the net translational force vanishes
    ///          from first principles and the projection is unnecessary. It
    ///          is kept behind this flag (default `false`) purely as a debug
    ///          / comparison aid; production code leaves it `false`.
    bool xc_grad_project_translation = false;

    /// @brief Regression override for the exact-exchange fraction (KS hybrids).
    /// @details When non-negative, this value replaces
    ///          `hybrid_exchange_fraction(xc_type)` as the effective
    ///          `alpha_hf` AND forces the hybrid Fock branch (K build) to
    ///          run even for a pure functional. Default `-1` means "use the
    ///          functional's own alpha_hf". Reserved for the ac-004
    ///          regression test (`alpha_hf=0` hybrid path must reproduce the
    ///          pure-DFT result bit-for-bit); production code leaves it `-1`.
    double hybrid_alpha_override = -1.0;

    /// @brief Optional externally supplied XC integration grid (KS only).
    /// @details When non-null, the KS Fock build uses this grid verbatim
    ///          instead of regenerating one from the current atom
    ///          positions. The finite-difference KS gradient oracle uses
    ///          it to hold the grid FIXED across the +/- displacements so
    ///          the FD result matches the fixed-grid analytic XC gradient
    ///          (Johnson-Gill-Pople convention; see gradient.cpp). The
    ///          pointee must outlive the scf_solve call. Not owned.
    const MolecularGrid* fixed_grid = nullptr;
};

/// @brief Optional density-matrix initial guess for SCF iterations.
struct SCFInitialGuess {
    xt::xtensor<double, 2> P_alpha; ///< Alpha or total density guess.
    xt::xtensor<double, 2>
        P_beta;             ///< Beta density guess for unrestricted runs.
    bool has_alpha = false; ///< Whether P_alpha should be used.
    bool has_beta = false;  ///< Whether P_beta should be used.
};

/// @brief Unified SCF result -- full energy decomposition and spin-resolved quantities.
struct SCFResult {
    // Core
    bool converged = false; ///< Whether SCF converged within `max_iter`.
    int iterations = 0;     ///< Number of iterations performed.
    double energy = 0;      ///< Total energy [Hartree].
    double residual = 0;    ///< Final convergence metric.

    // Energy decomposition (always populated)
    double one_electron_energy = 0; ///< sum_sigma Tr[Hcore * P^sigma].
    double coulomb_energy = 0;      ///< 0.5 * Tr[J * P_total].
    double exchange_energy = 0;     ///< HF exchange. Zero for pure KS.
    double xc_energy = 0;           ///< E_xc from grid. Zero for pure HF.
    double nuclear_repulsion = 0;   ///< V_nuc.

    // Spin info
    int n_alpha = 0;               ///< Number of alpha electrons.
    int n_beta = 0;                ///< Number of beta electrons.
    double spin_contamination = 0; ///< <S^2> - S(S+1). Zero for restricted.

    // Orbitals and densities (always spin-resolved)
    xt::xtensor<double, 1>
        eigenvalues_alpha; ///< Alpha orbital energies [Hartree].
    xt::xtensor<double, 1>
        eigenvalues_beta;           ///< Beta orbital energies [Hartree].
    xt::xtensor<double, 2> P_alpha; ///< Final alpha density matrix.
    xt::xtensor<double, 2> P_beta;  ///< Final beta density matrix.
    xt::xtensor<double, 2> C_alpha; ///< Final alpha MO coefficients.
    xt::xtensor<double, 2> C_beta;  ///< Final beta MO coefficients.

    // Diagnostics
    std::vector<IterationRecord> trace; ///< Per-iteration convergence trace.

    // Grid quality (KS only)
    double integrated_n_alpha = 0; ///< Grid-integrated alpha electron count.
    double integrated_n_beta = 0;  ///< Grid-integrated beta electron count.
};

/// @brief Convenience: RHF configuration.
/// @param ints Integral pack supplying electron count.
/// @param tol Density convergence tolerance.
/// @return `SCFConfig` initialized for restricted Hartree-Fock.
inline SCFConfig rhf_config(const IntegralPack& ints, double tol = 1e-9) {
    SCFConfig c;
    c.spin = SpinType::Restricted;
    c.method = MeanFieldType::HF;
    c.n_alpha = ints.n_occ();
    c.n_beta = ints.n_occ();
    c.tol = tol;
    c.damping_factor = 0.0; // MINAO guess lets us drop damping (PySCF default).
    c.init_guess = InitGuess::MINAO;
    return c;
}

/// @brief Convenience: UHF configuration.
/// @note DIIS disabled by default for UHF -- HOMO-LUMO mixing can destabilize DIIS
///       before the density is near stationary. Use damped SCF (matching original
///       uhf_solve behavior). Enable DIIS manually for advanced usage.
/// @param ints Integral pack retained for API symmetry with `rhf_config`.
/// @param n_alpha Number of alpha electrons.
/// @param n_beta Number of beta electrons.
/// @param tol Density convergence tolerance.
/// @return `SCFConfig` initialized for unrestricted Hartree-Fock.
inline SCFConfig uhf_config(const IntegralPack& ints, int n_alpha, int n_beta,
                            double tol = 1e-8) {
    (void)ints;
    SCFConfig c;
    c.spin = SpinType::Unrestricted;
    c.method = MeanFieldType::HF;
    c.n_alpha = n_alpha;
    c.n_beta = n_beta;
    c.tol = tol;
    c.use_diis = true;
    // UHF-specific: preserve legacy behaviour (unchanged in Phase 1 PySCF alignment).
    // NOTE: damping convention was flipped in 2026-04; legacy 0.3 (30% new) becomes 0.7 here.
    c.diis_max_vectors = 6;
    c.diis_start_iter = 5; // allow HOMO-LUMO mixing to settle before DIIS
    c.damping_factor =
        0.7; // heavy pre-DIIS damping (matches legacy 0.3 in old convention)
    c.init_guess = InitGuess::Hcore; // keep legacy UHF guess
    return c;
}

/// @brief Convenience: UKS configuration.
/// @param ints Integral pack retained for API symmetry with `rhf_config`.
/// @param n_alpha Number of alpha electrons.
/// @param n_beta Number of beta electrons.
/// @param xc Exchange-correlation functional.
/// @param tol Density convergence tolerance.
/// @return `SCFConfig` initialized for unrestricted Kohn-Sham DFT.
inline SCFConfig uks_config(const IntegralPack& ints, int n_alpha, int n_beta,
                            XCType xc = XCType::LSDA, double tol = 1e-8) {
    (void)ints;
    SCFConfig c;
    c.spin = SpinType::Unrestricted;
    c.method = MeanFieldType::KS;
    c.xc_type = xc;
    c.n_alpha = n_alpha;
    c.n_beta = n_beta;
    c.tol = tol;
    c.use_diis = true;
    // UKS-specific: preserve legacy behaviour (unchanged in Phase 1 PySCF alignment).
    c.diis_max_vectors = 6;
    c.diis_start_iter = 1;
    c.damping_factor = 0.5; // symmetric mix — invariant under convention flip
    c.init_guess = InitGuess::Hcore; // keep legacy UKS guess
    return c;
}

/// @brief Locked production default: closed-shell RKS / PBE.
///
/// Intended for the full benchmark / tutorial stack. Basis-set choice is
/// independent: the caller builds `IntegralPack` at whatever basis they
/// want (recommended `BasisSet::Def2SVP` for production; `BasisSet::Pople631G`
/// for quick tests until the d-shell XC-grid basis_eval ships).
///
/// Defaults:
///   - spin: Restricted (closed-shell). Throws-via-solver if n_electrons is odd.
///   - method: KS (Kohn-Sham)
///   - xc_type: PBE
///   - init_guess: MINAO (PySCF-matching, zero damping)
///   - grid: GridLevel::MatchPySCF -- the high-accuracy grid that
///     reproduces PySCF's default `level=3` (Treutler-Ahlrichs M4 radial
///     with per-element point counts, nwchem-style angular pruning). The
///     level-3 maximum Lebedev order is 302 for Z<=10 (pruned to
///     {50,86,266,302,266}) and 434 for Z=11..36 (pruned to
///     {50,86,350,434,350}). This brings the ATV-vs-PySCF def2-SVP energy
///     residue into the sub-microHartree regime; the lighter `SG1`
///     standard grid and the legacy uniform `SG0` grid remain selectable
///     via SCFConfig::grid_level.
///   - convergence: tol = 1e-9 Ha
///
/// @param ints Integral pack supplying the electron count (assumed even).
/// @return SCFConfig ready to feed `scf_solve`.
inline SCFConfig default_config(const IntegralPack& ints) {
    SCFConfig c;
    c.spin = SpinType::Restricted;
    c.method = MeanFieldType::KS;
    c.xc_type = XCType::PBE;
    c.n_alpha = ints.n_occ();
    c.n_beta = ints.n_occ();
    c.tol = 1e-9;
    c.use_diis = true;
    c.diis_max_vectors = 8;
    c.diis_start_iter = 1;
    c.damping_factor = 0.0;
    c.init_guess = InitGuess::MINAO;
    c.grid_level = GridLevel::MatchPySCF;
    return c;
}

/// @brief Convenience: closed-shell B3LYP global-hybrid KS configuration.
///
/// Clones `default_config` (Restricted, KS, MINAO guess,
/// `GridLevel::MatchPySCF`, `tol = 1e-9`, DIIS on) but selects the B3LYP
/// hybrid functional. B3LYP mixes `alpha_hf = 0.20` exact HF exchange with
/// the semi-local mix `0.08*Slater-X + 0.72*B88-X + 0.19*VWN_RPA-C +
/// 0.81*LYP-C` (libxc `XC_HYB_GGA_XC_B3LYP`).
///
/// @param ints Integral pack supplying the electron count (assumed even).
/// @return SCFConfig ready to feed `scf_solve` / `SCFSolver`.
inline SCFConfig b3lyp_config(const IntegralPack& ints) {
    SCFConfig c = default_config(ints);
    c.xc_type = XCType::B3LYP;
    return c;
}

/// @brief Convenience: closed-shell PBE0 global-hybrid KS configuration.
///
/// Clones `default_config` (Restricted, KS, MINAO guess,
/// `GridLevel::MatchPySCF`, `tol = 1e-9`, DIIS on) but selects the PBE0
/// hybrid functional. PBE0 mixes `alpha_hf = 0.25` exact HF exchange with
/// the semi-local mix `0.75*PBE-X + 1.00*PBE-C` (libxc
/// `XC_HYB_GGA_XC_PBEH`).
///
/// @param ints Integral pack supplying the electron count (assumed even).
/// @return SCFConfig ready to feed `scf_solve` / `SCFSolver`.
inline SCFConfig pbe0_config(const IntegralPack& ints) {
    SCFConfig c = default_config(ints);
    c.xc_type = XCType::PBE0;
    return c;
}

/// @brief Regression baseline: RHF / fast convergence, for STO-3G.
///
/// Intended for smoke tests and gradient-self-consistency regressions
/// where turnaround time beats fidelity. Same MINAO+no-damping pattern
/// as `rhf_config` but with a looser iteration budget.
///
/// @param ints Integral pack supplying the electron count.
/// @return SCFConfig for restricted Hartree-Fock.
inline SCFConfig baseline_config(const IntegralPack& ints) {
    SCFConfig c;
    c.spin = SpinType::Restricted;
    c.method = MeanFieldType::HF;
    c.n_alpha = ints.n_occ();
    c.n_beta = ints.n_occ();
    c.tol = 1e-9;
    c.use_diis = true;
    c.damping_factor = 0.0;
    c.init_guess = InitGuess::MINAO;
    c.max_iter = 50;
    return c;
}

/// @brief Run unified SCF solver.
///
/// Dispatches internally based on SCFConfig:
///   - Restricted HF:    F = Hcore + G(P)        [fused J-0.5K, DIIS]
///   - Unrestricted HF:  F^sigma = Hcore + J - K^sigma   [combined DIIS]
///   - Unrestricted KS:  F^sigma = Hcore + J + V^sigma_xc [grid XC, combined DIIS]
///
/// @param integrals Molecular integrals (S, Hcore, ERI, V_nuc, atoms, shells).
/// @param config    SCF configuration (method, spin, convergence, grid).
/// @param guess     Optional density-matrix initial guess.
/// @return SCFResult with full energy decomposition and spin-resolved quantities.
SCFResult scf_solve(const IntegralPack& integrals, const SCFConfig& config,
                    const SCFInitialGuess* guess = nullptr);

} // namespace scf
} // namespace cpu
} // namespace ATV
