/// @file cpu/scf/scf_types.h
/// @brief Configuration and result types for the unified CPU SCF solver.
///
/// Factorized design:
///   - SpinType: Restricted (1 channel) or Unrestricted (2 channels)
///   - MeanFieldType: HF (J - K) or KS (J + V_xc)
///   - SCFConfig: all parameters in one struct
///   - SCFResult: full energy decomposition, spin-resolved orbitals, diagnostics
///
/// @see cpu/linalg/{eigensolve,density}.h + cpu/scf/{scf_guess,spin_analysis}.h
/// @see cpu::scf::operators.h for J/K operator builders
/// @see cpu::scf::diis.h for DIIS accelerator
#pragma once

#include <stdexcept>
#include <string>
#include <vector>
#include <xtensor/containers/xtensor.hpp>

#include "atomiverse/cpu/scf/grid.h"
#include "atomiverse/cpu/scf/integrals.h"
#include "atomiverse/cpu/scf/xc_functional.h"

namespace ATV {
namespace cpu {
namespace scf {

/// @brief Per-iteration data for convergence analysis.
struct IterationRecord {
    int iter = 0;             ///< Iteration number (1-based).
    double energy = 0;        ///< Total energy at this iteration [Hartree].
    double residual = 0;      ///< Convergence residual.
    bool diis_active = false; ///< Whether DIIS was used this iteration.
};

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

    /// @brief Virtual-orbital level shift λ [Hartree] (Saunders–Hillier).
    ///        Applied inside the SCF loop with PySCF `get_fock` parity:
    ///        damping → DIIS extrapolation of the UNSHIFTED Fock → shift
    ///        F' = F + λ(S − S·D_occ·S) → diagonalize. Virtual eigenvalues
    ///        shift up by λ, occupied ones stay put — widening the
    ///        HOMO–LUMO gap to stabilize oscillating SCF. After
    ///        convergence one extra λ=0 diagonalization restores unshifted
    ///        reported orbital energies. |λ| <= 1e-4 is skipped entirely
    ///        (PySCF rule). Default 0.0 — behavior identical to before the
    ///        feature. Limitation: one scalar λ for both spin channels
    ///        (PySCF also allows a per-spin (λα, λβ) pair).
    double level_shift = 0.0;

    /// @brief Linear-dependence threshold on the raw overlap eigenvalues
    ///        (dimensionless). Eigen-directions of S with s_i below this
    ///        are dropped via canonical orthogonalization and the SCF
    ///        diagonalization switches to the two-step XᵀFX path
    ///        (n_mo = m <= n_ao). `0` disables filtering and keeps the
    ///        legacy Löwdin + dsygvd path bit-identical to prior behavior.
    ///        Healthy bases (nothing dropped) also stay on the fast path.
    ///        Note: PySCF >= 2.13 applies 1e-6 to the S eigenvalues by
    ///        default; Atomiverse uses 1e-8 on the raw (un-normalized)
    ///        eigenvalues.
    double lindep_threshold = 1e-8;

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
    ///          pointee must outlive the SCFSolver::run call. Not owned.
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

    // Orbitals and densities (always spin-resolved). When canonical
    // orthogonalization drops near-linear-dependent overlap directions
    // (SCFConfig::lindep_threshold), n_mo = m <= n_ao: the eigenvalue
    // vectors have length m and the C matrices are rectangular n_ao × m.
    // The density matrices stay n_ao × n_ao.
    xt::xtensor<double, 1>
        eigenvalues_alpha; ///< Alpha orbital energies [Hartree], length n_mo.
    xt::xtensor<double, 1>
        eigenvalues_beta; ///< Beta orbital energies [Hartree], length n_mo.
    xt::xtensor<double, 2>
        P_alpha; ///< Final alpha density matrix (n_ao × n_ao).
    xt::xtensor<double, 2> P_beta; ///< Final beta density matrix (n_ao × n_ao).
    xt::xtensor<double, 2>
        C_alpha; ///< Final alpha MO coefficients (n_ao × n_mo).
    xt::xtensor<double, 2>
        C_beta; ///< Final beta MO coefficients (n_ao × n_mo).

    // Diagnostics
    std::vector<IterationRecord> trace; ///< Per-iteration convergence trace.

    // Grid quality (KS only)
    double integrated_n_alpha = 0; ///< Grid-integrated alpha electron count.
    double integrated_n_beta = 0;  ///< Grid-integrated beta electron count.
};

/// @brief Occupation resolver: (electron count, spin multiplicity) →
///        (n_alpha, n_beta). Direct instantiation — no factory function.
///
/// Implements N_e = Σ_A Z_A − q (charge already folded into
/// `IntegralPack::n_electrons` by `build_integrals`), M = 2S + 1:
///   n_alpha = (N_e + M − 1) / 2,   n_beta = (N_e − M + 1) / 2.
///
/// Feed the resolved counts to the existing raw-count config builders:
/// @code
///   SpinResolution occ(ints.n_electrons, /*multiplicity=*/2);
///   auto cfg = uhf_config(ints, occ.n_alpha, occ.n_beta);
/// @endcode
///
/// Matches PySCF's convention contrast: PySCF's `mol.spin` is 2S =
/// n_alpha − n_beta, NOT the multiplicity M.
struct SpinResolution {
    int n_alpha; ///< Alpha electrons (dimensionless count).
    int n_beta;  ///< Beta electrons (dimensionless count).

    /// @brief Resolve and validate the occupation split.
    /// @param n_electrons Total electron count N_e (dimensionless, >= 0).
    /// @param multiplicity Spin multiplicity M = 2S + 1 (dimensionless, >= 1).
    /// @throws std::invalid_argument when n_electrons < 0, multiplicity < 1,
    ///         N_e and M − 1 differ in parity (n_alpha/n_beta would be
    ///         non-integer), or n_beta would be negative (M > N_e + 1).
    SpinResolution(int n_electrons, int multiplicity) {
        if (n_electrons < 0 || multiplicity < 1
            || (n_electrons + multiplicity - 1) % 2 != 0
            || n_electrons - multiplicity + 1 < 0) {
            throw std::invalid_argument(
                "SpinResolution: invalid electron count / multiplicity "
                "combination (n_electrons="
                + std::to_string(n_electrons)
                + ", multiplicity=" + std::to_string(multiplicity)
                + "): requires n_electrons >= 0, multiplicity >= 1, "
                  "matching parity (even N_e with odd M or odd N_e with "
                  "even M), and multiplicity <= n_electrons + 1");
        }
        n_alpha = (n_electrons + multiplicity - 1) / 2;
        n_beta = (n_electrons - multiplicity + 1) / 2;
    }
};

/// @brief Convenience: RHF configuration.
/// @note Odd-electron packs are rejected at the solve entry:
///       `SCFSolver::run` throws `std::invalid_argument` because
///       `n_occ() = n_electrons / 2` truncates, so
///       n_alpha + n_beta == n_electrons − 1 fails the occupation guard.
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
/// @note HOMO-LUMO mixing can destabilize DIIS before the density is near
///       stationary, so UHF uses delayed DIIS and heavy pre-DIIS damping.
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
    // UHF-specific delayed DIIS and heavy pre-DIIS damping.
    c.diis_max_vectors = 6;
    c.diis_start_iter = 5; // allow HOMO-LUMO mixing to settle before DIIS
    c.damping_factor = 0.7;
    c.init_guess = InitGuess::Hcore;
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
    // UKS-specific early DIIS and symmetric density mixing.
    c.diis_max_vectors = 6;
    c.diis_start_iter = 1;
    c.damping_factor = 0.5; // symmetric mix — invariant under convention flip
    c.init_guess = InitGuess::Hcore;
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
///   - spin: Restricted (closed-shell). Odd-electron packs are rejected at
///     the solve entry: `SCFSolver::run` throws `std::invalid_argument`
///     when n_alpha + n_beta != n_electrons (the truncated
///     `n_occ() = n_electrons / 2` sums to n_electrons − 1).
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
/// @param ints Integral pack supplying the electron count (must be even —
///        odd counts throw `std::invalid_argument` from `SCFSolver::run`).
/// @return SCFConfig ready to feed `SCFSolver`.
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
/// @return SCFConfig ready to feed `SCFSolver`.
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
/// @return SCFConfig ready to feed `SCFSolver`.
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

} // namespace scf
} // namespace cpu
} // namespace ATV
