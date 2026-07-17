/// @file cpu/scf/internal/scf_types.h
/// @brief **PRIVATE** SCFConfig + SpinResolution. Results live in public scf_result.h.
#pragma once

#include <stdexcept>
#include <string>

#include "atomiverse/cpu/scf/grid.h"
#include "atomiverse/cpu/scf/internal/integrals.h"
#include "atomiverse/cpu/scf/scf_result.h"
#include "atomiverse/cpu/scf/xc_functional.h"

namespace ATV {
namespace cpu {
namespace scf {

struct SCFConfig {
    /// @brief Restricted Hartree–Fock preset tag.
    struct RHF {};
    /// @brief Unrestricted Hartree–Fock preset tag.
    struct UHF {};
    /// @brief Unrestricted Kohn–Sham preset tag.
    struct UKS {};
    /// @brief Restricted Kohn–Sham / PBE production default preset tag.
    struct RKS {};
    /// @brief Closed-shell B3LYP hybrid preset tag.
    struct B3LYP {};
    /// @brief Closed-shell PBE0 hybrid preset tag.
    struct PBE0 {};
    /// @brief Fast RHF smoke/regression baseline preset tag.
    struct Baseline {};

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
    ///          by `SCFConfig(RKS{}, ...)`. When `grid_level` is not `SG0`,
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

    /// @brief Default config (field defaults only; not a method preset).
    SCFConfig() = default;

    /// @brief Restricted Hartree–Fock.
    /// @note Odd-electron packs are rejected at the solve entry:
    ///       `SCFSolver::run` throws because `n_occ() = n_electrons / 2`
    ///       truncates.
    /// @param ints Integral pack supplying electron count.
    /// @param tol  Density convergence tolerance [Hartree].
    SCFConfig(RHF, const IntegralPack& ints, double tol = 1e-9) {
        spin = SpinType::Restricted;
        method = MeanFieldType::HF;
        n_alpha = ints.n_occ();
        n_beta = ints.n_occ();
        this->tol = tol;
        damping_factor =
            0.0; // MINAO guess lets us drop damping (PySCF default).
        init_guess = InitGuess::MINAO;
    }

    /// @brief Unrestricted Hartree–Fock.
    /// @note HOMO-LUMO mixing can destabilize DIIS before the density is near
    ///       stationary, so UHF uses delayed DIIS and heavy pre-DIIS damping.
    /// @param n_alpha Number of alpha electrons.
    /// @param n_beta  Number of beta electrons.
    /// @param tol     Density convergence tolerance [Hartree].
    SCFConfig(UHF, int n_alpha, int n_beta, double tol = 1e-8) {
        spin = SpinType::Unrestricted;
        method = MeanFieldType::HF;
        this->n_alpha = n_alpha;
        this->n_beta = n_beta;
        this->tol = tol;
        use_diis = true;
        diis_max_vectors = 6;
        diis_start_iter = 5; // allow HOMO-LUMO mixing to settle before DIIS
        damping_factor = 0.7;
        init_guess = InitGuess::Hcore;
    }

    /// @brief Unrestricted Kohn–Sham DFT.
    /// @param n_alpha Number of alpha electrons.
    /// @param n_beta  Number of beta electrons.
    /// @param xc      Exchange-correlation functional.
    /// @param tol     Density convergence tolerance [Hartree].
    SCFConfig(UKS, int n_alpha, int n_beta, XCType xc = XCType::LSDA,
              double tol = 1e-8) {
        spin = SpinType::Unrestricted;
        method = MeanFieldType::KS;
        xc_type = xc;
        this->n_alpha = n_alpha;
        this->n_beta = n_beta;
        this->tol = tol;
        use_diis = true;
        diis_max_vectors = 6;
        diis_start_iter = 1;
        damping_factor = 0.5; // symmetric mix — invariant under convention flip
        init_guess = InitGuess::Hcore;
    }

    /// @brief Locked production default: closed-shell RKS / PBE.
    ///
    /// Defaults:
    ///   - spin: Restricted; method: KS; xc_type: PBE
    ///   - init_guess: MINAO; grid: GridLevel::MatchPySCF
    ///   - tol = 1e-9 Ha; DIIS on from iteration 1
    ///
    /// @param ints Integral pack supplying the electron count (must be even —
    ///        odd counts throw from `SCFSolver::run`).
    SCFConfig(RKS, const IntegralPack& ints) {
        spin = SpinType::Restricted;
        method = MeanFieldType::KS;
        xc_type = XCType::PBE;
        n_alpha = ints.n_occ();
        n_beta = ints.n_occ();
        tol = 1e-9;
        use_diis = true;
        diis_max_vectors = 8;
        diis_start_iter = 1;
        damping_factor = 0.0;
        init_guess = InitGuess::MINAO;
        grid_level = GridLevel::MatchPySCF;
    }

    /// @brief Closed-shell B3LYP global-hybrid KS (RKS defaults + B3LYP).
    /// @param ints Integral pack supplying the electron count (assumed even).
    SCFConfig(B3LYP, const IntegralPack& ints) : SCFConfig(RKS{}, ints) {
        xc_type = XCType::B3LYP;
    }

    /// @brief Closed-shell PBE0 global-hybrid KS (RKS defaults + PBE0).
    /// @param ints Integral pack supplying the electron count (assumed even).
    SCFConfig(PBE0, const IntegralPack& ints) : SCFConfig(RKS{}, ints) {
        xc_type = XCType::PBE0;
    }

    /// @brief Regression baseline: RHF / fast convergence, for STO-3G.
    /// @param ints Integral pack supplying the electron count.
    SCFConfig(Baseline, const IntegralPack& ints) {
        spin = SpinType::Restricted;
        method = MeanFieldType::HF;
        n_alpha = ints.n_occ();
        n_beta = ints.n_occ();
        tol = 1e-9;
        use_diis = true;
        damping_factor = 0.0;
        init_guess = InitGuess::MINAO;
        max_iter = 50;
    }
};

/// @brief Occupation resolver: (electron count, spin multiplicity) →
///        (n_alpha, n_beta). Direct instantiation — no factory function.
///
/// Implements N_e = Σ_A Z_A − q (charge already folded into
/// `IntegralPack::n_electrons` by its Frame constructor), M = 2S + 1:
///   n_alpha = (N_e + M − 1) / 2,   n_beta = (N_e − M + 1) / 2.
///
/// Feed the resolved counts to the UHF/UKS constructors:
/// @code
///   SpinResolution occ(ints.n_electrons, /*multiplicity=*/2);
///   SCFConfig cfg{SCFConfig::UHF{}, occ.n_alpha, occ.n_beta};
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

} // namespace scf
} // namespace cpu
} // namespace ATV
