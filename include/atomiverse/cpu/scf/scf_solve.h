/// @file cpu/scf/scf_solve.h
/// @brief Unified mean-field SCF solver — RHF, UHF, UKS through one engine.
///
/// Factorized design:
///   - SpinType: Restricted (1 channel) or Unrestricted (2 channels)
///   - MeanFieldType: HF (J - K) or KS (J + V_xc)
///   - SCFConfig: all parameters in one struct
///   - SCFResult: full energy decomposition, spin-resolved orbitals, diagnostics
///
/// The three existing solvers (rhf_solve, uhf_solve, uks_solve) become thin
/// wrappers around scf_solve() for backward compatibility.
///
/// @see cpu::scf::scf_utils.h for shared SCF utilities
/// @see cpu::scf::operators.h for J/K operator builders
/// @see cpu::scf::diis.h for DIIS accelerator
#pragma once

#include <vector>
#include <xtensor/xtensor.hpp>

#include "atomiverse/core/real.h"
#include "atomiverse/core/int.h"
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

/// @brief SCF configuration — all method/convergence/grid parameters.
struct SCFConfig {
    SpinType spin = SpinType::Unrestricted;   ///< Restricted or unrestricted spin treatment.
    MeanFieldType method = MeanFieldType::HF; ///< Hartree-Fock or Kohn-Sham mean field.
    XCType xc_type = XCType::LSDA;          ///< Only used when method == KS.
    int n_alpha = 0;                         ///< Alpha electrons.
    int n_beta = 0;                          ///< Beta electrons.

    // Convergence
    Real tol = Real(1e-8);                   ///< Convergence tolerance.
    int max_iter = 100;                      ///< Maximum SCF iterations.
    bool use_diis = true;                    ///< Enable DIIS acceleration.
    int diis_max_vectors = 6;                ///< DIIS history size.
    int diis_start_iter = 2;                 ///< Start DIIS after this many damped iterations.
    Real damping_factor = Real(0.5);         ///< Linear mixing for early/non-DIIS iterations.

    // Grid (KS only)
    int n_radial = 50;                       ///< Radial grid points per atom.
    int n_angular = 110;                     ///< Angular grid points per atom.
};

/// @brief Unified SCF result — full energy decomposition and spin-resolved quantities.
struct SCFResult {
    // Core
    bool converged = false;           ///< Whether SCF converged within `max_iter`.
    int iterations = 0;               ///< Number of iterations performed.
    Real energy = 0;                 ///< Total energy [Hartree].
    Real residual = 0;               ///< Final convergence metric.

    // Energy decomposition (always populated)
    Real one_electron_energy = 0;    ///< Σ_σ Tr[Hcore · P^σ].
    Real coulomb_energy = 0;         ///< 0.5 · Tr[J · P_total].
    Real exchange_energy = 0;        ///< HF exchange. Zero for pure KS.
    Real xc_energy = 0;             ///< E_xc from grid. Zero for pure HF.
    Real nuclear_repulsion = 0;      ///< V_nuc.

    // Spin info
    int n_alpha = 0;                  ///< Number of alpha electrons.
    int n_beta = 0;                   ///< Number of beta electrons.
    Real spin_contamination = 0;     ///< ⟨S²⟩ - S(S+1). Zero for restricted.

    // Orbitals and densities (always spin-resolved)
    xt::xtensor<Real, 1> eigenvalues_alpha; ///< Alpha orbital energies [Hartree].
    xt::xtensor<Real, 1> eigenvalues_beta;  ///< Beta orbital energies [Hartree].
    xt::xtensor<Real, 2> P_alpha;           ///< Final alpha density matrix.
    xt::xtensor<Real, 2> P_beta;            ///< Final beta density matrix.
    xt::xtensor<Real, 2> C_alpha;           ///< Final alpha MO coefficients.
    xt::xtensor<Real, 2> C_beta;            ///< Final beta MO coefficients.

    // Diagnostics
    std::vector<IterationRecord> trace;     ///< Per-iteration convergence trace.

    // Grid quality (KS only)
    Real integrated_n_alpha = 0;            ///< Grid-integrated alpha electron count.
    Real integrated_n_beta = 0;             ///< Grid-integrated beta electron count.
};

/// @brief Convenience: RHF configuration.
/// @param ints Integral pack supplying electron count.
/// @param tol Density convergence tolerance.
/// @return `SCFConfig` initialized for restricted Hartree-Fock.
inline SCFConfig rhf_config(const IntegralPack& ints, Real tol = Real(1e-8)) {
    SCFConfig c;
    c.spin = SpinType::Restricted;
    c.method = MeanFieldType::HF;
    c.n_alpha = ints.n_occ();
    c.n_beta = ints.n_occ();
    c.tol = tol;
    return c;
}

/// @brief Convenience: UHF configuration.
/// @note DIIS disabled by default for UHF — HOMO-LUMO mixing can destabilize DIIS
///       before the density is near stationary. Use damped SCF (matching original
///       uhf_solve behavior). Enable DIIS manually for advanced usage.
/// @param ints Integral pack retained for API symmetry with `rhf_config`.
/// @param n_alpha Number of alpha electrons.
/// @param n_beta Number of beta electrons.
/// @param tol Density convergence tolerance.
/// @return `SCFConfig` initialized for unrestricted Hartree-Fock.
inline SCFConfig uhf_config(const IntegralPack& ints, int n_alpha, int n_beta,
                            Real tol = Real(1e-8)) {
    (void)ints;
    SCFConfig c;
    c.spin = SpinType::Unrestricted;
    c.method = MeanFieldType::HF;
    c.n_alpha = n_alpha;
    c.n_beta = n_beta;
    c.tol = tol;
    c.use_diis = false;
    return c;
}

/// @brief Convenience: UKS configuration.
/// @note DIIS disabled by default for UKS — same reason as UHF (HOMO-LUMO mixing).
/// @param ints Integral pack retained for API symmetry with `rhf_config`.
/// @param n_alpha Number of alpha electrons.
/// @param n_beta Number of beta electrons.
/// @param xc Exchange-correlation functional.
/// @param tol Density convergence tolerance.
/// @return `SCFConfig` initialized for unrestricted Kohn-Sham DFT.
inline SCFConfig uks_config(const IntegralPack& ints, int n_alpha, int n_beta,
                            XCType xc = XCType::LSDA, Real tol = Real(1e-8)) {
    (void)ints;
    SCFConfig c;
    c.spin = SpinType::Unrestricted;
    c.method = MeanFieldType::KS;
    c.xc_type = xc;
    c.n_alpha = n_alpha;
    c.n_beta = n_beta;
    c.tol = tol;
    c.use_diis = false;
    return c;
}

/// @brief Run unified SCF solver.
///
/// Dispatches internally based on SCFConfig:
///   - Restricted HF:    F = Hcore + G(P)        [fused J-0.5K, DIIS]
///   - Unrestricted HF:  F^σ = Hcore + J - K^σ   [combined DIIS]
///   - Unrestricted KS:  F^σ = Hcore + J + V^σ_xc [grid XC, combined DIIS]
///
/// @param integrals Molecular integrals (S, Hcore, ERI, V_nuc, atoms, shells).
/// @param config    SCF configuration (method, spin, convergence, grid).
/// @return SCFResult with full energy decomposition and spin-resolved quantities.
SCFResult scf_solve(const IntegralPack& integrals, const SCFConfig& config);

} // namespace scf
} // namespace cpu
} // namespace ATV
