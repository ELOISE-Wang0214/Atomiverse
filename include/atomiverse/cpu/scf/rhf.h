/// @file cpu/scf/rhf.h
/// @brief Restricted Hartree–Fock (RHF) concrete mean-field method.
#pragma once

#include <string>

#include "atomiverse/cpu/scf/mean_field_method.h"

namespace ATV {
namespace cpu {
namespace scf {

/// @brief Result of the Seeger-Pople stability analysis of a converged RHF
///        solution (Seeger & Pople, J. Chem. Phys. 66, 3045 (1977),
///        DOI 10.1063/1.434318).
///
/// Internal = RHF→RHF (singlet channel of the orbital Hessian); external =
/// RHF→UHF spin-symmetry breaking (triplet channel). The solution is stable
/// in a channel iff the lowest eigenvalue of that channel's (A + B) block is
/// non-negative.
struct StabilityResult {
    /// @brief Lowest eigenvalue of the singlet (A + B) block [Hartree].
    double internal_eigenvalue = 0.0;
    /// @brief Lowest eigenvalue of the triplet (A + B) block [Hartree].
    double external_eigenvalue = 0.0;
    /// @brief Internal (RHF→RHF) stability: internal_eigenvalue ≥ 0 up to
    ///        the numerical-zero guard (dimensionless flag).
    bool internal_stable = true;
    /// @brief External (RHF→UHF triplet) stability: external_eigenvalue ≥ 0
    ///        up to the numerical-zero guard (dimensionless flag).
    bool external_stable = true;
    /// @brief MO coefficients rotated along the external instability mode
    ///        [n × n_mo, dimensionless] — a broken-symmetry α-channel
    ///        initial guess for UHF (bridge to `scf_guess.h`: build a
    ///        per-spin density P_α = C_occ C_occᵀ from these columns and
    ///        pass it via `SCFInitialGuess`). Equals the unrotated C when
    ///        the solution is externally stable.
    RealMatrix rotated_C;
};

/// @brief Restricted Hartree–Fock — a mean-field electronic-structure method
///        (a single Slater determinant in which every electron moves in the
///        averaged field of all the others) whose spin-up (α) and spin-down (β)
///        electrons share one common set of spatial orbitals: the closed-shell,
///        fully paired case.
class RHF : public MeanFieldMethod<RHF> {
  public:
    /// @brief Construct with the RHF defaults (restricted HF; MINAO initial
    ///        guess; no density damping; energy tolerance 1e-9 Ha).
    RHF();
    std::string name() const { return "rhf"; }

    /// @brief Seeger-Pople stability analysis of the most recent `run()`:
    ///        lowest eigenvalue of the singlet (internal, RHF→RHF) and
    ///        triplet (external, RHF→UHF) (A + B) orbital-Hessian blocks,
    ///        plus a broken-symmetry rotated coefficient matrix when the
    ///        solution is externally unstable. The dense working-dimension
    ///        ERI tensor is rebuilt from the private shell workspace.
    /// @return Eigenvalues [Hartree], stability flags, and `rotated_C`.
    /// @throws std::runtime_error when called before `run()`, or
    ///         std::invalid_argument when the system has no virtual
    ///         orbitals (e.g. He/STO-3G — no fallback).
    StabilityResult stability() const;

  protected:
    /// @brief Closed-shell occupation: n_alpha = n_beta = n_electrons / 2
    ///        (every occupied spatial orbital holds one α and one β electron).
    void resolve_occupation(int n_electrons) override;
};

} // namespace scf
} // namespace cpu
} // namespace ATV
