/// @file cpu/scf/stages/diagonalize_stage.h
/// @brief cpu::scf::DiagonalizeStage — the `scf::Phase::diagonalize` component.
///
/// Wraps the closed numeric kernel `internal::EigenSolver` (generalized
/// symmetric eigenproblem F C = S C ε) as an `scf::Phase`-bound `cpu::Component`,
/// with iteration state flowing through NAMED constructor references — no
/// config struct, no workspace bag, no Context slice. `launch()` reproduces the
/// shipping solver's diagonalize step byte-for-byte (`scf_solver.cpp:145-147`
/// path selection, `:285-291` restricted β reuse):
///   - `n_dropped > 0` → `EigenSolver::step_orthogonal(F, X)` (canonical
///     orthogonalization dropped near-linear-dependent columns);
///   - otherwise       → `EigenSolver::step(F, S)` (healthy-basis dsygvd);
///   - restricted reuses the α eigenpair for β (a copy, matching the reference).
///
/// The `EigenSolver` stays a closed kernel held PRIVATELY here; it does NOT
/// inherit `cpu::Component`.
#pragma once

#include <string>
#include <tuple>
#include <vector>

#include "atomiverse/cpu/component/component.h"
#include "atomiverse/cpu/scf/dense_matrix.h"
#include "atomiverse/cpu/scf/internal/eigen_solver.h"
#include "atomiverse/facade/stage.h"

namespace ATV {
namespace cpu {
namespace scf {

/// @brief Diagonalize stage: solves the generalized eigenproblem for the
///        orbitals, bound at `scf::Phase::diagonalize`.
///
/// State is wired by reference at construction: the read inputs (F/S/X and the
/// `n_dropped` path selector) as `const&`, the write outputs (C/eig) as `&`.
/// xtensor object identity is stable, so per-run resize never invalidates the
/// bound references. `role() == 0` — it writes only driver-owned matrices, no
/// frame field.
class DiagonalizeStage : public cpu::Component {
  public:
    /// @brief Wire the stage to the method driver's diagonalize state.
    /// @param restricted True for RHF/RKS (β reuses the α eigenpair); false for
    ///        UHF/UKS (β diagonalized independently).
    /// @param F_alpha α Fock matrix [Hartree], n × n.
    /// @param F_beta β Fock matrix [Hartree], n × n (ignored when restricted).
    /// @param S Overlap matrix (dimensionless), n × n — used on the healthy path.
    /// @param X Orthogonalization transform (dimensionless), n × m — used on the
    ///        column-dropped path.
    /// @param n_dropped Number of near-linear-dependent columns removed by
    ///        canonical orthogonalization; > 0 selects `step_orthogonal`.
    /// @param C_alpha [out] α MO coefficients (dimensionless), n × m.
    /// @param C_beta [out] β MO coefficients (dimensionless), n × m.
    /// @param eig_alpha [out] α orbital energies [Hartree], length m.
    /// @param eig_beta [out] β orbital energies [Hartree], length m.
    DiagonalizeStage(bool restricted, const RealMatrix& F_alpha,
                     const RealMatrix& F_beta, const RealMatrix& S,
                     const RealMatrix& X, const int& n_dropped,
                     RealMatrix& C_alpha, RealMatrix& C_beta,
                     RealVector& eig_alpha, RealVector& eig_beta)
        : restricted_(restricted), F_alpha_(F_alpha), F_beta_(F_beta), S_(S),
          X_(X), n_dropped_(n_dropped), C_alpha_(C_alpha), C_beta_(C_beta),
          eig_alpha_(eig_alpha), eig_beta_(eig_beta) {}

    /// @brief Component name.
    /// @return "scf_diagonalize".
    std::string name() const override { return "scf_diagonalize"; }

    /// @brief Default binding: `scf::Phase::diagonalize` @ 0.
    /// @return Single binding placing this at the diagonalize phase.
    std::vector<DefaultBinding> default_bindings() const override {
        return {ATV::make_binding(ATV::scf::Phase::diagonalize, 0)};
    }

    /// @brief One-time workspace pre-allocation for the closed eigensolver.
    ///        `launch()` allocates no fixed buffer of its own.
    /// @param ctx Build context (unused; state arrives by reference).
    void build(const cpu::BuildContext& ctx) override {
        (void)ctx;
        eigen_.prepare(static_cast<int>(S_.shape()[0]));
    }

    /// @brief Solve F C = S C ε into the wired C/eig references.
    ///
    /// Path selection by the wired `n_dropped` reference and restricted β reuse
    /// are bit-identical to `scf_solver.cpp` (`:145-147`, `:285-291`).
    /// @param ctx Launch context; self-guarded by `ctx.at(scf::Phase::diagonalize)`.
    void launch(const cpu::LaunchContext& ctx) override {
        if (!ctx.at(ATV::scf::Phase::diagonalize))
            return;
        auto diagonalize = [&](const RealMatrix& F) {
            return (n_dropped_ > 0) ? eigen_.step_orthogonal(F, X_)
                                    : eigen_.step(F, S_);
        };
        std::tie(eig_alpha_, C_alpha_) = diagonalize(F_alpha_);
        if (!restricted_) {
            std::tie(eig_beta_, C_beta_) = diagonalize(F_beta_);
        } else {
            C_beta_ = C_alpha_;
            eig_beta_ = eig_alpha_;
        }
    }

  private:
    bool restricted_;             ///< RHF/RKS β-reuse mode selector.
    const RealMatrix& F_alpha_;   ///< α Fock [Hartree] (read).
    const RealMatrix& F_beta_;    ///< β Fock [Hartree] (read).
    const RealMatrix& S_;         ///< Overlap (dimensionless, read).
    const RealMatrix& X_;         ///< Orthogonalization transform (read).
    const int& n_dropped_;        ///< Column-drop count → path selector (read).
    RealMatrix& C_alpha_;         ///< α MO coefficients (write).
    RealMatrix& C_beta_;          ///< β MO coefficients (write).
    RealVector& eig_alpha_;       ///< α orbital energies [Hartree] (write).
    RealVector& eig_beta_;        ///< β orbital energies [Hartree] (write).
    internal::EigenSolver eigen_; ///< Closed numeric kernel (NOT a Component).
};

} // namespace scf
} // namespace cpu
} // namespace ATV
