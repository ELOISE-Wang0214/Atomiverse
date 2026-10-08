/// @file cpu/scf/rhf.h
/// @brief Restricted Hartree–Fock (RHF) concrete mean-field method.
#pragma once

#include <string>
#include <vector>

#include "atomiverse/cpu/scf/mean_field_method.h"

namespace ATV {
namespace cpu {
namespace scf {

/// @brief Restricted Hartree–Fock — a mean-field electronic-structure method
///        (a single Slater determinant in which every electron moves in the
///        averaged field of all the others) whose spin-up (α) and spin-down (β)
///        electrons share one common set of spatial orbitals: the closed-shell,
///        fully paired case.
class RHF : public MeanFieldMethod<RHF> {
  public:
    /// @brief Construct with the RHF defaults (STO-3G; restricted HF; MINAO;
    ///        no density damping; energy tolerance 1e-9 Ha).
    ///
    /// Prefer the basis ctor: `RHF{BasisSet::STO3G}` — **do not** write
    /// `rhf.basis = …` after default construction.
    RHF();

    /// @brief Construct with an explicit basis (primary public form).
    /// @param basis_in Basis set (e.g. `BasisSet::STO3G`).
    explicit RHF(BasisSet basis_in);

    /// @brief Construct with basis and molecular charge.
    /// @param basis_in Basis set.
    /// @param charge_in Molecular charge [e] (electrons = ΣZ − charge).
    RHF(BasisSet basis_in, int charge_in);

    /// @brief Construct with basis, charge, SCF tolerance, and iteration cap.
    /// @param basis_in Basis set.
    /// @param charge_in Molecular charge [e].
    /// @param tol_in Energy / residual tolerance [Ha] (see base `tol`).
    /// @param max_iter_in Maximum SCF iterations [dimensionless].
    RHF(BasisSet basis_in, int charge_in, double tol_in, int max_iter_in = 100);

    std::string name() const { return "rhf"; }

    /// @brief Γ-point periodic RHF. Default false keeps the molecular
    ///        engine bit-identical. When true, Frame.box must be
    ///        orthorhombic; J, K, S, T, V and E_nn share one truncated
    ///        lattice-image window (no Ewald). A cell larger than the
    ///        Gaussian pair reach has only T = 0 and recovers molecular RHF.
    bool periodic = false;

  protected:
    /// @brief Closed-shell occupation: n_alpha = n_beta = n_electrons / 2
    ///        (every occupied spatial orbital holds one α and one β electron).
    void resolve_occupation(int n_electrons) override;
};

} // namespace scf
} // namespace cpu
} // namespace ATV
