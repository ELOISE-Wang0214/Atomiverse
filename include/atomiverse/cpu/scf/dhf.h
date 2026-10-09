/// @file cpu/scf/dhf.h
/// @brief Four-component Dirac–Hartree–Fock on j-adapted spinors.
///
/// Not a `MeanFieldMethod`: that driver stores real matrices. The SCF
/// is `pyscf.scf.dhf.DHF` (Dirac–Coulomb, complex 4-component). Analytic
/// nuclear gradients are not this product — `forces` throws
/// `dhf_nuclear_forces` (`pyscf.grad.dhf` is outside this product).
#pragma once

#include "atomiverse/cpu/component/hf.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/scf/basis_set.h"
#include "atomiverse/cpu/scf/gradient_result.h"

namespace ATV {
namespace cpu {
namespace scf {

/// @brief Dirac–Hartree–Fock in the kinetically balanced spinor basis.
///
/// Transplants `pyscf.scf.dhf.get_hcore`, `get_ovlp`, `get_occ`,
/// `energy_elec`, and the (LL|LL)+(SS|LL) Coulomb of `get_jk_coulomb`.
/// `run` writes `kEnergy` and `kEnergyElectronic` [eV]. `energy()` is
/// Hartree. Shells with l > 1 throw `dhf_ao_lgt1`. The (SS|SS) block
/// shifts H2 by about 1e-11 Ha and is not assembled. No Gaunt or Breit
/// term. No device path.
class DHF : public cpu::HF {
  public:
    /// @brief Basis set. Fixed at construction.
    const BasisSet basis;
    /// @brief Molecular charge [e]. Electrons = ΣZ − charge.
    const int charge;
    /// @brief |ΔE| convergence threshold [Hartree], two successive steps.
    const double tol;
    /// @brief Maximum SCF iterations.
    const int max_iter;

    /// @brief Construct Dirac–Hartree–Fock.
    /// @param basis_in Gaussian basis (l ≤ 1 on this path).
    /// @param charge_in Molecular charge [e].
    /// @param tol_in Energy change [Hartree].
    /// @param max_iter_in Iteration cap.
    DHF(BasisSet basis_in, int charge_in = 0, double tol_in = 1e-10,
        int max_iter_in = 80);

    /// @brief No analytic nuclear gradient (`dhf_nuclear_forces`).
    bool provides_force() const override { return false; }

    /// @brief Run SCF. Overwrites `kEnergy` and `kEnergyElectronic` [eV].
    /// @param frame Positions [Å], species in `kAtom_type`.
    void run(Frame& frame);

    /// @brief Run SCF without writing the Frame.
    /// @param frame Geometry (read only).
    void run(const Frame& frame);

    /// @brief Analytic forces are not implemented.
    /// @throws std::runtime_error naming `dhf_nuclear_forces`.
    GradientResult forces(Frame& frame);

    /// @brief Last total energy [Hartree].
    double energy() const { return energy_ha_; }
    /// @brief True after a converged `run`.
    bool converged() const { return converged_; }
    /// @brief SCF iterations of the last `run`.
    int iterations() const { return iterations_; }

  private:
    double energy_ha_ = 0.0;
    bool converged_ = false;
    int iterations_ = 0;
};

} // namespace scf
} // namespace cpu
} // namespace ATV
