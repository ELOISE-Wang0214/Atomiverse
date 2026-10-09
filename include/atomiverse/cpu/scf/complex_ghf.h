/// @file cpu/scf/complex_ghf.h
/// @brief Complex two-component Hartree–Fock (spin-orbital X2C-1e or
///        non-relativistic GHF).
///
/// Not a `MeanFieldMethod`: that driver stores real matrices. The SCF
/// here is the complex spinor product `pyscf.scf.GHF.x2c1e`
/// (`xuncontract=False`) or, with `x2c = false`, `pyscf.scf.ghf.GHF`.
/// Analytic nuclear gradients are not this product — `forces` throws
/// `complex_ghf_nuclear_forces` (`GHF.nuc_grad_method` /
/// `x2c.SCF.nuc_grad_method` raise NotImplementedError).
#pragma once

#include "atomiverse/cpu/component/hf.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/scf/basis_set.h"
#include "atomiverse/cpu/scf/gradient_result.h"

namespace ATV {
namespace cpu {
namespace scf {

/// @brief Complex generalized Hartree–Fock in the spin-orbital basis.
///
/// `x2c == true` (default) builds the one-electron Hamiltonian with
/// `SpinOrbitalX2CHelper.get_hcore` / `_x2c1e_get_hcore` (scalar
/// relativity and spin–orbit). `x2c == false` is non-relativistic
/// complex GHF (`get_hcore` = block-diagonal T+V). Two-electron
/// integrals stay non-relativistic (`ghf.get_jk` on the dense ERI).
/// `run` writes `kEnergy` and `kEnergyElectronic` [eV]. `energy()` is
/// Hartree. l > 3 on the X2C path throws `complex_ghf_ao_lgt3`.
class ComplexGHF : public cpu::HF {
  public:
    /// @brief Basis set. Fixed at construction.
    const BasisSet basis;
    /// @brief Molecular charge [e]. Electrons = ΣZ − charge.
    const int charge;
    /// @brief |ΔE| convergence threshold [Hartree], two successive steps.
    const double tol;
    /// @brief Maximum SCF iterations.
    const int max_iter;
    /// @brief Exact-two-component one-electron core when true.
    ///        False keeps the non-relativistic spin-orbital core.
    bool x2c = true;

    /// @brief Construct complex GHF.
    /// @param basis_in Gaussian basis.
    /// @param charge_in Molecular charge [e].
    /// @param tol_in Energy change [Hartree].
    /// @param max_iter_in Iteration cap.
    ComplexGHF(BasisSet basis_in, int charge_in = 0, double tol_in = 1e-10,
               int max_iter_in = 80);

    /// @brief No analytic nuclear gradient (`complex_ghf_nuclear_forces`).
    bool provides_force() const override { return false; }

    /// @brief Run SCF. Overwrites `kEnergy` and `kEnergyElectronic` [eV].
    /// @param frame Positions [Å], species in `kAtom_type`.
    void run(Frame& frame);

    /// @brief Run SCF without writing the Frame.
    /// @param frame Geometry (read only).
    void run(const Frame& frame);

    /// @brief Analytic forces are not implemented.
    /// @throws std::runtime_error naming `complex_ghf_nuclear_forces`.
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
