/// @file cpu/scf/krhf.h
/// @brief k-point restricted Hartree–Fock on an orthorhombic cell.
///
/// Public shape matches the other Gaussian mean-field methods:
/// `KRHF{basis, k_mesh}.run(frame)` writes `kEnergy` [eV] and
/// `energy()` returns Hartree. The SCF is the PySCF `KRHF` /
/// `FFTDF` kernel (`pbc.scf.khf`, `pbc.df.fft_jk`); it is not a
/// member of `MeanFieldMethod` (that class stays molecular / Γ).
/// Analytic nuclear forces are not this product — `forces` throws
/// naming `krhf_nuclear_forces`.
#pragma once

#include <array>

#include "atomiverse/cpu/component/hf.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/scf/basis_set.h"
#include "atomiverse/cpu/scf/gradient_result.h"

namespace ATV {
namespace cpu {
namespace scf {

/// @brief Restricted Hartree–Fock with Monkhorst–Pack k-point sampling.
///
/// All-electron FFT Coulomb and exchange (PySCF `FFTDF`, exxdiv
/// `'ewald'`) for s-type Gaussians on an orthorhombic `Frame.box`.
/// A shell with l > 0 throws (`krhf_ao_lgt0`).
class KRHF : public cpu::HF {
  public:
    /// @brief Basis set. Fixed at construction.
    const BasisSet basis;
    /// @brief Monkhorst–Pack divisions along a, b, c. Each entry ≥ 1.
    ///        `{1,1,1}` is the Γ point alone.
    const std::array<int, 3> k_mesh;
    /// @brief Cell charge [e]. Electrons per cell = ΣZ − charge;
    ///        the k-mesh total follows PySCF `tot_electrons(nk)`.
    const int charge;
    /// @brief |ΔE| convergence threshold [Hartree].
    const double tol;
    /// @brief Maximum SCF iterations.
    const int max_iter;

    /// @brief Construct k-point RHF.
    /// @param basis_in Gaussian basis (s functions only on this path).
    /// @param k_mesh_in Monkhorst–Pack mesh, each component ≥ 1.
    /// @param charge_in Cell charge [e].
    /// @param tol_in Energy change threshold [Hartree].
    /// @param max_iter_in Iteration cap.
    KRHF(BasisSet basis_in, std::array<int, 3> k_mesh_in, int charge_in = 0,
         double tol_in = 1e-9, int max_iter_in = 40);

    /// @brief This product has no analytic nuclear gradient.
    bool provides_force() const override { return false; }

    /// @brief Run SCF. Overwrites `kEnergy` and `kEnergyElectronic` [eV].
    /// @param frame Orthorhombic cell, positions [Å], species in `kAtom_type`.
    void run(Frame& frame);

    /// @brief Run SCF without writing the Frame.
    /// @param frame Orthorhombic cell (read only).
    void run(const Frame& frame);

    /// @brief Analytic forces are not implemented.
    /// @throws std::runtime_error naming `krhf_nuclear_forces`.
    GradientResult forces(Frame& frame);

    /// @brief Last converged total energy [Hartree].
    double energy() const { return energy_ha_; }

    /// @brief True after a converged `run`.
    bool converged() const { return converged_; }

    /// @brief SCF iterations of the last `run`.
    int iterations() const { return iterations_; }

  private:
    double energy_ha_ = 0.0;
    bool converged_ = false;
    int iterations_ = 0;

    void solve(const Frame& frame, bool write, Frame* out);
};

} // namespace scf
} // namespace cpu
} // namespace ATV
