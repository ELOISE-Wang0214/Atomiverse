/// @file cpu/scf/kuhf.h
/// @brief k-point unrestricted Hartree–Fock on an orthorhombic cell.
///
/// Public shape matches `KRHF`, with α/β counts:
/// `KUHF{basis, k_mesh, n_alpha, n_beta}.run(frame)` writes `kEnergy`
/// [eV] and `energy()` returns Hartree. The SCF is the PySCF `KUHF` /
/// `FFTDF` kernel (`pbc.scf.kuhf`, `pbc.df.fft_jk`). It is not a member
/// of `MeanFieldMethod` (that class stays molecular / Γ). Analytic
/// nuclear forces are not this product — `forces` throws naming
/// `kuhf_nuclear_forces`.
#pragma once

#include <array>

#include "atomiverse/cpu/component/hf.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/scf/basis_set.h"
#include "atomiverse/cpu/scf/gradient_result.h"

namespace ATV {
namespace cpu {
namespace scf {

/// @brief Unrestricted Hartree–Fock with Monkhorst–Pack k-point sampling.
///
/// All-electron FFT Coulomb and exchange (PySCF `FFTDF`, exxdiv
/// `'ewald'`) for s-type Gaussians on an orthorhombic `Frame.box`.
/// A shell with l > 0 throws (`krhf_ao_lgt0`). `n_alpha` and `n_beta`
/// are the electron counts across the whole mesh (`KUHF.nelec`), not a
/// per-k copy: PySCF sets `cell.spin = n_alpha - n_beta` and does not
/// scale that spin by the number of k-points.
class KUHF : public cpu::HF {
  public:
    /// @brief Basis set. Fixed at construction.
    const BasisSet basis;
    /// @brief Monkhorst–Pack divisions along a, b, c. Each entry ≥ 1.
    ///        `{1,1,1}` is the Γ point alone.
    const std::array<int, 3> k_mesh;
    /// @brief α electrons across the k-mesh (`KUHF.nelec` first entry).
    const int n_alpha;
    /// @brief β electrons across the k-mesh (`KUHF.nelec` second entry).
    const int n_beta;
    /// @brief Cell charge [e]. Electrons per mesh =
    ///        `(ΣZ)*nk − charge`, which `n_alpha + n_beta` must match.
    const int charge;
    /// @brief |ΔE| convergence threshold [Hartree].
    const double tol;
    /// @brief Maximum SCF iterations.
    const int max_iter;

    /// @brief Construct k-point UHF.
    /// @param basis_in Gaussian basis (s functions only on this path).
    /// @param k_mesh_in Monkhorst–Pack mesh, each component ≥ 1.
    /// @param n_alpha_in α electrons across the mesh.
    /// @param n_beta_in β electrons across the mesh.
    /// @param charge_in Cell charge [e].
    /// @param tol_in Energy change threshold [Hartree].
    /// @param max_iter_in Iteration cap.
    KUHF(BasisSet basis_in, std::array<int, 3> k_mesh_in, int n_alpha_in,
         int n_beta_in, int charge_in = 0, double tol_in = 1e-9,
         int max_iter_in = 40);

    /// @brief This product has no analytic nuclear gradient.
    bool provides_force() const override { return false; }

    /// @brief Run SCF. Overwrites `kEnergy` and `kEnergyElectronic` [eV].
    /// @param frame Orthorhombic cell, positions [Å], species in `kAtom_type`.
    void run(Frame& frame);

    /// @brief Run SCF without writing the Frame.
    /// @param frame Orthorhombic cell (read only).
    void run(const Frame& frame);

    /// @brief Analytic forces are not implemented.
    /// @throws std::runtime_error naming `kuhf_nuclear_forces`.
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
