/// @file cuda/scf/krhf.cuh
/// @brief cuda::scf::KRHF — k-point restricted Hartree–Fock on device.
///
/// Same public call as `cpu::scf::KRHF`: `run(Frame&)` writes `kEnergy`
/// [eV]. Setup packs geometry, the s-basis, the FFT mesh and the
/// core Hamiltonian. The SCF step (Coulomb, exchange, Fock, diagonalize,
/// density) runs on device. l > 0 throws `krhf_ao_lgt0`. Forces throw
/// `krhf_nuclear_forces`.
#pragma once

#include <array>

#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/scf/basis_set.h"
#include "atomiverse/cpu/scf/gradient_result.h"
#include "atomiverse/cuda/component/hf.cuh"

namespace ATV {
namespace cuda {
namespace scf {

/// @brief Device KRHF (FFTDF J/K, Ewald exchange divergence).
class KRHF : public cuda::HF {
  public:
    /// @brief Basis set. Fixed at construction.
    const cpu::scf::BasisSet basis;
    /// @brief Monkhorst–Pack divisions. Each entry ≥ 1.
    const std::array<int, 3> k_mesh;
    /// @brief Cell charge [e].
    const int charge;
    /// @brief |ΔE| threshold [Hartree].
    const double tol;
    /// @brief Iteration cap.
    const int max_iter;

    /// @brief Construct device k-point RHF.
    /// @param basis_in Gaussian basis (s functions on this path).
    /// @param k_mesh_in Monkhorst–Pack mesh.
    /// @param charge_in Cell charge [e].
    /// @param tol_in Energy change [Hartree].
    /// @param max_iter_in Iteration cap.
    KRHF(cpu::scf::BasisSet basis_in, std::array<int, 3> k_mesh_in,
         int charge_in = 0, double tol_in = 1e-9, int max_iter_in = 40);

    /// @brief No analytic nuclear gradient on this path.
    bool provides_force() const override { return false; }

    /// @brief SCF. Overwrites `kEnergy` and `kEnergyElectronic` [eV].
    /// @param frame Orthorhombic cell, positions [Å].
    void run(cpu::Frame& frame);

    /// @brief SCF without writing the Frame.
    /// @param frame Orthorhombic cell (read only).
    void run(const cpu::Frame& frame);

    /// @brief Throws `krhf_nuclear_forces`.
    cpu::scf::GradientResult forces(cpu::Frame& frame);

    /// @brief Last total energy [Hartree].
    double energy() const { return energy_ha_; }
    /// @brief True after a converged `run`.
    bool converged() const { return converged_; }
    /// @brief Iterations of the last `run`.
    int iterations() const { return iterations_; }

  private:
    double energy_ha_ = 0.0;
    bool converged_ = false;
    int iterations_ = 0;
    void solve(const cpu::Frame& frame, bool write, cpu::Frame* out);
};

} // namespace scf
} // namespace cuda
} // namespace ATV
