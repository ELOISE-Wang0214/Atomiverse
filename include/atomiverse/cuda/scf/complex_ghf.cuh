/// @file cuda/scf/complex_ghf.cuh
/// @brief cuda::scf::ComplexGHF — device complex spin-orbital Hartree–Fock.
///
/// Same public call as `cpu::scf::ComplexGHF`. The device Fock kernel
/// is not implemented: `run` throws `launch_complex_ghf_fock` and does
/// not build ERIs, the Fock matrix, or the density on the host.
#pragma once

#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/scf/basis_set.h"
#include "atomiverse/cpu/scf/gradient_result.h"
#include "atomiverse/cuda/component/hf.cuh"

namespace ATV {
namespace cuda {
namespace scf {

/// @brief Device complex GHF. The step kernel is absent.
class ComplexGHF : public cuda::HF {
  public:
    /// @brief Basis set. Fixed at construction.
    const cpu::scf::BasisSet basis;
    /// @brief Molecular charge [e].
    const int charge;
    /// @brief |ΔE| threshold [Hartree].
    const double tol;
    /// @brief Iteration cap.
    const int max_iter;
    /// @brief X2C-1e core when true. The device path throws either way.
    bool x2c = true;

    /// @brief Construct device complex GHF.
    /// @param basis_in Gaussian basis.
    /// @param charge_in Molecular charge [e].
    /// @param tol_in Energy change [Hartree].
    /// @param max_iter_in Iteration cap.
    ComplexGHF(cpu::scf::BasisSet basis_in, int charge_in = 0,
               double tol_in = 1e-10, int max_iter_in = 80);

    /// @brief No analytic nuclear gradient on this path.
    bool provides_force() const override { return false; }

    /// @brief Throws `launch_complex_ghf_fock`.
    /// @param frame Positions [Å], species in `kAtom_type`.
    void run(cpu::Frame& frame);

    /// @brief Throws `launch_complex_ghf_fock`.
    /// @param frame Geometry (read only).
    void run(const cpu::Frame& frame);

    /// @brief Throws `complex_ghf_nuclear_forces`.
    cpu::scf::GradientResult forces(cpu::Frame& frame);

    /// @brief Last total energy [Hartree]. Unset: `run` does not converge.
    double energy() const { return energy_ha_; }
    /// @brief Always false: the device kernel is absent.
    bool converged() const { return converged_; }
    /// @brief Iterations of the last `run` (zero: `run` throws).
    int iterations() const { return iterations_; }

  private:
    double energy_ha_ = 0.0;
    bool converged_ = false;
    int iterations_ = 0;
};

} // namespace scf
} // namespace cuda
} // namespace ATV
