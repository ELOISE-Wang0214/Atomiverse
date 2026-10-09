/// @file cuda/scf/complex_ghf.cu
/// @brief cuda::scf::ComplexGHF — refuse a host Fock fill.
///
/// The complex spin-orbital Fock, density, and diagonalization are not
/// on device. Naming the missing kernel is the contract: no call into
/// cpu:: native integrals and no host ERI upload.

#ifdef ATV_HAS_HF

#include "atomiverse/cuda/scf/complex_ghf.cuh"

#include <stdexcept>

namespace ATV {
namespace cuda {
namespace scf {

ComplexGHF::ComplexGHF(cpu::scf::BasisSet basis_in, int charge_in,
                       double tol_in, int max_iter_in)
    : basis(basis_in), charge(charge_in), tol(tol_in), max_iter(max_iter_in) {}

void ComplexGHF::run(const cpu::Frame& /*frame*/) {
    throw std::runtime_error(
        "cuda::scf::ComplexGHF::run: device complex spinor Fock is not "
        "implemented (missing kernel launch_complex_ghf_fock). No host "
        "ERI, Fock, density, or DIIS fill.");
}

void ComplexGHF::run(cpu::Frame& frame) {
    run(static_cast<const cpu::Frame&>(frame));
}

cpu::scf::GradientResult ComplexGHF::forces(cpu::Frame& /*frame*/) {
    throw std::runtime_error(
        "cuda::scf::ComplexGHF::forces: analytic complex-spinor nuclear "
        "forces are not implemented (missing kernel "
        "complex_ghf_nuclear_forces)");
}

} // namespace scf
} // namespace cuda
} // namespace ATV

#endif
