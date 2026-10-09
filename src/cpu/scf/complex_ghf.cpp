/// @file cpu/scf/complex_ghf.cpp
/// @brief Frame I/O for cpu::scf::ComplexGHF.

#include "atomiverse/cpu/scf/complex_ghf.h"

#include <stdexcept>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/core/units.h"
#include "atomiverse/cpu/scf/internal/complex_ghf.h"

namespace ATV {
namespace cpu {
namespace scf {

ComplexGHF::ComplexGHF(BasisSet basis_in, int charge_in, double tol_in,
                       int max_iter_in)
    : basis(basis_in), charge(charge_in), tol(tol_in), max_iter(max_iter_in) {}

void ComplexGHF::run(const Frame& frame) {
    const internal::ComplexGhfEnergy e =
        internal::complex_ghf_energy(frame, basis, charge, tol, max_iter, x2c);
    energy_ha_ = e.hartree;
    converged_ = e.converged;
    iterations_ = e.iterations;
}

void ComplexGHF::run(Frame& frame) {
    run(static_cast<const Frame&>(frame));
    const double e_ev = energy_ha_ * units::kHartreeInEv;
    frame.meta.get(kEnergy) = e_ev;
    frame.meta.get(kEnergyElectronic) = e_ev;
}

GradientResult ComplexGHF::forces(Frame& /*frame*/) {
    throw std::runtime_error(
        "ComplexGHF::forces: analytic complex-spinor nuclear forces are "
        "not implemented (missing kernel complex_ghf_nuclear_forces; "
        "PySCF GHF.nuc_grad_method and x2c.SCF.nuc_grad_method raise "
        "NotImplementedError)");
}

} // namespace scf
} // namespace cpu
} // namespace ATV
