/// @file cpu/scf/dhf.cpp
/// @brief Frame I/O for cpu::scf::DHF.

#include "atomiverse/cpu/scf/dhf.h"

#include <stdexcept>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/core/units.h"
#include "atomiverse/cpu/scf/internal/dhf.h"

namespace ATV {
namespace cpu {
namespace scf {

DHF::DHF(BasisSet basis_in, int charge_in, double tol_in, int max_iter_in)
    : basis(basis_in), charge(charge_in), tol(tol_in), max_iter(max_iter_in) {}

void DHF::run(const Frame& frame) {
    const internal::DhfEnergy e =
        internal::dhf_energy(frame, basis, charge, tol, max_iter);
    energy_ha_ = e.hartree;
    converged_ = e.converged;
    iterations_ = e.iterations;
}

void DHF::run(Frame& frame) {
    run(static_cast<const Frame&>(frame));
    const double e_ev = energy_ha_ * units::kHartreeInEv;
    frame.meta.get(kEnergy) = e_ev;
    frame.meta.get(kEnergyElectronic) = e_ev;
}

GradientResult DHF::forces(Frame& /*frame*/) {
    throw std::runtime_error(
        "DHF::forces: analytic Dirac–Hartree–Fock nuclear forces are not "
        "implemented (missing kernel dhf_nuclear_forces; PySCF entry is "
        "pyscf.grad.dhf.Gradients, outside this product)");
}

} // namespace scf
} // namespace cpu
} // namespace ATV
