/// @file cpu/scf/krhf.cpp
/// @brief cpu::scf::KRHF — Frame I/O around internal::krhf_energy.

#include "atomiverse/cpu/scf/krhf.h"

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/core/units.h"
#include "atomiverse/cpu/scf/internal/kpoint_hf.h"

#include <stdexcept>
#include <string>

namespace ATV {
namespace cpu {
namespace scf {

KRHF::KRHF(BasisSet basis_in, std::array<int, 3> k_mesh_in, int charge_in,
           double tol_in, int max_iter_in)
    : basis(basis_in), k_mesh(k_mesh_in), charge(charge_in), tol(tol_in),
      max_iter(max_iter_in) {}

void KRHF::solve(const Frame& frame, bool write, Frame* out) {
    const int mesh[3] = {k_mesh[0], k_mesh[1], k_mesh[2]};
    const internal::KrhfEnergy e =
        internal::krhf_energy(frame, basis, mesh, charge, tol, max_iter);
    energy_ha_ = e.hartree;
    converged_ = e.converged;
    iterations_ = e.iterations;
    if (write) {
        if (out == nullptr)
            throw std::logic_error(
                "KRHF::solve: write requested without frame");
        const double e_ev = energy_ha_ * units::kHartreeInEv;
        out->meta.get(kEnergy) = e_ev;
        out->meta.get(kEnergyElectronic) = e_ev;
    }
}

void KRHF::run(Frame& frame) { solve(frame, true, &frame); }

void KRHF::run(const Frame& frame) { solve(frame, false, nullptr); }

GradientResult KRHF::forces(Frame& /*frame*/) {
    throw std::runtime_error(
        "KRHF::forces: analytic k-point nuclear forces are not implemented "
        "(missing kernel krhf_nuclear_forces; PySCF entry is "
        "pyscf.pbc.grad.krhf.Gradients, outside this product)");
}

} // namespace scf
} // namespace cpu
} // namespace ATV
