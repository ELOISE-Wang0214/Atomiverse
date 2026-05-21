/// @file src/cpu/vib/hessian.cpp
/// @brief cpu::vib::Hessian — constructor wiring and run() implementation.

#include "atomiverse/cpu/vib/hessian.h"

#include "atomiverse/cpu/vib/eigen_solver.h"
#include "atomiverse/cpu/vib/mass_weighting.h"

namespace ATV {
namespace cpu {
namespace vib {

Hessian::Hessian() {
    // Auto-register the structural components, each bound to the
    // driver-owned shared workspace. Each carries a default binding
    // (setup / displace / accumulate / diagonalize), so plain add()
    // resolves the phase. On vib::Phase::diagonalize the ordering is
    // Symmetrize (-1) -> MassWeighting (0) -> EigenSolver (+1).
    add(ReferenceGeometry{workspace_.get()});
    add(DisplaceStep{workspace_.get()});
    add(HessianAccumulator{workspace_.get()});
    add(Symmetrize{workspace_.get()});
    add(MassWeighting{workspace_.get()});
    add(EigenSolver{workspace_.get()});
}

Hessian::Hessian(WithoutAccumulator /*tag*/) {
    // Test-only configuration: omit HessianAccumulator so the build-time
    // validation in ResponseDriver throws for the missing
    // Capability::Hessian component.
    add(ReferenceGeometry{workspace_.get()});
    add(DisplaceStep{workspace_.get()});
    add(Symmetrize{workspace_.get()});
}

HessianResult Hessian::run(cpu::Frame& frame) {
    // Size the shared workspace from the frame before the component loop
    // runs — no allocation happens inside any component launch().
    const int n_atoms = frame.atoms.size();
    workspace_->allocate(n_atoms);

    // Drive the 6N-displacement Response-category loop.
    ResponseDriver::run(frame);

    // Publish the assembled, symmetrized Hessian plus the vibrational
    // frequencies produced by MassWeighting + EigenSolver.
    last_result_.hessian      = workspace_->hessian;
    last_result_.n_atoms      = n_atoms;
    last_result_.units        = "eV/Angstrom^2";
    last_result_.frequencies  = workspace_->frequencies;
    last_result_.normal_modes = workspace_->normal_modes;
    last_result_.n_imaginary  = workspace_->n_imaginary;
    last_result_.is_linear    = workspace_->is_linear;
    return last_result_;
}

} // namespace vib
} // namespace cpu
} // namespace ATV
