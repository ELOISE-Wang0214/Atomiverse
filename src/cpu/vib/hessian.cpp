/// @file src/cpu/vib/hessian.cpp
/// @brief cpu::vib::Hessian — constructor wiring and run() implementation.

#include "atomiverse/cpu/vib/hessian.h"

namespace ATV {
namespace cpu {
namespace vib {

Hessian::Hessian() {
    // Auto-register the four structural components, each bound to the
    // driver-owned shared workspace. Each carries a default binding
    // (setup / displace / accumulate / diagonalize), so plain add()
    // resolves the phase.
    add(ReferenceGeometry{workspace_.get()});
    add(DisplaceStep{workspace_.get()});
    add(HessianAccumulator{workspace_.get()});
    add(Symmetrize{workspace_.get()});
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

    // Publish the assembled, symmetrized Hessian.
    last_result_.hessian = workspace_->hessian;
    last_result_.n_atoms = n_atoms;
    last_result_.units   = "eV/Angstrom^2";
    return last_result_;
}

} // namespace vib
} // namespace cpu
} // namespace ATV
