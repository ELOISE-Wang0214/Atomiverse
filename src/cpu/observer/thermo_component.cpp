/// @file cpu/observer/thermo_component.cpp
/// @brief CPU Thermo observer component implementation.

#include "atomiverse/cpu/observer/thermo_component.h"
#include "atomiverse/core/log.h"
#include "atomiverse/core/units.h"

namespace ATV {
namespace cpu {

Thermo::Thermo(int freq) : freq_(freq) {}

void Thermo::launch(const cpu::LaunchContext& ctx) {
    // MPI-aware: only rank 0 prints
    if (ctx.runtime.mpi_rank != 0) return;

    auto& f = ctx.frame;
    Int n = f.atoms.size();
    Int n_dof = 3 * n - 3;

    Real temperature = 0;
    if (n_dof > 0)
        temperature = Real(2) * f.kinetic / (static_cast<Real>(n_dof) * kBoltzmann);

    Real total_energy = f.energy + f.kinetic;

    if (ctx.step == 0) {
        ATV_INFO("  Step          PE              KE           Temp         TotalE");
        ATV_INFO("  ----  --------------- --------------- --------------- ---------------");
    }

    ATV_INFO("  {:4d}  {:15.8e} {:15.8e} {:15.8e} {:15.8e}",
             ctx.step, static_cast<double>(f.energy), static_cast<double>(f.kinetic),
             static_cast<double>(temperature), static_cast<double>(total_energy));
}

} // namespace cpu
} // namespace ATV
