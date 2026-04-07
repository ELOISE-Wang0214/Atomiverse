/// @file cpu/md/nose_hoover.cpp
/// @brief CPU Nose-Hoover thermostat implementation.

#include "atomiverse/cpu/md/nose_hoover.h"
#include "atomiverse/core/units.h"
#include "atomiverse/facade/stage.h"

#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace ATV {
namespace cpu {
namespace md {

NoseHoover::NoseHoover(double T_target, double tau, double dt)
    : T_target_(T_target), tau_(tau), dt_(dt) {}

void NoseHoover::build(const cpu::BuildContext& ctx) {
    int n = ctx.frame.atoms.size();
    int n_dof = 3 * n - 3;  // 3N - 3 (COM removal)
    Q_ = static_cast<double>(n_dof) * kBoltzmann * T_target_ * tau_ * tau_;
}

void NoseHoover::launch(const cpu::LaunchContext& ctx) {
    auto& f = ctx.frame;
    int n = f.atoms.size();
    int n_dof = 3 * n - 3;

    if (ctx.stage == static_cast<int>(ATV::md::Phase::initial_integrate) ||
        ctx.stage == static_cast<int>(ATV::md::Phase::final_integrate)) {

        // Compute current KE
        double ke = 0;
        #pragma omp parallel for schedule(static) reduction(+:ke)
        for (int i = 0; i < n; ++i) {
            double v2 = f.atoms.vx[i] * f.atoms.vx[i] + f.atoms.vy[i] * f.atoms.vy[i] + f.atoms.vz[i] * f.atoms.vz[i];
            ke += f.atoms.mass[i] * v2;
        }
        ke *= double(0.5) * kMvv2e;

        // Update thermostat variable
        xi_ += (dt_ * double(0.5)) * (double(2) * ke - static_cast<double>(n_dof) * kBoltzmann * T_target_) / Q_;

        // Rescale velocities
        double scale = std::exp(-dt_ * double(0.5) * xi_);
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < n; ++i) {
            f.atoms.vx[i] *= scale;
            f.atoms.vy[i] *= scale;
            f.atoms.vz[i] *= scale;
        }
    }
}

} // namespace md
} // namespace cpu
} // namespace ATV
