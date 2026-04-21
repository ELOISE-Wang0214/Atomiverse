/// @file cpu/md/virial_stress.cpp
/// @brief Kinetic-stress contribution.

#include "atomiverse/cpu/md/virial_stress.h"

#include "atomiverse/core/units.h"
#include "atomiverse/cpu/md/mpi/domain.h"
#include "atomiverse/cpu/md/owned_range.h"

namespace ATV {
namespace cpu {
namespace md {

void VirialStress::launch(const cpu::LaunchContext& ctx) {
    auto& f = ctx.frame;
    const int n = owned_count(ctx);

    double v_xx = 0, v_yy = 0, v_zz = 0, v_xy = 0, v_xz = 0, v_yz = 0;

    #pragma omp parallel for schedule(static) \
        reduction(+:v_xx) reduction(+:v_yy) reduction(+:v_zz) \
        reduction(+:v_xy) reduction(+:v_xz) reduction(+:v_yz)
    for (int i = 0; i < n; ++i) {
        const double m   = f.atoms.mass[i];
        const double vxi = f.atoms.vx[i];
        const double vyi = f.atoms.vy[i];
        const double vzi = f.atoms.vz[i];
        v_xx += m * vxi * vxi;
        v_yy += m * vyi * vyi;
        v_zz += m * vzi * vzi;
        v_xy += m * vxi * vyi;
        v_xz += m * vxi * vzi;
        v_yz += m * vyi * vzi;
    }

    // mass·v² is in amu·(Å/fs)², the kMvv2e factor brings it to eV.
    f.stress[0] += kMvv2e * v_xx;
    f.stress[1] += kMvv2e * v_yy;
    f.stress[2] += kMvv2e * v_zz;
    f.stress[3] += kMvv2e * v_xy;
    f.stress[4] += kMvv2e * v_xz;
    f.stress[5] += kMvv2e * v_yz;
}

} // namespace md
} // namespace cpu
} // namespace ATV
