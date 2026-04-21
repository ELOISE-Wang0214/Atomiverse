/// @file cpu/md/deform.cpp
/// @brief CPU uniaxial affine deformation.

#include "atomiverse/cpu/md/deform.h"

#include "atomiverse/core/mat3.h"
#include "atomiverse/cpu/md/mpi/domain.h"
#include "atomiverse/cpu/md/owned_range.h"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace ATV {
namespace cpu {
namespace md {

Deform::Deform(double strain_rate, double dt, int axis)
    : strain_rate_(strain_rate), dt_(dt), axis_(axis) {}

void Deform::launch(const cpu::LaunchContext& ctx) {
    auto& f = ctx.frame;
    const double scale = 1.0 + strain_rate_ * dt_;
    f.box(axis_, axis_) *= scale;

    const int n = owned_count(ctx);
    auto& coord = (axis_ == 0) ? f.atoms.x
                : (axis_ == 1) ? f.atoms.y
                                : f.atoms.z;
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i) {
        coord[i] *= scale;
    }
}

} // namespace md
} // namespace cpu
} // namespace ATV
