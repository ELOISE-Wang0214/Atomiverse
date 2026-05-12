/// @file cpu/pw/internal/crystal.cpp
/// @brief Si diamond crystal builder implementation.

#include "atomiverse/cpu/pw/internal/crystal.h"
#include "atomiverse/core/vec3.h"

namespace ATV {
namespace cpu {
namespace pw {
namespace internal {

cpu::Frame build_si_diamond(Real a) {
    cpu::Frame f(2);
    f.type_id[0] = 14; f.type_id[1] = 14;
    f.mass[0] = Real(28.0855); f.mass[1] = Real(28.0855);
    f.x[0] = 0;        f.y[0] = 0;        f.z[0] = 0;
    f.x[1] = a/Real(4); f.y[1] = a/Real(4); f.z[1] = a/Real(4);
    // FCC cell vectors
    f.box(0,0) = 0;          f.box(0,1) = a/Real(2); f.box(0,2) = a/Real(2);
    f.box(1,0) = a/Real(2);  f.box(1,1) = 0;          f.box(1,2) = a/Real(2);
    f.box(2,0) = a/Real(2);  f.box(2,1) = a/Real(2);  f.box(2,2) = 0;
    return f;
}

} // namespace internal
} // namespace pw
} // namespace cpu
} // namespace ATV
