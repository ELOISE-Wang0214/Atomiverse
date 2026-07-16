/// @file cpu/pw/internal/ewald.h
/// @brief 3D Ewald summation for periodic point charges.
#pragma once

#include "atomiverse/core/box.h"
#include "atomiverse/core/precision.h"
#include "atomiverse/core/vec3.h"
#include <vector>

namespace ATV {
namespace cpu {
namespace pw {
namespace internal {

class EwaldSummation {
    AccReal e_real_ = 0, e_recip_ = 0, e_self_ = 0;

  public:
    EwaldSummation(const std::vector<Vec3>& pos_ang,
                   const std::vector<AccReal>& charges, const ::ATV::Box& cell,
                   AccReal ecut_hartree);

    AccReal total_energy() const { return e_real_ + e_recip_ + e_self_; }
    AccReal real_energy() const { return e_real_; }
    AccReal recip_energy() const { return e_recip_; }
    AccReal self_energy() const { return e_self_; }
    AccReal alpha() const { return alpha_; }

  private:
    AccReal alpha_ = AccReal(1);
};

} // namespace internal
} // namespace pw
} // namespace cpu
} // namespace ATV
