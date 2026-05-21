/// @file cpu/pw/internal/ewald.h
/// @brief 3D Ewald summation for periodic point charges.
#pragma once

#include "atomiverse/core/box.h"
#include "atomiverse/core/real.h"
#include "atomiverse/core/vec3.h"
#include <vector>

namespace ATV {
namespace cpu {
namespace pw {
namespace internal {

class EwaldSummation {
    Real e_real_ = 0, e_recip_ = 0, e_self_ = 0;

  public:
    EwaldSummation(const std::vector<Vec3>& pos_ang,
                   const std::vector<Real>& charges, const ::ATV::Box& cell,
                   Real ecut_hartree);

    Real total_energy() const { return e_real_ + e_recip_ + e_self_; }
    Real real_energy() const { return e_real_; }
    Real recip_energy() const { return e_recip_; }
    Real self_energy() const { return e_self_; }
    Real alpha() const { return alpha_; }

  private:
    Real alpha_ = Real(1);
};

} // namespace internal
} // namespace pw
} // namespace cpu
} // namespace ATV
