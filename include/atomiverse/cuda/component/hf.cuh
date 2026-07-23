/// @file cuda/component/hf.cuh
/// @brief cuda::HF — pure polymorphic tag for Gaussian-basis mean-field methods.
///
/// Mirrors cpu::HF: a non-Component energy/force-source tag. Outer drivers
/// recognise an energy source via dynamic_cast<const HF*>(component_or_node).
#pragma once

namespace ATV {
namespace cuda {

/// @brief Gaussian-basis electronic structure (HF and hybrid KS) kind tag.
class HF {
  public:
    virtual ~HF() = default;
};

} // namespace cuda
} // namespace ATV
