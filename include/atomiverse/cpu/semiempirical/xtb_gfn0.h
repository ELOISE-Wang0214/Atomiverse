/// @file cpu/semiempirical/xtb_gfn0.h
/// @brief GFN0-xTB-style energy skeleton (H/C/O toy tables, CPU).
#pragma once

#include <string>

#include "atomiverse/cpu/component/semiempirical.h"

namespace ATV {
namespace cpu {
namespace semiempirical {

/// @brief Minimal GFN0-like tight-binding energy (second TB family vs DFTB2).
///
/// Toy element tables for H,C,O only — not production GFN parameters.
class XtbGfn0 : public cpu::Semiempirical {
  public:
    XtbGfn0() = default;
    std::string name() const override { return "xtb_gfn0"; }
    void build(const BuildContext& ctx) override;
    void launch(const LaunchContext& ctx) override;
    double last_energy() const { return energy_; }

  private:
    double energy_ = 0;
};

} // namespace semiempirical
} // namespace cpu
} // namespace ATV
