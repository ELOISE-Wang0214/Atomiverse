/// @file cpu/observer/thermo_component.h
/// @brief CPU Thermo observer component for MD simulations.
#pragma once

#include <string>

#include "atomiverse/cpu/component.h"

namespace ATV {
namespace cpu {

/// @brief Thermodynamic output component (CPU backend).
///
/// Prints step, PE, KE, temperature at the given frequency.
/// MPI-aware: only rank 0 prints.
class Thermo : public cpu::ComponentBase {
public:
    /// @param freq Output frequency (every freq steps). 0 = every step.
    explicit Thermo(int freq);

    std::string name() const override { return "thermo"; }
    int frequency() const override { return freq_; }
    void launch(const cpu::LaunchContext& ctx) override;

private:
    int freq_;
};

} // namespace cpu
} // namespace ATV
