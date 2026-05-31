/// @file cpu/potentials/bond.h
/// @brief LAMMPS-style abstract base for bonded 2-body styles (CPU).
///
/// Concrete bond styles (BondHarmonic, future BondMorse, BondFENE, ...)
/// derive from this class and implement compute().
#pragma once

#include "atomiverse/cpu/component.h"
#include "atomiverse/facade/component_role.h"

namespace ATV {
namespace cpu {

/// @brief Abstract base for bonded 2-body styles.
class Bond : public ComponentBase {
  public:
    ~Bond() override = default;

    ComponentRole role() const final { return ComponentRole::Evaluator; }
    Capability capabilities() const override {
        return Capability::Energy | Capability::Force | Capability::Stress;
    }

    void launch(const LaunchContext& ctx) final { compute(ctx); }

  protected:
    Bond() = default;

    /// @brief Hot-path override. No virtual call inside the per-bond loop.
    virtual void compute(const LaunchContext& ctx) = 0;
};

} // namespace cpu
} // namespace ATV
