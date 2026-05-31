/// @file cpu/potentials/angle.h
/// @brief LAMMPS-style abstract base for bonded 3-body styles (CPU).
#pragma once

#include "atomiverse/cpu/component.h"
#include "atomiverse/facade/component_role.h"

namespace ATV {
namespace cpu {

class Angle : public ComponentBase {
  public:
    ~Angle() override = default;

    ComponentRole role() const final { return ComponentRole::Evaluator; }
    Capability capabilities() const override {
        return Capability::Energy | Capability::Force;
    }

    void launch(const LaunchContext& ctx) final { compute(ctx); }

  protected:
    Angle() = default;

    virtual void compute(const LaunchContext& ctx) = 0;
};

} // namespace cpu
} // namespace ATV
