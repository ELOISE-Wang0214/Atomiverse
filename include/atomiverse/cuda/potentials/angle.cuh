/// @file cuda/potentials/angle.cuh
/// @brief LAMMPS-style abstract base for bonded 3-body styles (CUDA).
#pragma once

#include <cstdint>
#include <functional>
#include <vector>

#include "atomiverse/cuda/component.cuh"

namespace ATV {
namespace cuda {

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

    const std::vector<int>&
    select_owned_tuples(const LaunchContext& ctx, int n_tuples,
                        const std::function<std::int64_t(int)>& owner_tag_of);

  private:
    std::vector<int> owned_;
};

} // namespace cuda
} // namespace ATV
