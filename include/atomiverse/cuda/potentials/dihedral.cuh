/// @file cuda/potentials/dihedral.cuh
/// @brief LAMMPS-style abstract base for bonded 4-body styles (CUDA).
#pragma once

#include <cstdint>
#include <functional>
#include <vector>

#include "atomiverse/cuda/component.cuh"

namespace ATV {
namespace cuda {

class Dihedral : public ComponentBase {
  public:
    ~Dihedral() override = default;

    ComponentRole role() const final { return ComponentRole::Evaluator; }
    /// @brief Capabilities — Energy + Force + Stress. The kernel writes
    ///        per-tuple virial contributions into `frame.scalars->stress[]`
    ///        (see `cuda/potentials/dihedral_multi_harmonic.cu`).
    Capability capabilities() const override {
        return Capability::Energy | Capability::Force | Capability::Stress;
    }

    void launch(const LaunchContext& ctx) final { compute(ctx); }

  protected:
    Dihedral() = default;

    virtual void compute(const LaunchContext& ctx) = 0;

    const std::vector<int>&
    select_owned_tuples(const LaunchContext& ctx, int n_tuples,
                        const std::function<std::int64_t(int)>& owner_tag_of);

  private:
    std::vector<int> owned_;
};

} // namespace cuda
} // namespace ATV
