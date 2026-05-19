/// @file cpu/potentials/angle.h
/// @brief LAMMPS-style abstract base for bonded 3-body styles (CPU).
///
/// Same shape as `cpu::Bond`; differs only in the topology block walked
/// (`frame.angles`) and the LAMMPS owner-tag rule (centre atom owns).
#pragma once

#include <cstdint>
#include <functional>
#include <vector>

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

    /// @brief Multi-rank-filtered indices into `frame.angles`.
    ///
    /// LAMMPS convention: angle (i,j,k) owned by the rank holding `tag[j]`
    /// (the centre atom). Pass that rule via `owner_tag_of`. See
    /// `Bond::select_owned_tuples` for full semantics.
    const std::vector<int>& select_owned_tuples(
        const LaunchContext& ctx,
        int n_tuples,
        const std::function<std::int64_t(int)>& owner_tag_of);

private:
    std::vector<int> owned_;
    std::uint64_t cached_generation_ = ~std::uint64_t{0};
};

} // namespace cpu
} // namespace ATV
