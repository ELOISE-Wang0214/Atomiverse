/// @file cpu/potentials/dihedral.h
/// @brief LAMMPS-style abstract base for bonded 4-body styles (CPU).
///
/// Same shape as `cpu::Bond` / `cpu::Angle`; differs only in the topology
/// block walked (`frame.dihedrals`) and the LAMMPS owner-tag rule
/// (2nd-atom of the dihedral tuple owns).
#pragma once

#include <cstdint>
#include <functional>
#include <vector>

#include "atomiverse/cpu/component.h"
#include "atomiverse/facade/component_role.h"

namespace ATV {
namespace cpu {

class Dihedral : public ComponentBase {
public:
    ~Dihedral() override = default;

    ComponentRole role() const final { return ComponentRole::Evaluator; }
    Capability capabilities() const override {
        return Capability::Energy | Capability::Force;
    }

    void launch(const LaunchContext& ctx) final { compute(ctx); }

protected:
    Dihedral() = default;

    virtual void compute(const LaunchContext& ctx) = 0;

    /// @brief Multi-rank-filtered indices into `frame.dihedrals`.
    ///
    /// LAMMPS convention: dihedral (i,j,k,l) owned by the rank holding
    /// `tag[j]` (the 2nd atom). See `Bond::select_owned_tuples` for full
    /// semantics.
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
