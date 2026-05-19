/// @file cpu/potentials/bond.h
/// @brief LAMMPS-style abstract base for bonded 2-body styles (CPU).
///
/// Concrete bond styles (BondHarmonic, future BondMorse, BondFENE, ...)
/// derive from this class and implement compute(). The base absorbs the
/// multi-rank ownership cache (currently duplicated as `local_bonds_`
/// inside each bonded `launch()`).
#pragma once

#include <cstdint>
#include <functional>
#include <vector>

#include "atomiverse/cpu/component.h"
#include "atomiverse/facade/component_role.h"

namespace ATV {
namespace cpu {

/// @brief Abstract base for bonded 2-body styles.
///
/// Derived classes iterate `frame.bonds` (or the multi-rank-filtered
/// subset returned by `select_owned_tuples`) and compute energy + force
/// per bond. The base provides the ownership cache so each derived
/// `compute()` does not re-implement the LAMMPS owner-tag filter.
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

    /// @brief Multi-rank-filtered indices into `frame.bonds`.
    ///
    /// Returns the indices of bonds **owned by this rank** under the
    /// LAMMPS convention. The owner-tag rule is passed as a lambda so the
    /// same helper serves Bond (min of 2 tags), Angle (centre tag), and
    /// Dihedral (2nd-atom tag) via the parallel helpers on those bases.
    ///
    /// Single-rank: returns `0..n_bonds-1` (no filtering).
    /// Multi-rank: caches the index list against `dd->neighbor_generation`
    /// so the filter runs once per migration, not once per step.
    ///
    /// @param ctx        Launch context (provides frame + DD pointer).
    /// @param n_tuples   `frame.bonds.atomi.size()` — passed in to avoid
    ///                   depending on a specific xtensor field accessor here.
    /// @param owner_tag_of  Callback `idx → global tag of the owner atom`.
    /// @return Const reference to the cached index vector. Valid until the
    ///         next migration generation bump.
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
