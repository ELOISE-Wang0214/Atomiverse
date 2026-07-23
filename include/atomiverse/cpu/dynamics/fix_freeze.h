/// @file cpu/dynamics/fix_freeze.h
/// @brief Freeze / immobilize atoms by type — LAMMPS-class `fix setforce 0`.
///
/// E1 electrode walls and any immobile subset: zero forces after physical
/// force evaluation, restore build-time positions and zero velocities after
/// each integrate half-step. Reports \(3 N_{\mathrm{frozen}}\) via
/// `ConstraintProvider` so thermostats use the correct \(N_{\mathrm{dof}}\).
#pragma once

#include <string>
#include <vector>

#include "atomiverse/core/binding.h"
#include "atomiverse/core/precision.h"
#include "atomiverse/core/stage.h"
#include "atomiverse/cpu/component/constraint.h"
#include "atomiverse/cpu/component/constraint_protocol.h"

namespace ATV {
namespace cpu {
namespace md {

/// @brief Immobilize all atoms of a given type (electrode / wall subset).
///
/// Selection is by atom `type` (one type id). Mirrors LAMMPS
/// `group wall type N` + `fix setforce 0 0 0` + `velocity wall set 0 0 0`.
///
/// Default bindings:
/// - `post_force @ +900` — zero forces after physical (+ Langevin) forces
/// - `initial_integrate @ +200` — after VV (+ NH): restore \(x\), zero \(v\)
/// - `final_integrate @ +200` — after VV (+ NH): zero \(v\) (and restore \(x\))
class FixFreeze : public cpu::Constraint, public cpu::ConstraintProvider {
  public:
    /// @param atom_type Atom type id to freeze (e.g. wall type 2).
    explicit FixFreeze(int atom_type);

    std::string name() const override { return "fix_freeze"; }

    const cpu::ConstraintProvider* as_constraint() const override {
        return this;
    }

    /// @brief Global scalar constraints: \(3 N_{\mathrm{frozen}}\) (fix \(x,y,z\)).
    int n_constraints(const cpu::Context& ctx,
                      const cpu::Frame& frame) const override;

    std::vector<DefaultBinding> default_bindings() const override {
        return {ATV::make_binding(ATV::md::Phase::post_force, +900),
                ATV::make_binding(ATV::md::Phase::initial_integrate, +200),
                ATV::make_binding(ATV::md::Phase::final_integrate, +200)};
    }

    void build(const cpu::BuildContext& ctx) override;
    void launch(const cpu::LaunchContext& ctx) override;

    /// @return Frozen atom type id.
    int atom_type() const noexcept { return type_; }

  private:
    int type_;
    /// Build-time positions of every owned atom of `type_` (local rows).
    std::vector<int> local_idx_;
    std::vector<StorageReal> x0_, y0_, z0_;
    /// Global frozen count (for DOF; single-rank = local size).
    int n_frozen_global_ = 0;
};

} // namespace md
} // namespace cpu
} // namespace ATV
