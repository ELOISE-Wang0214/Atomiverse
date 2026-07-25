/// @file cpu/dynamics/fix_settle.h
/// @brief Analytic SETTLE constraint for rigid 3-site water (Miyamoto–Kollman).
///
/// Specialized alternative to general RATTLE for TIP3P-like triangles.
/// Ref: Miyamoto & Kollman, J. Comput. Chem. 13, 952 (1992).
#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "atomiverse/core/precision.h"
#include "atomiverse/core/stage.h"
#include "atomiverse/cpu/component/constraint.h"
#include "atomiverse/cpu/component/constraint_protocol.h"

namespace ATV {
namespace cpu {
namespace md {

/// @brief SETTLE: analytic rigid water constraint (one molecule or list).
///
/// Constructor takes OH and HH target distances [Å]. Add waters by O/H tags
/// via `add_water(tag_O, tag_H1, tag_H2)`.
class FixSettle : public Constraint, public ConstraintProvider {
  public:
    /// @param dt Timestep [fs] (for velocity stage).
    /// @param d_oh O–H distance [Å] (TIP3P default 0.9572).
    /// @param d_hh H–H distance [Å] (TIP3P default 1.5139).
    FixSettle(double dt, double d_oh = 0.9572, double d_hh = 1.5139);

    /// @brief Register one water by stable tags.
    void add_water(TagInt tag_O, TagInt tag_H1, TagInt tag_H2);

    std::string name() const override { return "settle"; }

    std::vector<DefaultBinding> default_bindings() const override {
        return {ATV::make_binding(ATV::md::Phase::initial_integrate, -500),
                ATV::make_binding(ATV::md::Phase::post_integrate, 0),
                ATV::make_binding(ATV::md::Phase::final_integrate, +200)};
    }

    const ConstraintProvider* as_constraint() const override { return this; }

    int n_constraints(const cpu::Context& /*ctx*/,
                      const cpu::Frame& /*frame*/) const override {
        return 3 * static_cast<int>(waters_.size());
    }

    void build(const BuildContext& ctx) override;
    void launch(const LaunchContext& ctx) override;

  private:
    struct Water {
        TagInt o, h1, h2;
        int io = -1, ih1 = -1, ih2 = -1;
    };

    void resolve_rows_(const Frame& f);
    void settle_positions_(Frame& f) const;
    void settle_velocities_(Frame& f) const;
    void snapshot_(const Frame& f);

    double dt_;
    double d_oh_, d_hh_;
    std::vector<Water> waters_;
    // Reference positions for velocity stage (pre-move).
    std::vector<double> x0_, y0_, z0_;
    int n_atoms_ = 0;
};

} // namespace md
} // namespace cpu
} // namespace ATV
