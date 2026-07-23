/// @file cpu/minimize/anisotropic_cell_relax.h
/// @brief Orthorhombic anisotropic cell relax (independent s_x,s_y,s_z).
#pragma once

#include <string>
#include <vector>

#include "atomiverse/core/stage.h"
#include "atomiverse/cpu/component/optimizer.h"

namespace ATV {
namespace cpu {
namespace minimizer {

/// @brief Diagonal-cell steepest-descent Optimizer (anisotropic skeleton).
///
/// Reads Voigt stress on `kStress`; independently scales axes toward
/// target pressure \(P_{\mathrm{target}}\) (same on each axis for v1).
class AnisotropicCellRelax : public Optimizer {
  public:
    double P_target = 0.0;   ///< Target pressure [eV/Å³].
    double alpha0 = 0.05;    ///< Scale gain (dimensionless).
    double max_scale = 0.05; ///< Cap on |s−1| per axis.
    double bulk_mod = 1.0;   ///< Soft bulk modulus [eV/Å³].

    AnisotropicCellRelax() = default;
    explicit AnisotropicCellRelax(double P_target_ev_a3, double alpha = 0.05);

    std::string name() const override { return "anisotropic_cell_relax"; }

    std::vector<DefaultBinding> default_bindings() const override {
        return {ATV::make_binding(ATV::minimizer::Phase::search_direction, 0),
                ATV::make_binding(ATV::minimizer::Phase::line_search, 0),
                ATV::make_binding(ATV::minimizer::Phase::convergence_check, 0)};
    }

    void build(const BuildContext& ctx) override;
    void launch(const LaunchContext& ctx) override;

    double pressure_residual() const { return p_residual_; }
    double pressure() const { return pressure_; }

  private:
    void apply_scales_(Frame& f, double sx, double sy, double sz) const;

    double pressure_ = 0;
    double p_residual_ = 0;
    double sx_ = 1, sy_ = 1, sz_ = 1;
    double energy_prev_ = 0;
    int step_ = 0;
    bool has_prev_ = false;
};

} // namespace minimizer
} // namespace cpu
} // namespace ATV
