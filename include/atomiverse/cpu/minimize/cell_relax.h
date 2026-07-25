/// @file cpu/minimize/cell_relax.h
/// @brief Isotropic variable-cell relaxation Optimizer (classical virial).
///
/// Uses PairLJCut-style Voigt stress on `kStress` (xx,yy,zz,…) to form the
/// scalar pressure P = (σ_xx+σ_yy+σ_zz)/(3V) [eV/Å³] and isotropically scales
/// the Box (+ affine atom coordinates) toward a target pressure.
///
/// Termination is **not** ForceConverged by default — use
/// `cpu::PressureConverged` (criteria.h). Atom positions only scale affinely
/// with the cell; there is no independent atomic L-BFGS in this component.
///
/// Breadth-first: isotropic classical only. Anisotropic / joint atom+cell /
/// PW stress later.
#pragma once

#include <string>
#include <vector>

#include "atomiverse/core/stage.h"
#include "atomiverse/cpu/component/optimizer.h"

namespace ATV {
namespace cpu {
namespace minimizer {

/// @brief Isotropic cell steepest-descent Optimizer.
class IsotropicCellRelax : public Optimizer {
  public:
    double P_target =
        0.0;              ///< Target pressure [eV/Å³] (0 = vacuum/zero stress).
    double alpha0 = 0.05; ///< Relative scale step gain (dimensionless).
    double max_scale = 0.05; ///< Cap on |s−1| per step (dimensionless).
    double bulk_mod = 1.0;   ///< Soft bulk modulus scale [eV/Å³] for step size.

    /// @brief Construct with zero target pressure and default step gains.
    IsotropicCellRelax() = default;

    /// @param P_target_ev_a3 Target pressure [eV/Å³].
    /// @param alpha Relative scale gain (dimensionless).
    explicit IsotropicCellRelax(double P_target_ev_a3, double alpha = 0.05);

    std::string name() const override { return "isotropic_cell_relax"; }

    std::vector<DefaultBinding> default_bindings() const override {
        return {ATV::make_binding(ATV::minimizer::Phase::search_direction, 0),
                ATV::make_binding(ATV::minimizer::Phase::line_search, 0),
                ATV::make_binding(ATV::minimizer::Phase::convergence_check, 0)};
    }

    void build(const BuildContext& ctx) override;
    void launch(const LaunchContext& ctx) override;

    /// @brief Latest pressure residual |P − P_target| [eV/Å³].
    double pressure_residual() const { return p_residual_; }

    /// @brief Latest pressure [eV/Å³].
    double pressure() const { return pressure_; }

  private:
    static double pressure_from_frame_(const Frame& f);
    void apply_scale_(Frame& f, double s) const;

    double pressure_ = 0;
    double p_residual_ = 0;
    double scale_step_ = 1.0; ///< Multiplicative scale to apply at line_search.
    double energy_prev_ = 0;
    int step_ = 0;
    bool has_prev_ = false;
};

} // namespace minimizer
} // namespace cpu
} // namespace ATV
