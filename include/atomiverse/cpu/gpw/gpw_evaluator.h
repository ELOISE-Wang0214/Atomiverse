/// @file cpu/gpw/gpw_evaluator.h
/// @brief GPWEvaluator — GPW DFT energy component for CPU drivers.
///
/// Mirrors `cpu::pw::PWEvaluator` / `cpu::scf::SCFEvaluator`: a self-contained
/// Evaluator that wraps the inner `gpw_scf` solver. `build()` records the
/// per-element basis/pseudopotential templates; each `launch()` assembles the
/// atoms at the frame's current geometry (Angstrom -> Bohr), runs the GPW SCF,
/// and writes the total energy to `frame.scalars.energy` [eV].
///
/// `role() == Evaluator`, `capabilities() == Energy`. Plugs into the Evaluate /
/// MD / Minimize drivers exactly like the other ES evaluators (forces are a
/// later capability).
#pragma once

#include <string>
#include <vector>

#include "atomiverse/cpu/component/component.h"
#include "atomiverse/cpu/gpw/gpw_solver.h"
#include "atomiverse/facade/component_role.h"

namespace ATV {
namespace cpu {
namespace gpw {

/// @brief Per-element basis + pseudopotential, keyed by atomic number.
struct GpwElement {
    int z = 0;            ///< atomic number, matched against frame.atoms.type.
    GaussianBasis basis;  ///< contracted Gaussian basis for this element.
    GthPotential pot;     ///< GTH pseudopotential for this element.
};

/// @brief Configuration for the GPW evaluator (grid/XC/electrostatics + elements).
struct GpwEvaluatorConfig {
    double cutoff_ry = 300.0;
    int max_iter = 200;
    double tol = 1e-7;
    double mix = 0.3;
    XcKind xc = XcKind::PADE;
    bool nonperiodic = false;
    bool compute_forces = false;      ///< also compute finite-difference forces.
    double fd_step = 1e-3;            ///< FD displacement [Bohr] for forces.
    std::vector<GpwElement> elements; ///< one per distinct element in the system.
};

/// @brief GPW total-energy evaluator component.
class GPWEvaluator : public cpu::ComponentBase {
  public:
    explicit GPWEvaluator(GpwEvaluatorConfig config) : config_(std::move(config)) {}

    std::string name() const override { return "gpw_evaluator"; }
    ComponentRole role() const override { return ComponentRole::Evaluator; }
    Capability capabilities() const override {
        return config_.compute_forces ? (Capability::Energy | Capability::Force)
                                      : Capability::Energy;
    }

    void build(const cpu::BuildContext& ctx) override;
    void launch(const cpu::LaunchContext& ctx) override;

    /// @brief Last total energy [Hartree] (frame carries the eV value).
    double last_energy_hartree() const { return last_.e_total; }
    const GpwResult& last_result() const { return last_; }

  private:
    GpwEvaluatorConfig config_;
    GpwResult last_{};
};

} // namespace gpw
} // namespace cpu
} // namespace ATV
