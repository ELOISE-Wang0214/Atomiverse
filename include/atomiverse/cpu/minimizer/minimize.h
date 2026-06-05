/// @file cpu/minimizer/minimize.h
/// @brief cpu::minimizer::Minimize — concrete energy-minimizer driver.
///
/// `Minimize` is the user-facing CPU outer-driver under the
/// `OptimizeDriver` category. It iterates an optimizer component
/// (LBFGS / CG / future FIRE) until the registered convergence-check
/// component reports converged, or until `max_iter` is reached.
///
/// Exposes the same LAMMPS-shaped categorical setter API as
/// `cpu::md::MD` (minus the time-evolution setters): set_pair,
/// set_bond, set_angle, set_dihedral, set_neighbor, set_optimizer,
/// set_observer, set_record, set_thermo — all except set_optimizer
/// are inherited from CPUDriver.
///
/// `add()` and `on(phase).add()` from the underlying `Driver<>`
/// remain public for user-defined hooks.
#pragma once

#include <string>
#include <type_traits>

#include "atomiverse/cpu/md/init.h"
#include "atomiverse/cpu/optimize_driver.h"
#include "atomiverse/cpu/pair_evaluator.h"
#include "atomiverse/cpu/snapshot_publisher.h"
#include "atomiverse/facade/stage.h"

namespace ATV {
namespace cpu {
namespace minimizer {

/// @brief Energy-minimization driver — convergence-aware Optimize loop.
///
/// Auto-adds ZeroForces and ZeroEnergy at force phase. Wires the shared
/// `cpu::Publish` into `end_of_step @ +10` so async observers
/// see snapshots.
///
/// Validation rules inherited from `OptimizeDriver`: requires ≥ 1
/// `Optimizer` component (e.g. `LBFGS`, `CG`) and ≥ 1 `Evaluator`
/// with `Capability::Energy` or `Capability::Force` (e.g. `SCFEvaluator`).
///
/// Inheritance: Minimize → OptimizeDriver → CPUDriver<minimizer::Phase> →
/// Driver<minimizer::Phase, cpu::device_tag>.
/// CPUDriver provides the shared categorical setters and slot variables.
class Minimize : public OptimizeDriver {
  public:
    Minimize() {
        // Phase-agnostic Zero* from cpu::md, bound to minimizer phases.
        this->on(ATV::minimizer::Phase::force)
            .add(cpu::md::ZeroForces{}, -1000)
            .add(cpu::md::ZeroEnergy{}, -999);
        this->on(ATV::minimizer::Phase::end_of_step)
            .add(cpu::Publish{}, 10);
    }

    // ─── Minimizer-specific categorical setters ─────────────────────
    //
    // set_pair, set_bond, set_angle, set_dihedral, set_neighbor,
    // set_observer, set_record, set_thermo are inherited from CPUDriver.

    /// @brief Configure the optimization algorithm (LAMMPS `min_style`).
    template <typename O>
    Minimize& set_optimizer(O&& o) {
        this->set_singular_default_(std::forward<O>(o), optimizer_set_,
                                    "set_optimizer");
        return *this;
    }

  protected:
    void install_defaults() override {
        const auto s = this->scan_pending_infra_();
        this->throw_if_missing_neighbor_list_(
            s, "Minimize",
            "Add `opt.set_neighbor(cpu::md::NeighborList<>{cutoff, skin})`.");
    }

  private:
    std::optional<std::string> optimizer_set_;
};

} // namespace minimizer
} // namespace cpu
} // namespace ATV
