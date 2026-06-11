/// @file cpu/evaluate.h
/// @brief cpu::Evaluate — concrete one-pass evaluator driver.
///
/// `cpu::Evaluate` is the user-facing CPU outer-driver under the
/// `EvaluateDriver` category. It runs the registered evaluator
/// components exactly once on the current `cpu::Frame` and returns —
/// no time advance, no minimization loop.
///
/// Single-point energy / force / property workflows go here. Exposes
/// the same LAMMPS-shaped categorical setter API as `cpu::md::MD`
/// (minus the time-evolution and optimizer setters): set_pair,
/// set_bond, set_angle, set_dihedral, set_neighbor, set_observer,
/// set_record — all inherited from CPUDriver.
///
/// Auto-adds ZeroForces and ZeroEnergy at `evaluate::Phase::prepare`
/// (orders -1000 / -999) so frame fields start clean before evaluators
/// write to them.
///
/// `add()` and `on(phase).add()` from the underlying `Driver<>`
/// remain public for user-defined hooks (e.g. SCF evaluators
/// registered via `add(SCFEvaluator{...})`).
///
/// Inheritance: Evaluate → EvaluateDriver → CPUDriver<evaluate::Phase> →
/// Driver<evaluate::Phase, cpu::device_tag>.
/// CPUDriver provides the shared categorical setters and slot variables.
///
/// @code{.cpp}
///   // SCF energy + forces (default)
///   cpu::Evaluate eval;
///   eval.add(cpu::scf::SCFEvaluator{cpu::scf::default_config(integrals)});
///   eval.build();
///   eval.run(frame);
///
///   // Force-field force evaluation (LAMMPS-shaped)
///   cpu::Evaluate eval2;
///   eval2.set_pair(cpu::md::LJ{eps, sigma, cutoff});
///   eval2.set_neighbor(cpu::md::NeighborList<>{cutoff, skin});
///   eval2.build();
///   eval2.run(frame);
/// @endcode
#pragma once

#include "atomiverse/cpu/evaluate_driver.h"
#include "atomiverse/cpu/pair_evaluator.h"
#include "atomiverse/cpu/zero.h"
#include "atomiverse/facade/stage.h"

namespace ATV {
namespace cpu {

/// @brief One-pass evaluator driver.
///
/// Categorical setters (set_pair, set_bond, set_angle, set_dihedral,
/// set_neighbor, set_observer, set_record) are inherited from CPUDriver.
class Evaluate : public EvaluateDriver {
  public:
    Evaluate() {
        this->on(ATV::evaluate::Phase::prepare)
            .add(cpu::ZeroForces{}, -1000);
        this->on(ATV::evaluate::Phase::prepare)
            .add(cpu::ZeroEnergy{}, -999);
    }

  protected:
    void install_defaults() override {
        const auto s = this->scan_pending_infra_();
        this->throw_if_missing_neighbor_list_(
            s, "Evaluate",
            "Add `eval.set_neighbor(cpu::md::NeighborList<>{cutoff, skin})`.");
    }
};

} // namespace cpu
} // namespace ATV
