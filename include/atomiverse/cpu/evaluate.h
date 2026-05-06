/// @file cpu/evaluate.h
/// @brief cpu::Evaluate — concrete one-pass evaluator driver.
///
/// `cpu::Evaluate` is the user-facing CPU outer-driver under the
/// `EvaluateDriver` category. It runs the registered evaluator
/// components exactly once on the current `cpu::Frame` and returns —
/// no time advance, no minimization loop.
///
/// Single-point energy / force / property workflows go here:
///
/// @code
///   // SCF energy + forces (default)
///   cpu::Evaluate eval;
///   cpu::scf::SCFEvaluator scf{cpu::scf::default_config(integrals)};
///   eval.add(scf);
///   eval.build();
///   eval.run(frame);
///   // frame.energy + frame.atoms.f* populated.
///
///   // SCF energy only (no gradient)
///   cpu::Evaluate eval2;
///   cpu::scf::SCFEvaluator scf2{cfg};
///   scf2.with_gradient(false);
///   eval2.add(scf2);
///   eval2.build();
///   eval2.run(frame);
///   // frame.energy populated; forces unchanged.
///
///   // Force-field force evaluation
///   cpu::Evaluate eval3;
///   eval3.add(cpu::md::LJ{eps, sigma, cutoff});
///   eval3.add(cpu::md::NeighborList<>{cutoff});
///   eval3.build();
///   eval3.run(frame);
/// @endcode
///
/// Auto-adds `cpu::md::ZeroForces` and `cpu::md::ZeroEnergy` at
/// `evaluate::Phase::prepare` (orders -1000 / -999) so frame fields
/// start clean before evaluators write to them. The MD-flavoured Zero*
/// components are reused via explicit `on(phase).add()` binding —
/// their launch bodies are phase-agnostic (they just zero the frame
/// fields).
#pragma once

#include "atomiverse/cpu/evaluate_driver.h"
#include "atomiverse/cpu/md/init.h"
#include "atomiverse/facade/stage.h"

namespace ATV {
namespace cpu {

/// @brief One-pass evaluator driver.
class Evaluate : public EvaluateDriver {
public:
    Evaluate() {
        // Override the MD-flavoured default_bindings (which target
        // md::Phase::pre_force) by binding explicitly to evaluate::Phase.
        this->on(ATV::evaluate::Phase::prepare).add(cpu::md::ZeroForces{}, -1000);
        this->on(ATV::evaluate::Phase::prepare).add(cpu::md::ZeroEnergy{},  -999);
    }
};

} // namespace cpu
} // namespace ATV
