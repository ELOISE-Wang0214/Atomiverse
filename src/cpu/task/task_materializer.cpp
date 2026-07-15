/// @file cpu/task/task_materializer.cpp
/// @brief TaskSpec registration replay onto MD / Evaluate pipelines.

#include "atomiverse/cpu/task/task_materializer.h"

#include <stdexcept>
#include <type_traits>
#include <utility>

#include "atomiverse/cpu/md/init.h"
#include "atomiverse/cpu/md/nose_hoover.h"
#include "atomiverse/cpu/md/velocity_verlet.h"
#include "atomiverse/cpu/neighbor/neighbor_list.h"
#include "atomiverse/cpu/potentials/angle_harmonic.h"
#include "atomiverse/cpu/potentials/bond_harmonic.h"
#include "atomiverse/cpu/potentials/dihedral_multi_harmonic.h"
#include "atomiverse/cpu/potentials/pair_lj_cut.h"
#include "atomiverse/cpu/recorder/thermo.h"
#include "atomiverse/cpu/semiempirical/am1_evaluator.h"

namespace ATV {
namespace cpu {
namespace task {

namespace {

/// Replay one entry onto the MD driver exactly as the bindings registered
/// it: default-binding `add()` when `binding` is empty, explicit
/// `on(phase).add(comp, order)` otherwise.
void register_entry(md::MD& md, const ComponentEntry& entry) {
    const auto reg = [&](auto comp) {
        if (entry.binding)
            md.on(entry.binding->phase)
                .add(std::move(comp), entry.binding->order);
        else
            md.add(std::move(comp));
    };
    std::visit(
        [&](const auto& p) {
            using P = std::decay_t<decltype(p)>;
            if constexpr (std::is_same_v<P, PairLJCutParams>)
                reg(cpu::PairLJCut{p.epsilon, p.sigma, p.cutoff});
            else if constexpr (std::is_same_v<P, NeighborListParams>)
                reg(cpu::md::NeighborList<>{p.cutoff, p.skin});
            else if constexpr (std::is_same_v<P, BondHarmonicParams>)
                reg(cpu::BondHarmonic{p.K, p.r0});
            else if constexpr (std::is_same_v<P, AngleHarmonicParams>)
                reg(cpu::AngleHarmonic{p.K, p.theta0});
            else if constexpr (std::is_same_v<P, DihedralMultiHarmonicParams>)
                reg(cpu::DihedralMultiHarmonic{p.Cn});
            else if constexpr (std::is_same_v<P, VelocityVerletParams>)
                reg(cpu::md::VelocityVerlet{p.dt});
            else if constexpr (std::is_same_v<P, InitVelocityParams>)
                reg(cpu::md::InitVelocity{p.temperature, p.seed});
            else if constexpr (std::is_same_v<P, NoseHooverParams>)
                reg(cpu::md::NoseHoover{p.T_start, p.tau_T, p.dt});
            else if constexpr (std::is_same_v<P, ThermoParams>)
                reg(cpu::Thermo{p.freq});
            else if constexpr (std::is_same_v<P, AM1EvaluatorParams>)
                // Unreachable: the MD TaskSpec ctor rejects AM1 entries.
                throw std::logic_error(
                    "task_materializer: AM1Evaluator entry inside an MD "
                    "task (TaskSpec invariant breached)");
        },
        entry.params);
}

} // namespace

TaskMaterializer::TaskMaterializer(const TaskSpec& spec)
    : kind_(spec.kind()), n_steps_(spec.n_steps()) {
    if (kind_ == TaskKind::MD) {
        md_.emplace();
        for (const auto& entry : spec.entries())
            register_entry(*md_, entry);
        md_->build();
        return;
    }
    // SinglePoint (v1: AM1). The TaskSpec ctor guarantees exactly one
    // AM1 entry with no explicit binding; the evaluator's phase is pinned
    // here (evaluate::Phase::evaluate), not on the wire.
    const auto& p = std::get<AM1EvaluatorParams>(spec.entries().front().params);
    semiempirical::AM1Config cfg;
    cfg.charge = p.charge;
    cfg.multiplicity = p.multiplicity;
    cfg.max_iter = p.max_iter;
    cfg.density_tol = p.density_tol;
    cfg.energy_tol = p.energy_tol;
    eval_.emplace();
    eval_->on(ATV::evaluate::Phase::evaluate)
        .add(semiempirical::AM1Evaluator{cfg});
    eval_->build();
}

void TaskMaterializer::run(cpu::Frame& frame) {
    if (kind_ == TaskKind::MD)
        md_->run(frame, static_cast<int>(n_steps_));
    else
        eval_->run(frame);
}

} // namespace task
} // namespace cpu
} // namespace ATV
