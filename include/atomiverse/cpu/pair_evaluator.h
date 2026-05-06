/// @file cpu/pair_evaluator.h
/// @brief PairEvaluator — backend-root base for cutoff-pair force fields.
///
/// Lives at `cpu/` root (not `cpu/md/`) because cutoff-pair forces are
/// shared infrastructure across multiple domains:
///   - MD (`cpu::md::MD`) — pair forces in the per-step force kernel.
///   - Minimizer (`cpu::minimizer::Minimize`) — same evaluators reused
///     in geometry optimisation.
///   - Evaluator (`cpu::Evaluate`) — single-point pair energy / force.
///   - Future Monte Carlo / WIDOM / replica-exchange drivers — same
///     pair-cutoff abstraction, different outer loop.
///
/// Concrete subclasses (today: `cpu::md::LJ`) live in their domain
/// directory; only the abstraction is at root.
///
/// `needs_neighbor_list()` is forced to `true` here so subclasses
/// don't have to remember it — every cutoff-pair force field needs
/// a `NeighborList` upstream. The CPU category drivers throw at
/// `build()` time when this returns true and no `set_neighbor()`
/// was made.
#pragma once

#include "atomiverse/core/real.h"
#include "atomiverse/cpu/component.h"

namespace ATV {
namespace cpu {

/// @brief Backend-root base for cutoff-pair force fields.
///
/// Subclasses must override `pair_cutoff()` to report the interaction
/// cutoff in Ångströms. `needs_neighbor_list()` is always `true` for
/// pair-cutoff evaluators.
class PairEvaluator : public cpu::ComponentBase {
public:
    /// @brief Pair interaction cutoff [Å]. Used by category drivers
    /// to size MPI halo exchange and surfaced in error messages when
    /// a `NeighborList` is missing.
    virtual Real pair_cutoff() const = 0;

    /// @brief Pair-cutoff evaluators always need a neighbor list.
    bool needs_neighbor_list() const override { return true; }
};

} // namespace cpu
} // namespace ATV
