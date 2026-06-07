/// @file cpu/kspace/kspace.h
/// @brief Abstract base for reciprocal-space (long-range) electrostatics — CPU.
///
/// A first-class force *family*, sibling to `cpu::Pair` / `cpu::Bond` /
/// `cpu::Angle` / `cpu::Dihedral`: the force engine composes one KSpace
/// alongside them (LAMMPS `Force::kspace`, `src/force.h:95`). Concrete styles
/// (`KSpaceEwald`, `KSpacePME`) implement only the reciprocal-space
/// `compute()`; the short-range real-space coulomb is a separate evaluator
/// (`cpu::CoulLong`). The two share the Ewald splitting parameter α.
///
/// The reciprocal sum needs a *global* structure factor
/// `S(k) = Σ_i q_i e^{ik·r_i}`, so the multi-rank reduction (MPI_Allreduce over
/// partial S(k)) is a KSpace-specific communication distinct from the per-atom
/// halo (force-infra-04 formalizes that hook on this base).
#pragma once

#include "atomiverse/cpu/component.h"

namespace ATV {
namespace cpu {

/// @brief Abstract base for reciprocal-space electrostatics styles (CPU).
class KSpace : public cpu::ComponentBase {
  public:
    ~KSpace() override = default;

    ComponentRole role() const final { return ComponentRole::Evaluator; }

    /// @brief Framework dispatch — `final`, calls `compute()` (mirrors
    ///        `cpu::Pair` / `cpu::Bond`).
    void launch(const LaunchContext& ctx) final { compute(ctx); }

  protected:
    KSpace() = default;

    /// @brief Hot-path override — derived class runs its reciprocal-space sum
    ///        (direct k-sum for Ewald, FFT for PME).
    virtual void compute(const LaunchContext& ctx) = 0;
};

} // namespace cpu
} // namespace ATV
