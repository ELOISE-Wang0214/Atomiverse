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

#include <cstdint>

#include "atomiverse/cpu/component/long_range_protocol.h"
#include "atomiverse/cpu/component/potential.h"

namespace ATV {
namespace cpu {

/// @brief Observable cold-event ledger for reciprocal-space components.
struct KSpaceLifecycleStats {
    std::uint64_t builds = 0;
    std::uint64_t cell_refreshes = 0;
    std::uint64_t capacity_rebuilds = 0;
};

/// @brief Abstract base for reciprocal-space electrostatics styles (CPU).
class KSpace : public cpu::Potential, public cpu::LongRangeParticipant {
  public:
    ~KSpace() override = default;

    /* role via Potential/kind base */
    const LongRangeParticipant* as_long_range_participant() const final {
        return this;
    }

    /// @return Cumulative lifecycle events for this component instance.
    const KSpaceLifecycleStats& lifecycle_stats() const noexcept {
        return lifecycle_stats_;
    }

    /// @brief Framework dispatch — `final`, calls `compute()` (mirrors
    ///        `cpu::Pair` / `cpu::Bond`).
    void launch(const LaunchContext& ctx) final { compute(ctx); }

  protected:
    KSpace() = default;

    void note_build() noexcept { ++lifecycle_stats_.builds; }
    void note_cell_refresh() noexcept { ++lifecycle_stats_.cell_refreshes; }
    void note_capacity_rebuild() noexcept {
        ++lifecycle_stats_.capacity_rebuilds;
    }

    /// @brief Hot-path override — derived class runs its reciprocal-space sum
    ///        (direct k-sum for Ewald, FFT for PME).
    virtual void compute(const LaunchContext& ctx) = 0;

  private:
    KSpaceLifecycleStats lifecycle_stats_{};
};

} // namespace cpu
} // namespace ATV
