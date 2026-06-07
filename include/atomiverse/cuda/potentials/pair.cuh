/// @file cuda/potentials/pair.cuh
/// @brief LAMMPS-style abstract base for non-bonded pair styles (CUDA).
///
/// Mirror of `cpu::Pair`. Independent hierarchy — no shared ancestor
/// between CPU and CUDA potential layers; only the names match.
#pragma once

#include "atomiverse/cuda/component.cuh"

namespace ATV {
namespace cuda {

/// @brief Abstract base for non-bonded pair styles (CUDA backend).
class Pair : public ComponentBase {
  public:
    ~Pair() override = default;

    ComponentRole role() const final { return ComponentRole::Evaluator; }
    Capability capabilities() const override {
        return Capability::Energy | Capability::Force | Capability::Stress;
    }
    /// @brief Pair styles depend on a neighbor list — used by
    ///        `cuda::md::MD::install_infra_` to detect a missing
    ///        `cuda::NeighborList` registration without `dynamic_cast`.
    bool needs_neighbor_list() const final { return true; }

    /// @brief Framework dispatch — `final`, calls `compute()`.
    void launch(const LaunchContext& ctx) final { compute(ctx); }

    /// @brief Per-type pair cutoff [Å]. Reports the value stored at
    /// construction; derived classes can override to compute from
    /// per-type tables.
    double pair_cutoff() const { return cutoff_; }

    /// @brief Engine-facing hook: enable full 1-2/1-3/1-4 special-bonds
    ///        exclusions (s12=s13=s14=0). Default no-op so the force engine
    ///        (`cuda::md::MultiGpu`) can request exclusions on any `Pair`
    ///        without knowing its concrete type. A distinct name from the
    ///        derived ref-qualified `set_special_bonds(…)` setters (which a
    ///        non-ref-qualified base virtual could not override).
    virtual void enable_full_exclusions() {}

  protected:
    explicit Pair(double cutoff) : cutoff_(cutoff) {}

    /// @brief Hot-path override — derived class launches its CUDA kernel.
    virtual void compute(const LaunchContext& ctx) = 0;

    double cutoff_ = 0.0;
};

} // namespace cuda
} // namespace ATV
