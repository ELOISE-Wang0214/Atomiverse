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

  protected:
    explicit Pair(double cutoff) : cutoff_(cutoff) {}

    /// @brief Hot-path override — derived class launches its CUDA kernel.
    virtual void compute(const LaunchContext& ctx) = 0;

    double cutoff_ = 0.0;
};

} // namespace cuda
} // namespace ATV
