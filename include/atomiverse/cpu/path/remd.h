/// @file cpu/path/remd.h
/// @brief REMD — replica-exchange MD preset over IdentityReplication.
///
/// Per-replica owned MD trees (no Component::clone). Budget-terminating.
/// Exchange on even/odd adjacent pairs; collapse to lowest-T slot.
#pragma once

#include <cstdint>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/dynamics/md.h"
#include "atomiverse/cpu/replica/replica_expansion.h"
#include "atomiverse/cpu/replica/replica_set.h"

namespace ATV {
namespace cpu {
namespace path {

/// @brief Simple REMD driver: outer exchange loop + per-replica MD programs.
///
/// Not a new Driver category — a method preset that owns N `cpu::md::MD`
/// instances and a `ReplicaSet`. Parameters are constructor arguments /
/// public members — no `RemdConfig` bag.
class Remd {
  public:
    /// @brief Strictly increasing temperature ladder [K], size ≥ 2.
    std::vector<double> temperatures;
    /// @brief Nested MD timestep [fs].
    double dt = 1.0;
    /// @brief Nested MD steps between exchange rounds.
    int steps_per_exchange = 10;
    /// @brief Master seed (deterministic streams).
    std::uint64_t seed = 1;

    /// @param temperatures Strictly increasing T [K], size ≥ 2.
    /// @param dt Nested MD timestep [fs].
    /// @param steps_per_exchange Nested MD steps between exchanges (≥ 1).
    /// @param seed Master RNG seed.
    explicit Remd(std::vector<double> temperatures, double dt = 1.0,
                  int steps_per_exchange = 10, std::uint64_t seed = 1);

    /// @brief Install the force program builder: called once per replica at setup
    /// with the replica temperature. Must register pair/neighbor/etc. on `md`.
    using ProgramBuilder = void (*)(cpu::md::MD& md, double temperature,
                                    void* user);
    void set_program_builder(ProgramBuilder builder, void* user = nullptr);

    /// @brief Expand, build N MD programs, run `n_exchanges` rounds, collapse to T0.
    void run(cpu::Frame& frame, int n_exchanges);

    int exchange_attempts() const { return attempts_; }
    int exchange_accepts() const { return accepts_; }

  private:
    void ensure_built_(cpu::Frame& frame);
    void propagate_();
    void exchange_round_(bool even_edges);

    ProgramBuilder builder_ = nullptr;
    void* builder_user_ = nullptr;
    std::unique_ptr<IdentityReplication> expansion_;
    std::unique_ptr<ReplicaSet> band_;
    std::vector<std::unique_ptr<cpu::md::MD>> programs_;
    std::mt19937_64 exchange_rng_;
    bool built_ = false;
    int attempts_ = 0;
    int accepts_ = 0;
};

} // namespace path
} // namespace cpu
} // namespace ATV
