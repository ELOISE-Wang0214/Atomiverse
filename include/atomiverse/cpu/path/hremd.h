/// @file cpu/path/hremd.h
/// @brief Hamiltonian REMD — replica exchange on a λ ladder at fixed T.
#pragma once

#include <cstdint>
#include <memory>
#include <random>
#include <vector>

#include "atomiverse/cpu/dynamics/md.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/replica/replica_expansion.h"
#include "atomiverse/cpu/replica/replica_set.h"

namespace ATV {
namespace cpu {
namespace path {

/// @brief Configuration for Hamiltonian REMD (fixed temperature, λ ladder).
struct HRemdConfig {
    std::vector<double> lambdas; ///< Distinct λ values, size ≥ 2.
    double temperature = 300.0;  ///< Fixed T [K] for all replicas.
    double dt = 1.0;             ///< Nested MD timestep [fs].
    int steps_per_exchange = 10; ///< Nested MD steps between exchange rounds.
    std::uint64_t seed = 1;      ///< Master seed.
};

/// @brief Hamiltonian REMD: outer exchange loop + per-replica MD programs.
///
/// Exchange uses the fixed-T Hamiltonian criterion with dual energy
/// evaluation U_i(x_j), U_j(x_i) on adjacent λ pairs. Sibling of `Remd`
/// (temperature ladder); does not replace it.
class HRemd {
  public:
    explicit HRemd(HRemdConfig cfg);

    /// @brief Program builder: register forces for one λ (e.g. PairLJSoft).
    using ProgramBuilder = void (*)(cpu::md::MD& md, double lambda, void* user);
    void set_program_builder(ProgramBuilder builder, void* user = nullptr);

    /// @brief Optional energy evaluator for dual U(λ) on a foreign frame.
    /// If unset, uses each program's last `kEnergy` only (same-λ energies —
    /// valid only when soft potentials expose λ on a shared component).
    /// Preferred: call `set_energy_fn` to compute U at arbitrary λ.
    using EnergyFn = double (*)(const cpu::Frame& frame, double lambda,
                                void* user);
    void set_energy_fn(EnergyFn fn, void* user = nullptr);

    void run(cpu::Frame& frame, int n_exchanges);

    int exchange_attempts() const { return attempts_; }
    int exchange_accepts() const { return accepts_; }

  private:
    void ensure_built_(cpu::Frame& frame);
    void propagate_();
    void exchange_round_(bool even_edges);
    double energy_at_(const cpu::Frame& frame, double lambda);

    HRemdConfig cfg_;
    ProgramBuilder builder_ = nullptr;
    void* builder_user_ = nullptr;
    EnergyFn energy_fn_ = nullptr;
    void* energy_user_ = nullptr;
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
