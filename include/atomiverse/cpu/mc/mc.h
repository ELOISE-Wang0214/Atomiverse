/// @file cpu/mc/mc.h
/// @brief cpu::mc::MC — concrete canonical (NVT) Metropolis Monte Carlo driver.
#pragma once

#include <cstdint>
#include <optional>
#include <string>
#include <utility>

#include "atomiverse/core/precision.h"
#include "atomiverse/cpu/component/zero.h"
#include "atomiverse/cpu/driver/sampling_driver.h"
#include "atomiverse/facade/stage.h"

namespace ATV {
namespace cpu {
namespace mc {

/// @brief CPU canonical (NVT) Metropolis Monte Carlo driver.
///
/// The concrete Sampling-category driver: it runs a fixed budget of trial
/// moves, each proposed by a `Mover` (v1: `TranslateMove`), scored by any
/// registered `Capability::Energy` evaluator (e.g. `PairLJCut`), and accepted
/// or rejected by the inline Metropolis rule in `SamplingDriver::run()`. It
/// reuses the MD component set unchanged — energy evaluator, neighbor list,
/// frame — and mirrors the shape of `cpu::md::MD`.
///
/// The constructor binds `cpu::ZeroEnergy` at the head of the `energy` phase
/// (order -999) so the evaluator accumulates onto a freshly zeroed
/// `frame.scalars.energy` each trial. Categorical setters `set_pair` /
/// `set_neighbor` / `set_observer` / `set_record` / `set_thermo` are inherited
/// from `CPUDriver`; `set_mover` is added here (mirroring MD's
/// `set_integrator`). The RNG, proposal changeset, and cached accepted energy
/// live on the `SamplingDriver` base.
///
/// @note For MC the neighbor list should rebuild whenever an accepted move
///       displaces an atom past the skin — construct it with `delay = 0`
///       (`NeighborList<>{cutoff, skin, 1, 0, true}`) so the displacement
///       check is never deferred.
///
/// Usage:
/// @code{.cpp}
///   cpu::mc::MC mc{/*T [K]=*/300.0, /*seed=*/12345u};
///   mc.set_pair(cpu::PairLJCut{eps, sigma, cutoff});
///   mc.set_neighbor(cpu::md::NeighborList<>{cutoff, skin, 1, 0, true});
///   mc.set_mover(cpu::mc::TranslateMove{/*delta [Å]=*/0.2});
///   mc.set_thermo(100);
///   mc.build();
///   mc.run(frame, 100000);            // 100k trial moves
///   double e = mc.current_energy();   // [eV]; mc.acceptance_ratio() in [0,1]
/// @endcode
class MC : public SamplingDriver {
  public:
    /// @brief Construct the NVT sampler.
    /// @param temperature_K Ensemble temperature [K] (must be > 0).
    /// @param seed Seed for the single reproducible RNG stream.
    /// @throws std::invalid_argument if `temperature_K <= 0`.
    MC(AccReal temperature_K, std::uint64_t seed)
        : SamplingDriver(temperature_K, seed) {
        // Zero the energy accumulator at the head of the trial-energy phase so
        // the evaluator accumulates onto a clean slate each move (mirrors MD's
        // ZeroEnergy at pre_force -999).
        this->on(ATV::sampling::Phase::energy).add(cpu::ZeroEnergy{}, -999);
    }

    /// @brief Configure the trial-move generator (mirrors MD `set_integrator`).
    /// Singular: re-call with the same Mover is idempotent, a different Mover
    /// throws. Binds the Mover at its `default_bindings()` (`propose` phase).
    /// @param m Mover component (e.g. `TranslateMove`).
    /// @return `*this` for chaining.
    template <typename M>
    MC& set_mover(M&& m) {
        this->set_singular_default_(std::forward<M>(m), mover_set_, "set_mover");
        return *this;
    }

  protected:
    /// @brief Enforce that a pair evaluator that needs a neighbor list has one.
    void install_defaults() override {
        const auto s = this->scan_pending_infra_();
        this->throw_if_missing_neighbor_list_(
            s, "MC",
            "Add `mc.set_neighbor(cpu::md::NeighborList<>{cutoff, skin, 1, 0, "
            "true})`.");
    }

  private:
    std::optional<std::string> mover_set_; ///< set_mover singular slot.
};

} // namespace mc
} // namespace cpu
} // namespace ATV
