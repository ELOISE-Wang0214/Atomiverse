/// @file cpu/potentials/pair_dpd.h
/// @brief Groot–Warren DPD: soft conservative + pairwise D/R (CPU).
#pragma once

#include <cstdint>
#include <random>
#include <string>

#include "atomiverse/cpu/potentials/pair.h"

namespace ATV {
namespace cpu {

/// @brief Dissipative particle dynamics pair force (Groot–Warren 1997).
///
/// Conservative: \(F^C = a\,w(r)\hat r\), \(w=1-r/r_c\).
/// Dissipative / random with \(\sigma^2 = 2\gamma k_B T\) and weight \(w\).
/// Random stream is deterministic from `seed` and step via `set_step`.
///
/// Energy accumulates only the conservative pair potential
/// \(U=\frac a2\sum_{ij} w(r)^2\) (integral of F^C).
class PairDpd : public Pair {
  public:
    /// @param a Soft repulsion amplitude [eV/Å].
    /// @param gamma Dissipation coefficient [eV·fs/Å²] (mass·time⁻¹ scale).
    /// @param temperature Target T [K] for fluctuation–dissipation.
    /// @param cutoff Interaction range r_c [Å].
    /// @param dt Timestep used for random force scaling [fs].
    /// @param seed RNG seed (dimensionless).
    PairDpd(double a, double gamma, double temperature, double cutoff,
            double dt, std::uint64_t seed = 1);

    std::string name() const override { return "dpd"; }

    /// @brief Publish current MD step for reproducible random forces.
    void set_step(int step) { step_ = step; }

    const NeighborConsumer* as_neighbor_consumer() const final {
        return nullptr;
    }

  protected:
    void compute(const LaunchContext& ctx) override;

  private:
    double a_ = 0;
    double gamma_ = 0;
    double temperature_ = 0;
    double dt_ = 1;
    double cutoff_sq_ = 0;
    double sigma_ = 0; ///< √(2 γ k_B T)
    std::uint64_t seed_ = 1;
    int step_ = 0;
};

} // namespace cpu
} // namespace ATV
