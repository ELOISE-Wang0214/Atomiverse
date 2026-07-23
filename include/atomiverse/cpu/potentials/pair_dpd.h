/// @file cpu/potentials/pair_dpd.h
/// @brief Groot–Warren DPD: soft conservative + pairwise D/R (CPU).
#pragma once

#include <cstdint>
#include <random>
#include <string>
#include <vector>

#include "atomiverse/cpu/potentials/pair.h"

namespace ATV {
namespace cpu {

/// @brief Dissipative particle dynamics pair force (Groot–Warren 1997).
///
/// Conservative: \(F^C = a\,w(r)\hat r\), \(w=1-r/r_c\).
/// Multi-type: optional `a_matrix[type_i][type_j]` (symmetric).
/// Dissipative / random with \(\sigma^2 = 2\gamma k_B T\) and weight \(w\).
class PairDpd : public Pair {
  public:
    /// @param a Soft repulsion amplitude [eV/Å] (single-type).
    PairDpd(double a, double gamma, double temperature, double cutoff,
            double dt, std::uint64_t seed = 1);

    /// @brief Multi-component: square a_matrix[t_i][t_j] [eV/Å], size ntypes².
    PairDpd(std::vector<std::vector<double>> a_matrix, double gamma,
            double temperature, double cutoff, double dt,
            std::uint64_t seed = 1);

    std::string name() const override { return "dpd"; }

    void set_step(int step) { step_ = step; }

    const NeighborConsumer* as_neighbor_consumer() const final {
        return nullptr;
    }

  protected:
    void compute(const LaunchContext& ctx) override;

  private:
    double a_default_ = 0;
    std::vector<std::vector<double>> a_mat_;
    bool multi_ = false;
    double gamma_ = 0;
    double temperature_ = 0;
    double dt_ = 1;
    double cutoff_sq_ = 0;
    double sigma_ = 0;
    std::uint64_t seed_ = 1;
    int step_ = 0;

    double a_of_(int ti, int tj) const;
};

} // namespace cpu
} // namespace ATV
