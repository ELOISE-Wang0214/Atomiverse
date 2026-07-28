/// @file cpu/potentials/pair_gaussian_coul.h
/// @brief Real-space Coulomb with Gaussian electrode charges (CPU).
///
/// Evaluates the Reed / MetalWalls Gaussian charge energy and forces over
/// the published pair list (min-image) with a hard cutoff. Per-type width
/// \(\eta\) [Å⁻¹]: \(\eta>0\) ⇒ Gaussian electrode site; \(\eta\le 0\) ⇒ point
/// charge (electrolyte). Optional Gaussian self-energy (default on).
///
/// A `NeighborConsumer` through `Pair`: the driver requires a `NeighborList`
/// rather than letting the style drop, unannounced, onto the O(N²) sweep it
/// keeps for the no-list case.
#pragma once

#include <string>
#include <vector>

#include "atomiverse/cpu/potentials/pair.h"

namespace ATV {
namespace cpu {

/// @brief Gaussian / point Coulomb pair potential for electrode MD footholds.
///
/// Parameters are constructor arguments — no public `*Config` bag.
class PairGaussianCoul : public Pair {
  public:
    /// @brief Construct with a hard cutoff and per-type Gaussian widths.
    ///
    /// @param cutoff      Pair cutoff [Å]; must be > 0.
    /// @param eta_by_type Width \(\eta\) [Å⁻¹] indexed by `kAtom_type`.
    ///                    Entry \(\le 0\) means point charge for that type.
    ///                    Empty vector ⇒ every atom is a point charge.
    /// @param include_self If true (default), add Gaussian self-energy
    ///                     \(\frac12\texttt{coulomb}\,q_i^2 G_{ii}\) for
    ///                     Gaussian types.
    PairGaussianCoul(double cutoff, std::vector<double> eta_by_type,
                     bool include_self = true);

    std::string name() const override { return "gaussian/coul"; }

    // NeighborConsumer via Pair — production depth uses ctx.md.pair_list when
    // published; falls back to O(N²) if no list (small textbook cells).

    /// @brief Whether Gaussian self-energy is included.
    bool include_self() const { return include_self_; }

    /// @brief Width table (copy) [Å⁻¹].
    const std::vector<double>& eta_by_type() const { return eta_by_type_; }

  protected:
    void compute(const LaunchContext& ctx) override;

  private:
    double cutoff_sq_ = 0;
    std::vector<double> eta_by_type_;
    bool include_self_ = true;

    double eta_of(int type) const {
        if (type < 0 || static_cast<std::size_t>(type) >= eta_by_type_.size())
            return 0.0; // unknown type ⇒ point
        return eta_by_type_[static_cast<std::size_t>(type)];
    }
};

} // namespace cpu
} // namespace ATV
