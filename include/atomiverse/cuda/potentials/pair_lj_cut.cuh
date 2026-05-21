/// @file cuda/potentials/pair_lj_cut.cuh
/// @brief LAMMPS-style PairLJCut — 12-6 with shifted potential (CUDA).
///
/// Was `ATV::LJ` before the LAMMPS-naming refactor. The legacy header
/// `lj.cuh` provides `using LJ = PairLJCut;` so existing user code
/// compiles unchanged.
#pragma once

#include <string>

#include "atomiverse/cuda/potentials/pair.cuh"

namespace ATV {
namespace cuda {

/// @brief Lennard-Jones 12-6 pair potential evaluator (CUDA).
///
/// Must be explicitly bound to a force phase (cross-driver component).
/// Requires a populated PairList in DeviceFrame (from NeighborList).
class PairLJCut : public Pair {
  public:
    /// @brief Construct with LJ parameters.
    /// @param epsilon Well depth [eV].
    /// @param sigma   Collision diameter [Å].
    /// @param cutoff  Interaction cutoff [Å].
    PairLJCut(double epsilon, double sigma, double cutoff);

    std::string name() const override { return "lj"; }

    double epsilon() const { return epsilon_; }
    double sigma() const { return sigma_; }
    double cutoff() const {
        return cutoff_;
    } // shadows base for ergonomic access

  protected:
    void compute(const LaunchContext& ctx) override;

  private:
    double epsilon_; ///< Well depth [eV].
    double sigma_;   ///< Collision diameter [Å].
    double v_shift_; ///< Potential shift at cutoff [eV].
};

} // namespace cuda
} // namespace ATV
