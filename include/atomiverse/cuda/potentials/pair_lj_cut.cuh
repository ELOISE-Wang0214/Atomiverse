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
/// Requires a populated PairList published on `cuda::Context::pair_list`
/// by a NeighborList component.
class PairLJCut : public Pair {
  public:
    /// @brief Construct with LJ parameters.
    /// @param epsilon Well depth [eV].
    /// @param sigma   Collision diameter [Å].
    /// @param cutoff  Interaction cutoff [Å].
    PairLJCut(double epsilon, double sigma, double cutoff);

    std::string name() const override { return "lj"; }

    /// @brief Enable LAMMPS-style `special_bonds` 1-2 / 1-3 / 1-4
    ///        pair exclusions.
    ///
    /// Mirrors LAMMPS `special_bonds lj/coul s12 s13 s14` semantics.
    /// The CUDA arm supports two combinations in this revision (per
    /// `msu-pe-cuda-components-05-cuda-lj-exclusions` design):
    ///   - `(1.0, 1.0, 1.0)` — no exclusions, kernel takes fast path.
    ///   - All three equal (typically `(0.0, 0.0, 0.0)` for Dreiding-UA
    ///     polyethylene): every 1-2/1-3/1-4 pair is multiplied by the
    ///     shared scale `s`. `s == 0` ⇒ completely skipped.
    /// Other combinations throw `std::invalid_argument`. The supported
    /// pair (0,0,0) is what MSU PE protocol uses.
    ///
    /// Exclusion partner lists are enumerated from `frame.bonds /
    /// angles / dihedrals` at `build()` time, sorted+uniqued per atom,
    /// and uploaded to device as a CSR. The kernel does a per-pair
    /// binary search against the CSR.
    ///
    /// @param s12 Scale for 1-2 (bonded) pairs.
    /// @param s13 Scale for 1-3 (across-angle) pairs.
    /// @param s14 Scale for 1-4 (across-dihedral) pairs.
    /// @return `*this`, chainable.
    PairLJCut& set_special_bonds(double s12, double s13, double s14);

    void build(const BuildContext& ctx) override;
    void destroy(const DestroyContext& ctx) override;

    double epsilon() const { return epsilon_; }
    double sigma() const { return sigma_; }
    double cutoff() const { return cutoff_; }

  protected:
    void compute(const LaunchContext& ctx) override;

  private:
    double epsilon_; ///< Well depth [eV].
    double sigma_;   ///< Collision diameter [Å].
    double v_shift_; ///< Potential shift at cutoff [eV].

    // Special-bonds configuration (host-side).
    double s12_ = 1.0, s13_ = 1.0, s14_ = 1.0;
    double s_excl_ = 1.0; ///< Single scale used by kernel when collapsed.
    bool has_exclusions_ = false;

    // Device CSR — d_excl_offsets_[n_atoms+1], d_excl_partners_[n_excl_total_].
    // Both `nullptr` when has_exclusions_ is false (kernel fast path).
    int* d_excl_offsets_ = nullptr;
    int* d_excl_partners_ = nullptr;
    int n_excl_total_ = 0;
};

} // namespace cuda
} // namespace ATV
