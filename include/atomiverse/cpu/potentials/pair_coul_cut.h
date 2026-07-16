/// @file cpu/potentials/pair_coul_cut.h
/// @brief Real-space Coulomb with hard cutoff (LAMMPS `pair_style coul/cut`) — CPU.
#pragma once

#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

#include "atomiverse/core/precision.h"
#include "atomiverse/cpu/component/pair_evaluator.h"
#include "atomiverse/cpu/neighbor/pair_list.h"
#include "atomiverse/facade/component_role.h"

namespace ATV {
namespace cpu {

/// @brief Truncated Coulomb pair style (no Ewald/PME partner required).
///
/// Optional special-bonds exclusions (LAMMPS convention): 1-2 / 1-3 / 1-4 pairs
/// with scale 0 are dropped via a filtered neighbor CSR, rebuilt whenever the
/// upstream PairList version advances (same contract as `PairLJCut`).
class PairCoulCut : public cpu::PairPotential {
  public:
    /// @param cutoff Real-space cutoff [Å].
    explicit PairCoulCut(double cutoff);

    std::string name() const override { return "coul/cut"; }
    /* role via Potential/kind base */
    /// @brief Enable LAMMPS-style special-bonds exclusions for Coulomb.
    /// @param s12 Scale for 1-2 pairs (0 = exclude fully; nonzero scales not
    ///            yet applied — only 0/1 gates are supported, matching LJ).
    /// @param s13 Scale for 1-3 pairs.
    /// @param s14 Scale for 1-4 pairs.
    PairCoulCut& set_special_bonds(double s12, double s13, double s14) &;
    PairCoulCut set_special_bonds(double s12, double s13, double s14) &&;

    void build(const BuildContext& ctx) override;
    void destroy(const DestroyContext& ctx) override;

  protected:
    void compute(const LaunchContext& ctx) override;

  private:
    void refresh_exclusion_map(const Frame& frame, int n_owners);
    void rebuild_filtered_list(const Frame& frame, const md::PairList& pairs);

    double cutoff_sq_;
    double s12_ = 1.0, s13_ = 1.0, s14_ = 1.0;
    bool any_excl_ = false;

    std::unordered_map<std::int64_t, std::vector<std::int64_t>>
        excluded_by_tag_;

    std::vector<int> filtered_neighbors_;
    std::vector<int> filtered_numneigh_;
    std::vector<int> filtered_firstneigh_;
    std::vector<double> filtered_shift_x_;
    std::vector<double> filtered_shift_y_;
    std::vector<double> filtered_shift_z_;
    bool filtered_has_shifts_ = false;
    int filtered_inum_ = 0;
    unsigned long long filtered_version_ = 0;

    int n_threads_ = 1;
    std::vector<std::vector<double>> fx_thr_, fy_thr_, fz_thr_;
};

} // namespace cpu
} // namespace ATV
