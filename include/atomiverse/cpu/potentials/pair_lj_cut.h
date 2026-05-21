/// @file cpu/potentials/pair_lj_cut.h
/// @brief LAMMPS-style PairLJCut — 12-6 with shifted potential (CPU).
///
/// Was `ATV::cpu::md::LJ` before the LAMMPS-naming refactor. The legacy
/// header `cpu/md/lj.h` provides `using LJ = PairLJCut;` so user code
/// (`sim.set_pair(cpu::md::LJ{...})`) compiles unchanged.
#pragma once

#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

#include "atomiverse/cpu/md/pair_list.h"
#include "atomiverse/cpu/potentials/pair.h"

namespace ATV {
namespace cpu {

/// @brief Lennard-Jones 12-6 with shifted potential.
///
/// V(r) = 4ε[(σ/r)¹² − (σ/r)⁶] − V_shift
/// F(r) = (24ε/r)·[2·(σ/r)¹² − (σ/r)⁶]·r̂
///
/// Optional special-bonds exclusions (LAMMPS convention): pairs that are
/// 1-2 (bonded), 1-3 (across an angle vertex) or 1-4 (across a dihedral)
/// are skipped during the pair loop. Enable via `set_special_bonds()`.
/// Lists are derived from `frame.bonds`/`frame.angles`/`frame.dihedrals`
/// at `build()` time and consulted with a per-atom binary search at
/// `compute()` time. Default: no exclusions.
///
/// Thread-parallel with per-thread force reduction buffers
/// (Newton's 3rd law). xsimd ket-loop in the inner pair walk (Track B.2,
/// 1.21× over scalar on NEON64).
class PairLJCut : public Pair {
  public:
    /// @param epsilon Well depth [eV].
    /// @param sigma   Size parameter [Å].
    /// @param cutoff  Cutoff distance [Å].
    PairLJCut(double epsilon, double sigma, double cutoff);

    std::string name() const override { return "lj"; }

    /// @brief Enable/disable LAMMPS-style "special bonds" exclusions.
    /// See pre-refactor `cpu::md::LJ::set_special_bonds` for full
    /// semantics — unchanged by the rename.
    PairLJCut& set_special_bonds(double s12, double s13, double s14) &;
    PairLJCut set_special_bonds(double s12, double s13, double s14) &&;

    void build(const BuildContext& ctx) override;
    void destroy(const DestroyContext& ctx) override;

  protected:
    void compute(const LaunchContext& ctx) override;

  private:
    /// @brief Rebuild `filtered_{numneigh,firstneigh,neighbors}` from
    /// `pairs` by walking each owner i's neighbors, resolving partner
    /// tags, and dropping any pair (i,j) whose j_tag appears in i's
    /// excluded set. Called from `compute()` only when `pairs.version`
    /// has advanced past the last snapshot we filtered.
    void rebuild_filtered_list(const Frame& frame, const md::PairList& pairs);

    double epsilon_, sigma_;
    double sigma6_, cutoff_sq_, v_shift_;

    // Special-bonds: when any of these is < 1, the corresponding bonded
    // pair list is built at build() time and consulted at launch() time.
    double s12_ = 1.0, s13_ = 1.0, s14_ = 1.0;

    // Tag-indexed exclusion lists. See pre-refactor LJ for full docs;
    // unchanged by the rename.
    std::unordered_map<std::int64_t, std::vector<std::int64_t>>
        excluded_by_tag_;

    // Per-thread force reduction buffers
    int n_threads_ = 1;
    std::vector<std::vector<double>> fx_thr_, fy_thr_, fz_thr_;

    // Exclusion-filtered neighbor CSR — see pre-refactor LJ for full docs.
    std::vector<int> filtered_neighbors_;
    std::vector<int> filtered_numneigh_;
    std::vector<int> filtered_firstneigh_;
    int filtered_inum_ = 0;
    unsigned long long filtered_version_ = 0;
};

} // namespace cpu
} // namespace ATV
