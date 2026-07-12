/// @file cpu/potentials/pair_lj_cut.h
/// @brief LAMMPS-style PairLJCut — 12-6 with shifted potential (CPU).
///
/// Was `ATV::cpu::md::LJ` before the LAMMPS-naming refactor.
#pragma once

#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

#include "atomiverse/cpu/neighbor/pair_list.h"
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
/// Lists are derived from `frame.atoms (per-atom bonded topology)`/`frame.atoms (per-atom angle topology)`/`frame.atoms (per-atom dihedral topology)`
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

    /// @brief Per-type mode: the coefficients come from `set_type_coeffs`.
    ///
    /// Use this constructor when the system has more than one element. It takes NO
    /// epsilon/sigma, because there is no such thing as a system-wide one — supplying
    /// dummy values that the type table then overrides would be a fallback value, and
    /// this engine does not have those. `build()` throws if no coefficients were set.
    ///
    /// @code
    ///   // TIP3P: the LJ site is on oxygen only.
    ///   sim.set_pair(cpu::PairLJCut{cutoff}
    ///                    .set_type_coeffs(8, 8, eps_OO, sigma_OO));
    ///   // Hydrogen has no diagonal, so (H,H) and (O,H) mix to eps = 0 and never
    ///   // interact — which is exactly right, and is checked, not assumed.
    /// @endcode
    ///
    /// @param cutoff Interaction cutoff [Angstrom].
    explicit PairLJCut(double cutoff);

    /// @brief Set the LJ coefficients for one type pair (LAMMPS `pair_coeff i j`).
    ///
    /// Without this, `PairLJCut` applies ONE epsilon/sigma to every atom in the system.
    /// That is fine for argon and useless for anything with more than one element: in
    /// TIP3P the Lennard-Jones site sits on oxygen only, so a single-sigma run forces
    /// the hydrogens to wear the oxygen's sigma = 3.15 A, and at liquid density
    /// neighbouring waters overlap so violently that no constraint solver can hold a
    /// molecule together. (That is not a hypothetical either — it is why
    /// regressions/md/rigid_water has to run dilute.)
    ///
    /// Types are the frame's `atoms.type[i]` values. Calling this at all switches the
    /// component to the per-type table; the off-diagonal entries you do not set are
    /// filled by **Lorentz-Berthelot** mixing from the diagonal:
    ///
    ///     sigma_ij = (sigma_ii + sigma_jj) / 2       eps_ij = sqrt(eps_ii * eps_jj)
    ///
    /// so setting only the like-pairs (i,i) is usually enough. An explicit (i,j) call
    /// overrides the mixed value.
    ///
    /// **Performance note.** The single-epsilon path keeps the hoisted xsimd ket loop —
    /// that is the kernel benchmarked against LAMMPS, and it is untouched. A per-type
    /// table cannot hoist its coefficients out of the inner loop, so the multi-type
    /// path is SCALAR. It is correct and it is not fast; vectorizing it needs a SIMD
    /// gather over the coefficient table, which nobody has needed yet.
    ///
    /// @param type_i First atom type.
    /// @param type_j Second atom type.
    /// @param epsilon Well depth [eV].
    /// @param sigma Length scale [Angstrom].
    PairLJCut& set_type_coeffs(int type_i, int type_j, double epsilon, double sigma) &;
    PairLJCut set_type_coeffs(int type_i, int type_j, double epsilon,
                              double sigma) &&;

    /// @brief True when a per-type coefficient table is active (scalar path).
    bool is_multi_type() const { return n_types_ > 0; }

  private:
    /// @brief Complete the table by Lorentz-Berthelot mixing and precompute
    ///        sigma^6 + the per-pair energy shift. Called once from build().
    void finalize_type_table_();

  public:

    std::string name() const override { return "lj"; }

    /// @brief Enable/disable LAMMPS-style "special bonds" exclusions.
    /// See pre-refactor `cpu::PairLJCut::set_special_bonds` for full
    /// semantics — unchanged by the rename.
    PairLJCut& set_special_bonds(double s12, double s13, double s14) &;
    PairLJCut set_special_bonds(double s12, double s13, double s14) &&;

    /// @brief Enable/disable per-pair virial (stress) accumulation.
    ///
    /// The pair virial (6 Voigt components) feeds the pressure / barostat
    /// path only; a plain NVE run never consumes it. Computing it costs
    /// ~20-25% of the inner-loop arithmetic (6 extra batch products + a
    /// horizontal reduce per atom). Disabling it when no stress consumer is
    /// active mirrors LAMMPS's per-step `vflag` gating. Forces and energy
    /// are unaffected — bit-identical either way. Default: enabled.
    PairLJCut& set_compute_stress(bool on) &;
    PairLJCut set_compute_stress(bool on) &&;

    void build(const BuildContext& ctx) override;
    void destroy(const DestroyContext& ctx) override;

    /// @brief Resident heap bytes: per-thread force accumulators (OMP>1 only)
    /// plus the exclusion map + filtered neighbor list (special-bonds runs).
    /// Classified under the "Pair" category by the driver memory summary.
    std::size_t memory_bytes() const override;

  protected:
    void compute(const LaunchContext& ctx) override;

  private:
    /// @brief Rebuild `filtered_{numneigh,firstneigh,neighbors}` from
    /// `pairs` by walking each owner i's neighbors, resolving partner
    /// tags, and dropping any pair (i,j) whose j_tag appears in i's
    /// excluded set. Called from `compute()` only when `pairs.version`
    /// has advanced past the last snapshot we filtered.
    void rebuild_filtered_list(const Frame& frame, const md::PairList& pairs);

    /// @brief (Re)build `excluded_by_tag_` from the CURRENT per-atom
    /// `nspecial`/`special` topology of the first `n_owners` rows.
    ///
    /// Must run whenever the owned-atom set changes (multi-rank
    /// `Migrate`): the per-atom special lists travel WITH the atom
    /// (pack/unpack/move_row), but a map cached once at `build()` goes
    /// stale the moment an atom migrates into a rank whose build-time map
    /// never held that tag — its 1-2/1-3/1-4 partners then escape the
    /// exclusion filter and the LJ kernel sees a ~1.5 Å bonded pair,
    /// producing a catastrophic repulsive spike (MSU-PE multi-rank NPT
    /// collapse). Called from `rebuild_filtered_list`, which is itself
    /// gated on `pairs.version` so the refresh tracks every NL rebuild.
    void refresh_exclusion_map(const Frame& frame, int n_owners);

    double epsilon_, sigma_;

    /// Per-type coefficient tables, flattened `n_types_ x n_types_`. Empty = the
    /// single-epsilon fast path. Filled lazily by `set_type_coeffs`, then completed by
    /// Lorentz-Berthelot mixing at `build()`.
    bool per_type_mode_ = false; ///< Constructed with the cutoff-only ctor.
    int n_types_ = 0;
    std::vector<double> eps_tab_;   ///< eps_ij [eV]
    std::vector<double> sigma_tab_; ///< sigma_ij [A]
    std::vector<double> sig6_tab_;  ///< sigma_ij^6 [A^6]
    std::vector<double> shift_tab_; ///< energy shift at the cutoff, per pair [eV]
    std::vector<char> set_tab_;     ///< 1 = explicitly set by the user
    double sigma6_, cutoff_sq_, v_shift_;

    // Special-bonds: when any of these is < 1, the corresponding bonded
    // pair list is built at build() time and consulted at launch() time.
    double s12_ = 1.0, s13_ = 1.0, s14_ = 1.0;

    // When false, the inner loop skips the 6-component pair virial and the
    // final stress write (LAMMPS-style vflag gating). Forces/energy are
    // unchanged. Default true to preserve behavior for pressure/NPT runs.
    bool compute_stress_ = true;

    // True when any special-bonds tier is fully excluded (scale 0). Gates
    // the exclusion-filtered CSR path independently of whether the current
    // exclusion map happens to be non-empty, so the filter (and its
    // per-rebuild map refresh) stays active across migrations even on a
    // step where this rank momentarily owns no excluded pair.
    bool any_excl_ = false;

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
    // Per-pair PBC shifts parallel to filtered_neighbors_ (copied from the
    // upstream PairList when it carries shifts) so the kernel can skip the
    // per-pair min-image. Empty when the upstream list has no shifts.
    std::vector<double> filtered_shift_x_;
    std::vector<double> filtered_shift_y_;
    std::vector<double> filtered_shift_z_;
    bool filtered_has_shifts_ = false;
    int filtered_inum_ = 0;
    unsigned long long filtered_version_ = 0;
};

} // namespace cpu
} // namespace ATV
