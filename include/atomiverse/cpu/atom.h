/// @file cpu/atom.h
/// @brief cpu::Atom -- LAMMPS-shaped per-atom storage with capacity-on-demand growth.
///
/// `cpu::Atom` is the CPU backend's canonical per-atom container. It owns:
///   - SoA fields for positions, velocities, forces, mass, charge, type, tag, image
///   - Per-atom bonded topology (`num_bond[i]`, `bond_atom[i][m]`, ...)
///   - `nlocal` / `nghost` population counts and `nmax` allocation capacity
///   - The `tag → local index` map (`Atom::map`), populated by the comm
///     components after `borders()`
///
/// **Naming follows LAMMPS Atom.** `nlocal`, `nghost`, `nmax`, `map`,
/// `num_bond`, `bond_atom`, `bond_type`, `nspecial`, `special` -- all
/// lowercase names match `lammps/src/atom.h` so cross-reference is
/// 1:1 (modulo C-style `double **x` vs xtensor `xt::xtensor<double, 1> x`).
///
/// **Growth semantics.** Capacity (`nmax`) is dynamic. The `grow()` method
/// doubles `nmax` (LAMMPS DELTA growth) and resizes every per-atom array.
/// `grow()` is only ever called from `cpu::comm::Exchange::launch` or
/// `cpu::comm::Borders::launch` -- both at the `pre_force` phase -- so
/// no force kernel observes an invalidated pointer mid-step. Cached
/// `.data()` pointers in force-evaluator inner loops are safe for the
/// duration of a single `launch()` call.
///
/// @see lammps/src/atom.h for the upstream model
/// @see cpu/comm/comm.h for the components that drive nlocal/nghost
/// @see cpu/comm/atom_pack.h for pack/unpack helpers
#pragma once

#include <cstdint>
#include <unordered_map>

#include <xtensor/containers/xtensor.hpp>
#include <xtensor/generators/xbuilder.hpp>

namespace ATV {
namespace cpu {

/// @brief CPU backend per-atom container with dynamic capacity and LAMMPS-style
///        per-atom topology.
///
/// Replaces the deleted `cpu::Frame::Atoms` + `cpu::Frame::Bonds/Angles/Dihedrals`
/// sub-structs. Hot-path discipline: stores only fundamental physical
/// quantities + the `tag → local` map (which is intrinsically tied to the
/// per-atom array indexing and has nowhere else sensible to live).
class Atom {
  public:
    // ── Per-atom topology bounds (compile-time) ────────────────────────────

    /// @brief Maximum bonds per atom. PE/Dreiding: 2; SAFE/CL polymers: 4.
    /// A runtime assertion in the materialization scatter loop catches
    /// over-valent atoms before they corrupt the per-atom storage at
    /// migration time.
    static constexpr int kMaxBondPerAtom = 6;

    /// @brief Maximum angles per atom (center-atom-owned).
    static constexpr int kMaxAnglePerAtom = 6;

    /// @brief Maximum dihedrals per atom (atom-j-owned, i.e. atom 2 of 4).
    static constexpr int kMaxDihedralPerAtom = 8;

    /// @brief Maximum 1-2 + 1-3 + 1-4 special-neighbor tags per atom.
    /// Used by pair evaluators for non-bonded exclusion.
    static constexpr int kMaxSpecial = 16;

    // ── Construction ───────────────────────────────────────────────────────

    /// @brief Construct an Atom with a starting capacity.
    /// @param initial_nmax Initial backing-store capacity. Per-atom arrays
    ///        are allocated and zero-initialized to this size. Pass the
    ///        global atom count (or a sufficient upper bound) at
    ///        materialization time; subsequent ranks will shrink their
    ///        owned count after the first `Comm::exchange`.
    explicit Atom(int initial_nmax = 0);

    /// @brief Grow capacity to at least `target_nmax`, doubling from the
    ///        current `nmax` (LAMMPS DELTA growth). No-op when
    ///        `target_nmax <= nmax`. Every per-atom array is resized in
    ///        place; trailing rows are zero-initialized.
    ///
    /// @warning Invalidates any cached `.data()` pointer into a per-atom
    ///          array. Only called from `pre_force` phase (inside
    ///          `cpu::comm::Exchange::launch` or `Borders::launch`); force
    ///          kernels at the `force` phase never observe a stale pointer
    ///          mid-`launch()`.
    void grow(int target_nmax);

    // ── Population counts (LAMMPS Atom::nlocal/nghost/nmax) ────────────────

    /// @brief Number of locally owned atoms. Owners occupy `[0, nlocal)`.
    int nlocal = 0;

    /// @brief Number of ghost atoms received from neighboring ranks.
    /// Ghosts occupy `[nlocal, nlocal + nghost)`.
    int nghost = 0;

    /// @brief Allocated capacity. Every per-atom array has size `nmax`.
    /// `nlocal + nghost <= nmax` is the invariant maintained by `grow()`.
    int nmax = 0;

    /// @brief Convenience accessor for total live row count.
    /// @return `nlocal + nghost`.
    int n_total() const noexcept { return nlocal + nghost; }

    /// @brief Owner count (alias for `nlocal`). Lets `for (int i=0; i<atoms.size(); ++i)`
    /// continue to work; reading owners-only is the dominant idiom across
    /// the CPU codebase. LAMMPS reads `atom->nlocal` directly; this is the
    /// C++ container-API analogue.
    int size() const noexcept { return nlocal; }

    // ── SoA per-atom fields (all sized to `nmax`) ──────────────────────────

    xt::xtensor<double, 1> x;       ///< x positions [Angstrom].
    xt::xtensor<double, 1> y;       ///< y positions [Angstrom].
    xt::xtensor<double, 1> z;       ///< z positions [Angstrom].
    xt::xtensor<double, 1> vx;      ///< x velocities [Angstrom/fs].
    xt::xtensor<double, 1> vy;      ///< y velocities [Angstrom/fs].
    xt::xtensor<double, 1> vz;      ///< z velocities [Angstrom/fs].
    xt::xtensor<double, 1> fx;      ///< x forces [eV/Angstrom].
    xt::xtensor<double, 1> fy;      ///< y forces [eV/Angstrom].
    xt::xtensor<double, 1> fz;      ///< z forces [eV/Angstrom].
    xt::xtensor<double, 1> mass;    ///< Masses [amu].
    xt::xtensor<double, 1> charge;  ///< Partial charges [elementary charge].
    xt::xtensor<int, 1>    type;    ///< Atom type indices.

    /// @brief Global atom identity, persistent across migration and restart.
    /// Default-initialized to identity (`tag(i) == i`) on `grow()`.
    xt::xtensor<std::int64_t, 1> tag;

    /// @brief Per-atom PBC crossing counts, shape `[nmax, 3]`.
    /// `image(i, d) += ±1` each time atom `i` crosses the box face along
    /// dimension `d`. Used to reconstruct unwrapped trajectories (MSD,
    /// diffusion).
    xt::xtensor<std::int8_t, 2> image;

    // ── Per-atom bonded topology (LAMMPS Atom::num_bond / bond_atom / ...) ─

    /// @brief Number of bonds atom `i` participates in (as the owner).
    /// LAMMPS owner convention with Newton bond off: owner is the atom
    /// with the smaller tag.
    xt::xtensor<int, 1> num_bond;

    /// @brief Partner global tag for each bond, shape `[nmax, kMaxBondPerAtom]`.
    /// `bond_atom(i, m)` is the global tag of the partner of atom `i`'s
    /// `m`-th bond. The partner's local row index is resolved per-call
    /// via `map[bond_atom(i, m)]`.
    xt::xtensor<std::int64_t, 2> bond_atom;

    /// @brief Bond type index, shape `[nmax, kMaxBondPerAtom]`.
    xt::xtensor<int, 2> bond_type;

    /// @brief Number of angles atom `i` is the center of (owner convention:
    /// center atom owns).
    xt::xtensor<int, 1> num_angle;

    /// @brief Non-center partner tags for each angle, shape
    /// `[nmax, 2 * kMaxAnglePerAtom]`. For atom `i`'s `m`-th angle:
    /// `angle_atom(i, 2*m + 0)` is the i-partner tag (atom 0 of the
    /// `i-j-k` triple), `angle_atom(i, 2*m + 1)` is the k-partner tag.
    /// Atom `i` itself is the center `j`.
    xt::xtensor<std::int64_t, 2> angle_atom;

    /// @brief Angle type index, shape `[nmax, kMaxAnglePerAtom]`.
    xt::xtensor<int, 2> angle_type;

    /// @brief Number of dihedrals atom `i` is the j (second) atom of.
    /// LAMMPS owner convention: atom-j owns for an `i-j-k-l` quadruple.
    xt::xtensor<int, 1> num_dihedral;

    /// @brief Non-j partner tags for each dihedral, shape
    /// `[nmax, 3 * kMaxDihedralPerAtom]`. For atom `i`'s `m`-th dihedral
    /// (where `i` is atom `j`): `dihedral_atom(i, 3*m + 0/1/2)` are the
    /// `i`, `k`, `l` partner tags (the dihedral's atom 0, 2, 3).
    xt::xtensor<std::int64_t, 2> dihedral_atom;

    /// @brief Dihedral type index, shape `[nmax, kMaxDihedralPerAtom]`.
    xt::xtensor<int, 2> dihedral_type;

    /// @brief 1-2 / 1-3 / 1-4 neighbor counts, shape `[nmax, 3]`.
    /// `nspecial(i, 0)` is the number of 1-2 (direct bond) neighbors;
    /// `nspecial(i, 1)` includes 1-3 (through angle); `nspecial(i, 2)`
    /// includes 1-4 (through dihedral). Cumulative counts.
    xt::xtensor<int, 2> nspecial;

    /// @brief Special-neighbor tag list, shape `[nmax, kMaxSpecial]`.
    /// `special(i, 0 .. nspecial(i, 0))` are 1-2 partners; the next
    /// `nspecial(i, 1) - nspecial(i, 0)` are 1-3; and so on. Used by
    /// pair evaluators for non-bonded exclusion (`special_bonds`).
    xt::xtensor<std::int64_t, 2> special;

    // ── Global tag → local row index map (LAMMPS Atom::map) ────────────────

    /// @brief Tag-to-local-row map. Owners-only after `Comm::exchange`;
    /// owners + ghosts after `Comm::borders`. Maintained by the comm
    /// components -- never written by force kernels.
    std::unordered_map<std::int64_t, int> map;

    /// @brief Empty the map. LAMMPS `Atom::map_clear` equivalent.
    void map_clear() { map.clear(); }

    /// @brief Rebuild the map for rows `[begin, end)` from current
    /// `tag` values. LAMMPS `Atom::map_set` equivalent. Inserts (not
    /// replaces) -- callers typically call `map_clear()` first.
    /// @param begin First row to insert.
    /// @param end Past-the-end row.
    void map_set(int begin, int end);

    // ── Per-atom topology population (LAMMPS Atom::data_atoms equivalents) ─
    //
    // Used by `cpu::materialize` and the LAMMPS-data reader to scatter a
    // flat (global) topology specification into the per-atom layout
    // before any MPI partition. Tags `ti/tj/tk/tl` are global atom tags;
    // each helper looks the tag up in `map` to find the local row and
    // appends to that row's bond / angle / dihedral list.
    //
    // **Pre-requisites for the caller:**
    //   - `nlocal` is set to the total atom count
    //   - `tag[i]` is set for `i in [0, nlocal)`
    //   - `map_set(0, nlocal)` has been called so tag → local works

    /// @brief Add a bond (ti, tj, type). LAMMPS Newton-bond-ON ownership:
    /// the smaller-tag atom owns the bond; the bond appears only in that
    /// atom's `bond_atom` list. Partner force is gathered via reverse_comm
    /// when the partner is a ghost.
    void add_bond(std::int64_t ti, std::int64_t tj, int type);

    /// @brief Add an angle (ti, tj, tk, type). Center atom `tj` owns the
    /// angle; the angle appears in `tj`'s `angle_atom` list with the two
    /// non-center partners (ti, tk) stored as a pair.
    void add_angle(std::int64_t ti, std::int64_t tj, std::int64_t tk, int type);

    /// @brief Add a dihedral (ti, tj, tk, tl, type). Atom-j (`tj`, the
    /// 2nd atom) owns the dihedral; the dihedral appears in `tj`'s
    /// `dihedral_atom` list with (ti, tk, tl) stored as a triple.
    void add_dihedral(std::int64_t ti, std::int64_t tj, std::int64_t tk,
                      std::int64_t tl, int type);

    /// @brief Build the 1-2 / 1-3 / 1-4 `nspecial` + `special` exclusion
    /// list by BFS over the per-atom bond graph.
    void build_special_list();

    /// @brief Reorder per-atom topology rows to match the current `tag`
    ///        permutation.
    ///
    /// Use case: a checkpoint loader (Zarr) overwrites `tag[i]` with a
    /// permutation of `0..n-1` without touching the per-atom topology
    /// arrays. After the load, row `i` holds atom-`tag[i]`'s state but
    /// `num_bond[i] / bond_atom[i][m] / ...` still reflect what was
    /// originally at row `i` (probably what was at atom-`i` if the
    /// materializer assumed identity tags). This helper moves each
    /// row's topology to wherever its tag now lives, so row `i`'s
    /// topology again describes atom-`tag[i]`.
    ///
    /// Assumes `tag[0..nlocal)` is a permutation of `[0, nlocal)` —
    /// i.e. tags are dense, 0-indexed, no gaps. Throws otherwise.
    void reshuffle_topology_to_tags();
};

} // namespace cpu
} // namespace ATV
