/// @file cpu/frame/atom.h
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
/// @see cpu/comm/halo_exchange.h, cpu/comm/atom_migration.h for the
///      components that drive nlocal/nghost
/// @see cpu/comm/atom_pack.h for pack/unpack helpers
#pragma once

#include <cstdint>
#include <unordered_map>

#include <xtensor/containers/xtensor.hpp>
#include <xtensor/generators/xbuilder.hpp>

#include "atomiverse/core/precision.h"

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
    // ── Per-atom topology strides (runtime, LAMMPS Atom::bond_per_atom) ────
    //
    // Column widths of the wide bonded arrays below. All start at 0 and grow
    // on demand as topology is loaded (`add_bond` / `add_angle` /
    // `add_dihedral` / `build_special_list`), so each system pays exactly its
    // own valence: a linear PE chain sizes to {2, 1, 1, 6} instead of a
    // compile-time worst case — ~3.5x less topology memory per atom.
    // Stride growth re-strides the whole array (O(nmax) copy) and only ever
    // happens at topology-load time; the strides are then fixed for the run.
    // MPI note: every rank loads the full global topology before partitioning
    // (add_* self-filters via `map`), so strides converge to identical values
    // on all ranks; `unpack_exchange` additionally grows defensively if a
    // migrating atom arrives with more entries than the local stride.

    /// @brief Columns of `bond_atom` / `bond_type` (max bonds owned per atom).
    int bond_per_atom = 0;

    /// @brief Max angles owned per atom; `angle_atom` has `2 *` this columns.
    int angle_per_atom = 0;

    /// @brief Max dihedrals owned per atom; `dihedral_atom` has `3 *` this.
    int dihedral_per_atom = 0;

    /// @brief Columns of `special` (max 1-2 + 1-3 + 1-4 partners per atom).
    int maxspecial = 0;

    /// @brief Grow `bond_per_atom` to at least `target` (no-op if already
    /// wide enough), re-striding `bond_atom` / `bond_type` in place.
    void grow_bond_per_atom(int target);

    /// @brief Grow `angle_per_atom` to at least `target`, re-striding
    /// `angle_atom` / `angle_type`.
    void grow_angle_per_atom(int target);

    /// @brief Grow `dihedral_per_atom` to at least `target`, re-striding
    /// `dihedral_atom` / `dihedral_type`.
    void grow_dihedral_per_atom(int target);

    /// @brief Grow `maxspecial` to at least `target`, re-striding `special`.
    void grow_maxspecial(int target);

    /// @brief Grow all four topology strides to at least the given values and
    /// mark the atom_style molecular. Use when the per-atom valence is known
    /// up front (avoids the incremental re-strides) or when populating the
    /// wide arrays directly instead of via `add_bond` / `add_angle` / ....
    void reserve_topology(int bond, int angle, int dihedral, int special_n);

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

    /// @brief Reallocate every per-atom array to exactly `new_nmax` rows,
    ///        preserving rows `[0, min(nmax, new_nmax))` by explicit copy and
    ///        zero-filling any new tail (tag tail gets the identity default).
    ///        Supports BOTH growth and shrink. Unlike `grow()` (which trusts
    ///        `xt::xtensor::resize` to keep the prefix), this copies into fresh
    ///        storage, so it is safe to shrink — used to release the global-N
    ///        backing store each rank loads before partitioning down to its
    ///        owned subset, and as the grow-on-overflow path for halo/migration
    ///        when a rank starts under-provisioned. Caller updates any external
    ///        capacity mirror (`DomainDecomposition::capacity`) and the
    ///        tag→row `map` (entries past `new_nmax` are dropped on shrink).
    /// @warning Invalidates cached `.data()` pointers (same contract as grow).
    void set_capacity(int new_nmax);

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

    /// @brief Resident heap bytes of all per-atom SoA arrays (sized to `nmax`),
    /// the lazily-allocated bonded topology arrays, and the `tag → local` map.
    /// Consumed by the built-in driver memory summary ("Atoms" category).
    std::size_t memory_bytes() const;

    // ── SoA per-atom fields (all sized to `nmax`) ──────────────────────────
    //
    // Precision tiers (core/precision.h): storage fields (positions, velocities,
    // mass, charge) are `StorageReal` — lowerable to float via
    // -DATV_REAL_STORAGE=float for memory-bound runs. Forces are `AccReal` (they
    // are `+=` accumulators — never economized, LAMMPS double-accumulation rule).
    // Kernels read StorageReal scalars and widen to double (WorkReal) for the
    // math; see .claude/specs/precision-migration.md.

    xt::xtensor<StorageReal, 1> x;      ///< x positions [Angstrom].
    xt::xtensor<StorageReal, 1> y;      ///< y positions [Angstrom].
    xt::xtensor<StorageReal, 1> z;      ///< z positions [Angstrom].
    xt::xtensor<StorageReal, 1> vx;     ///< x velocities [Angstrom/fs].
    xt::xtensor<StorageReal, 1> vy;     ///< y velocities [Angstrom/fs].
    xt::xtensor<StorageReal, 1> vz;     ///< z velocities [Angstrom/fs].
    xt::xtensor<AccReal, 1> fx;         ///< x forces [eV/Angstrom]. AccReal (accumulator).
    xt::xtensor<AccReal, 1> fy;         ///< y forces [eV/Angstrom]. AccReal (accumulator).
    xt::xtensor<AccReal, 1> fz;         ///< z forces [eV/Angstrom]. AccReal (accumulator).
    xt::xtensor<StorageReal, 1> mass;   ///< Masses [amu].
    xt::xtensor<StorageReal, 1> charge; ///< Partial charges [elementary charge].
    xt::xtensor<int, 1>    type;        ///< Atom type indices. (LocalInt, Phase 3)

    /// @brief Global atom identity, persistent across migration and restart.
    /// Default-initialized to identity (`tag(i) == i`) on `grow()`.
    xt::xtensor<TagInt, 1> tag;

    /// @brief Allocate a fresh, never-used global tag (LAMMPS `tag_extend`).
    ///
    /// Required by any path that CREATES atoms (GCMC insertion): after a
    /// delete-compaction the identity seeding of `grow()`/`set_capacity()`
    /// would hand a new row a tag some surviving atom already carries.
    /// Lazily seeds a monotonic counter at `max(live tags) + 1` on first
    /// call; deleted tags are never reused, so uniqueness holds across any
    /// interleaving of inserts and deletes. Single-rank only (no MPI
    /// allreduce over the max) — matching MC's v1 single-rank contract.
    /// @return A tag no live (or previously live) atom carries.
    TagInt acquire_tag();

    /// @brief Next tag `acquire_tag()` hands out; `< 0` = not yet seeded.
    TagInt next_tag_ = -1;

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

    /// @brief Partner global tag for each bond, shape `[nmax, bond_per_atom]`.
    /// `bond_atom(i, m)` is the global tag of the partner of atom `i`'s
    /// `m`-th bond. The partner's local row index is resolved per-call
    /// via `map[bond_atom(i, m)]`.
    xt::xtensor<TagInt, 2> bond_atom;

    /// @brief Bond type index, shape `[nmax, bond_per_atom]`.
    xt::xtensor<int, 2> bond_type;

    /// @brief Number of angles atom `i` is the center of (owner convention:
    /// center atom owns).
    xt::xtensor<int, 1> num_angle;

    /// @brief Non-center partner tags for each angle, shape
    /// `[nmax, 2 * angle_per_atom]`. For atom `i`'s `m`-th angle:
    /// `angle_atom(i, 2*m + 0)` is the i-partner tag (atom 0 of the
    /// `i-j-k` triple), `angle_atom(i, 2*m + 1)` is the k-partner tag.
    /// Atom `i` itself is the center `j`.
    xt::xtensor<TagInt, 2> angle_atom;

    /// @brief Angle type index, shape `[nmax, angle_per_atom]`.
    xt::xtensor<int, 2> angle_type;

    /// @brief Number of dihedrals atom `i` is the j (second) atom of.
    /// LAMMPS owner convention: atom-j owns for an `i-j-k-l` quadruple.
    xt::xtensor<int, 1> num_dihedral;

    /// @brief Non-j partner tags for each dihedral, shape
    /// `[nmax, 3 * dihedral_per_atom]`. For atom `i`'s `m`-th dihedral
    /// (where `i` is atom `j`): `dihedral_atom(i, 3*m + 0/1/2)` are the
    /// `i`, `k`, `l` partner tags (the dihedral's atom 0, 2, 3).
    xt::xtensor<TagInt, 2> dihedral_atom;

    /// @brief Dihedral type index, shape `[nmax, dihedral_per_atom]`.
    xt::xtensor<int, 2> dihedral_type;

    /// @brief 1-2 / 1-3 / 1-4 neighbor counts, shape `[nmax, 3]`.
    /// `nspecial(i, 0)` is the number of 1-2 (direct bond) neighbors;
    /// `nspecial(i, 1)` includes 1-3 (through angle); `nspecial(i, 2)`
    /// includes 1-4 (through dihedral). Cumulative counts.
    xt::xtensor<int, 2> nspecial;

    /// @brief Special-neighbor tag list, shape `[nmax, maxspecial]`.
    /// `special(i, 0 .. nspecial(i, 0))` are 1-2 partners; the next
    /// `nspecial(i, 1) - nspecial(i, 0)` are 1-3; and so on. Used by
    /// pair evaluators for non-bonded exclusion (`special_bonds`).
    xt::xtensor<TagInt, 2> special;

    // ── Lazy bonded storage (LAMMPS atom_style atomic vs molecular) ────────
    //
    // The wide per-bond/angle/dihedral/special partner+type arrays above
    // (bond_atom..special) are allocated ONLY for systems that actually carry
    // bonded topology, at the runtime strides above. A pure-LJ system (no
    // bonds) keeps them empty (all strides 0) and pays nothing — matching
    // LAMMPS atom_style atomic. The small COUNT arrays (num_bond/num_angle/
    // num_dihedral/nspecial) stay allocated-and-zeroed always, so every reader
    // that loops `for m in [0, num_bond[i])` is a no-op on an LJ system and
    // needs no guard; the wide arrays are touched only after a bond is added.
    bool has_bonded_topology_ = false;

    /// @brief Whether the wide bonded arrays are allocated (topology present).
    bool has_bonded() const noexcept { return has_bonded_topology_; }

    /// @brief Mark the atom_style molecular. Idempotent. Called automatically
    /// by the first `add_bond/add_angle/add_dihedral`. The wide arrays
    /// themselves are sized by the stride growers (`grow_bond_per_atom`/...
    /// or `reserve_topology`); `grow`/`set_capacity` then keep them sized.
    void ensure_bonded_storage();

    // ── Global tag → local row index map (LAMMPS Atom::map) ────────────────

    /// @brief Tag-to-local-row map. Owners-only after `Comm::exchange`;
    /// owners + ghosts after `Comm::borders`. Maintained by the comm
    /// components -- never written by force kernels.
    std::unordered_map<TagInt, int> map;

    /// @brief Empty the map. LAMMPS `Atom::map_clear` equivalent.
    void map_clear() { map.clear(); }

    /// @brief Empty the map AND free its bucket storage. Used at run
    /// setup for systems with no bonded topology: the map exists solely
    /// so the `add_bond/add_angle/add_dihedral` loaders and bonded force
    /// kernels can resolve global tags, and costs ~50 B/atom of hash-node
    /// heap when populated (e.g. ~400 MB dead weight on an 8M-atom LJ
    /// run). `map_clear()` alone keeps the bucket array alive.
    void map_release() { std::unordered_map<TagInt, int>{}.swap(map); }

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
    void add_bond(TagInt ti, TagInt tj, int type);

    /// @brief Add an angle (ti, tj, tk, type). Center atom `tj` owns the
    /// angle; the angle appears in `tj`'s `angle_atom` list with the two
    /// non-center partners (ti, tk) stored as a pair.
    void add_angle(TagInt ti, TagInt tj, TagInt tk, int type);

    /// @brief Add a dihedral (ti, tj, tk, tl, type). Atom-j (`tj`, the
    /// 2nd atom) owns the dihedral; the dihedral appears in `tj`'s
    /// `dihedral_atom` list with (ti, tk, tl) stored as a triple.
    void add_dihedral(TagInt ti, TagInt tj, TagInt tk,
                      TagInt tl, int type);

    /// @brief Build the 1-2 / 1-3 / 1-4 `nspecial` + `special` exclusion
    /// list by BFS over the per-atom bond graph.
    void build_special_list();

    /// @brief Copy one full per-atom row — every SoA field plus the per-atom
    ///        bonded topology — from row `src` into row `dst` (no-op when
    ///        `dst == src`).
    ///
    /// The single row-move primitive: MPI migration compacts survivors with
    /// it, and Monte Carlo insert/delete reuses it (a GCMC deletion is
    /// `copy_row(i, nlocal - 1); --nlocal;`). Does NOT touch `nlocal` or the
    /// tag→row `map` — callers own the population count and rebuild the map.
    /// @param dst Destination row index (must be `< nmax`).
    /// @param src Source row index (must be `< nmax`).
    void copy_row(int dst, int src);

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
