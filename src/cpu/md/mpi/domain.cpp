/// @file src/cpu/md/mpi/domain.cpp
/// @brief `DomainDecomposition` implementation: single-rank `bind`, plus
/// multi-rank `init_topology`, `compute_sublo_subhi`, and `partition_atoms`
/// (Phase 2).

#include "atomiverse/cpu/md/mpi/domain.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <unordered_set>

#include "atomiverse/core/mat3.h"
#include "atomiverse/cpu/frame.h"

namespace ATV {
namespace cpu {
namespace md {
namespace mpi {

void DomainDecomposition::bind(const Frame& frame) {
    const int n = frame.atoms.size();
    n_local = n;
    n_ghost = 0;
    capacity = n;

    tag_to_local.clear();
    tag_to_local.reserve(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        tag_to_local.emplace(frame.atoms.tag(i), i);
    }
}

void DomainDecomposition::compute_sublo_subhi(const ATV::cpu::Box& box) {
    for (int d = 0; d < 3; ++d) {
        const double L = box.length(d);
        if (static_cast<int>(cuts[d].size()) == dims[d] + 1) {
            // Non-uniform cut planes (set by `RCBBalancer`). cuts[d] is in
            // box-length fractions, so the partition survives NPT/Deform box
            // rescaling without recomputation.
            sublo[d] = cuts[d][static_cast<std::size_t>(coords[d])] * L;
            subhi[d] = cuts[d][static_cast<std::size_t>(coords[d] + 1)] * L;
        } else {
            const double slab = L / static_cast<double>(dims[d]);
            sublo[d] = slab * static_cast<double>(coords[d]);
            subhi[d] = sublo[d] + slab;
        }
    }
}

#ifdef ATV_HAS_MPI

void DomainDecomposition::release_comm(bool /*mpi_is_finalizing*/) {
    // We do not own a Cartesian sub-comm (we use MPI_COMM_WORLD directly --
    // see init_topology). Just clear the handle.
    cart_comm = MPI_COMM_NULL;
}

void DomainDecomposition::init_topology(MPI_Comm world) {
    int world_size = 0;
    int my_rank = 0;
    MPI_Comm_size(world, &world_size);
    MPI_Comm_rank(world, &my_rank);

    if (world_size <= 1) {
        // Single rank: keep defaults (dims={1,1,1}, no cart_comm).
        return;
    }

    // Do NOT use MPI_Cart_create -- on OpenMPI 5.0.x (macOS) the resulting
    // sub-communicator triggers a segfault inside MPI_Finalize when any
    // user MPI_Comm_free has been issued or when the same process repeatedly
    // creates Cartesian comms (e.g. across GoogleTest cases). Compute the
    // brick decomposition directly on top of the world communicator using
    // simple rank arithmetic; this is functionally equivalent and avoids
    // the OpenMPI bug entirely.
    cart_comm = world;

    // `MPI_Dims_create` is safe -- it's a pure algorithm with no state.
    // Optional env-var override `ATV_FORCE_DIMS=Nx,Ny,Nz` skips the
    // automatic 3D decomposition and uses the given factorization. This
    // is a diagnostic/workaround: the face-only halo exchange misses
    // edge/corner ghosts in true 3D decompositions, so for
    // bonded-heavy systems (e.g. PE chains in 2x2x2) a 1D `Nx,1,1`
    // override produces correct multi-rank behavior at the cost of
    // less-uniform slab geometry. See
    // `.claude/notes/cpu_mpi_halo_edge_corner_bug_2026-05-27.md`.
    dims[0] = dims[1] = dims[2] = 0;
    if (const char* env = std::getenv("ATV_FORCE_DIMS")) {
        int dx, dy, dz;
        if (std::sscanf(env, "%d,%d,%d", &dx, &dy, &dz) == 3
            && dx * dy * dz == world_size) {
            dims[0] = dx;
            dims[1] = dy;
            dims[2] = dz;
        }
    }
    if (dims[0] == 0) {
        MPI_Dims_create(world_size, 3, dims);
    }

    // Convention: rank index = ix * (dy*dz) + iy * dz + iz.
    coords[0] = my_rank / (dims[1] * dims[2]);
    coords[1] = (my_rank / dims[2]) % dims[1];
    coords[2] = my_rank % dims[2];

    auto rank_at = [&](int cx, int cy, int cz) {
        // Periodic wrap on all axes.
        cx = ((cx % dims[0]) + dims[0]) % dims[0];
        cy = ((cy % dims[1]) + dims[1]) % dims[1];
        cz = ((cz % dims[2]) + dims[2]) % dims[2];
        return cx * (dims[1] * dims[2]) + cy * dims[2] + cz;
    };

    // Face order: -x, +x, -y, +y, -z, +z.
    face_neighbors[0] = rank_at(coords[0] - 1, coords[1], coords[2]);
    face_neighbors[1] = rank_at(coords[0] + 1, coords[1], coords[2]);
    face_neighbors[2] = rank_at(coords[0], coords[1] - 1, coords[2]);
    face_neighbors[3] = rank_at(coords[0], coords[1] + 1, coords[2]);
    face_neighbors[4] = rank_at(coords[0], coords[1], coords[2] - 1);
    face_neighbors[5] = rank_at(coords[0], coords[1], coords[2] + 1);
}

namespace {

// Orthorhombic PBC wrap: bring `p` into [0, L).
inline double wrap_to_cell(double p, double L) {
    if (L <= 0.0)
        return p;
    double q = p - std::floor(p / L) * L;
    if (q < 0.0)
        q += L;
    if (q >= L)
        q -= L;
    return q;
}

// Rank index (along one dim) that owns wrapped position `q` with `dims_d`
// equal-width slabs over `[0, L)`. Clamps to [0, dims_d - 1] against FP
// rounding at the upper edge.
inline int rank_along_dim(double q, double L, int dims_d) {
    if (dims_d <= 1)
        return 0;
    const int r = static_cast<int>(q * static_cast<double>(dims_d) / L);
    return std::max(0, std::min(dims_d - 1, r));
}

} // namespace

void DomainDecomposition::partition_atoms(ATV::cpu::Frame& frame) {
    const int world_size = dims[0] * dims[1] * dims[2];
    if (world_size <= 1 || cart_comm == MPI_COMM_NULL) {
        // Single-rank fallback: identical to `bind`.
        bind(frame);
        return;
    }

    compute_sublo_subhi(frame.box);

    const double Lx = frame.box.length(0);
    const double Ly = frame.box.length(1);
    const double Lz = frame.box.length(2);

    // First-run vs subsequent-run scoping. On the first call (n_local == 0)
    // we scan the entire user-supplied atom array because we don't yet know
    // which atoms belong to this rank. On subsequent calls -- e.g. when a
    // user invokes `Driver::run()` twice in a row -- the array's tail rows
    // hold stale ghost positions and old owner copies left over by halo and
    // by the previous compaction. Re-scanning those would let the same
    // physical atom be picked again as an extra owner (its stale FCC
    // position is unchanged, its post-step compacted copy is at a different
    // slot), and the resulting two near-identical "owner" rows blow up LJ.
    // Phase 2 has no atom migration, so previously-owned atoms still belong
    // to this rank; we only need to re-scan that range.
    const bool first_partition = (n_local == 0);
    const int n_to_iterate = first_partition ? frame.atoms.size() : n_local;
    std::vector<int> owned_indices;
    owned_indices.reserve(static_cast<std::size_t>(n_to_iterate));

    // First-partition deduplication. When the caller passed an
    // already-partitioned frame (e.g. polyethylene_tension's per-stage
    // driver rebuild: stage 1 partitioned the frame into [0..n_local) +
    // stale tail; stage 2's fresh DD enters here with n_local==0 so the
    // tail rows are also scanned), the same tag can appear at both the
    // compacted slot and the stale slot. Without dedup both copies are
    // claimed, LJ then finds two atoms at identical positions, and the
    // run blows up at stage 2 step 0. Tag-keyed `first-wins` keeps the
    // live (compacted) copy and drops the stale one — physics is
    // identical to a clean global-frame first-partition.
    std::unordered_set<std::int64_t> seen_tags;
    if (first_partition)
        seen_tags.reserve(static_cast<std::size_t>(n_to_iterate));

    for (int i = 0; i < n_to_iterate; ++i) {
        if (first_partition) {
            const std::int64_t t = frame.atoms.tag(i);
            if (!seen_tags.insert(t).second)
                continue; // duplicate tag from stale tail — keep first
        }
        const double xw = wrap_to_cell(frame.atoms.x(i), Lx);
        const double yw = wrap_to_cell(frame.atoms.y(i), Ly);
        const double zw = wrap_to_cell(frame.atoms.z(i), Lz);
        const int cx = rank_along_dim(xw, Lx, dims[0]);
        const int cy = rank_along_dim(yw, Ly, dims[1]);
        const int cz = rank_along_dim(zw, Lz, dims[2]);
        if (cx == coords[0] && cy == coords[1] && cz == coords[2]) {
            owned_indices.push_back(i);
        }
    }

    n_local = static_cast<int>(owned_indices.size());
    n_ghost = 0;
    capacity =
        frame.atoms.size(); // keep the whole allocation as slack for ghosts

    // Compact owned atoms to the front of every per-atom array. Iterating
    // in increasing `new_i` guarantees we read `old_i >= new_i` so the
    // destination has not yet been overwritten.
    auto compact_1d_d = [&](xt::xtensor<double, 1>& arr) {
        for (int ni = 0; ni < n_local; ++ni) {
            const int oi = owned_indices[static_cast<std::size_t>(ni)];
            if (oi != ni)
                arr(ni) = arr(oi);
        }
    };
    auto compact_1d_i = [&](xt::xtensor<int, 1>& arr) {
        for (int ni = 0; ni < n_local; ++ni) {
            const int oi = owned_indices[static_cast<std::size_t>(ni)];
            if (oi != ni)
                arr(ni) = arr(oi);
        }
    };
    auto compact_1d_i64 = [&](xt::xtensor<std::int64_t, 1>& arr) {
        for (int ni = 0; ni < n_local; ++ni) {
            const int oi = owned_indices[static_cast<std::size_t>(ni)];
            if (oi != ni)
                arr(ni) = arr(oi);
        }
    };

    compact_1d_d(frame.atoms.x);
    compact_1d_d(frame.atoms.y);
    compact_1d_d(frame.atoms.z);
    compact_1d_d(frame.atoms.vx);
    compact_1d_d(frame.atoms.vy);
    compact_1d_d(frame.atoms.vz);
    compact_1d_d(frame.atoms.fx);
    compact_1d_d(frame.atoms.fy);
    compact_1d_d(frame.atoms.fz);
    compact_1d_d(frame.atoms.mass);
    compact_1d_i(frame.atoms.type);
    compact_1d_i64(frame.atoms.tag);
    for (int ni = 0; ni < n_local; ++ni) {
        const int oi = owned_indices[static_cast<std::size_t>(ni)];
        if (oi != ni) {
            frame.atoms.image(ni, 0) = frame.atoms.image(oi, 0);
            frame.atoms.image(ni, 1) = frame.atoms.image(oi, 1);
            frame.atoms.image(ni, 2) = frame.atoms.image(oi, 2);
        }
    }

    // Per-atom topology must compact along with the SoA fields, else
    // num_bond[i] / bond_atom[i][m] become disconnected from tag[i].
    // (Bug surfaced as bonded_extension = ~60 Å instead of ~3 Å — the
    // resolver was looking up the right partner tag but reading the
    // wrong owner row's bond list.)
    auto& a = frame.atoms;

    for (int ni = 0; ni < n_local; ++ni) {
        const int oi = owned_indices[static_cast<std::size_t>(ni)];
        if (oi == ni)
            continue;
        a.num_bond(ni) = a.num_bond(oi);
        for (int m = 0; m < cpu::Atom::kMaxBondPerAtom; ++m) {
            a.bond_atom(ni, m) = a.bond_atom(oi, m);
            a.bond_type(ni, m) = a.bond_type(oi, m);
        }
        a.num_angle(ni) = a.num_angle(oi);
        for (int m = 0; m < 2 * cpu::Atom::kMaxAnglePerAtom; ++m)
            a.angle_atom(ni, m) = a.angle_atom(oi, m);
        for (int m = 0; m < cpu::Atom::kMaxAnglePerAtom; ++m)
            a.angle_type(ni, m) = a.angle_type(oi, m);
        a.num_dihedral(ni) = a.num_dihedral(oi);
        for (int m = 0; m < 3 * cpu::Atom::kMaxDihedralPerAtom; ++m)
            a.dihedral_atom(ni, m) = a.dihedral_atom(oi, m);
        for (int m = 0; m < cpu::Atom::kMaxDihedralPerAtom; ++m)
            a.dihedral_type(ni, m) = a.dihedral_type(oi, m);
        a.nspecial(ni, 0) = a.nspecial(oi, 0);
        a.nspecial(ni, 1) = a.nspecial(oi, 1);
        a.nspecial(ni, 2) = a.nspecial(oi, 2);
        for (int m = 0; m < cpu::Atom::kMaxSpecial; ++m)
            a.special(ni, m) = a.special(oi, m);
    }


    tag_to_local.clear();
    tag_to_local.reserve(static_cast<std::size_t>(n_local));
    for (int i = 0; i < n_local; ++i) {
        tag_to_local.emplace(frame.atoms.tag(i), i);
    }

    // No topology filter at partition. Each rank keeps the FULL frame.bonds
    // / angles / dihedrals so that:
    //   * bonded force kernels can iterate all tuples and filter by their
    //     owner-tag (LAMMPS-style ownership rule) at force-eval time;
    //   * LJ::set_special_bonds() can build a complete tag-indexed
    //     exclusion map without missing bonds whose owner lives on another
    //     rank -- otherwise an owner-ghost LJ pair that is bonded but
    //     ownership-routed away would not be excluded and would blow up.
    //
    // Trade-off: each rank scans all O(N) tuples per force call instead
    // of just its O(N/ranks) owned subset. For < 1M-atom systems this is
    // a ~1 ms tax per step. A future optimisation can split: keep full
    // topology for the LJ exclusion build, filter for the force pass.
    ++neighbor_generation;
}

#endif // ATV_HAS_MPI

} // namespace mpi
} // namespace md
} // namespace cpu
} // namespace ATV
