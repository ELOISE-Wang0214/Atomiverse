/// @file src/cpu/md/mpi/domain.cpp
/// @brief `DomainDecomposition` implementation: single-rank `bind`, plus
/// multi-rank `init_topology`, `compute_sublo_subhi`, and `partition_atoms`
/// (Phase 2).

#include "atomiverse/cpu/md/mpi/domain.h"

#include <algorithm>
#include <cmath>

#include "atomiverse/cpu/frame.h"
#include "atomiverse/core/mat3.h"

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

void DomainDecomposition::compute_sublo_subhi(const ATV::Mat3& box) {
    for (int d = 0; d < 3; ++d) {
        const double L = box(d, d);
        const double slab = L / static_cast<double>(dims[d]);
        sublo[d] = slab * static_cast<double>(coords[d]);
        subhi[d] = sublo[d] + slab;
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
    dims[0] = dims[1] = dims[2] = 0;
    MPI_Dims_create(world_size, 3, dims);

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
    face_neighbors[0] = rank_at(coords[0] - 1, coords[1],     coords[2]    );
    face_neighbors[1] = rank_at(coords[0] + 1, coords[1],     coords[2]    );
    face_neighbors[2] = rank_at(coords[0],     coords[1] - 1, coords[2]    );
    face_neighbors[3] = rank_at(coords[0],     coords[1] + 1, coords[2]    );
    face_neighbors[4] = rank_at(coords[0],     coords[1],     coords[2] - 1);
    face_neighbors[5] = rank_at(coords[0],     coords[1],     coords[2] + 1);
}

namespace {

// Orthorhombic PBC wrap: bring `p` into [0, L).
inline double wrap_to_cell(double p, double L) {
    if (L <= 0.0) return p;
    double q = p - std::floor(p / L) * L;
    if (q < 0.0) q += L;
    if (q >= L) q -= L;
    return q;
}

// Rank index (along one dim) that owns wrapped position `q` with `dims_d`
// equal-width slabs over `[0, L)`. Clamps to [0, dims_d - 1] against FP
// rounding at the upper edge.
inline int rank_along_dim(double q, double L, int dims_d) {
    if (dims_d <= 1) return 0;
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

    const double Lx = frame.box(0, 0);
    const double Ly = frame.box(1, 1);
    const double Lz = frame.box(2, 2);

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

    for (int i = 0; i < n_to_iterate; ++i) {
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

    n_local  = static_cast<int>(owned_indices.size());
    n_ghost  = 0;
    capacity = frame.atoms.size();  // keep the whole allocation as slack for ghosts

    // Compact owned atoms to the front of every per-atom array. Iterating
    // in increasing `new_i` guarantees we read `old_i >= new_i` so the
    // destination has not yet been overwritten.
    auto compact_1d_d = [&](xt::xtensor<double, 1>& arr) {
        for (int ni = 0; ni < n_local; ++ni) {
            const int oi = owned_indices[static_cast<std::size_t>(ni)];
            if (oi != ni) arr(ni) = arr(oi);
        }
    };
    auto compact_1d_i = [&](xt::xtensor<int, 1>& arr) {
        for (int ni = 0; ni < n_local; ++ni) {
            const int oi = owned_indices[static_cast<std::size_t>(ni)];
            if (oi != ni) arr(ni) = arr(oi);
        }
    };
    auto compact_1d_i64 = [&](xt::xtensor<std::int64_t, 1>& arr) {
        for (int ni = 0; ni < n_local; ++ni) {
            const int oi = owned_indices[static_cast<std::size_t>(ni)];
            if (oi != ni) arr(ni) = arr(oi);
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

    tag_to_local.clear();
    tag_to_local.reserve(static_cast<std::size_t>(n_local));
    for (int i = 0; i < n_local; ++i) {
        tag_to_local.emplace(frame.atoms.tag(i), i);
    }

    // Topology filter (Phase 0.3 final): on the first partition, keep the
    // tuples whose OWNER atom (per LAMMPS-style rule) is owned by this
    // rank. The other atoms in the tuple may live on neighbour ranks --
    // the bonded force kernels resolve those via `tag_to_local` after
    // `HaloExchangeForward` has populated ghost slots and added their
    // tags to the map.
    //
    //   bond  (i,j):    owner = rank holding atom with min(tag[i], tag[j])
    //   angle (i,j,k):  owner = rank holding atom j (centre)
    //   dihedral i-j-k-l: owner = rank holding atom j (2nd atom)
    //
    // Atom indices stored in `Frame::bonds/angles/dihedrals` carry GLOBAL
    // TAGS, never local row indices. The pre-partition Frame has
    // `tag(i) == i` (default), so the user's existing `atomi = i` setup
    // is correct as-is. Force kernels resolve through `tag_to_local`.
    if (first_partition) {
        auto owns_tag = [&](std::int64_t tag) {
            return tag_to_local.find(tag) != tag_to_local.end();
        };

        // Bonds: owner = rank of lower-tag atom.
        {
            const int n_b = frame.bonds.size();
            std::vector<int> ki, kj, kt;
            ki.reserve(static_cast<std::size_t>(n_b));
            kj.reserve(static_cast<std::size_t>(n_b));
            kt.reserve(static_cast<std::size_t>(n_b));
            for (int b = 0; b < n_b; ++b) {
                const std::int64_t ti = frame.bonds.atomi(b);
                const std::int64_t tj = frame.bonds.atomj(b);
                const std::int64_t owner_tag = (ti < tj) ? ti : tj;
                if (!owns_tag(owner_tag)) continue;
                ki.push_back(static_cast<int>(ti));
                kj.push_back(static_cast<int>(tj));
                kt.push_back(frame.bonds.type(b));
            }
            const int new_n = static_cast<int>(ki.size());
            frame.bonds.atomi = xt::zeros<int>({new_n});
            frame.bonds.atomj = xt::zeros<int>({new_n});
            frame.bonds.type  = xt::zeros<int>({new_n});
            for (int b = 0; b < new_n; ++b) {
                frame.bonds.atomi[b] = ki[static_cast<std::size_t>(b)];
                frame.bonds.atomj[b] = kj[static_cast<std::size_t>(b)];
                frame.bonds.type[b]  = kt[static_cast<std::size_t>(b)];
            }
        }

        // Angles: owner = rank of centre atom (atomj).
        {
            const int n_a = frame.angles.size();
            std::vector<int> ki, kj, kk, kt;
            ki.reserve(static_cast<std::size_t>(n_a));
            kj.reserve(static_cast<std::size_t>(n_a));
            kk.reserve(static_cast<std::size_t>(n_a));
            kt.reserve(static_cast<std::size_t>(n_a));
            for (int a = 0; a < n_a; ++a) {
                if (!owns_tag(frame.angles.atomj(a))) continue;
                ki.push_back(frame.angles.atomi(a));
                kj.push_back(frame.angles.atomj(a));
                kk.push_back(frame.angles.atomk(a));
                kt.push_back(frame.angles.type(a));
            }
            const int new_n = static_cast<int>(ki.size());
            frame.angles.atomi = xt::zeros<int>({new_n});
            frame.angles.atomj = xt::zeros<int>({new_n});
            frame.angles.atomk = xt::zeros<int>({new_n});
            frame.angles.type  = xt::zeros<int>({new_n});
            for (int a = 0; a < new_n; ++a) {
                frame.angles.atomi[a] = ki[static_cast<std::size_t>(a)];
                frame.angles.atomj[a] = kj[static_cast<std::size_t>(a)];
                frame.angles.atomk[a] = kk[static_cast<std::size_t>(a)];
                frame.angles.type[a]  = kt[static_cast<std::size_t>(a)];
            }
        }

        // Dihedrals: owner = rank of 2nd atom (atomj).
        {
            const int n_d = frame.dihedrals.size();
            std::vector<int> ki, kj, kk, kl, kt;
            ki.reserve(static_cast<std::size_t>(n_d));
            kj.reserve(static_cast<std::size_t>(n_d));
            kk.reserve(static_cast<std::size_t>(n_d));
            kl.reserve(static_cast<std::size_t>(n_d));
            kt.reserve(static_cast<std::size_t>(n_d));
            for (int d = 0; d < n_d; ++d) {
                if (!owns_tag(frame.dihedrals.atomj(d))) continue;
                ki.push_back(frame.dihedrals.atomi(d));
                kj.push_back(frame.dihedrals.atomj(d));
                kk.push_back(frame.dihedrals.atomk(d));
                kl.push_back(frame.dihedrals.atoml(d));
                kt.push_back(frame.dihedrals.type(d));
            }
            const int new_n = static_cast<int>(ki.size());
            frame.dihedrals.atomi = xt::zeros<int>({new_n});
            frame.dihedrals.atomj = xt::zeros<int>({new_n});
            frame.dihedrals.atomk = xt::zeros<int>({new_n});
            frame.dihedrals.atoml = xt::zeros<int>({new_n});
            frame.dihedrals.type  = xt::zeros<int>({new_n});
            for (int d = 0; d < new_n; ++d) {
                frame.dihedrals.atomi[d] = ki[static_cast<std::size_t>(d)];
                frame.dihedrals.atomj[d] = kj[static_cast<std::size_t>(d)];
                frame.dihedrals.atomk[d] = kk[static_cast<std::size_t>(d)];
                frame.dihedrals.atoml[d] = kl[static_cast<std::size_t>(d)];
                frame.dihedrals.type[d]  = kt[static_cast<std::size_t>(d)];
            }
        }
    }

    ++neighbor_generation;
}

#endif // ATV_HAS_MPI

} // namespace mpi
} // namespace md
} // namespace cpu
} // namespace ATV
