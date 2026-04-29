/// @file cpu/md/mpi/domain.h
/// @brief `DomainDecomposition` -- MD-specific spatial domain decomposition state.
///
/// Owned by `cpu::Backend`; non-owning pointer published through `cpu::Context::dd`.
/// In single-rank mode the topology is `dims={1,1,1}` and `n_ghost==0`; the
/// struct is still allocated so components can read it without null-checks beyond
/// the usual `mpi_size==1` shortcut.
///
/// @see .claude/specs/parallel-md-engine.md  § A.1
#pragma once

#include <array>
#include <cstdint>
#include <unordered_map>
#include <vector>

#ifdef ATV_HAS_MPI
#include <mpi.h>
#endif

namespace ATV {

struct Mat3;

namespace cpu {

class Frame;

namespace md {
namespace mpi {

/// @brief MD spatial domain decomposition state for a single rank.
///
/// Holds (a) MPI Cartesian topology, (b) this rank's owned spatial slab,
/// (c) authoritative atom counts (`n_local`, `n_ghost`, `capacity`),
/// (d) a generation counter that consumers (neighbor lists, halo lists)
/// rebuild against, and (e) derived bookkeeping (send lists, tag→local map).
///
/// `Frame::atoms.size()` is the *allocated capacity* in multi-rank mode.
/// In Phase 1 single-rank mode this still equals `n_local`; the loop-bound
/// audit (replacing `atoms.size()` with `n_local` in force kernels) is
/// deferred to Phase 2 where it becomes observable.
struct DomainDecomposition {
#ifdef ATV_HAS_MPI
    /// @brief Cartesian sub-communicator created by `MPI_Cart_create`.
    MPI_Comm cart_comm = MPI_COMM_NULL;
#endif

    /// @brief Cartesian process grid dimensions. Single-rank: {1,1,1}.
    int dims[3] = {1, 1, 1};

    /// @brief This rank's coords in the Cartesian grid. Single-rank: {0,0,0}.
    int coords[3] = {0, 0, 0};

    /// @brief Periodic flags per dimension (1=periodic, 0=non-periodic).
    int periods[3] = {1, 1, 1};

    /// @brief Face neighbor ranks: order is -x, +x, -y, +y, -z, +z.
    /// Single-rank: all zeros (this rank is its own neighbor along periodic axes).
    std::array<int, 6> face_neighbors{};

    /// @brief This rank's owned spatial slab, low corner [Angstrom].
    double sublo[3]{};

    /// @brief This rank's owned spatial slab, high corner [Angstrom].
    double subhi[3]{};

    /// @brief Number of locally owned atoms (authoritative count).
    int n_local = 0;

    /// @brief Number of ghost copies received from neighbors.
    int n_ghost = 0;

    /// @brief Allocated atom-array capacity. `n_local + n_ghost <= capacity`.
    int capacity = 0;

    /// @brief Generation counter, bumped on migration. Neighbor-list / bonded /
    /// halo-list consumers rebuild when they observe a bump.
    std::uint64_t neighbor_generation = 0;

    /// @brief Per-face send lists: owned-atom indices shadowed to each face neighbor.
    /// Order matches `face_neighbors`. Empty in single-rank mode. Rebuilt by
    /// `HaloExchangeForward::launch()` each step; consumed by both
    /// `HaloExchangeForward` (to pack outgoing positions) and
    /// `HaloExchangeReverse` (to add incoming ghost forces back to owners).
    std::array<std::vector<int>, 6> send_lists{};

    /// @brief Per-face ghost counts received in the most recent forward halo.
    /// `recv_counts[f]` ghosts from face `f` occupy frame slots
    /// `[n_local + sum(recv_counts[0..f]), n_local + sum(recv_counts[0..f+1]))`.
    /// Used by `HaloExchangeReverse` to know how many ghost forces to ship
    /// back per face. Reset to zero on partition; set by Forward.
    std::array<int, 6> recv_counts{};

    /// @brief Map from global atom tag to local row index.
    ///
    /// Includes both owned atoms (slots `[0, n_local)`) and current-step
    /// ghost atoms (slots `[n_local, n_local + n_ghost)`), so bonded force
    /// kernels can resolve a tuple's atom-tags to local rows regardless of
    /// whether the atom is an owner or a halo ghost.
    ///
    /// Owner entries are rebuilt on every `partition_atoms()` /
    /// `AtomMigration`. Ghost entries are appended by `HaloExchangeForward`
    /// and erased back out at the start of the next halo via the
    /// companion `ghost_tags` log.
    std::unordered_map<std::int64_t, int> tag_to_local;

    /// @brief Tags of the ghost atoms that `HaloExchangeForward` last
    /// inserted into `tag_to_local`, in slot order. Used to roll back
    /// ghost entries before the next halo write so a different ghost set
    /// does not leave stale tag→slot bindings behind.
    std::vector<std::int64_t> ghost_tags;

    /// @brief Per-axis cut planes as **fractions** of the box length, in `[0, 1]`.
    ///
    /// `cuts[d]` has either size `0` (uniform partition along axis `d`, default)
    /// or size `dims[d] + 1` with `cuts[d].front() == 0.0` and `cuts[d].back() == 1.0`.
    /// When non-empty, `compute_sublo_subhi()` reads from these instead of the
    /// uniform `coord[d] / dims[d]` formula, and `AtomMigration::dest_rank_for`
    /// binary-searches the cuts to decide each owner's destination rank.
    ///
    /// Stored as fractions so the partition survives NPT box rescaling without
    /// needing to be recomputed every step. Populated by `RCBBalancer` (Phase 4)
    /// or any future load-balance component; left empty in default Phase 1–3
    /// runs (uniform-brick decomposition).
    ///
    /// All ranks must agree on these values for migration to be consistent.
    /// `RCBBalancer` derives them from a global atom histogram via
    /// `MPI_Allreduce`, which guarantees consistency by construction.
    std::array<std::vector<double>, 3> cuts{};

    /// @brief Bind to a frame: populate `n_local`, `capacity`, `tag_to_local`.
    ///
    /// Single-rank initialization. Sets `n_local = capacity = frame.atoms.size()`,
    /// `n_ghost = 0`, and rebuilds `tag_to_local` from `frame.atoms.tag`.
    /// Does not touch MPI topology fields.
    ///
    /// @param frame Host-side frame whose atoms this DD describes.
    void bind(const Frame& frame);

#ifdef ATV_HAS_MPI
    /// @brief Release the Cartesian sub-communicator if still alive.
    ///
    /// Called by `Backend::shutdown()` before `MPI_Finalize` so OpenMPI
    /// does not segfault trying to clean up a user comm that outlived the
    /// runtime. Safe to call when `cart_comm == MPI_COMM_NULL`.
    /// @param mpi_is_finalizing If true, the MPI runtime is shutting
    ///        down and we skip the `MPI_Comm_free` (avoid double-free
    ///        against an already-dead runtime).
    void release_comm(bool mpi_is_finalizing = false);
#endif

    /// @brief Compute and store this rank's spatial slab from the box.
    ///
    /// Phase 2 (orthorhombic): `sublo[d] = (L_d / dims[d]) * coords[d]`,
    /// `subhi[d] = sublo[d] + L_d / dims[d]`. Idempotent. Single-rank
    /// (`dims == {1,1,1}`): the slab is the full box. Not MPI-gated: the
    /// slab record is useful for any spatial component, e.g. the neighbor
    /// list and future observer shell analysis.
    /// @param box Simulation box (orthorhombic in Phase 2).
    void compute_sublo_subhi(const ATV::Mat3& box);

#ifdef ATV_HAS_MPI
    /// @brief Build the 3D Cartesian topology for `world` and discover the
    /// 6 face neighbors.
    ///
    /// `MPI_Dims_create` chooses a balanced `dims` factorization of the
    /// world size; `MPI_Cart_create` produces `cart_comm` with periodic
    /// boundaries on all axes; `MPI_Cart_coords` populates `coords`;
    /// `MPI_Cart_shift` populates `face_neighbors` in the order
    /// `{-x, +x, -y, +y, -z, +z}`.
    ///
    /// Single-rank world (`world_size == 1`): the call is a no-op and
    /// `dims` stays `{1,1,1}`.
    ///
    /// Calling this twice is allowed; the existing `cart_comm` (if any)
    /// is freed before the new one is created.
    /// @param world Parent communicator (typically `MPI_COMM_WORLD`).
    void init_topology(MPI_Comm world);

    /// @brief Partition atoms across ranks by spatial position.
    ///
    /// Each rank inspects all atoms in `frame.atoms[0..size())`, picks the
    /// subset whose (PBC-wrapped) positions fall in this rank's slab,
    /// compacts those atoms to the front of every per-atom array, and
    /// updates `n_local`, `n_ghost = 0`, `capacity = atoms.size()`. Stale
    /// rows in `[n_local, capacity)` are reused later by `HaloExchangeForward`
    /// as ghost slots.
    ///
    /// **Phase 2 contract:** every rank starts with the same `frame`.
    /// Migration after integration is **not** performed (Phase 3). Setups
    /// where atoms drift beyond their initial slab within the run produce
    /// wrong physics under Phase 2; the spec mandates short runs at low
    /// temperature for the exit-gate test.
    ///
    /// Single-rank: behaves like `bind(frame)`.
    /// @param frame Frame whose atoms this rank's slab claims a subset of.
    void partition_atoms(ATV::cpu::Frame& frame);
#endif
};

} // namespace mpi
} // namespace md
} // namespace cpu
} // namespace ATV
