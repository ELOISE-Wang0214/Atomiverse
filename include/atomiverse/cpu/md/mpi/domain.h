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
    /// Order matches `face_neighbors`. Empty in single-rank mode.
    std::array<std::vector<int>, 6> send_lists{};

    /// @brief Map from global atom tag to local row index. Rebuilt on migration.
    std::unordered_map<std::int64_t, int> tag_to_local;

    /// @brief Bind to a frame: populate `n_local`, `capacity`, `tag_to_local`.
    ///
    /// Single-rank initialization. Sets `n_local = capacity = frame.atoms.size()`,
    /// `n_ghost = 0`, and rebuilds `tag_to_local` from `frame.atoms.tag`.
    /// Does not touch MPI topology fields.
    ///
    /// @param frame Host-side frame whose atoms this DD describes.
    void bind(const Frame& frame);
};

} // namespace mpi
} // namespace md
} // namespace cpu
} // namespace ATV
