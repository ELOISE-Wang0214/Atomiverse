/// @file src/cpu/md/mpi/halo_exchange.cpp
/// @brief `HaloExchangeForward` / `HaloExchangeReverse` real implementation.
///
/// Phase 2: orthorhombic, 6 face neighbors, non-blocking `MPI_Isend/Irecv`.
/// Ghost slots are appended in face order `{-x, +x, -y, +y, -z, +z}` so
/// `Reverse` can mirror the layout via `dd->recv_counts[]`.

#include "atomiverse/cpu/md/mpi/halo_exchange.h"

#include <array>
#include <cstring>
#include <numeric>
#include <vector>

#include "atomiverse/cpu/frame.h"
#include "atomiverse/cpu/md/mpi/domain.h"

#ifdef ATV_HAS_MPI
#include <mpi.h>
#endif

namespace ATV {
namespace cpu {
namespace md {
namespace mpi {

#ifdef ATV_HAS_MPI

namespace {

/// @brief Per-face MPI buffers. Owned per `HaloExchangeForward` /
/// `HaloExchangeReverse` instance; allocated in `build()`, sized to the
/// initial `dd->capacity`.
struct FaceBufs {
    std::array<std::vector<double>, 6> send_pos;   // packed positions per face
    std::array<std::vector<double>, 6> recv_pos;   // received positions per face
    std::array<std::vector<double>, 6> send_force; // packed forces (Reverse)
    std::array<std::vector<double>, 6> recv_force; // received forces (Reverse)
};

/// @brief Helper: opposite face index (toggle low bit of {-x, +x, -y, +y, -z, +z}).
inline int opposite_face(int f) { return f ^ 1; }

/// @brief Map face index to {dim, sign}. dim = f / 2, sign = (f & 1) ? +1 : -1.
inline int face_dim(int f)  { return f >> 1; }
inline int face_sign(int f) { return (f & 1) ? +1 : -1; }

/// @brief PBC shift along `dim` for atoms going out face `f` from this rank.
///
/// If this rank sits on the periodic edge in the +dim direction
/// (`coords[dim] == dims[dim]-1`) and we are sending out the +dim face,
/// the receiving rank is at the opposite edge of the global box;
/// pre-subtract `L_dim` from the packed coord so the receiver sees the
/// atom in its own coordinate frame. Symmetric for the -dim edge.
inline double pbc_shift_for_face(const DomainDecomposition& dd,
                                  int f, double L_dim) {
    const int dim = face_dim(f);
    const int sign = face_sign(f);
    if (sign > 0 && dd.coords[dim] == dd.dims[dim] - 1) return -L_dim;
    if (sign < 0 && dd.coords[dim] == 0)               return +L_dim;
    return 0.0;
}

/// @brief Decide which owned atoms qualify as ghosts for each face.
///
/// An atom qualifies for face f along dim/sign when its position is within
/// `cutoff` of the slab boundary in that direction. Wraps `dd->send_lists[f]`
/// in place; `vector<int>` capacity is preserved across rebuilds.
void rebuild_send_lists(const Frame& frame, DomainDecomposition& dd, double cutoff) {
    for (int f = 0; f < 6; ++f) dd.send_lists[f].clear();

    const int n_local = dd.n_local;
    const double* x = frame.atoms.x.data();
    const double* y = frame.atoms.y.data();
    const double* z = frame.atoms.z.data();

    for (int i = 0; i < n_local; ++i) {
        // -x: atom near sublo[0]
        if (x[i] - dd.sublo[0] < cutoff) dd.send_lists[0].push_back(i);
        // +x: atom near subhi[0]
        if (dd.subhi[0] - x[i] < cutoff) dd.send_lists[1].push_back(i);
        // -y / +y
        if (y[i] - dd.sublo[1] < cutoff) dd.send_lists[2].push_back(i);
        if (dd.subhi[1] - y[i] < cutoff) dd.send_lists[3].push_back(i);
        // -z / +z
        if (z[i] - dd.sublo[2] < cutoff) dd.send_lists[4].push_back(i);
        if (dd.subhi[2] - z[i] < cutoff) dd.send_lists[5].push_back(i);
    }
}

/// @brief Pack packed positions for face `f` into a contiguous buffer.
///
/// Layout: `[x0, y0, z0, x1, y1, z1, ...]`. Applies the per-face PBC shift
/// to the packed coordinate for atoms whose owner sits on the periodic
/// edge along the send direction.
void pack_positions_for_face(const Frame& frame, const DomainDecomposition& dd,
                              int f, std::vector<double>& buf) {
    const auto& list = dd.send_lists[f];
    const int n = static_cast<int>(list.size());
    buf.resize(static_cast<std::size_t>(n) * 3);

    const int dim = face_dim(f);
    const double shift = pbc_shift_for_face(dd, f, frame.box(dim, dim));

    const double* x = frame.atoms.x.data();
    const double* y = frame.atoms.y.data();
    const double* z = frame.atoms.z.data();

    for (int k = 0; k < n; ++k) {
        const int i = list[static_cast<std::size_t>(k)];
        double xi = x[i], yi = y[i], zi = z[i];
        if (dim == 0) xi += shift;
        else if (dim == 1) yi += shift;
        else if (dim == 2) zi += shift;
        buf[static_cast<std::size_t>(3 * k + 0)] = xi;
        buf[static_cast<std::size_t>(3 * k + 1)] = yi;
        buf[static_cast<std::size_t>(3 * k + 2)] = zi;
    }
}

/// @brief Per-face buffer storage owned by the calling component.
/// Static thread-local? No -- we use a per-component instance via a small
/// helper handle so that Forward and Reverse do not share buffer state
/// (they share only `dd->send_lists` and `dd->recv_counts`).
FaceBufs& bufs_for(HaloExchangeForward* /*comp*/) {
    static thread_local FaceBufs s; return s;
}
FaceBufs& bufs_for(HaloExchangeReverse* /*comp*/) {
    static thread_local FaceBufs s; return s;
}

void prealloc_bufs(FaceBufs& b, int initial_cap) {
    for (int f = 0; f < 6; ++f) {
        // Conservative: each face could ship up to `initial_cap` atoms, each
        // 3 doubles for positions, 3 doubles for forces. Bounded by
        // 6 * initial_cap * 6 * 8 B = ~288 B per atom of capacity, i.e.
        // ~290 KB at 1000 atoms -- trivial.
        b.send_pos[f].reserve(static_cast<std::size_t>(initial_cap) * 3);
        b.recv_pos[f].reserve(static_cast<std::size_t>(initial_cap) * 3);
        b.send_force[f].reserve(static_cast<std::size_t>(initial_cap) * 3);
        b.recv_force[f].reserve(static_cast<std::size_t>(initial_cap) * 3);
    }
}

} // namespace

void HaloExchangeForward::build(const cpu::BuildContext& ctx) {
    if (cutoff_ <= 0.0) return;
    auto* dd = ctx.runtime.dd;
    if (!dd || ctx.runtime.mpi_size <= 1) return;
    prealloc_bufs(bufs_for(this), dd->capacity);
    for (int f = 0; f < 6; ++f) dd->send_lists[f].reserve(static_cast<std::size_t>(dd->capacity));
}

void HaloExchangeForward::launch(const cpu::LaunchContext& ctx) {
    if (ctx.runtime.mpi_size <= 1 || cutoff_ <= 0.0) return;
    auto* dd = ctx.runtime.dd;
    if (!dd || dd->cart_comm == MPI_COMM_NULL) return;

    Frame& frame = ctx.frame;
    FaceBufs& bufs = bufs_for(this);

    rebuild_send_lists(frame, *dd, cutoff_);

    // Phase 1: exchange counts on each face so the receivers know how big
    // each incoming buffer will be.
    int send_counts[6];
    int recv_counts[6];
    MPI_Request count_reqs[12];
    int n_count_reqs = 0;
    for (int f = 0; f < 6; ++f) {
        send_counts[f] = static_cast<int>(dd->send_lists[f].size());
        recv_counts[f] = 0;
    }
    for (int f = 0; f < 6; ++f) {
        const int neighbor = dd->face_neighbors[static_cast<std::size_t>(f)];
        if (neighbor == MPI_PROC_NULL) continue;
        // Send my count for face f to neighbor; tag = f.
        MPI_Isend(&send_counts[f], 1, MPI_INT, neighbor, f,
                  dd->cart_comm, &count_reqs[n_count_reqs++]);
        // Receive their count from the opposite-face direction; the
        // neighbor's view of "my face f" is their face f^1, so tag = f^1.
        MPI_Irecv(&recv_counts[f], 1, MPI_INT, neighbor, opposite_face(f),
                  dd->cart_comm, &count_reqs[n_count_reqs++]);
    }
    if (n_count_reqs > 0)
        MPI_Waitall(n_count_reqs, count_reqs, MPI_STATUSES_IGNORE);

    // Phase 2: pack positions on each face, then exchange them.
    for (int f = 0; f < 6; ++f) {
        if (dd->face_neighbors[static_cast<std::size_t>(f)] == MPI_PROC_NULL) {
            bufs.send_pos[f].clear();
        } else {
            pack_positions_for_face(frame, *dd, f, bufs.send_pos[f]);
        }
        bufs.recv_pos[f].assign(static_cast<std::size_t>(recv_counts[f]) * 3, 0.0);
    }

    MPI_Request pos_reqs[12];
    int n_pos_reqs = 0;
    for (int f = 0; f < 6; ++f) {
        const int neighbor = dd->face_neighbors[static_cast<std::size_t>(f)];
        if (neighbor == MPI_PROC_NULL) continue;
        if (send_counts[f] > 0) {
            MPI_Isend(bufs.send_pos[f].data(),
                      static_cast<int>(bufs.send_pos[f].size()),
                      MPI_DOUBLE, neighbor, 100 + f,
                      dd->cart_comm, &pos_reqs[n_pos_reqs++]);
        }
        if (recv_counts[f] > 0) {
            MPI_Irecv(bufs.recv_pos[f].data(),
                      static_cast<int>(bufs.recv_pos[f].size()),
                      MPI_DOUBLE, neighbor, 100 + opposite_face(f),
                      dd->cart_comm, &pos_reqs[n_pos_reqs++]);
        }
    }
    if (n_pos_reqs > 0)
        MPI_Waitall(n_pos_reqs, pos_reqs, MPI_STATUSES_IGNORE);

    // Phase 3: unpack received positions into ghost slots, in face order.
    int slot = dd->n_local;
    for (int f = 0; f < 6; ++f) {
        const int n = recv_counts[f];
        for (int k = 0; k < n; ++k) {
            frame.atoms.x[slot] = bufs.recv_pos[f][static_cast<std::size_t>(3 * k + 0)];
            frame.atoms.y[slot] = bufs.recv_pos[f][static_cast<std::size_t>(3 * k + 1)];
            frame.atoms.z[slot] = bufs.recv_pos[f][static_cast<std::size_t>(3 * k + 2)];
            ++slot;
        }
        dd->recv_counts[static_cast<std::size_t>(f)] = n;
    }
    dd->n_ghost = slot - dd->n_local;
}

void HaloExchangeReverse::build(const cpu::BuildContext& ctx) {
    auto* dd = ctx.runtime.dd;
    if (!dd || ctx.runtime.mpi_size <= 1) return;
    prealloc_bufs(bufs_for(this), dd->capacity);
}

void HaloExchangeReverse::launch(const cpu::LaunchContext& ctx) {
    if (ctx.runtime.mpi_size <= 1) return;
    auto* dd = ctx.runtime.dd;
    if (!dd || dd->cart_comm == MPI_COMM_NULL) return;
    if (dd->n_ghost == 0) return;

    Frame& frame = ctx.frame;
    FaceBufs& bufs = bufs_for(this);

    // Pack ghost forces by face (reverse direction sends ghost forces back
    // to their owners). Layout mirrors the recv layout established by Forward.
    int slot = dd->n_local;
    int send_force_counts[6] = {0, 0, 0, 0, 0, 0};
    int recv_force_counts[6] = {0, 0, 0, 0, 0, 0};
    for (int f = 0; f < 6; ++f) {
        const int n = dd->recv_counts[static_cast<std::size_t>(f)];
        send_force_counts[f] = n;
        bufs.send_force[f].resize(static_cast<std::size_t>(n) * 3);
        for (int k = 0; k < n; ++k) {
            bufs.send_force[f][static_cast<std::size_t>(3 * k + 0)] = frame.atoms.fx[slot];
            bufs.send_force[f][static_cast<std::size_t>(3 * k + 1)] = frame.atoms.fy[slot];
            bufs.send_force[f][static_cast<std::size_t>(3 * k + 2)] = frame.atoms.fz[slot];
            ++slot;
        }
        // The corresponding inbound count is `send_lists[f].size() * 3`:
        // the neighbor will return forces for the atoms I sent it in Forward,
        // and those atoms are at indices `dd->send_lists[f]` in my local frame.
        recv_force_counts[f] = static_cast<int>(dd->send_lists[static_cast<std::size_t>(f)].size());
        bufs.recv_force[f].assign(static_cast<std::size_t>(recv_force_counts[f]) * 3, 0.0);
    }

    MPI_Request reqs[12];
    int n_reqs = 0;
    for (int f = 0; f < 6; ++f) {
        const int neighbor = dd->face_neighbors[static_cast<std::size_t>(f)];
        if (neighbor == MPI_PROC_NULL) continue;
        if (send_force_counts[f] > 0) {
            MPI_Isend(bufs.send_force[f].data(),
                      static_cast<int>(bufs.send_force[f].size()),
                      MPI_DOUBLE, neighbor, 200 + f,
                      dd->cart_comm, &reqs[n_reqs++]);
        }
        if (recv_force_counts[f] > 0) {
            MPI_Irecv(bufs.recv_force[f].data(),
                      static_cast<int>(bufs.recv_force[f].size()),
                      MPI_DOUBLE, neighbor, 200 + opposite_face(f),
                      dd->cart_comm, &reqs[n_reqs++]);
        }
    }
    if (n_reqs > 0)
        MPI_Waitall(n_reqs, reqs, MPI_STATUSES_IGNORE);

    // Add received ghost forces back into local atoms. The k-th atom in
    // `dd->send_lists[f]` corresponds to the k-th triple in `recv_force[f]`.
    for (int f = 0; f < 6; ++f) {
        const auto& list = dd->send_lists[static_cast<std::size_t>(f)];
        const int n = static_cast<int>(list.size());
        for (int k = 0; k < n; ++k) {
            const int i = list[static_cast<std::size_t>(k)];
            frame.atoms.fx[i] += bufs.recv_force[f][static_cast<std::size_t>(3 * k + 0)];
            frame.atoms.fy[i] += bufs.recv_force[f][static_cast<std::size_t>(3 * k + 1)];
            frame.atoms.fz[i] += bufs.recv_force[f][static_cast<std::size_t>(3 * k + 2)];
        }
    }
}

#else // !ATV_HAS_MPI

void HaloExchangeForward::build(const cpu::BuildContext& /*ctx*/) {}
void HaloExchangeForward::launch(const cpu::LaunchContext& /*ctx*/) {}
void HaloExchangeReverse::build(const cpu::BuildContext& /*ctx*/) {}
void HaloExchangeReverse::launch(const cpu::LaunchContext& /*ctx*/) {}

#endif // ATV_HAS_MPI

} // namespace mpi
} // namespace md
} // namespace cpu
} // namespace ATV
