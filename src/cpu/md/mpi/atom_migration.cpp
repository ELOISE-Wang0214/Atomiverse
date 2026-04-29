/// @file src/cpu/md/mpi/atom_migration.cpp
/// @brief `AtomMigration` -- Phase 0.2 / P3.1 implementation.
///
/// For each owner whose post-integration position lies outside this rank's
/// spatial slab, the atom (and all its per-atom fields) is shipped to the
/// rank whose slab contains the new position. Implemented as a single
/// `MPI_Alltoallv` over byte buffers so every world rank participates and
/// either packs/unpacks zero bytes or some multiple of `kAtomPayloadBytes`.
///
/// Bonded-tuple migration is **not** handled here; that is tracked under
/// Phase 3.2 of msu-pe-validation.md.

#include "atomiverse/cpu/md/mpi/atom_migration.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
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

// Per-atom payload layout (must stay in sync between pack and unpack).
//   doubles: x, y, z, vx, vy, vz, fx, fy, fz, mass        (10)
//   int32:   type                                          (1)
//   int64:   tag                                           (1)
//   int8:    image[3]                                      (3)
constexpr int kPayloadDoubles = 10;
constexpr int kPayloadInt32   = 1;
constexpr int kPayloadInt64   = 1;
constexpr int kPayloadInt8    = 3;
constexpr int kAtomPayloadBytes =
    sizeof(double)  * kPayloadDoubles
  + sizeof(std::int32_t) * kPayloadInt32
  + sizeof(std::int64_t) * kPayloadInt64
  + sizeof(std::int8_t)  * kPayloadInt8;

// Wrap `q` into the half-open interval `[0, L)`. Identical to the helper
// in `domain.cpp`; duplicated here to keep the migration pack self-contained.
inline double wrap_to_box(double q, double L) {
    if (L <= 0.0) return q;
    double w = q - std::floor(q / L) * L;
    if (w < 0.0)  w += L;
    if (w >= L)   w -= L;
    return w;
}

inline int rank_along_dim(double q, double L, int dims_d) {
    if (dims_d <= 1) return 0;
    const int r = static_cast<int>(q * static_cast<double>(dims_d) / L);
    return std::max(0, std::min(dims_d - 1, r));
}

// Locate the slab index `k` (in [0, dims_d)) whose [cuts[k], cuts[k+1])
// contains the box-length fraction `frac` (already wrapped into [0, 1)).
// `cuts` has size `dims_d + 1` with `cuts.front() == 0` and `cuts.back() == 1`.
inline int slab_along_dim_from_cuts(double frac,
                                    const std::vector<double>& cuts,
                                    int dims_d) {
    if (dims_d <= 1) return 0;
    // upper_bound gives the first cuts[i] strictly greater than frac;
    // its predecessor is the slab index. Clamp against FP rounding at the
    // top edge.
    auto it = std::upper_bound(cuts.begin(), cuts.end(), frac);
    int k = static_cast<int>(std::distance(cuts.begin(), it)) - 1;
    if (k < 0) k = 0;
    if (k > dims_d - 1) k = dims_d - 1;
    return k;
}

// Compute the rank that owns position (x, y, z) under the current DD.
//
// When `dd.cuts[d]` is empty along axis d, the uniform `coord = floor(q*dims/L)`
// formula applies. When non-empty (set by `RCBBalancer`), the destination is
// derived from the box-length fraction via binary search through the cut
// planes -- which guarantees the same rank assignment that
// `DomainDecomposition::compute_sublo_subhi` would produce on the destination
// side.
inline int dest_rank_for(const DomainDecomposition& dd, double x, double y, double z,
                          double Lx, double Ly, double Lz) {
    auto axis_coord = [&](int d, double q, double L) {
        const double wrapped = wrap_to_box(q, L);
        if (static_cast<int>(dd.cuts[d].size()) == dd.dims[d] + 1) {
            const double frac = (L > 0.0) ? (wrapped / L) : 0.0;
            return slab_along_dim_from_cuts(frac, dd.cuts[d], dd.dims[d]);
        }
        return rank_along_dim(wrapped, L, dd.dims[d]);
    };
    const int cx = axis_coord(0, x, Lx);
    const int cy = axis_coord(1, y, Ly);
    const int cz = axis_coord(2, z, Lz);
    return cx * (dd.dims[1] * dd.dims[2]) + cy * dd.dims[2] + cz;
}

// Pack one atom's per-atom fields into `buf` at `offset`. Returns the new
// offset after writing `kAtomPayloadBytes` bytes.
std::size_t pack_atom(std::vector<char>& buf, std::size_t offset,
                       const Frame& frame, int i) {
    char* p = buf.data() + offset;
    const double dvals[kPayloadDoubles] = {
        frame.atoms.x[i],  frame.atoms.y[i],  frame.atoms.z[i],
        frame.atoms.vx[i], frame.atoms.vy[i], frame.atoms.vz[i],
        frame.atoms.fx[i], frame.atoms.fy[i], frame.atoms.fz[i],
        frame.atoms.mass[i],
    };
    std::memcpy(p, dvals, sizeof(dvals));
    p += sizeof(dvals);
    const std::int32_t typ = static_cast<std::int32_t>(frame.atoms.type[i]);
    std::memcpy(p, &typ, sizeof(typ));
    p += sizeof(typ);
    const std::int64_t tag = frame.atoms.tag[i];
    std::memcpy(p, &tag, sizeof(tag));
    p += sizeof(tag);
    const std::int8_t img[3] = {
        frame.atoms.image(i, 0), frame.atoms.image(i, 1), frame.atoms.image(i, 2),
    };
    std::memcpy(p, img, sizeof(img));
    return offset + kAtomPayloadBytes;
}

// Unpack one atom's per-atom fields into local frame slot `slot`. Returns
// the source offset after reading `kAtomPayloadBytes` bytes.
std::size_t unpack_atom(const std::vector<char>& buf, std::size_t offset,
                         Frame& frame, int slot) {
    const char* p = buf.data() + offset;
    double dvals[kPayloadDoubles];
    std::memcpy(dvals, p, sizeof(dvals));
    p += sizeof(dvals);
    frame.atoms.x[slot]    = dvals[0];
    frame.atoms.y[slot]    = dvals[1];
    frame.atoms.z[slot]    = dvals[2];
    frame.atoms.vx[slot]   = dvals[3];
    frame.atoms.vy[slot]   = dvals[4];
    frame.atoms.vz[slot]   = dvals[5];
    frame.atoms.fx[slot]   = dvals[6];
    frame.atoms.fy[slot]   = dvals[7];
    frame.atoms.fz[slot]   = dvals[8];
    frame.atoms.mass[slot] = dvals[9];
    std::int32_t typ;
    std::memcpy(&typ, p, sizeof(typ));
    p += sizeof(typ);
    frame.atoms.type[slot] = static_cast<int>(typ);
    std::int64_t tag;
    std::memcpy(&tag, p, sizeof(tag));
    p += sizeof(tag);
    frame.atoms.tag[slot] = tag;
    std::int8_t img[3];
    std::memcpy(img, p, sizeof(img));
    frame.atoms.image(slot, 0) = img[0];
    frame.atoms.image(slot, 1) = img[1];
    frame.atoms.image(slot, 2) = img[2];
    return offset + kAtomPayloadBytes;
}

// Move owner row `oi` into slot `ni` (no-op when `oi == ni`).
void move_row(Frame& frame, int ni, int oi) {
    if (ni == oi) return;
    frame.atoms.x[ni]    = frame.atoms.x[oi];
    frame.atoms.y[ni]    = frame.atoms.y[oi];
    frame.atoms.z[ni]    = frame.atoms.z[oi];
    frame.atoms.vx[ni]   = frame.atoms.vx[oi];
    frame.atoms.vy[ni]   = frame.atoms.vy[oi];
    frame.atoms.vz[ni]   = frame.atoms.vz[oi];
    frame.atoms.fx[ni]   = frame.atoms.fx[oi];
    frame.atoms.fy[ni]   = frame.atoms.fy[oi];
    frame.atoms.fz[ni]   = frame.atoms.fz[oi];
    frame.atoms.mass[ni] = frame.atoms.mass[oi];
    frame.atoms.type[ni] = frame.atoms.type[oi];
    frame.atoms.tag[ni]  = frame.atoms.tag[oi];
    frame.atoms.image(ni, 0) = frame.atoms.image(oi, 0);
    frame.atoms.image(ni, 1) = frame.atoms.image(oi, 1);
    frame.atoms.image(ni, 2) = frame.atoms.image(oi, 2);
}

} // namespace

void AtomMigration::launch(const cpu::LaunchContext& ctx) {
    if (ctx.runtime.mpi_size <= 1) return;
    auto* dd = ctx.runtime.dd;
    if (!dd || dd->cart_comm == MPI_COMM_NULL) return;

    Frame& frame = ctx.frame;
    const double Lx = frame.box(0, 0);
    const double Ly = frame.box(1, 1);
    const double Lz = frame.box(2, 2);

    int my_rank = 0;
    MPI_Comm_rank(dd->cart_comm, &my_rank);
    const int world_size = dd->dims[0] * dd->dims[1] * dd->dims[2];

    const bool first_launch =
        last_observed_generation_ == ~std::uint64_t{0};
    const bool generation_changed =
        !first_launch && dd->neighbor_generation != last_observed_generation_;
    if ((ctx.step % every_) != 0 && !generation_changed) {
        return;
    }

    // Classify each owner. `keep_indices` holds the local rows to compact
    // down; `migrate_indices_by_rank[r]` holds rows we'll ship to rank `r`.
    std::vector<int> keep_indices;
    keep_indices.reserve(static_cast<std::size_t>(dd->n_local));
    std::vector<std::vector<int>> migrate_indices_by_rank(static_cast<std::size_t>(world_size));

    int n_to_migrate = 0;
    for (int i = 0; i < dd->n_local; ++i) {
        const int dest = dest_rank_for(*dd, frame.atoms.x[i], frame.atoms.y[i],
                                       frame.atoms.z[i], Lx, Ly, Lz);
        if (dest == my_rank) {
            keep_indices.push_back(i);
        } else {
            migrate_indices_by_rank[static_cast<std::size_t>(dest)].push_back(i);
            ++n_to_migrate;
        }
    }

    // Early exit: if NO rank has anything to migrate, skip all collective
    // work. We still need a global agreement on this -- a single rank with
    // a migration would deadlock everyone if others bailed out early.
    int any_local = (n_to_migrate > 0) ? 1 : 0;
    int any_global = 0;
    MPI_Allreduce(&any_local, &any_global, 1, MPI_INT, MPI_LOR,
                  dd->cart_comm);
    if (any_global == 0) {
        last_observed_generation_ = dd->neighbor_generation;
        return;
    }

    // Pack per-rank send buffers into a single contiguous byte buffer for
    // MPI_Alltoallv (which needs flat arrays + per-rank counts/displs).
    std::vector<int> send_counts(static_cast<std::size_t>(world_size), 0);
    std::vector<int> send_displs(static_cast<std::size_t>(world_size), 0);
    std::vector<int> recv_counts(static_cast<std::size_t>(world_size), 0);
    std::vector<int> recv_displs(static_cast<std::size_t>(world_size), 0);

    std::size_t total_send_bytes = 0;
    for (int r = 0; r < world_size; ++r) {
        const std::size_t n = migrate_indices_by_rank[static_cast<std::size_t>(r)].size();
        send_displs[static_cast<std::size_t>(r)] = static_cast<int>(total_send_bytes);
        send_counts[static_cast<std::size_t>(r)] = static_cast<int>(n * kAtomPayloadBytes);
        total_send_bytes += n * kAtomPayloadBytes;
    }

    std::vector<char> send_buf(total_send_bytes);
    std::size_t off = 0;
    for (int r = 0; r < world_size; ++r) {
        for (int idx : migrate_indices_by_rank[static_cast<std::size_t>(r)]) {
            off = pack_atom(send_buf, off, frame, idx);
        }
    }

    // Exchange counts (in bytes) so every rank knows how much to receive
    // from every other rank.
    MPI_Alltoall(send_counts.data(), 1, MPI_INT,
                 recv_counts.data(), 1, MPI_INT, dd->cart_comm);

    std::size_t total_recv_bytes = 0;
    for (int r = 0; r < world_size; ++r) {
        recv_displs[static_cast<std::size_t>(r)] = static_cast<int>(total_recv_bytes);
        total_recv_bytes += static_cast<std::size_t>(recv_counts[static_cast<std::size_t>(r)]);
    }

    std::vector<char> recv_buf(total_recv_bytes);

    MPI_Alltoallv(send_buf.data(), send_counts.data(), send_displs.data(), MPI_BYTE,
                  recv_buf.data(), recv_counts.data(), recv_displs.data(), MPI_BYTE,
                  dd->cart_comm);

    // Compact kept owners to the front of every per-atom array.
    const int new_n_kept = static_cast<int>(keep_indices.size());
    for (int ni = 0; ni < new_n_kept; ++ni) {
        const int oi = keep_indices[static_cast<std::size_t>(ni)];
        if (oi != ni) move_row(frame, ni, oi);
    }

    // Append received atoms to the local owner range.
    const int n_received = static_cast<int>(total_recv_bytes / kAtomPayloadBytes);
    const int new_total  = new_n_kept + n_received;
    if (new_total > dd->capacity) {
        // Phase 0.2: never grow allocation in launch(). Capacity is fixed at
        // partition time. If a real workload trips this, the user must
        // pre-size the frame with enough slack (Phase 3 will move this to
        // the driver-level rebuild_frame() hook).
        // Fall through and clamp to capacity to avoid OOB; downstream tests
        // will report the dropped atoms as a global count mismatch.
    }
    const int n_appendable = std::min(n_received, dd->capacity - new_n_kept);
    std::size_t r_off = 0;
    for (int k = 0; k < n_appendable; ++k) {
        r_off = unpack_atom(recv_buf, r_off, frame, new_n_kept + k);
    }

    dd->n_local = new_n_kept + n_appendable;
    dd->n_ghost = 0;  // any ghost slots from the previous step are now stale;
                      // HaloExchangeForward will rewrite them this step.
    ++dd->neighbor_generation;

    // Rebuild tag → local index map. Bonded indices that key off `tag_to_local`
    // re-resolve through this map on the next neighbor / bonded build.
    dd->tag_to_local.clear();
    dd->tag_to_local.reserve(static_cast<std::size_t>(dd->n_local));
    for (int i = 0; i < dd->n_local; ++i) {
        dd->tag_to_local.emplace(frame.atoms.tag[i], i);
    }

    last_observed_generation_ = dd->neighbor_generation;
}

#else // !ATV_HAS_MPI

void AtomMigration::launch(const cpu::LaunchContext& /*ctx*/) {}

#endif // ATV_HAS_MPI

} // namespace mpi
} // namespace md
} // namespace cpu
} // namespace ATV
