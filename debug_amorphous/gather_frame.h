/// @file debug_amorphous/gather_frame.h
/// @brief Inter-stage all-gather for multi-rank Frame state.
///
/// After a multi-rank `cpu::md::MD::run()`, every rank holds only the
/// atoms in its slab inside `frame.atoms[0 : n_local)`; the remainder of
/// the array is stale (mostly old ghost slots). When the user then
/// constructs a SECOND driver and calls `run()` again, that driver's
/// `partition_atoms()` sees `n_local == 0` and rescans the whole array
/// (including the stale tail) -- with garbage going into the partition.
///
/// `gather_full_frame_inplace()` fixes the handoff: it `MPI_Allgatherv`s
/// every rank's owned atoms back into `[0, n_total)` on every rank, in
/// canonical (sorted-by-global-tag) order, and resets the `dd->n_local`
/// to the global count so the next partition starts from a clean,
/// fully-replicated state.
///
/// Single-rank: no-op (every atom is already on this rank).
#pragma once

#ifdef ATV_HAS_MPI
#include <mpi.h>
#endif

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <vector>

#include "atomiverse/cpu/frame.h"
#include "atomiverse/cpu/md/mpi/domain.h"

namespace ATV {
namespace debug_amorphous {

/// @brief Capture-n_local helper: a Component that records `dd->n_local` to
/// a user-supplied int at every end_of_step, leaving the int holding the
/// rank's owned-count after the driver returns.
class CaptureNLocal : public cpu::ComponentBase {
public:
    explicit CaptureNLocal(int* sink) : sink_(sink) {}
    std::string name() const override { return "capture_n_local"; }
    std::vector<DefaultBinding> default_bindings() const override {
        return {{static_cast<int>(ATV::md::Phase::end_of_step), 30}};
    }
    void launch(const cpu::LaunchContext& ctx) override {
        const auto* dd = ctx.runtime.dd;
        *sink_ = (dd && ctx.runtime.mpi_size > 1)
                     ? dd->n_local : ctx.frame.atoms.size();
    }
private:
    int* sink_;
};

inline void gather_full_frame_inplace(cpu::Frame& frame,
                                      int n_local,
                                      int world_size,
                                      MPI_Comm comm) {
#ifdef ATV_HAS_MPI
    if (world_size <= 1) return;

    // Per-atom payload: x, y, z, vx, vy, vz, fx, fy, fz, mass (10 doubles)
    // + tag (i64 in a double slot via reinterpret) + type (i32 in a double).
    // Use MPI_DOUBLE so the reduction is a single Allgatherv.
    constexpr int kPerAtomDoubles = 13;

    const int capacity = frame.atoms.size();

    // Pack our owned atoms.
    std::vector<double> sendbuf(static_cast<std::size_t>(n_local) * kPerAtomDoubles);
    for (int i = 0; i < n_local; ++i) {
        double* p = sendbuf.data() + static_cast<std::size_t>(i) * kPerAtomDoubles;
        p[0]  = frame.atoms.x[i];
        p[1]  = frame.atoms.y[i];
        p[2]  = frame.atoms.z[i];
        p[3]  = frame.atoms.vx[i];
        p[4]  = frame.atoms.vy[i];
        p[5]  = frame.atoms.vz[i];
        p[6]  = frame.atoms.fx[i];
        p[7]  = frame.atoms.fy[i];
        p[8]  = frame.atoms.fz[i];
        p[9]  = frame.atoms.mass[i];
        // Stash tag (i64) in a double slot via memcpy (preserves bits).
        const std::int64_t tag = frame.atoms.tag[i];
        std::memcpy(&p[10], &tag, sizeof(double));
        // Stash type (i32) and image (3 i8) packed into one double.
        std::int32_t typ_and_image[2] = {
            static_cast<std::int32_t>(frame.atoms.type[i]),
            (static_cast<std::int32_t>(frame.atoms.image(i, 0)) & 0xFF)
            | ((static_cast<std::int32_t>(frame.atoms.image(i, 1)) & 0xFF) << 8)
            | ((static_cast<std::int32_t>(frame.atoms.image(i, 2)) & 0xFF) << 16)
        };
        std::memcpy(&p[11], typ_and_image, sizeof(typ_and_image));
        // Last slot holds nothing right now; reserve for future extension.
        p[12] = 0.0;
    }

    // Discover per-rank send sizes.
    std::vector<int> recv_counts(static_cast<std::size_t>(world_size));
    int my_count = n_local * kPerAtomDoubles;
    MPI_Allgather(&my_count, 1, MPI_INT,
                  recv_counts.data(), 1, MPI_INT,
                  comm);
    std::vector<int> recv_displs(static_cast<std::size_t>(world_size), 0);
    for (int r = 1; r < world_size; ++r)
        recv_displs[r] = recv_displs[r - 1] + recv_counts[r - 1];

    int total_doubles = 0;
    for (int r = 0; r < world_size; ++r) total_doubles += recv_counts[r];

    std::vector<double> recvbuf(static_cast<std::size_t>(total_doubles));
    MPI_Allgatherv(sendbuf.data(), my_count, MPI_DOUBLE,
                   recvbuf.data(), recv_counts.data(), recv_displs.data(),
                   MPI_DOUBLE, comm);

    const int n_total = total_doubles / kPerAtomDoubles;
    if (n_total > capacity) {
        // Should never happen unless the frame was resized between stages.
        std::fprintf(stderr,
            "[gather_full_frame] n_total=%d > capacity=%d -- aborting.\n",
            n_total, capacity);
        return;
    }

    // Sort entries by global tag so the gathered layout is deterministic
    // across runs (and across rank counts).
    std::vector<int> order(static_cast<std::size_t>(n_total));
    for (int i = 0; i < n_total; ++i) order[i] = i;
    auto tag_of = [&](int slot) {
        std::int64_t t;
        std::memcpy(&t, &recvbuf[slot * kPerAtomDoubles + 10], sizeof(t));
        return t;
    };
    std::sort(order.begin(), order.end(),
              [&](int a, int b) { return tag_of(a) < tag_of(b); });

    for (int dst = 0; dst < n_total; ++dst) {
        const int src = order[dst];
        const double* p =
            recvbuf.data() + static_cast<std::size_t>(src) * kPerAtomDoubles;
        frame.atoms.x[dst]    = p[0];
        frame.atoms.y[dst]    = p[1];
        frame.atoms.z[dst]    = p[2];
        frame.atoms.vx[dst]   = p[3];
        frame.atoms.vy[dst]   = p[4];
        frame.atoms.vz[dst]   = p[5];
        frame.atoms.fx[dst]   = p[6];
        frame.atoms.fy[dst]   = p[7];
        frame.atoms.fz[dst]   = p[8];
        frame.atoms.mass[dst] = p[9];
        std::int64_t tag;
        std::memcpy(&tag, &p[10], sizeof(tag));
        frame.atoms.tag[dst]  = tag;
        std::int32_t typ_and_image[2];
        std::memcpy(typ_and_image, &p[11], sizeof(typ_and_image));
        frame.atoms.type[dst] = typ_and_image[0];
        frame.atoms.image(dst, 0) = static_cast<std::int8_t>(
            (typ_and_image[1] & 0xFF) << 24 >> 24);
        frame.atoms.image(dst, 1) = static_cast<std::int8_t>(
            ((typ_and_image[1] >> 8) & 0xFF) << 24 >> 24);
        frame.atoms.image(dst, 2) = static_cast<std::int8_t>(
            ((typ_and_image[1] >> 16) & 0xFF) << 24 >> 24);
    }
#else
    (void)frame; (void)n_local; (void)world_size; (void)comm;
#endif
}

} // namespace debug_amorphous
} // namespace ATV
