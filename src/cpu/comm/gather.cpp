/// @file src/cpu/md/mpi/gather.cpp
/// @brief `gather_to_root` -- collect distributed frame owners to one rank.
///
/// Implementation uses one `MPI_Gather` for per-rank owner counts plus a
/// batch of `MPI_Gatherv` calls for the per-atom fields. After gather the
/// root sorts atoms by `tag` so the I/O order is decomposition-independent.
///
/// @see .claude/specs/parallel-md-engine.md  § A.9

#include "atomiverse/cpu/comm/gather.h"

#include <algorithm>
#include <cstdint>
#include <numeric>
#include <vector>

#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/comm/domain.h"

#ifdef ATV_HAS_MPI
#include <mpi.h>
#endif

namespace ATV {
namespace cpu {
namespace md {
namespace mpi {

namespace {

void copy_owners_into(GatheredFrame& out, const Frame& frame, int n_local) {
    out.n_total = n_local;
    out.box = frame.box.cell();
    // Hoisted once at this launch-boundary function's entry (single-rank
    // fallback path — no loop over atoms here, one assign() per field).
    const StorageReal* __restrict sp_x = frame.span(kAtom_x);
    const StorageReal* __restrict sp_y = frame.span(kAtom_y);
    const StorageReal* __restrict sp_z = frame.span(kAtom_z);
    const StorageReal* __restrict sp_vx = frame.span(kAtom_vx);
    const StorageReal* __restrict sp_vy = frame.span(kAtom_vy);
    const StorageReal* __restrict sp_vz = frame.span(kAtom_vz);
    const AccReal* __restrict sp_fx = frame.span(kAtom_fx);
    const AccReal* __restrict sp_fy = frame.span(kAtom_fy);
    const AccReal* __restrict sp_fz = frame.span(kAtom_fz);
    const StorageReal* __restrict sp_mass = frame.span(kAtom_mass);
    const int* __restrict sp_type = frame.span(kAtom_type);
    const TagInt* __restrict sp_tag = frame.span(kAtom_tag);
    const std::int8_t* __restrict sp_image = frame.span(kAtom_image);
    out.x.assign(sp_x, sp_x + n_local);
    out.y.assign(sp_y, sp_y + n_local);
    out.z.assign(sp_z, sp_z + n_local);
    out.vx.assign(sp_vx, sp_vx + n_local);
    out.vy.assign(sp_vy, sp_vy + n_local);
    out.vz.assign(sp_vz, sp_vz + n_local);
    out.fx.assign(sp_fx, sp_fx + n_local);
    out.fy.assign(sp_fy, sp_fy + n_local);
    out.fz.assign(sp_fz, sp_fz + n_local);
    out.mass.assign(sp_mass, sp_mass + n_local);
    out.type.assign(sp_type, sp_type + n_local);
    out.tag.assign(sp_tag, sp_tag + n_local);
    out.image.assign(sp_image, sp_image + 3 * n_local);
}

// Reorder all per-atom fields in `g` so that `g.tag` is ascending. Stable so
// equal tags (which would indicate a bug upstream) preserve their order.
void sort_by_tag(GatheredFrame& g) {
    const int n = g.n_total;
    if (n <= 1)
        return;
    std::vector<int> perm(static_cast<std::size_t>(n));
    std::iota(perm.begin(), perm.end(), 0);
    std::stable_sort(perm.begin(), perm.end(),
                     [&](int a, int b) { return g.tag[a] < g.tag[b]; });

    auto apply_perm_double = [&](std::vector<double>& v) {
        std::vector<double> tmp(v.size());
        for (int i = 0; i < n; ++i)
            tmp[i] = v[perm[i]];
        v.swap(tmp);
    };
    auto apply_perm_int = [&](std::vector<int>& v) {
        std::vector<int> tmp(v.size());
        for (int i = 0; i < n; ++i)
            tmp[i] = v[perm[i]];
        v.swap(tmp);
    };
    auto apply_perm_int64 = [&](std::vector<std::int64_t>& v) {
        std::vector<std::int64_t> tmp(v.size());
        for (int i = 0; i < n; ++i)
            tmp[i] = v[perm[i]];
        v.swap(tmp);
    };

    apply_perm_double(g.x);
    apply_perm_double(g.y);
    apply_perm_double(g.z);
    apply_perm_double(g.vx);
    apply_perm_double(g.vy);
    apply_perm_double(g.vz);
    apply_perm_double(g.fx);
    apply_perm_double(g.fy);
    apply_perm_double(g.fz);
    apply_perm_double(g.mass);
    apply_perm_int(g.type);
    apply_perm_int64(g.tag);

    // image: 3 entries per atom in flat row-major layout.
    std::vector<std::int8_t> tmp_img(g.image.size());
    for (int i = 0; i < n; ++i) {
        const int src = perm[i];
        tmp_img[3 * i + 0] = g.image[3 * src + 0];
        tmp_img[3 * i + 1] = g.image[3 * src + 1];
        tmp_img[3 * i + 2] = g.image[3 * src + 2];
    }
    g.image.swap(tmp_img);
}

} // namespace

#ifdef ATV_HAS_MPI

GatheredFrame gather_to_root(const cpu::Frame& frame,
                             const DomainDecomposition& dd, int root) {
    GatheredFrame out;
    out.root_rank = root;
    out.box = frame.box.cell();

    if (dd.cart_comm == MPI_COMM_NULL) {
        // Single-rank fallback: direct copy of owners.
        copy_owners_into(out, frame, dd.n_local);
        return out;
    }

    int world_size = 0;
    int my_rank = 0;
    MPI_Comm_size(dd.cart_comm, &world_size);
    MPI_Comm_rank(dd.cart_comm, &my_rank);

    if (world_size <= 1) {
        copy_owners_into(out, frame, dd.n_local);
        return out;
    }

    const int n_local = dd.n_local;

    // Step 1: every rank ships its owner count to root.
    std::vector<int> recv_counts;
    if (my_rank == root)
        recv_counts.assign(static_cast<std::size_t>(world_size), 0);
    MPI_Gather(&n_local, 1, MPI_INT, recv_counts.data(), 1, MPI_INT, root,
               dd.cart_comm);

    // Compute total + per-rank displacements on root.
    std::vector<int> displs;
    int n_total = 0;
    if (my_rank == root) {
        displs.assign(static_cast<std::size_t>(world_size), 0);
        for (int r = 0; r < world_size; ++r) {
            displs[static_cast<std::size_t>(r)] = n_total;
            n_total += recv_counts[static_cast<std::size_t>(r)];
        }
    }

    // Broadcast n_total so non-root ranks know the global atom count
    // (useful for ergonomic checks at call sites without re-reducing).
    MPI_Bcast(&n_total, 1, MPI_INT, root, dd.cart_comm);
    out.n_total = n_total;

    // Allocate root buffers.
    if (my_rank == root) {
        out.x.assign(static_cast<std::size_t>(n_total), 0.0);
        out.y.assign(static_cast<std::size_t>(n_total), 0.0);
        out.z.assign(static_cast<std::size_t>(n_total), 0.0);
        out.vx.assign(static_cast<std::size_t>(n_total), 0.0);
        out.vy.assign(static_cast<std::size_t>(n_total), 0.0);
        out.vz.assign(static_cast<std::size_t>(n_total), 0.0);
        out.fx.assign(static_cast<std::size_t>(n_total), 0.0);
        out.fy.assign(static_cast<std::size_t>(n_total), 0.0);
        out.fz.assign(static_cast<std::size_t>(n_total), 0.0);
        out.mass.assign(static_cast<std::size_t>(n_total), 0.0);
        out.type.assign(static_cast<std::size_t>(n_total), 0);
        out.tag.assign(static_cast<std::size_t>(n_total), 0);
        out.image.assign(static_cast<std::size_t>(3 * n_total), 0);
    }

    // Helper: gather one per-atom real field. `sendbuf` may be single-precision
    // storage (StorageReal); widen to a double staging buffer so the gathered
    // output is double regardless of the storage tier. Gather is a low-frequency
    // output path (rank-0 collection), so the copy is negligible.
    auto gather_double = [&](const auto* sendbuf, double* recvbuf) {
        std::vector<double> send(sendbuf, sendbuf + n_local);
        MPI_Gatherv(send.data(), n_local, MPI_DOUBLE, recvbuf,
                    my_rank == root ? recv_counts.data() : nullptr,
                    my_rank == root ? displs.data() : nullptr, MPI_DOUBLE, root,
                    dd.cart_comm);
    };

    // Hoisted once at the launch boundary; every gather call below reads
    // through these instead of re-resolving per field.
    const StorageReal* __restrict sp_x = frame.span(kAtom_x);
    const StorageReal* __restrict sp_y = frame.span(kAtom_y);
    const StorageReal* __restrict sp_z = frame.span(kAtom_z);
    const StorageReal* __restrict sp_vx = frame.span(kAtom_vx);
    const StorageReal* __restrict sp_vy = frame.span(kAtom_vy);
    const StorageReal* __restrict sp_vz = frame.span(kAtom_vz);
    const AccReal* __restrict sp_fx = frame.span(kAtom_fx);
    const AccReal* __restrict sp_fy = frame.span(kAtom_fy);
    const AccReal* __restrict sp_fz = frame.span(kAtom_fz);
    const StorageReal* __restrict sp_mass = frame.span(kAtom_mass);
    const int* __restrict sp_type = frame.span(kAtom_type);
    const TagInt* __restrict sp_tag = frame.span(kAtom_tag);
    const std::int8_t* __restrict sp_image = frame.span(kAtom_image);

    gather_double(sp_x, my_rank == root ? out.x.data() : nullptr);
    gather_double(sp_y, my_rank == root ? out.y.data() : nullptr);
    gather_double(sp_z, my_rank == root ? out.z.data() : nullptr);
    gather_double(sp_vx, my_rank == root ? out.vx.data() : nullptr);
    gather_double(sp_vy, my_rank == root ? out.vy.data() : nullptr);
    gather_double(sp_vz, my_rank == root ? out.vz.data() : nullptr);
    gather_double(sp_fx, my_rank == root ? out.fx.data() : nullptr);
    gather_double(sp_fy, my_rank == root ? out.fy.data() : nullptr);
    gather_double(sp_fz, my_rank == root ? out.fz.data() : nullptr);
    gather_double(sp_mass, my_rank == root ? out.mass.data() : nullptr);

    MPI_Gatherv(sp_type, n_local, MPI_INT,
                my_rank == root ? out.type.data() : nullptr,
                my_rank == root ? recv_counts.data() : nullptr,
                my_rank == root ? displs.data() : nullptr, MPI_INT, root,
                dd.cart_comm);

    MPI_Gatherv(sp_tag, n_local, MPI_INT64_T,
                my_rank == root ? out.tag.data() : nullptr,
                my_rank == root ? recv_counts.data() : nullptr,
                my_rank == root ? displs.data() : nullptr, MPI_INT64_T, root,
                dd.cart_comm);

    // image is shape [n, 3] -- ship 3 bytes per atom.
    std::vector<int> img_counts;
    std::vector<int> img_displs;
    if (my_rank == root) {
        img_counts.assign(static_cast<std::size_t>(world_size), 0);
        img_displs.assign(static_cast<std::size_t>(world_size), 0);
        for (int r = 0; r < world_size; ++r) {
            img_counts[static_cast<std::size_t>(r)] =
                3 * recv_counts[static_cast<std::size_t>(r)];
            img_displs[static_cast<std::size_t>(r)] =
                3 * displs[static_cast<std::size_t>(r)];
        }
    }
    MPI_Gatherv(sp_image, 3 * n_local, MPI_INT8_T,
                my_rank == root ? out.image.data() : nullptr,
                my_rank == root ? img_counts.data() : nullptr,
                my_rank == root ? img_displs.data() : nullptr, MPI_INT8_T, root,
                dd.cart_comm);

    if (my_rank == root) {
        sort_by_tag(out);
    }
    return out;
}

#else // !ATV_HAS_MPI

GatheredFrame gather_to_root(const cpu::Frame& frame,
                             const DomainDecomposition& dd, int root) {
    GatheredFrame out;
    out.root_rank = root;
    copy_owners_into(out, frame, dd.n_local);
    return out;
}

#endif // ATV_HAS_MPI

cpu::Frame to_frame(const GatheredFrame& g) {
    const int n = static_cast<int>(g.x.size());
    cpu::Frame f(n);
    f.box.set_cell(g.box);
    // Hoisted after `f`'s construction (dominates the loop below), never
    // before -- f is constructed here, not passed in.
    StorageReal* __restrict sp_x = f.span(kAtom_x);
    StorageReal* __restrict sp_y = f.span(kAtom_y);
    StorageReal* __restrict sp_z = f.span(kAtom_z);
    StorageReal* __restrict sp_vx = f.span(kAtom_vx);
    StorageReal* __restrict sp_vy = f.span(kAtom_vy);
    StorageReal* __restrict sp_vz = f.span(kAtom_vz);
    AccReal* __restrict sp_fx = f.span(kAtom_fx);
    AccReal* __restrict sp_fy = f.span(kAtom_fy);
    AccReal* __restrict sp_fz = f.span(kAtom_fz);
    StorageReal* __restrict sp_mass = f.span(kAtom_mass);
    int* __restrict sp_type = f.span(kAtom_type);
    TagInt* __restrict sp_tag = f.span(kAtom_tag);
    for (int i = 0; i < n; ++i) {
        sp_x[i] = static_cast<StorageReal>(g.x[i]);
        sp_y[i] = static_cast<StorageReal>(g.y[i]);
        sp_z[i] = static_cast<StorageReal>(g.z[i]);
        sp_vx[i] = static_cast<StorageReal>(g.vx[i]);
        sp_vy[i] = static_cast<StorageReal>(g.vy[i]);
        sp_vz[i] = static_cast<StorageReal>(g.vz[i]);
        sp_fx[i] = static_cast<AccReal>(g.fx[i]);
        sp_fy[i] = static_cast<AccReal>(g.fy[i]);
        sp_fz[i] = static_cast<AccReal>(g.fz[i]);
        sp_mass[i] = static_cast<StorageReal>(g.mass[i]);
        sp_type[i] = g.type[i];
        sp_tag[i] = static_cast<TagInt>(g.tag[i]);
    }
    return f;
}

} // namespace mpi
} // namespace md
} // namespace cpu
} // namespace ATV
