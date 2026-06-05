/// @file src/cpu/comm/atom_migration.cpp
/// @brief `AtomMigration` -- MPI_Alltoallv owner migration implementation.
///
/// For each owner whose post-integration position lies outside this rank's
/// spatial slab, the atom (and all its per-atom fields) is shipped to the
/// rank whose slab contains the new position. Implemented as a single
/// `MPI_Alltoallv` over byte buffers so every world rank participates and
/// either packs/unpacks zero bytes or some multiple of `kAtomPayloadBytes`.
///
/// Bonded-tuple migration is **not** handled here; that is tracked under
/// See msu-pe-validation.md (MSU-PE multi-rank validation).

#include "atomiverse/cpu/comm/atom_migration.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <numeric>
#include <vector>

#include "atomiverse/cpu/comm/atom_pack.h"
#include "atomiverse/cpu/frame.h"
#include "atomiverse/cpu/md/mpi/domain.h"

#ifdef ATV_HAS_MPI
#include <mpi.h>
#endif

namespace ATV {
namespace cpu {
namespace comm {

#ifdef ATV_HAS_MPI

namespace {

// Per-atom payload layout: worst-case fixed-size dump from
// `comm::pack_exchange`. Carries all SoA fields plus per-atom topology
// (num_bond + bond_atom + bond_type + angles + dihedrals + nspecial +
// special). Atoms with fewer bonds than the max simply leave the
// trailing slots zero-padded — receiver reads `num_bond` first and
// only consumes that many partner entries.
constexpr int kAtomPayloadDoubles = comm::kMaxExchangeDoubles;
constexpr int kAtomPayloadBytes = sizeof(double) * kAtomPayloadDoubles;

inline int rank_along_dim(double q, double L, int dims_d) {
    if (dims_d <= 1)
        return 0;
    const int r = static_cast<int>(q * static_cast<double>(dims_d) / L);
    return std::max(0, std::min(dims_d - 1, r));
}

// Locate the slab index `k` (in [0, dims_d)) whose [cuts[k], cuts[k+1])
// contains the box-length fraction `frac` (already wrapped into [0, 1)).
// `cuts` has size `dims_d + 1` with `cuts.front() == 0` and `cuts.back() == 1`.
inline int slab_along_dim_from_cuts(double frac,
                                    const std::vector<double>& cuts,
                                    int dims_d) {
    if (dims_d <= 1)
        return 0;
    // upper_bound gives the first cuts[i] strictly greater than frac;
    // its predecessor is the slab index. Clamp against FP rounding at the
    // top edge.
    auto it = std::upper_bound(cuts.begin(), cuts.end(), frac);
    int k = static_cast<int>(std::distance(cuts.begin(), it)) - 1;
    if (k < 0)
        k = 0;
    if (k > dims_d - 1)
        k = dims_d - 1;
    return k;
}

// Compute the rank that owns position (x, y, z) under the current DD.
//
// When `dd.cuts[d]` is empty along axis d, the uniform `coord = floor(q*dims/L)`
// formula applies. When non-empty (set by `RCBBalancer`), the destination is
// derived from the box-length fraction via binary search through the cut
// planes -- which guarantees the same rank assignment that
// `md::mpi::DomainDecomposition::compute_sublo_subhi` would produce on the destination
// side.
inline int dest_rank_for(const md::mpi::DomainDecomposition& dd, double x, double y,
                         double z, const cpu::Box& box) {
    auto axis_coord = [&](int d, double q) {
        const double L = box.length(d);
        const double wrapped = box.wrap_axis(q, d);
        if (static_cast<int>(dd.cuts[d].size()) == dd.dims[d] + 1) {
            const double frac = (L > 0.0) ? (wrapped / L) : 0.0;
            return slab_along_dim_from_cuts(frac, dd.cuts[d], dd.dims[d]);
        }
        return rank_along_dim(wrapped, L, dd.dims[d]);
    };
    const int cx = axis_coord(0, x);
    const int cy = axis_coord(1, y);
    const int cz = axis_coord(2, z);
    return cx * (dd.dims[1] * dd.dims[2]) + cy * dd.dims[2] + cz;
}

// Pack one atom (full SoA + per-atom topology) into `buf` at `offset`.
// Delegates to `comm::pack_exchange` for the payload layout; pads with
// zeros up to the worst-case `kAtomPayloadDoubles` so every atom
// occupies the same byte count on the wire (atom_migration's
// `MPI_Alltoallv` needs uniform per-atom strides).
std::size_t pack_atom(std::vector<char>& buf, std::size_t offset,
                      const Frame& frame, int i) {
    // pack_exchange writes between min and kMaxExchangeDoubles doubles.
    double* dst = reinterpret_cast<double*>(buf.data() + offset);
    // Zero-pad the slot first so the trailing topology entries don't
    // carry stale data.
    std::memset(dst, 0, kAtomPayloadBytes);
    comm::pack_exchange(frame.atoms, i, dst);
    return offset + kAtomPayloadBytes;
}

std::size_t unpack_atom(const std::vector<char>& buf, std::size_t offset,
                        Frame& frame, int slot) {
    const double* src =
        reinterpret_cast<const double*>(buf.data() + offset);
    comm::unpack_exchange(frame.atoms, slot, src);
    return offset + kAtomPayloadBytes;
}

// Move owner row `oi` into slot `ni` (no-op when `oi == ni`). Handles
// both SoA fields and per-atom topology.
void move_row(Frame& frame, int ni, int oi) {
    if (ni == oi)
        return;
    auto& a = frame.atoms;
    a.x[ni] = a.x[oi];
    a.y[ni] = a.y[oi];
    a.z[ni] = a.z[oi];
    a.vx[ni] = a.vx[oi];
    a.vy[ni] = a.vy[oi];
    a.vz[ni] = a.vz[oi];
    a.fx[ni] = a.fx[oi];
    a.fy[ni] = a.fy[oi];
    a.fz[ni] = a.fz[oi];
    a.mass[ni] = a.mass[oi];
    a.charge[ni] = a.charge[oi];
    a.type[ni] = a.type[oi];
    a.tag[ni] = a.tag[oi];
    a.image(ni, 0) = a.image(oi, 0);
    a.image(ni, 1) = a.image(oi, 1);
    a.image(ni, 2) = a.image(oi, 2);

    // Per-atom topology — must follow the row.
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

} // namespace

void AtomMigration::launch(const cpu::LaunchContext& ctx) {
    if (ctx.runtime.mpi_size <= 1)
        return;
    auto* dd = ctx.runtime.dd;
    if (!dd || dd->cart_comm == MPI_COMM_NULL)
        return;

    Frame& frame = ctx.frame;

    int my_rank = 0;
    MPI_Comm_rank(dd->cart_comm, &my_rank);
    const int world_size = dd->dims[0] * dd->dims[1] * dd->dims[2];

    const bool first_launch = last_observed_generation_ == ~std::uint64_t{0};
    const bool generation_changed =
        !first_launch && dd->neighbor_generation != last_observed_generation_;
    if ((ctx.step % every_) != 0 && !generation_changed) {
        return;
    }

    // Classify each owner. `keep_indices` holds the local rows to compact
    // down; `migrate_indices_by_rank[r]` holds rows we'll ship to rank `r`.
    std::vector<int> keep_indices;
    keep_indices.reserve(static_cast<std::size_t>(dd->n_local));
    std::vector<std::vector<int>> migrate_indices_by_rank(
        static_cast<std::size_t>(world_size));

    int n_to_migrate = 0;
    for (int i = 0; i < dd->n_local; ++i) {
        const int dest = dest_rank_for(*dd, frame.atoms.x[i], frame.atoms.y[i],
                                       frame.atoms.z[i], frame.box);
        if (dest == my_rank) {
            keep_indices.push_back(i);
        } else {
            migrate_indices_by_rank[static_cast<std::size_t>(dest)].push_back(
                i);
            ++n_to_migrate;
        }
    }

    int any_local = (n_to_migrate > 0) ? 1 : 0;
    int any_global = 0;
    MPI_Allreduce(&any_local, &any_global, 1, MPI_INT, MPI_LOR, dd->cart_comm);
    if (any_global == 0) {
        last_observed_generation_ = dd->neighbor_generation;
        return;
    }

    std::vector<int> send_counts(static_cast<std::size_t>(world_size), 0);
    std::vector<int> send_displs(static_cast<std::size_t>(world_size), 0);
    std::vector<int> recv_counts(static_cast<std::size_t>(world_size), 0);
    std::vector<int> recv_displs(static_cast<std::size_t>(world_size), 0);

    std::size_t total_send_bytes = 0;
    for (int r = 0; r < world_size; ++r) {
        const std::size_t n =
            migrate_indices_by_rank[static_cast<std::size_t>(r)].size();
        send_displs[static_cast<std::size_t>(r)] =
            static_cast<int>(total_send_bytes);
        send_counts[static_cast<std::size_t>(r)] =
            static_cast<int>(n * kAtomPayloadBytes);
        total_send_bytes += n * kAtomPayloadBytes;
    }

    std::vector<char> send_buf(total_send_bytes);
    std::size_t off = 0;
    for (int r = 0; r < world_size; ++r) {
        for (int idx : migrate_indices_by_rank[static_cast<std::size_t>(r)]) {
            off = pack_atom(send_buf, off, frame, idx);
        }
    }

    MPI_Alltoall(send_counts.data(), 1, MPI_INT, recv_counts.data(), 1, MPI_INT,
                 dd->cart_comm);

    std::size_t total_recv_bytes = 0;
    for (int r = 0; r < world_size; ++r) {
        recv_displs[static_cast<std::size_t>(r)] =
            static_cast<int>(total_recv_bytes);
        total_recv_bytes +=
            static_cast<std::size_t>(recv_counts[static_cast<std::size_t>(r)]);
    }

    std::vector<char> recv_buf(total_recv_bytes);

    MPI_Alltoallv(send_buf.data(), send_counts.data(), send_displs.data(),
                  MPI_BYTE, recv_buf.data(), recv_counts.data(),
                  recv_displs.data(), MPI_BYTE, dd->cart_comm);

    const int new_n_kept = static_cast<int>(keep_indices.size());
    for (int ni = 0; ni < new_n_kept; ++ni) {
        const int oi = keep_indices[static_cast<std::size_t>(ni)];
        if (oi != ni)
            move_row(frame, ni, oi);
    }

    const int n_received =
        static_cast<int>(total_recv_bytes / kAtomPayloadBytes);
    // Grow on demand: ranks no longer hold the global-N backing store, so an
    // influx that exceeds capacity must enlarge the arrays — never silently
    // drop incoming owners (that would lose atoms and corrupt conservation).
    const int n_needed = new_n_kept + n_received;
    if (n_needed > frame.atoms.nmax) {
        frame.atoms.set_capacity(n_needed + n_needed / 2 + 1024);
        dd->capacity = frame.atoms.nmax;
    }
    std::size_t r_off = 0;
    for (int k = 0; k < n_received; ++k) {
        r_off = unpack_atom(recv_buf, r_off, frame, new_n_kept + k);
    }

    dd->n_local = new_n_kept + n_received;
    dd->n_ghost = 0;
    ++dd->neighbor_generation;

    // Rebuild tag → local index map. Bonded indices that key off `tag_to_local`
    // re-resolve through this map on the next neighbor / bonded build.
    dd->tag_to_local.clear();
    dd->tag_to_local.reserve(static_cast<std::size_t>(dd->n_local));
    for (int i = 0; i < dd->n_local; ++i) {
        dd->tag_to_local.emplace(frame.atoms.tag[i], i);
    }

    // Also clear `ghost_tags`. The previous step's ghost-tag log is now
    // stale, and worse, leaving it populated lets `HaloExchangeForward`'s
    // full-rebuild path `erase` those tags from the FRESH `tag_to_local`
    // we just rebuilt. Any atom that was a ghost on this rank
    // pre-migration AND just migrated IN as a now-owned atom would have
    // its owner entry silently removed, causing the bonded
    // `select_owned_tuples` filter to drop every bond owned by that atom
    // on every rank (no rank sees it as `owned`, slot < n_local). For
    // MSU PE that was ~15 bonds dropped per migration step, ~5 eV PE
    // loss per halo cycle, accumulating to a catastrophic ε̇ runaway by
    // step ~1000.
    dd->ghost_tags.clear();

    last_observed_generation_ = dd->neighbor_generation;
}

#else // !ATV_HAS_MPI

void AtomMigration::launch(const cpu::LaunchContext& /*ctx*/) {}

#endif // ATV_HAS_MPI

} // namespace comm
} // namespace cpu
} // namespace ATV
