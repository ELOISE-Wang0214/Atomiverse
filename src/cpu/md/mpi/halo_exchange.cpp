/// @file src/cpu/md/mpi/halo_exchange.cpp
/// @brief `HaloExchangeForward` / `HaloExchangeReverse` real implementation.
///
/// Phase 2: orthorhombic, 6 face neighbors, non-blocking `MPI_Isend/Irecv`.
/// Ghost slots are appended in face order `{-x, +x, -y, +y, -z, +z}` so
/// `Reverse` can mirror the layout via `dd->recv_counts[]`.

#include "atomiverse/cpu/md/mpi/halo_exchange.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <memory>
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

// ── Per-instance buffer storage (PIMPL) ──────────────────────────────
//
// `FaceBufs` holds only `std::vector<double>` arrays, so the type is
// usable regardless of whether MPI is linked in; the send-/recv-force
// arrays simply stay empty on the !MPI path. Keeping the full definition
// here (rather than under `#ifdef ATV_HAS_MPI`) lets the non-MPI stubs
// and the destructors in both classes see a complete type.

namespace detail {

struct FaceBufs {
    std::array<std::vector<double>, 6> send_pos;   // packed positions per face
    std::array<std::vector<double>, 6> recv_pos;   // received positions per face
    std::array<std::vector<double>, 6> send_force; // packed forces (Reverse)
    std::array<std::vector<double>, 6> recv_force; // received forces (Reverse)
};

} // namespace detail

// ── Special members (defined here so the PIMPL dtor sees complete FaceBufs) ──

HaloExchangeForward::HaloExchangeForward() = default;
HaloExchangeForward::HaloExchangeForward(double cutoff) : cutoff_(cutoff) {}
HaloExchangeForward::~HaloExchangeForward() = default;
HaloExchangeForward::HaloExchangeForward(HaloExchangeForward&&) noexcept = default;
HaloExchangeForward& HaloExchangeForward::operator=(HaloExchangeForward&&) noexcept = default;

HaloExchangeReverse::HaloExchangeReverse() = default;
HaloExchangeReverse::~HaloExchangeReverse() = default;
HaloExchangeReverse::HaloExchangeReverse(HaloExchangeReverse&&) noexcept = default;
HaloExchangeReverse& HaloExchangeReverse::operator=(HaloExchangeReverse&&) noexcept = default;

std::size_t HaloExchangeForward::send_buf_capacity_for_test() const {
    return bufs_ ? bufs_->send_pos[0].capacity() : 0;
}

// ── Topology-scan helpers (available with or without MPI) ────────────

namespace {

inline double pair_dist_sq(const Frame& f, int i, int j) {
    double dx = f.atoms.x[i] - f.atoms.x[j];
    double dy = f.atoms.y[i] - f.atoms.y[j];
    double dz = f.atoms.z[i] - f.atoms.z[j];
    // Orthorhombic minimum-image. LAMMPS-format data files wrap atoms into
    // `[0, L)`, so a bond that crosses the periodic boundary appears as a
    // raw span of ~L instead of its physical length. Without the wrap the
    // scanned bonded extension can exceed the box (e.g. 110 Å in an 80 Å
    // PE cell), which in turn inflates the halo cutoff past every slab
    // width and causes HaloExchangeForward to tag every owned atom as a
    // ghost on every face — blowing n_ghost past `dd->capacity`.
    const double lx = f.box(0, 0);
    const double ly = f.box(1, 1);
    const double lz = f.box(2, 2);
    if (lx > 0.0) dx -= lx * std::round(dx / lx);
    if (ly > 0.0) dy -= ly * std::round(dy / ly);
    if (lz > 0.0) dz -= lz * std::round(dz / lz);
    return dx * dx + dy * dy + dz * dz;
}

/// @brief Try to resolve a global atom tag to a local row. Returns -1 if
/// the atom is not present on this rank (either as owner or ghost). When
/// `dd == nullptr` (pre-partition / single-rank path) the tag is assumed
/// to equal the row index, matching the default `tag[i] = i` layout set
/// up by `Frame::Atoms`'s constructor.
inline int resolve_tag(const DomainDecomposition* dd, int atom_tag) {
    if (dd == nullptr) return atom_tag;
    auto it = dd->tag_to_local.find(static_cast<std::int64_t>(atom_tag));
    return (it == dd->tag_to_local.end()) ? -1 : it->second;
}

/// @brief Max intra-tuple span across every tuple whose atoms all resolve
/// (locally) on this rank. Each rank's local max is post-combined with an
/// `MPI_Allreduce(MAX)` in `compute_bonded_extension_mpi` so the global
/// bonded span is visible to every caller even when every cross-rank
/// tuple is partial on every rank.
double compute_bonded_extension_local(const Frame& frame,
                                       const DomainDecomposition* dd) {
    double max_sq = 0.0;

    const int nb = frame.bonds.size();
    for (int k = 0; k < nb; ++k) {
        const int i = resolve_tag(dd, frame.bonds.atomi[k]);
        const int j = resolve_tag(dd, frame.bonds.atomj[k]);
        if (i < 0 || j < 0) continue;
        max_sq = std::max(max_sq, pair_dist_sq(frame, i, j));
    }

    const int na = frame.angles.size();
    for (int k = 0; k < na; ++k) {
        const int i = resolve_tag(dd, frame.angles.atomi[k]);
        const int j = resolve_tag(dd, frame.angles.atomj[k]);
        const int m = resolve_tag(dd, frame.angles.atomk[k]);
        if (i < 0 || j < 0 || m < 0) continue;
        max_sq = std::max(max_sq, pair_dist_sq(frame, i, j));
        max_sq = std::max(max_sq, pair_dist_sq(frame, j, m));
        max_sq = std::max(max_sq, pair_dist_sq(frame, i, m));
    }

    const int nd = frame.dihedrals.size();
    for (int k = 0; k < nd; ++k) {
        const int i = resolve_tag(dd, frame.dihedrals.atomi[k]);
        const int j = resolve_tag(dd, frame.dihedrals.atomj[k]);
        const int m = resolve_tag(dd, frame.dihedrals.atomk[k]);
        const int l = resolve_tag(dd, frame.dihedrals.atoml[k]);
        if (i < 0 || j < 0 || m < 0 || l < 0) continue;
        max_sq = std::max(max_sq, pair_dist_sq(frame, i, j));
        max_sq = std::max(max_sq, pair_dist_sq(frame, i, m));
        max_sq = std::max(max_sq, pair_dist_sq(frame, i, l));
        max_sq = std::max(max_sq, pair_dist_sq(frame, j, m));
        max_sq = std::max(max_sq, pair_dist_sq(frame, j, l));
        max_sq = std::max(max_sq, pair_dist_sq(frame, m, l));
    }

    return (max_sq > 0.0) ? std::sqrt(max_sq) : 0.0;
}

} // namespace

double HaloExchangeForward::compute_bonded_extension(const cpu::Frame& frame) {
    // Public, no-DD-provided variant: assumes `frame.atoms.tag(i) == i`
    // (the default layout produced by `Frame::Atoms`) so that tag values
    // in `frame.bonds/angles/dihedrals` can be used as direct row indices
    // into `frame.atoms.x/y/z`. This holds for single-rank callers, for
    // pre-partition callers, and for all current unit tests. Rank-local
    // multi-rank scans must go through `compute_bonded_extension_local`
    // with the rank's `DomainDecomposition` to avoid reading wrong rows
    // in a post-`partition_atoms` compacted frame.
    return compute_bonded_extension_local(frame, /*dd=*/nullptr);
}

#ifdef ATV_HAS_MPI

namespace {

using detail::FaceBufs;

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

/// @brief Pack packed positions + tag for face `f` into a contiguous buffer.
///
/// Layout: `[x0, y0, z0, tag0, x1, y1, z1, tag1, ...]` -- 4 doubles per
/// atom. Tag is cast to double; atom counts up to 2^53 are exact. Applies
/// the per-face PBC shift to the packed coordinate for atoms whose owner
/// sits on the periodic edge along the send direction.
constexpr int kPosTagDoublesPerAtom = 4;

void pack_positions_for_face(const Frame& frame, const DomainDecomposition& dd,
                              int f, std::vector<double>& buf) {
    const auto& list = dd.send_lists[f];
    const int n = static_cast<int>(list.size());
    buf.resize(static_cast<std::size_t>(n) * kPosTagDoublesPerAtom);

    const int dim = face_dim(f);
    const double shift = pbc_shift_for_face(dd, f, frame.box(dim, dim));

    const double* x = frame.atoms.x.data();
    const double* y = frame.atoms.y.data();
    const double* z = frame.atoms.z.data();
    const std::int64_t* tag = frame.atoms.tag.data();

    for (int k = 0; k < n; ++k) {
        const int i = list[static_cast<std::size_t>(k)];
        double xi = x[i], yi = y[i], zi = z[i];
        if (dim == 0) xi += shift;
        else if (dim == 1) yi += shift;
        else if (dim == 2) zi += shift;
        const std::size_t off = static_cast<std::size_t>(k) * kPosTagDoublesPerAtom;
        buf[off + 0] = xi;
        buf[off + 1] = yi;
        buf[off + 2] = zi;
        buf[off + 3] = static_cast<double>(tag[i]);
    }
}

void prealloc_bufs(FaceBufs& b, int initial_cap) {
    for (int f = 0; f < 6; ++f) {
        // Conservative: each face could ship up to `initial_cap` atoms.
        // Forward sends 4 doubles per atom (xyz + tag). Reverse sends 3
        // doubles per atom (force). ~0.3 KB per ghost slot -- trivial.
        b.send_pos[f].reserve(static_cast<std::size_t>(initial_cap) * kPosTagDoublesPerAtom);
        b.recv_pos[f].reserve(static_cast<std::size_t>(initial_cap) * kPosTagDoublesPerAtom);
        b.send_force[f].reserve(static_cast<std::size_t>(initial_cap) * 3);
        b.recv_force[f].reserve(static_cast<std::size_t>(initial_cap) * 3);
    }
}

} // namespace

void HaloExchangeForward::build(const cpu::BuildContext& ctx) {
    // Compute the bonded skin extension from the rank-local topology.
    // Multi-rank: scan uses `dd->tag_to_local` to map topology atom-tags
    // to the post-partition row indices; tuples whose atoms are not all
    // resolvable locally contribute nothing on this rank, and the global
    // max is reconciled with an `MPI_Allreduce(MAX)` below. Single-rank:
    // skip the allreduce; the local scan alone is the global answer.
    auto* dd = ctx.runtime.dd;
    bonded_extension_ = compute_bonded_extension_local(ctx.frame, dd);

    if (cutoff_ <= 0.0) return;
    if (!dd || ctx.runtime.mpi_size <= 1) return;

#ifdef ATV_HAS_MPI
    if (dd->cart_comm != MPI_COMM_NULL) {
        double local = bonded_extension_;
        double global = 0.0;
        MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_MAX, dd->cart_comm);
        bonded_extension_ = global;
    }
#endif

    if (!bufs_) bufs_ = std::make_unique<detail::FaceBufs>();
    prealloc_bufs(*bufs_, dd->capacity);
    for (int f = 0; f < 6; ++f) dd->send_lists[f].reserve(static_cast<std::size_t>(dd->capacity));
}

void HaloExchangeForward::launch(const cpu::LaunchContext& ctx) {
    if (ctx.runtime.mpi_size <= 1 || cutoff_ <= 0.0) return;
    auto* dd = ctx.runtime.dd;
    if (!dd || dd->cart_comm == MPI_COMM_NULL) return;

    Frame& frame = ctx.frame;
    FaceBufs& bufs = *bufs_;

    // Roll back the previous step's ghost tag→slot entries before we
    // overwrite the ghost rows. Owners are not touched here.
    for (auto t : dd->ghost_tags) {
        dd->tag_to_local.erase(t);
    }
    dd->ghost_tags.clear();

    // Effective halo = max(pair cutoff, bonded skin): the wider of the
    // two is sufficient for both pair partners and bonded-tuple members
    // near the boundary.
    const double halo_cut = (cutoff_ > bonded_extension_)
                                ? cutoff_ : bonded_extension_;
    rebuild_send_lists(frame, *dd, halo_cut);

    // When `dims[dim] == 1` the face neighbors along that dim are *this
    // rank itself* (periodic self-wrap). OpenMPI 5.0.x on macOS occasionally
    // traps on self-Isend/Irecv patterns; mathematically the self-exchange
    // is a copy-from-myself anyway, so we short-circuit it: the atoms
    // already live in our slab. Short-circuit means `recv_counts[f] == 0`
    // and no MPI requests posted for that face.
    int my_rank = 0;
    MPI_Comm_rank(dd->cart_comm, &my_rank);

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
        if (neighbor == my_rank) continue;  // self-face: owner has atoms already
        MPI_Isend(&send_counts[f], 1, MPI_INT, neighbor, f,
                  dd->cart_comm, &count_reqs[n_count_reqs++]);
        MPI_Irecv(&recv_counts[f], 1, MPI_INT, neighbor, opposite_face(f),
                  dd->cart_comm, &count_reqs[n_count_reqs++]);
    }
    if (n_count_reqs > 0)
        MPI_Waitall(n_count_reqs, count_reqs, MPI_STATUSES_IGNORE);

    // Phase 2: pack positions+tag on each face, then exchange them. Each
    // packed atom is 4 doubles (x, y, z, tag-as-double). Send/recv counts
    // are inflated accordingly when posting MPI Isend/Irecv.
    for (int f = 0; f < 6; ++f) {
        if (dd->face_neighbors[static_cast<std::size_t>(f)] == MPI_PROC_NULL) {
            bufs.send_pos[f].clear();
        } else {
            pack_positions_for_face(frame, *dd, f, bufs.send_pos[f]);
        }
        bufs.recv_pos[f].assign(static_cast<std::size_t>(recv_counts[f])
                                * kPosTagDoublesPerAtom, 0.0);
    }

    MPI_Request pos_reqs[12];
    int n_pos_reqs = 0;
    for (int f = 0; f < 6; ++f) {
        const int neighbor = dd->face_neighbors[static_cast<std::size_t>(f)];
        if (neighbor == MPI_PROC_NULL) continue;
        if (neighbor == my_rank) continue;  // self-face short-circuit
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

    // Phase 3: unpack received positions+tag into ghost slots, in face
    // order. Each ghost's tag is also pushed into `dd->tag_to_local` so
    // bonded force kernels on this rank can resolve cross-rank tuple
    // partners as ghost rows.
    int slot = dd->n_local;
    for (int f = 0; f < 6; ++f) {
        const int n = recv_counts[f];
        for (int k = 0; k < n; ++k) {
            const std::size_t off = static_cast<std::size_t>(k) * kPosTagDoublesPerAtom;
            frame.atoms.x[slot] = bufs.recv_pos[f][off + 0];
            frame.atoms.y[slot] = bufs.recv_pos[f][off + 1];
            frame.atoms.z[slot] = bufs.recv_pos[f][off + 2];
            const std::int64_t tag =
                static_cast<std::int64_t>(bufs.recv_pos[f][off + 3]);
            frame.atoms.tag[slot] = tag;
            dd->tag_to_local[tag] = slot;
            dd->ghost_tags.push_back(tag);
            ++slot;
        }
        dd->recv_counts[static_cast<std::size_t>(f)] = n;
    }
    dd->n_ghost = slot - dd->n_local;
}

void HaloExchangeReverse::build(const cpu::BuildContext& ctx) {
    auto* dd = ctx.runtime.dd;
    if (!dd || ctx.runtime.mpi_size <= 1) return;
    if (!bufs_) bufs_ = std::make_unique<detail::FaceBufs>();
    prealloc_bufs(*bufs_, dd->capacity);
}

void HaloExchangeReverse::launch(const cpu::LaunchContext& ctx) {
    if (ctx.runtime.mpi_size <= 1) return;
    auto* dd = ctx.runtime.dd;
    if (!dd || dd->cart_comm == MPI_COMM_NULL) return;
    if (dd->n_ghost == 0) return;

    Frame& frame = ctx.frame;
    FaceBufs& bufs = *bufs_;

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

    int my_rank = 0;
    MPI_Comm_rank(dd->cart_comm, &my_rank);

    MPI_Request reqs[12];
    int n_reqs = 0;
    for (int f = 0; f < 6; ++f) {
        const int neighbor = dd->face_neighbors[static_cast<std::size_t>(f)];
        if (neighbor == MPI_PROC_NULL) continue;
        if (neighbor == my_rank) continue;  // self-face short-circuit
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
