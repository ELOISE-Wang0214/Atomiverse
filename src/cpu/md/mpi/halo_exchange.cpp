/// @file src/cpu/md/mpi/halo_exchange.cpp
/// @brief `HaloExchangeForward` / `HaloExchangeReverse` real implementation.
///
/// Orthorhombic, LAMMPS-style multi-stage forward comm: 3 sequential axis
/// stages (`x`, then `y`, then `z`), within each stage the ±-direction send
/// pair is dispatched concurrently. The key invariant is that each axis
/// stage's send list scans `[0, n_local + n_ghost_accumulated_so_far)` —
/// so ghosts received from an earlier axis are eligible to be **forwarded**
/// to a later-axis neighbor. After three stages every rank has the
/// owned-atom plus 1-, 2- and 3-axis-shifted ghost copies it needs, which
/// reaches all 26 cartesian neighbors via cascade (6 face + 12 edge + 8
/// corner). The previous single-round 6-face exchange dropped every
/// edge/corner ghost which silently dropped bonded tuples whose partner
/// lived diagonally — see `.claude/notes/cpu_mpi_halo_edge_corner_bug_2026-05-27.md`.
///
/// Ghost slots are appended in face order `{-x, +x, -y, +y, -z, +z}` so
/// `Reverse` can mirror the layout via `dd->recv_counts[]`; `Reverse`
/// processes the axes in reverse order so the cascade plays backwards.

#include "atomiverse/cpu/md/mpi/halo_exchange.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <memory>
#include <numeric>
#include <vector>

#include "atomiverse/core/log.h"
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
    std::array<std::vector<double>, 6> send_pos; // packed positions per face
    std::array<std::vector<double>, 6> recv_pos; // received positions per face
    std::array<std::vector<double>, 6> send_force; // packed forces (Reverse)
    std::array<std::vector<double>, 6> recv_force; // received forces (Reverse)
};

} // namespace detail

// ── Special members (defined here so the PIMPL dtor sees complete FaceBufs) ──

HaloExchangeForward::HaloExchangeForward() = default;
HaloExchangeForward::HaloExchangeForward(double cutoff, int every)
    : cutoff_(cutoff), every_(every < 1 ? 1 : every) {}
HaloExchangeForward::~HaloExchangeForward() = default;
HaloExchangeForward::HaloExchangeForward(HaloExchangeForward&&) noexcept =
    default;
HaloExchangeForward&
HaloExchangeForward::operator=(HaloExchangeForward&&) noexcept = default;

HaloExchangeReverse::HaloExchangeReverse() = default;
HaloExchangeReverse::~HaloExchangeReverse() = default;
HaloExchangeReverse::HaloExchangeReverse(HaloExchangeReverse&&) noexcept =
    default;
HaloExchangeReverse&
HaloExchangeReverse::operator=(HaloExchangeReverse&&) noexcept = default;

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
    const double lx = f.box.length(0);
    const double ly = f.box.length(1);
    const double lz = f.box.length(2);
    if (lx > 0.0)
        dx -= lx * std::round(dx / lx);
    if (ly > 0.0)
        dy -= ly * std::round(dy / ly);
    if (lz > 0.0)
        dz -= lz * std::round(dz / lz);
    return dx * dx + dy * dy + dz * dz;
}

/// @brief Try to resolve a global atom tag to a local row. Returns -1 if
/// the atom is not present on this rank (either as owner or ghost). When
/// `dd == nullptr` (pre-partition / single-rank path) the tag is assumed
/// to equal the row index, matching the default `tag[i] = i` layout set
/// up by `Frame::Atoms`'s constructor.
inline int resolve_tag(const DomainDecomposition* dd, int atom_tag) {
    if (dd == nullptr)
        return atom_tag;
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
        if (i < 0 || j < 0)
            continue;
        max_sq = std::max(max_sq, pair_dist_sq(frame, i, j));
    }

    const int na = frame.angles.size();
    for (int k = 0; k < na; ++k) {
        const int i = resolve_tag(dd, frame.angles.atomi[k]);
        const int j = resolve_tag(dd, frame.angles.atomj[k]);
        const int m = resolve_tag(dd, frame.angles.atomk[k]);
        if (i < 0 || j < 0 || m < 0)
            continue;
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
        if (i < 0 || j < 0 || m < 0 || l < 0)
            continue;
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
inline int face_dim(int f) { return f >> 1; }
inline int face_sign(int f) { return (f & 1) ? +1 : -1; }

/// @brief PBC shift along `dim` for atoms going out face `f` from this rank.
///
/// If this rank sits on the periodic edge in the +dim direction
/// (`coords[dim] == dims[dim]-1`) and we are sending out the +dim face,
/// the receiving rank is at the opposite edge of the global box;
/// pre-subtract `L_dim` from the packed coord so the receiver sees the
/// atom in its own coordinate frame. Symmetric for the -dim edge.
inline double pbc_shift_for_face(const DomainDecomposition& dd, int f,
                                 double L_dim) {
    const int dim = face_dim(f);
    const int sign = face_sign(f);
    if (sign > 0 && dd.coords[dim] == dd.dims[dim] - 1)
        return -L_dim;
    if (sign < 0 && dd.coords[dim] == 0)
        return +L_dim;
    return 0.0;
}

/// @brief Build the ±dim send pair, scanning frame indices `[0, sources_end)`.
///
/// Multi-stage forward comm scans owned atoms in the first axis stage, then
/// owned + axis-0 ghosts in the second stage, etc. The eligibility test is
/// the same proximity check used before (within `cutoff` of the slab
/// boundary along the send direction); the only change is the scan upper
/// bound, which is the caller's responsibility.
///
/// `dd.send_lists[2*dim]` (= face `−dim`) and `dd.send_lists[2*dim+1]` (=
/// face `+dim`) are cleared and refilled. Capacity is preserved.
void rebuild_send_lists_for_axis(const Frame& frame, DomainDecomposition& dd,
                                 int dim, int sources_end, double cutoff) {
    const int f_lo = 2 * dim;
    const int f_hi = 2 * dim + 1;
    dd.send_lists[static_cast<std::size_t>(f_lo)].clear();
    dd.send_lists[static_cast<std::size_t>(f_hi)].clear();

    const double* q = (dim == 0) ? frame.atoms.x.data()
                      : (dim == 1) ? frame.atoms.y.data()
                                   : frame.atoms.z.data();
    const double lo = dd.sublo[dim];
    const double hi = dd.subhi[dim];

    for (int i = 0; i < sources_end; ++i) {
        if (q[i] - lo < cutoff)
            dd.send_lists[static_cast<std::size_t>(f_lo)].push_back(i);
        if (hi - q[i] < cutoff)
            dd.send_lists[static_cast<std::size_t>(f_hi)].push_back(i);
    }
}

/// @brief Pack packed positions + tag for face `f` into a contiguous buffer.
///
/// Layout: `[x0, y0, z0, tag0, x1, y1, z1, tag1, ...]` -- 4 doubles per
/// atom. Tag is cast to double; atom counts up to 2^53 are exact. Applies
/// the per-face PBC shift to the packed coordinate for atoms whose owner
/// sits on the periodic edge along the send direction.
constexpr int kPosTagDoublesPerAtom = 4;
constexpr int kPosOnlyDoublesPerAtom = 3;

void pack_positions_for_face(const Frame& frame, const DomainDecomposition& dd,
                             int f, std::vector<double>& buf) {
    const auto& list = dd.send_lists[f];
    const int n = static_cast<int>(list.size());
    buf.resize(static_cast<std::size_t>(n) * kPosTagDoublesPerAtom);

    const int dim = face_dim(f);
    const double shift = pbc_shift_for_face(dd, f, frame.box.length(dim));

    const double* x = frame.atoms.x.data();
    const double* y = frame.atoms.y.data();
    const double* z = frame.atoms.z.data();
    const std::int64_t* tag = frame.atoms.tag.data();

    for (int k = 0; k < n; ++k) {
        const int i = list[static_cast<std::size_t>(k)];
        double xi = x[i], yi = y[i], zi = z[i];
        if (dim == 0)
            xi += shift;
        else if (dim == 1)
            yi += shift;
        else if (dim == 2)
            zi += shift;
        const std::size_t off =
            static_cast<std::size_t>(k) * kPosTagDoublesPerAtom;
        buf[off + 0] = xi;
        buf[off + 1] = yi;
        buf[off + 2] = zi;
        buf[off + 3] = static_cast<double>(tag[i]);
    }
}

/// @brief Pack only `(x,y,z)` for face `f` (position-only refresh between
/// reneighbor events). Ghost slot indices on the receiver stay fixed, so
/// we do not re-send tags.
void pack_positions_only_for_face(const Frame& frame,
                                  const DomainDecomposition& dd, int f,
                                  std::vector<double>& buf) {
    const auto& list = dd.send_lists[f];
    const int n = static_cast<int>(list.size());
    buf.resize(static_cast<std::size_t>(n) * kPosOnlyDoublesPerAtom);

    const int dim = face_dim(f);
    const double shift = pbc_shift_for_face(dd, f, frame.box.length(dim));

    const double* x = frame.atoms.x.data();
    const double* y = frame.atoms.y.data();
    const double* z = frame.atoms.z.data();

    for (int k = 0; k < n; ++k) {
        const int i = list[static_cast<std::size_t>(k)];
        double xi = x[i], yi = y[i], zi = z[i];
        if (dim == 0)
            xi += shift;
        else if (dim == 1)
            yi += shift;
        else if (dim == 2)
            zi += shift;
        const std::size_t off =
            static_cast<std::size_t>(k) * kPosOnlyDoublesPerAtom;
        buf[off + 0] = xi;
        buf[off + 1] = yi;
        buf[off + 2] = zi;
    }
}

void prealloc_bufs(FaceBufs& b, int initial_cap) {
    for (int f = 0; f < 6; ++f) {
        // Conservative: each face could ship up to `initial_cap` atoms.
        // Forward sends 4 doubles per atom (xyz + tag). Reverse sends 3
        // doubles per atom (force). ~0.3 KB per ghost slot -- trivial.
        b.send_pos[f].reserve(static_cast<std::size_t>(initial_cap)
                              * kPosTagDoublesPerAtom);
        b.recv_pos[f].reserve(static_cast<std::size_t>(initial_cap)
                              * kPosTagDoublesPerAtom);
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

    if (cutoff_ <= 0.0)
        return;
    if (!dd || ctx.runtime.mpi_size <= 1)
        return;

#ifdef ATV_HAS_MPI
    if (dd->cart_comm != MPI_COMM_NULL) {
        double local = bonded_extension_;
        double global = 0.0;
        MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_MAX, dd->cart_comm);
        bonded_extension_ = global;
    }
#endif

    if (!bufs_)
        bufs_ = std::make_unique<detail::FaceBufs>();
    prealloc_bufs(*bufs_, dd->capacity);
    for (int f = 0; f < 6; ++f)
        dd->send_lists[f].reserve(static_cast<std::size_t>(dd->capacity));
}

// Per-axis sub-stage helper for the forward path. Both `full_rebuild` and
// the position-only refresh path walk the same three axis stages; what
// differs is whether send lists are rebuilt and whether the tag column is
// included in the packet. Encapsulating one axis as a single function
// keeps the orchestrator readable.
namespace {

/// @brief Drive one axis-`dim` sub-stage of the cascading forward halo.
///
/// `slot_in_out` is the running ghost-slot write cursor: on entry it is
/// the next row to write (`n_local` before axis 0; after every axis it
/// advances by the total ghosts received in that axis). `frame_cap` is
/// the per-atom array capacity used for the overflow guard. When
/// `full_rebuild`, the function rescans `[0, sources_end)` for ghost
/// candidates, packs positions+tag, registers each received tag in
/// `tag_to_local` / `ghost_tags`, and stores per-face `recv_counts`. When
/// `!full_rebuild`, send lists and `recv_counts` are reused from the last
/// rebuild and only positions get refreshed.
///
/// Returns `true` on success, `false` on capacity overflow (the orchestrator
/// is expected to have already triggered `MPI_Abort` in that case via its
/// own logging, but the caller may want to short-circuit cleanup paths).
inline void forward_axis_stage(Frame& frame, DomainDecomposition& dd,
                               FaceBufs& bufs, int dim, int sources_end,
                               double cutoff, bool full_rebuild,
                               int frame_cap, int my_rank,
                               double bonded_extension_for_log,
                               int& slot_in_out) {
    const int f_lo = 2 * dim;
    const int f_hi = 2 * dim + 1;
    const std::size_t fz_lo = static_cast<std::size_t>(f_lo);
    const std::size_t fz_hi = static_cast<std::size_t>(f_hi);

    if (full_rebuild) {
        rebuild_send_lists_for_axis(frame, dd, dim, sources_end, cutoff);
    }

    const int neighbor_lo = dd.face_neighbors[fz_lo];
    const int neighbor_hi = dd.face_neighbors[fz_hi];

    int send_counts[2] = {static_cast<int>(dd.send_lists[fz_lo].size()),
                          static_cast<int>(dd.send_lists[fz_hi].size())};
    int recv_counts[2];

    if (full_rebuild) {
        recv_counts[0] = 0;
        recv_counts[1] = 0;
        MPI_Request count_reqs[4];
        int n_creqs = 0;
        auto post_count = [&](int s) {
            const int face_s = (s == 0) ? f_lo : f_hi;
            const int neighbor = (s == 0) ? neighbor_lo : neighbor_hi;
            if (neighbor == MPI_PROC_NULL || neighbor == my_rank)
                return;
            MPI_Isend(&send_counts[s], 1, MPI_INT, neighbor, face_s,
                      dd.cart_comm, &count_reqs[n_creqs++]);
            MPI_Irecv(&recv_counts[s], 1, MPI_INT, neighbor,
                      opposite_face(face_s), dd.cart_comm,
                      &count_reqs[n_creqs++]);
        };
        post_count(0);
        post_count(1);
        if (n_creqs > 0)
            MPI_Waitall(n_creqs, count_reqs, MPI_STATUSES_IGNORE);
    } else {
        recv_counts[0] = dd.recv_counts[fz_lo];
        recv_counts[1] = dd.recv_counts[fz_hi];
    }

    const int per_atom_doubles =
        full_rebuild ? kPosTagDoublesPerAtom : kPosOnlyDoublesPerAtom;

    for (int s = 0; s < 2; ++s) {
        const int face_s = (s == 0) ? f_lo : f_hi;
        const int neighbor = (s == 0) ? neighbor_lo : neighbor_hi;
        auto& sbuf = bufs.send_pos[static_cast<std::size_t>(face_s)];
        auto& rbuf = bufs.recv_pos[static_cast<std::size_t>(face_s)];
        if (neighbor == MPI_PROC_NULL) {
            sbuf.clear();
        } else if (full_rebuild) {
            pack_positions_for_face(frame, dd, face_s, sbuf);
        } else {
            pack_positions_only_for_face(frame, dd, face_s, sbuf);
        }
        rbuf.assign(static_cast<std::size_t>(recv_counts[s])
                        * static_cast<std::size_t>(per_atom_doubles),
                    0.0);
    }

    MPI_Request pos_reqs[4];
    int n_preqs = 0;
    for (int s = 0; s < 2; ++s) {
        const int face_s = (s == 0) ? f_lo : f_hi;
        const int neighbor = (s == 0) ? neighbor_lo : neighbor_hi;
        if (neighbor == MPI_PROC_NULL || neighbor == my_rank)
            continue;
        if (send_counts[s] > 0) {
            MPI_Isend(bufs.send_pos[static_cast<std::size_t>(face_s)].data(),
                      static_cast<int>(
                          bufs.send_pos[static_cast<std::size_t>(face_s)]
                              .size()),
                      MPI_DOUBLE, neighbor, 100 + face_s, dd.cart_comm,
                      &pos_reqs[n_preqs++]);
        }
        if (recv_counts[s] > 0) {
            MPI_Irecv(bufs.recv_pos[static_cast<std::size_t>(face_s)].data(),
                      static_cast<int>(
                          bufs.recv_pos[static_cast<std::size_t>(face_s)]
                              .size()),
                      MPI_DOUBLE, neighbor, 100 + opposite_face(face_s),
                      dd.cart_comm, &pos_reqs[n_preqs++]);
        }
    }
    if (n_preqs > 0)
        MPI_Waitall(n_preqs, pos_reqs, MPI_STATUSES_IGNORE);

    if (full_rebuild) {
        // Capacity guard for this axis. The per-atom xtensor arrays are
        // sized to `frame.atoms.size()` (== original global atom count);
        // overflowing here silently corrupts memory and eventually trips
        // LJ / bonded kernels.
        const int new_ghosts = recv_counts[0] + recv_counts[1];
        if (slot_in_out + new_ghosts > frame_cap) {
            ATV_CRITICAL(
                "halo_exchange: rank {} ghost capacity overflow at axis {}: "
                "slot={} + new={} > frame.atoms.size()={}. "
                "send=[{},{}] recv=[{},{}] cutoff={} bonded_ext={}",
                my_rank, dim, slot_in_out, new_ghosts, frame_cap,
                send_counts[0], send_counts[1], recv_counts[0], recv_counts[1],
                cutoff, bonded_extension_for_log);
            MPI_Abort(dd.cart_comm, 1);
        }
    }

    // Unpack into ghost slots — for full rebuilds register tags too. The
    // unpack order ALSO defines the slot layout that `Reverse` mirrors
    // via `dd->recv_counts`.
    int slot = slot_in_out;
    for (int s = 0; s < 2; ++s) {
        const int face_s = (s == 0) ? f_lo : f_hi;
        const int n = recv_counts[s];
        const auto& rbuf = bufs.recv_pos[static_cast<std::size_t>(face_s)];
        if (full_rebuild) {
            for (int k = 0; k < n; ++k) {
                const std::size_t off =
                    static_cast<std::size_t>(k) * kPosTagDoublesPerAtom;
                frame.atoms.x[slot] = rbuf[off + 0];
                frame.atoms.y[slot] = rbuf[off + 1];
                frame.atoms.z[slot] = rbuf[off + 2];
                const std::int64_t tag =
                    static_cast<std::int64_t>(rbuf[off + 3]);
                frame.atoms.tag[slot] = tag;
                dd.tag_to_local[tag] = slot;
                dd.ghost_tags.push_back(tag);
                ++slot;
            }
            dd.recv_counts[static_cast<std::size_t>(face_s)] = n;
        } else {
            for (int k = 0; k < n; ++k) {
                const std::size_t off =
                    static_cast<std::size_t>(k) * kPosOnlyDoublesPerAtom;
                frame.atoms.x[slot] = rbuf[off + 0];
                frame.atoms.y[slot] = rbuf[off + 1];
                frame.atoms.z[slot] = rbuf[off + 2];
                ++slot;
            }
        }
    }
    slot_in_out = slot;
}

} // namespace

void HaloExchangeForward::launch(const cpu::LaunchContext& ctx) {
    if (ctx.runtime.mpi_size <= 1 || cutoff_ <= 0.0)
        return;
    auto* dd = ctx.runtime.dd;
    if (!dd || dd->cart_comm == MPI_COMM_NULL)
        return;

    Frame& frame = ctx.frame;
    FaceBufs& bufs = *bufs_;

    // Decide whether this step needs a full rebuild of per-face send lists
    // and ghost-slot assignments. LAMMPS-style cadence: rebuild every
    // `every_` steps, plus any step where migration has bumped the
    // neighbor generation (owner indices changed → any saved send_list is
    // stale). In between, only ghost positions are refreshed; slot indices
    // and the tag→slot map carry through unchanged, which lets downstream
    // `NeighborList` reuse its CSR from the last reneighbor step.
    ++ago_;
    const bool generation_changed =
        dd->neighbor_generation != last_rebuild_generation_;
    const bool full_rebuild = generation_changed || (ago_ >= every_);

    int my_rank = 0;
    MPI_Comm_rank(dd->cart_comm, &my_rank);

    const double halo_cut =
        (cutoff_ > bonded_extension_) ? cutoff_ : bonded_extension_;
    const int frame_cap = static_cast<int>(frame.atoms.tag.size());

    if (full_rebuild) {
        // Effective halo = max(pair cutoff, bonded skin): the wider of
        // the two is sufficient for both pair partners and bonded-tuple
        // members near the boundary.
        for (auto t : dd->ghost_tags)
            dd->tag_to_local.erase(t);
        dd->ghost_tags.clear();

        // 3 axis stages; each stage's send list scans [0, slot) so atoms
        // received in earlier axes get *forwarded* to later-axis
        // neighbors. After axis 2 (z), each rank has all 26 cartesian
        // neighbor ghosts via cascade.
        int slot = dd->n_local;
        for (int dim = 0; dim < 3; ++dim) {
            forward_axis_stage(frame, *dd, bufs, dim, slot, halo_cut,
                               /*full_rebuild=*/true, frame_cap, my_rank,
                               bonded_extension_, slot);
        }
        dd->n_ghost = slot - dd->n_local;

        // Signal to ghost-index-caching consumers (the NeighborList CSR
        // is the load-bearing one) that ghost slot assignments may have
        // changed. Without this bump the NL pair list keeps `j` indices
        // pointing at slots whose atoms have been replaced by the
        // current rebuild, producing a per-rebuild PE / virial
        // discontinuity that drifts the NPT extended Hamiltonian
        // (~1-3 eV per `every_=10` rebuild → catastrophic at step
        // ~1000). See `.claude/notes/msu_stage2_conserved_drift_2026-05-27.md`
        // for why we use a halo-local counter instead of bumping the
        // shared `neighbor_generation`.
        ++dd->halo_generation;

        ago_ = 0;
        last_rebuild_generation_ = dd->neighbor_generation;
        return;
    }

    // Position-only refresh path (no reneighbor). Same three axis stages
    // play out; per-axis ghost positions must be refreshed in order so
    // a later axis sees the up-to-date forwarded ghosts (a later-axis
    // send list may reference an earlier-axis ghost slot).
    int slot = dd->n_local;
    for (int dim = 0; dim < 3; ++dim) {
        forward_axis_stage(frame, *dd, bufs, dim, /*sources_end=*/0, halo_cut,
                           /*full_rebuild=*/false, frame_cap, my_rank,
                           bonded_extension_, slot);
    }
}

void HaloExchangeReverse::build(const cpu::BuildContext& ctx) {
    auto* dd = ctx.runtime.dd;
    if (!dd || ctx.runtime.mpi_size <= 1)
        return;
    if (!bufs_)
        bufs_ = std::make_unique<detail::FaceBufs>();
    prealloc_bufs(*bufs_, dd->capacity);
}

void HaloExchangeReverse::launch(const cpu::LaunchContext& ctx) {
    if (ctx.runtime.mpi_size <= 1)
        return;
    auto* dd = ctx.runtime.dd;
    if (!dd || dd->cart_comm == MPI_COMM_NULL)
        return;
    if (dd->n_ghost == 0)
        return;

    Frame& frame = ctx.frame;
    FaceBufs& bufs = *bufs_;

    int my_rank = 0;
    MPI_Comm_rank(dd->cart_comm, &my_rank);

    // Per-face ghost slot start positions, fixed by the forward halo's
    // unpack order: slot_start[f] = n_local + Σ_{j<f} recv_counts[j].
    int slot_start[6];
    {
        int s = dd->n_local;
        for (int f = 0; f < 6; ++f) {
            slot_start[f] = s;
            s += dd->recv_counts[static_cast<std::size_t>(f)];
        }
    }

    // Reverse the axis cascade: z first, then y, then x. Forces accumulated
    // on a swap-K ghost row by a later-axis reverse stage are then sent
    // back further down the cascade in earlier-axis stages. Because the
    // forward axis-K send list referenced rows in [0, n_local + Σ_{j<K}
    // recv_counts[j]), the reverse accumulates correctly into rows that
    // are either owned (axis 0) or earlier-axis ghosts (axes 1/2). By the
    // time axis 0 (x) finishes, all force has flowed back to owners.
    for (int dim = 2; dim >= 0; --dim) {
        const int f_lo = 2 * dim;
        const int f_hi = 2 * dim + 1;
        const std::size_t fz_lo = static_cast<std::size_t>(f_lo);
        const std::size_t fz_hi = static_cast<std::size_t>(f_hi);

        const int neighbor_lo = dd->face_neighbors[fz_lo];
        const int neighbor_hi = dd->face_neighbors[fz_hi];

        // For face K's reverse: we send the K-ghost slot forces back to the
        // rank that owns them (= face_neighbors[K], the same rank we
        // exchanged with in forward), and we receive forces for the atoms
        // we sent forward (those live at indices in dd->send_lists[K] on
        // this rank).
        int send_force_counts[2] = {dd->recv_counts[fz_lo],
                                    dd->recv_counts[fz_hi]};
        int recv_force_counts[2] = {
            static_cast<int>(dd->send_lists[fz_lo].size()),
            static_cast<int>(dd->send_lists[fz_hi].size())};

        for (int s = 0; s < 2; ++s) {
            const int face_s = (s == 0) ? f_lo : f_hi;
            const int n = send_force_counts[s];
            auto& sbuf = bufs.send_force[static_cast<std::size_t>(face_s)];
            sbuf.resize(static_cast<std::size_t>(n) * 3);
            int slot = slot_start[face_s];
            for (int k = 0; k < n; ++k) {
                sbuf[static_cast<std::size_t>(3 * k + 0)] = frame.atoms.fx[slot];
                sbuf[static_cast<std::size_t>(3 * k + 1)] = frame.atoms.fy[slot];
                sbuf[static_cast<std::size_t>(3 * k + 2)] = frame.atoms.fz[slot];
                ++slot;
            }
            bufs.recv_force[static_cast<std::size_t>(face_s)].assign(
                static_cast<std::size_t>(recv_force_counts[s]) * 3, 0.0);
        }

        MPI_Request reqs[4];
        int n_reqs = 0;
        for (int s = 0; s < 2; ++s) {
            const int face_s = (s == 0) ? f_lo : f_hi;
            const int neighbor = (s == 0) ? neighbor_lo : neighbor_hi;
            if (neighbor == MPI_PROC_NULL || neighbor == my_rank)
                continue;
            if (send_force_counts[s] > 0) {
                MPI_Isend(
                    bufs.send_force[static_cast<std::size_t>(face_s)].data(),
                    static_cast<int>(
                        bufs.send_force[static_cast<std::size_t>(face_s)]
                            .size()),
                    MPI_DOUBLE, neighbor, 200 + face_s, dd->cart_comm,
                    &reqs[n_reqs++]);
            }
            if (recv_force_counts[s] > 0) {
                MPI_Irecv(
                    bufs.recv_force[static_cast<std::size_t>(face_s)].data(),
                    static_cast<int>(
                        bufs.recv_force[static_cast<std::size_t>(face_s)]
                            .size()),
                    MPI_DOUBLE, neighbor, 200 + opposite_face(face_s),
                    dd->cart_comm, &reqs[n_reqs++]);
            }
        }
        if (n_reqs > 0)
            MPI_Waitall(n_reqs, reqs, MPI_STATUSES_IGNORE);

        // Accumulate received forces into send_list indices. Critical: for
        // dim > 0 these indices may include ghost-row indices from earlier
        // axes — those forces will be re-packed and sent further along the
        // cascade in the next (lower-dim) iteration of this loop.
        for (int s = 0; s < 2; ++s) {
            const int face_s = (s == 0) ? f_lo : f_hi;
            const auto& list = dd->send_lists[static_cast<std::size_t>(face_s)];
            const auto& rbuf = bufs.recv_force[static_cast<std::size_t>(face_s)];
            const int n = static_cast<int>(list.size());
            for (int k = 0; k < n; ++k) {
                const int i = list[static_cast<std::size_t>(k)];
                frame.atoms.fx[i] += rbuf[static_cast<std::size_t>(3 * k + 0)];
                frame.atoms.fy[i] += rbuf[static_cast<std::size_t>(3 * k + 1)];
                frame.atoms.fz[i] += rbuf[static_cast<std::size_t>(3 * k + 2)];
            }
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
