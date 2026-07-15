/// @file src/cpu/comm/halo_exchange.cpp
/// @brief `HaloSend` / `HaloRecv` real implementation.
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

#include "atomiverse/cpu/comm/halo_exchange.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <memory>
#include <numeric>
#include <vector>

#include "atomiverse/core/log.h"
#include "atomiverse/cpu/comm/domain.h"
#include "atomiverse/cpu/frame/frame.h"

#ifdef ATV_HAS_MPI
#include <mpi.h>
#endif

namespace ATV {
namespace cpu {
namespace comm {

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

namespace {
/// @brief Flatten the declared fields of every bound participant into
/// `out`, skipping null/degenerate descriptors. One-time (first launch).
void append_declared_fields(
    const std::vector<const HaloParticipant*>& participants, bool forward,
    std::vector<HaloParticipant::CommField>& out) {
    for (const auto* p : participants) {
        const auto fs =
            forward ? p->comm_forward_fields() : p->comm_reverse_fields();
        for (const auto& f : fs)
            if (f.data && f.width > 0)
                out.push_back(f);
    }
}
} // namespace

// ── Special members (defined here so the PIMPL dtor sees complete FaceBufs) ──

HaloSend::HaloSend() = default;
HaloSend::HaloSend(double cutoff, int every)
    : cutoff_(cutoff), every_(every < 1 ? 1 : every) {}
HaloSend::~HaloSend() = default;
HaloSend::HaloSend(HaloSend&&) noexcept = default;
HaloSend& HaloSend::operator=(HaloSend&&) noexcept = default;

HaloRecv::HaloRecv() = default;
HaloRecv::~HaloRecv() = default;
HaloRecv::HaloRecv(HaloRecv&&) noexcept = default;
HaloRecv& HaloRecv::operator=(HaloRecv&&) noexcept = default;

std::size_t HaloSend::send_buf_capacity_for_test() const {
    return bufs_ ? bufs_->send_pos[0].capacity() : 0;
}

namespace {
// Resident bytes of a FaceBufs (4 arrays × 6 per-face double vectors).
std::size_t facebufs_bytes(const detail::FaceBufs* b) {
    if (!b)
        return 0;
    // size() (live elements) ≈ resident; reserved slack is lazily mapped.
    std::size_t n = 0;
    for (int f = 0; f < 6; ++f)
        n += (b->send_pos[f].size() + b->recv_pos[f].size()
              + b->send_force[f].size() + b->recv_force[f].size())
             * sizeof(double);
    return n;
}
} // namespace

std::size_t HaloSend::memory_bytes() const {
    return facebufs_bytes(bufs_.get());
}
std::size_t HaloRecv::memory_bytes() const {
    return facebufs_bytes(bufs_.get());
}

// ── Topology-scan helpers (available with or without MPI) ────────────

namespace {

/// @brief Try to resolve a global atom tag to a local row. Returns -1 if
/// the atom is not present on this rank (either as owner or ghost). When
/// `dd == nullptr` (pre-partition / single-rank path) the tag is assumed
/// to equal the row index, matching the default `tag[i] = i` layout set
/// up by `Frame::Atoms`'s constructor.
inline int resolve_tag(const md::mpi::DomainDecomposition* dd, int atom_tag) {
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
                                      const md::mpi::DomainDecomposition* dd) {
    // Reads per-atom bonded topology from `frame.atoms` (was: flat
    // bonds/angles/dihedrals blocks on Frame; those are deleted).
    // Single-rank: `dd` may be null; we treat all atoms as locally
    // resolvable. Multi-rank: `dd->tag_to_local` (still maintained by
    // the old halo / migration path until full switchover) is the
    // resolver — same semantics as before, just sourced from a
    // different layout.
    double max_sq = 0.0;
    const auto& a = frame.atoms;

    // Hoisted ONCE at the launch boundary; `pair_d2` below only ever reads
    // through these (topology arrays num_bond/bond_atom/... stay direct
    // member access -- they are not migrated per the frame-block-02 scope).
    const StorageReal* __restrict sp_x = frame.span(kAtom_x);
    const StorageReal* __restrict sp_y = frame.span(kAtom_y);
    const StorageReal* __restrict sp_z = frame.span(kAtom_z);

    auto pair_d2 = [&](int i, int j) {
        double dx = sp_x[i] - sp_x[j];
        double dy = sp_y[i] - sp_y[j];
        double dz = sp_z[i] - sp_z[j];
        frame.box.min_image(dx, dy, dz);
        return dx * dx + dy * dy + dz * dz;
    };

    auto resolve_via_dd = [&](std::int64_t tag) -> int {
        if (dd == nullptr)
            return static_cast<int>(tag);
        auto it = dd->tag_to_local.find(tag);
        return (it == dd->tag_to_local.end()) ? -1 : it->second;
    };

    // Iterate the *compacted* owner range (`dd->n_local` in multi-rank;
    // `atom.nlocal` single-rank). Beyond that range, rows hold stale
    // pre-partition topology that would resolve to current rows whose
    // positions belong to entirely different atoms — giving spurious
    // 50+ Å bond distances and an astronomically wrong bonded extent.
    const TagInt* __restrict sp_tag = frame.span(kAtom_tag);
    const int nlocal_iter = (dd != nullptr) ? dd->n_local : a.nlocal;
    for (int i_owner = 0; i_owner < nlocal_iter; ++i_owner) {
        const std::int64_t ti = sp_tag[i_owner];
        // Bonds: owner is the smaller-tag atom; partner stored in
        // bond_atom[i][m].
        const int nb = a.num_bond[i_owner];
        for (int m = 0; m < nb; ++m) {
            const int i = resolve_via_dd(ti);
            const int j = resolve_via_dd(a.bond_atom(i_owner, m));
            if (i < 0 || j < 0)
                continue;
            max_sq = std::max(max_sq, pair_d2(i, j));
        }
        // Angles: center atom owns; non-center partners are
        // (angle_atom(i, 2m), angle_atom(i, 2m+1)).
        const int na = a.num_angle[i_owner];
        for (int m = 0; m < na; ++m) {
            const int j = resolve_via_dd(ti); // center
            const int i = resolve_via_dd(a.angle_atom(i_owner, 2 * m + 0));
            const int k = resolve_via_dd(a.angle_atom(i_owner, 2 * m + 1));
            if (i < 0 || j < 0 || k < 0)
                continue;
            max_sq = std::max(max_sq, pair_d2(i, j));
            max_sq = std::max(max_sq, pair_d2(j, k));
            max_sq = std::max(max_sq, pair_d2(i, k));
        }
        // Dihedrals: atom-j owns; non-j partners are (i, k, l) at offsets
        // 3m+0/1/2.
        const int nd = a.num_dihedral[i_owner];
        for (int m = 0; m < nd; ++m) {
            const int j = resolve_via_dd(ti);
            const int i = resolve_via_dd(a.dihedral_atom(i_owner, 3 * m + 0));
            const int k = resolve_via_dd(a.dihedral_atom(i_owner, 3 * m + 1));
            const int l = resolve_via_dd(a.dihedral_atom(i_owner, 3 * m + 2));
            if (i < 0 || j < 0 || k < 0 || l < 0)
                continue;
            max_sq = std::max(max_sq, pair_d2(i, j));
            max_sq = std::max(max_sq, pair_d2(i, k));
            max_sq = std::max(max_sq, pair_d2(i, l));
            max_sq = std::max(max_sq, pair_d2(j, k));
            max_sq = std::max(max_sq, pair_d2(j, l));
            max_sq = std::max(max_sq, pair_d2(k, l));
        }
    }
    return (max_sq > 0.0) ? std::sqrt(max_sq) : 0.0;
}

} // namespace

double HaloSend::compute_bonded_extension(const cpu::Frame& frame) {
    // Public, no-DD-provided variant: assumes `frame.atoms.tag(i) == i`
    // (the default layout produced by `Frame::Atoms`) so that tag values
    // in `frame.atoms per-atom topology` can be used as direct row indices
    // into `frame.atoms.x/y/z`. This holds for single-rank callers, for
    // pre-partition callers, and for all current unit tests. Rank-local
    // multi-rank scans must go through `compute_bonded_extension_local`
    // with the rank's `md::mpi::DomainDecomposition` to avoid reading wrong rows
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
inline double pbc_shift_for_face(const md::mpi::DomainDecomposition& dd, int f,
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
void rebuild_send_lists_for_axis(const Frame& frame,
                                 md::mpi::DomainDecomposition& dd, int dim,
                                 int sources_end, double cutoff) {
    const int f_lo = 2 * dim;
    const int f_hi = 2 * dim + 1;
    dd.send_lists[static_cast<std::size_t>(f_lo)].clear();
    dd.send_lists[static_cast<std::size_t>(f_hi)].clear();

    const StorageReal* __restrict q = (dim == 0)   ? frame.span(kAtom_x)
                                      : (dim == 1) ? frame.span(kAtom_y)
                                                   : frame.span(kAtom_z);
    const double lo = dd.sublo[dim];
    const double hi = dd.subhi[dim];

    for (int i = 0; i < sources_end; ++i) {
        if (q[i] - lo < cutoff)
            dd.send_lists[static_cast<std::size_t>(f_lo)].push_back(i);
        if (hi - q[i] < cutoff)
            dd.send_lists[static_cast<std::size_t>(f_hi)].push_back(i);
    }
}

/// @brief Pack the border message (positions + tag + charge) for face `f`.
///
/// Layout: `[x0, y0, z0, tag0, q0, <extra fields>, x1, ...]` --
/// `kBorderDoublesPerAtom + Σ declared forward widths` doubles per atom. Sent only at
/// full rebuilds (reneighbor), mirroring LAMMPS `borders()`/AtomVec carrying `q`
/// to ghosts: charge is static per run, so ghosts get it once and the per-step
/// forward refresh stays positions-only. Tag is cast to double (exact to 2^53).
/// Applies the per-face PBC shift to the packed coordinate for atoms whose owner
/// sits on the periodic edge along the send direction.
constexpr int kBorderDoublesPerAtom = 5; // x, y, z, tag, q
constexpr int kPosOnlyDoublesPerAtom = 3;
constexpr int kForceDoublesPerAtom = 3; // fx, fy, fz

using CommField = HaloParticipant::CommField;

/// @brief Total per-atom doubles contributed by declared extra fields.
int extra_width(const std::vector<CommField>& fields) {
    int w = 0;
    for (const auto& f : fields)
        w += f.width;
    return w;
}

/// @brief Append the declared extra fields for send-list atom `i` at
/// `buf[off]` (contiguous, participant order). Returns doubles written.
std::size_t pack_extra_fields(const std::vector<CommField>& fields, int i,
                              std::vector<double>& buf, std::size_t off) {
    for (const auto& f : fields) {
        const std::size_t base =
            static_cast<std::size_t>(i) * static_cast<std::size_t>(f.width);
        for (int k = 0; k < f.width; ++k)
            buf[off++] = static_cast<double>(f.data[base + k]);
    }
    return off;
}

/// @brief Write the extra-field payload at `rbuf[off]` into ghost row
/// `slot` (forward: overwrite). Returns doubles consumed.
std::size_t unpack_extra_fields(const std::vector<CommField>& fields, int slot,
                                const std::vector<double>& rbuf,
                                std::size_t off) {
    for (const auto& f : fields) {
        const std::size_t base =
            static_cast<std::size_t>(slot) * static_cast<std::size_t>(f.width);
        for (int k = 0; k < f.width; ++k)
            f.data[base + k] = static_cast<StorageReal>(rbuf[off++]);
    }
    return off;
}

/// @brief ACCUMULATE the extra-field payload at `rbuf[off]` onto row `i`
/// (reverse: exactly-once partial sums, mirroring the intrinsic forces).
std::size_t accumulate_extra_fields(const std::vector<CommField>& fields, int i,
                                    const std::vector<double>& rbuf,
                                    std::size_t off) {
    for (const auto& f : fields) {
        const std::size_t base =
            static_cast<std::size_t>(i) * static_cast<std::size_t>(f.width);
        for (int k = 0; k < f.width; ++k)
            f.data[base + k] += static_cast<StorageReal>(rbuf[off++]);
    }
    return off;
}

void pack_positions_for_face(const Frame& frame,
                             const md::mpi::DomainDecomposition& dd, int f,
                             std::vector<double>& buf,
                             const std::vector<CommField>& xfields) {
    const auto& list = dd.send_lists[f];
    const int n = static_cast<int>(list.size());
    const int per_atom = kBorderDoublesPerAtom + extra_width(xfields);
    buf.resize(static_cast<std::size_t>(n) * per_atom);

    const int dim = face_dim(f);
    const double shift = pbc_shift_for_face(dd, f, frame.box.length(dim));

    const StorageReal* __restrict x = frame.span(kAtom_x);
    const StorageReal* __restrict y = frame.span(kAtom_y);
    const StorageReal* __restrict z = frame.span(kAtom_z);
    const TagInt* __restrict tag = frame.span(kAtom_tag);
    const StorageReal* __restrict q = frame.span(kAtom_charge);

    for (int k = 0; k < n; ++k) {
        const int i = list[static_cast<std::size_t>(k)];
        double xi = x[i], yi = y[i], zi = z[i];
        if (dim == 0)
            xi += shift;
        else if (dim == 1)
            yi += shift;
        else if (dim == 2)
            zi += shift;
        std::size_t off = static_cast<std::size_t>(k) * per_atom;
        buf[off + 0] = xi;
        buf[off + 1] = yi;
        buf[off + 2] = zi;
        buf[off + 3] = static_cast<double>(tag[i]);
        buf[off + 4] = q[i]; // static per run; carried to ghosts at reneighbor
        pack_extra_fields(xfields, i, buf, off + kBorderDoublesPerAtom);
    }
}

/// @brief Pack `(x,y,z)` plus declared extra fields for face `f`
/// (position-only refresh between reneighbor events). Ghost slot indices
/// on the receiver stay fixed, so we do not re-send tags; extra forward
/// fields refresh every step alongside positions.
void pack_positions_only_for_face(const Frame& frame,
                                  const md::mpi::DomainDecomposition& dd, int f,
                                  std::vector<double>& buf,
                                  const std::vector<CommField>& xfields) {
    const auto& list = dd.send_lists[f];
    const int n = static_cast<int>(list.size());
    const int per_atom = kPosOnlyDoublesPerAtom + extra_width(xfields);
    buf.resize(static_cast<std::size_t>(n) * per_atom);

    const int dim = face_dim(f);
    const double shift = pbc_shift_for_face(dd, f, frame.box.length(dim));

    const StorageReal* __restrict x = frame.span(kAtom_x);
    const StorageReal* __restrict y = frame.span(kAtom_y);
    const StorageReal* __restrict z = frame.span(kAtom_z);

    for (int k = 0; k < n; ++k) {
        const int i = list[static_cast<std::size_t>(k)];
        double xi = x[i], yi = y[i], zi = z[i];
        if (dim == 0)
            xi += shift;
        else if (dim == 1)
            yi += shift;
        else if (dim == 2)
            zi += shift;
        std::size_t off = static_cast<std::size_t>(k) * per_atom;
        buf[off + 0] = xi;
        buf[off + 1] = yi;
        buf[off + 2] = zi;
        pack_extra_fields(xfields, i, buf, off + kPosOnlyDoublesPerAtom);
    }
}

void prealloc_bufs(FaceBufs& b, int initial_cap) {
    for (int f = 0; f < 6; ++f) {
        // Conservative: each face could ship up to `initial_cap` atoms.
        // Border (full-rebuild) sends xyz + tag + q (+ declared forward
        // fields); the per-step forward sends xyz (+ declared forward
        // fields). Reverse sends fx,fy,fz (+ declared reverse fields).
        b.send_pos[f].reserve(static_cast<std::size_t>(initial_cap)
                              * kBorderDoublesPerAtom);
        b.recv_pos[f].reserve(static_cast<std::size_t>(initial_cap)
                              * kBorderDoublesPerAtom);
        b.send_force[f].reserve(static_cast<std::size_t>(initial_cap) * 3);
        b.recv_force[f].reserve(static_cast<std::size_t>(initial_cap) * 3);
    }
}

} // namespace

void HaloSend::build(const cpu::BuildContext& ctx) {
    // Compute the bonded skin extension from the rank-local topology.
    // Multi-rank: scan uses `dd->tag_to_local` to map topology atom-tags
    // to the post-partition row indices; tuples whose atoms are not all
    // resolvable locally contribute nothing on this rank, and the global
    // max is reconciled with an `MPI_Allreduce(MAX)` below. Single-rank:
    // skip the allreduce; the local scan alone is the global answer.
    auto* dd = ctx.runtime.md.dd;
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

    // The ghost shell cannot be wider than the subdomain it is grown from: the halo
    // exchanges with IMMEDIATE neighbour ranks only, so a shell deeper than one
    // subdomain would need atoms from a rank two hops away, and those simply never
    // arrive. Until now that failed SILENTLY — the ghost region was quietly truncated
    // and the forces came out wrong. It surfaced the first time a component asked for
    // a ghost-of-ghost halo (EAM, halo = 2*cutoff): on a box that was big enough for
    // a normal pair halo but not for a doubled one, the decomposed forces differed
    // from a brute-force reference by 1.6 eV/A, with no diagnostic of any kind.
    //
    // Refuse instead. A user who hits this needs a bigger box, fewer ranks, or a
    // smaller cutoff — and must be told so, not handed wrong physics.
    const double need = halo_cutoff();
    for (int d = 0; d < 3; ++d) {
        if (dd->dims[d] <= 1) continue; // not decomposed along this axis
        const double sub = dd->subhi[d] - dd->sublo[d];
        if (need > sub) {
            throw std::runtime_error(
                "HaloSend: the requested ghost shell (" + std::to_string(need)
                + " A) is wider than this rank's subdomain along axis "
                + std::to_string(d) + " (" + std::to_string(sub)
                + " A). The halo only exchanges with immediate neighbour ranks, so a "
                  "deeper shell cannot be filled and the ghost region would be "
                  "silently truncated. Use a larger box, fewer ranks along that axis, "
                  "or a smaller cutoff. (A many-body style such as EAM asks for a "
                  "ghost-of-ghost shell of 2*cutoff.)");
        }
    }

    if (!bufs_)
        bufs_ = std::make_unique<detail::FaceBufs>();
    // Each face ships only the boundary layer (cutoff-deep slab, typically
    // 5-15% of owners; the cascade re-forwarding adds a little). A quarter
    // of capacity is generous; std::vector growth covers outliers.
    const int face_cap = dd->capacity / 4 + 1024;
    prealloc_bufs(*bufs_, face_cap);
    for (int f = 0; f < 6; ++f)
        dd->send_lists[f].reserve(static_cast<std::size_t>(face_cap));
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
/// `xfields` is the flattened list of declared extra forward fields; their
/// payload rides both message shapes (border and position-only refresh)
/// after the intrinsic columns.
inline void forward_axis_stage(Frame& frame, md::mpi::DomainDecomposition& dd,
                               FaceBufs& bufs, int dim, int sources_end,
                               double cutoff, bool full_rebuild, int& frame_cap,
                               int my_rank, double bonded_extension_for_log,
                               int& slot_in_out,
                               const std::vector<CommField>& xfields) {
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
        (full_rebuild ? kBorderDoublesPerAtom : kPosOnlyDoublesPerAtom)
        + extra_width(xfields);

    for (int s = 0; s < 2; ++s) {
        const int face_s = (s == 0) ? f_lo : f_hi;
        const int neighbor = (s == 0) ? neighbor_lo : neighbor_hi;
        auto& sbuf = bufs.send_pos[static_cast<std::size_t>(face_s)];
        auto& rbuf = bufs.recv_pos[static_cast<std::size_t>(face_s)];
        if (neighbor == MPI_PROC_NULL) {
            sbuf.clear();
        } else if (full_rebuild) {
            pack_positions_for_face(frame, dd, face_s, sbuf, xfields);
        } else {
            pack_positions_only_for_face(frame, dd, face_s, sbuf, xfields);
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
            MPI_Isend(
                bufs.send_pos[static_cast<std::size_t>(face_s)].data(),
                static_cast<int>(
                    bufs.send_pos[static_cast<std::size_t>(face_s)].size()),
                MPI_DOUBLE, neighbor, 100 + face_s, dd.cart_comm,
                &pos_reqs[n_preqs++]);
        }
        if (recv_counts[s] > 0) {
            MPI_Irecv(
                bufs.recv_pos[static_cast<std::size_t>(face_s)].data(),
                static_cast<int>(
                    bufs.recv_pos[static_cast<std::size_t>(face_s)].size()),
                MPI_DOUBLE, neighbor, 100 + opposite_face(face_s), dd.cart_comm,
                &pos_reqs[n_preqs++]);
        }
    }
    if (n_preqs > 0)
        MPI_Waitall(n_preqs, pos_reqs, MPI_STATUSES_IGNORE);

    if (full_rebuild) {
        // Capacity guard for this axis. Ranks no longer hold the global-N
        // backing store (partition_atoms shrinks to owned + headroom to avoid
        // per-rank memory scaling with the total system size), so a thin
        // subdomain can legitimately need more ghost rows than the headroom.
        // Grow on demand (explicit-copy realloc preserving the live rows)
        // instead of aborting; update the shared capacity mirrors and the
        // by-ref `frame_cap` so the remaining axis stages see the new size.
        const int new_ghosts = recv_counts[0] + recv_counts[1];
        if (slot_in_out + new_ghosts > frame_cap) {
            const int needed = slot_in_out + new_ghosts;
            // 50% slack so a cascade of small overflows doesn't realloc each axis.
            frame.atoms.set_capacity(needed + needed / 2 + 1024);
            frame_cap = frame.atoms.nmax;
            dd.capacity = frame.atoms.nmax;
            ATV_DEBUG("halo_exchange: rank {} grew ghost capacity to {} "
                      "(needed {}) at axis {}",
                      my_rank, frame_cap, needed, dim);
        }
    }

    // Unpack into ghost slots — for full rebuilds register tags too. The
    // unpack order ALSO defines the slot layout that `Reverse` mirrors
    // via `dd->recv_counts`.
    //
    // Hoisted AFTER the capacity guard above (a `set_capacity` there would
    // invalidate any span resolved before it) and dominating both the
    // full_rebuild and refresh unpack loops below.
    StorageReal* __restrict sp_x = frame.span(kAtom_x);
    StorageReal* __restrict sp_y = frame.span(kAtom_y);
    StorageReal* __restrict sp_z = frame.span(kAtom_z);
    TagInt* __restrict sp_tag = frame.span(kAtom_tag);
    StorageReal* __restrict sp_charge = frame.span(kAtom_charge);
    int slot = slot_in_out;
    // Ghost tag→slot tracking exists solely for bonded-tuple resolution
    // (ResolveView / select_owned_tuples). Pure pairwise systems never
    // query ghost tags, so skip the per-ghost hash inserts and the
    // companion ghost_tags log entirely — the map stays empty (~50 B per
    // ghost saved, plus the per-rebuild insert/erase churn).
    const bool track_tags = frame.atoms.has_bonded();
    for (int s = 0; s < 2; ++s) {
        const int face_s = (s == 0) ? f_lo : f_hi;
        const int n = recv_counts[s];
        const auto& rbuf = bufs.recv_pos[static_cast<std::size_t>(face_s)];
        if (full_rebuild) {
            for (int k = 0; k < n; ++k) {
                const std::size_t off =
                    static_cast<std::size_t>(k)
                    * static_cast<std::size_t>(per_atom_doubles);
                sp_x[slot] = static_cast<StorageReal>(rbuf[off + 0]);
                sp_y[slot] = static_cast<StorageReal>(rbuf[off + 1]);
                sp_z[slot] = static_cast<StorageReal>(rbuf[off + 2]);
                const std::int64_t tag =
                    static_cast<std::int64_t>(rbuf[off + 3]);
                sp_tag[slot] = static_cast<TagInt>(tag);
                sp_charge[slot] =
                    static_cast<StorageReal>(rbuf[off + 4]); // ghost charge (border)
                unpack_extra_fields(xfields, slot, rbuf,
                                    off + kBorderDoublesPerAtom);
                if (track_tags) {
                    dd.tag_to_local[tag] = slot;
                    dd.ghost_tags.push_back(tag);
                }
                ++slot;
            }
            dd.recv_counts[static_cast<std::size_t>(face_s)] = n;
        } else {
            for (int k = 0; k < n; ++k) {
                const std::size_t off =
                    static_cast<std::size_t>(k)
                    * static_cast<std::size_t>(per_atom_doubles);
                sp_x[slot] = static_cast<StorageReal>(rbuf[off + 0]);
                sp_y[slot] = static_cast<StorageReal>(rbuf[off + 1]);
                sp_z[slot] = static_cast<StorageReal>(rbuf[off + 2]);
                unpack_extra_fields(xfields, slot, rbuf,
                                    off + kPosOnlyDoublesPerAtom);
                ++slot;
            }
        }
    }
    slot_in_out = slot;
}

// Zero forces on ghost rows [n_local, n_local + n_ghost). ZeroForces only
// clears owned rows (zeroing the full backing store is ~global-N work per
// rank); the force kernels accumulate into ghost rows, and HaloRecv
// ships those contributions back to owners — so ghost forces must start each
// step at zero. Called by HaloSend after the current n_ghost is set.
void zero_ghost_forces(Frame& frame, const md::mpi::DomainDecomposition& dd,
                       const std::vector<CommField>& reverse_fields) {
    const int lo = dd.n_local;
    const int hi = dd.n_local + dd.n_ghost;
    AccReal* __restrict sp_fx = frame.span(kAtom_fx);
    AccReal* __restrict sp_fy = frame.span(kAtom_fy);
    AccReal* __restrict sp_fz = frame.span(kAtom_fz);
    for (int i = lo; i < hi; ++i) {
        sp_fx[i] = 0.0;
        sp_fy[i] = 0.0;
        sp_fz[i] = 0.0;
    }
    // Reverse extra fields share the force semantics: ghost rows are
    // partial-sum accumulators and must start each step at zero.
    for (const auto& f : reverse_fields) {
        for (int i = lo; i < hi; ++i) {
            const std::size_t base =
                static_cast<std::size_t>(i) * static_cast<std::size_t>(f.width);
            for (int k = 0; k < f.width; ++k)
                f.data[base + k] = StorageReal(0);
        }
    }
}

} // namespace

void HaloSend::launch(const cpu::LaunchContext& ctx) {
    if (ctx.runtime.mpi_size <= 1 || cutoff_ <= 0.0)
        return;
    auto* dd = ctx.runtime.md.dd;
    if (!dd || dd->cart_comm == MPI_COMM_NULL)
        return;

    Frame& frame = ctx.frame;
    FaceBufs& bufs = *bufs_;

    // One-time flatten of bound participants' field declarations. Deferred
    // to first launch (not build) because participant arrays are allocated
    // in the participants' own build(), whose plan order may follow ours.
    if (!fields_cached_) {
        append_declared_fields(participants_, /*forward=*/true, fwd_fields_);
        append_declared_fields(participants_, /*forward=*/false, rev_fields_);
        fields_cached_ = true;
    }

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
    // Mutable: forward_axis_stage grows the frame (and bumps frame_cap) when a
    // shrunk rank needs more ghost rows than its headroom.
    int frame_cap = static_cast<int>(frame.atoms.tag.size());

    if (full_rebuild) {
        // Re-check the build-time geometry guards: an NPT/deform run can
        // shrink the subdomain AFTER build() validated it (both the
        // shell-width and the dims==2 self-overlap condition are
        // box-dependent). Same fail-loud contract as build().
        for (int d = 0; d < 3; ++d) {
            if (dd->dims[d] <= 1)
                continue;
            const double sub = dd->subhi[d] - dd->sublo[d];
            if (halo_cut > sub)
                throw std::runtime_error(
                    "HaloSend: the box shrank until the ghost shell (" +
                    std::to_string(halo_cut) + " A) no longer fits axis " +
                    std::to_string(d) + "'s subdomain (" +
                    std::to_string(sub) + " A). See HaloSend::build for "
                    "the geometry contract.");
        }

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
                               bonded_extension_, slot, fwd_fields_);
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
        zero_ghost_forces(frame, *dd, rev_fields_);
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
                           bonded_extension_, slot, fwd_fields_);
    }
    zero_ghost_forces(frame, *dd, rev_fields_);
}

// ── HaloRefresh — a forward field refresh at an arbitrary phase ───────

void HaloRefresh::build(const cpu::BuildContext& /*ctx*/) {}

void HaloRefresh::launch(const cpu::LaunchContext& ctx) {
    auto* dd = ctx.runtime.md.dd;
    if (!dd || ctx.runtime.mpi_size <= 1 || dd->cart_comm == MPI_COMM_NULL)
        return;
    // NO early return on `n_ghost == 0`. Every rank must reach the Sendrecv below: a
    // rank with no ghosts of its own is still a NEIGHBOUR that some other rank is
    // waiting to hear from, and bailing out here deadlocks the whole run. (It did.)
    // Ranks with nothing to send simply exchange zero-length messages.

    Frame& frame = ctx.frame;
    const bool pos = (field_ == Field::positions);

    // Source arrays. Positions carry the periodic shift; velocities must NOT — a
    // ghost across the boundary lives at x - L, but it moves at the same v.
    StorageReal* __restrict sx = pos ? frame.span(kAtom_x) : frame.span(kAtom_vx);
    StorageReal* __restrict sy = pos ? frame.span(kAtom_y) : frame.span(kAtom_vy);
    StorageReal* __restrict sz = pos ? frame.span(kAtom_z) : frame.span(kAtom_vz);

    // Replay the cascade EXACTLY as the last full rebuild laid it out: axis 0,1,2;
    // within an axis, the -face then the +face; ghost slots appended from n_local.
    // Reusing dd->send_lists and dd->recv_counts is what makes this a refresh rather
    // than a second, conflicting halo.
    // Both faces of an axis go out CONCURRENTLY (Isend/Irecv + Waitall), never as two
    // sequential Sendrecv calls. With two ranks the lo and hi neighbour of an axis are
    // the SAME rank, and a sequential pair crosses its tags — each side blocks sending
    // face 0 while waiting to receive face 1. It deadlocks, and it deadlocked. The
    // rebuild path already does it this way, for exactly this reason.
    std::array<std::vector<double>, 2> sbuf, rbuf;
    int slot = dd->n_local;
    const int me = ctx.runtime.mpi_rank;

    for (int dim = 0; dim < 3; ++dim) {
        const int face[2] = {2 * dim, 2 * dim + 1};
        const double L = frame.box.length(dim);

        MPI_Request reqs[4];
        int n_req = 0;

        for (int s = 0; s < 2; ++s) {
            const int f = face[s];
            const int nbr = dd->face_neighbors[static_cast<std::size_t>(f)];
            const auto& list = dd->send_lists[static_cast<std::size_t>(f)];
            const int n_send = static_cast<int>(list.size());
            const int n_recv = dd->recv_counts[static_cast<std::size_t>(f)];

            const double shift = pos ? pbc_shift_for_face(*dd, f, L) : 0.0;
            auto& sb = sbuf[static_cast<std::size_t>(s)];
            auto& rb = rbuf[static_cast<std::size_t>(s)];
            sb.resize(static_cast<std::size_t>(n_send) * 3);
            for (int k = 0; k < n_send; ++k) {
                const int i = list[static_cast<std::size_t>(k)];
                sb[static_cast<std::size_t>(3 * k) + 0] = sx[i] + shift;
                sb[static_cast<std::size_t>(3 * k) + 1] = sy[i];
                sb[static_cast<std::size_t>(3 * k) + 2] = sz[i];
            }
            rb.assign(static_cast<std::size_t>(n_recv) * 3, 0.0);

            if (nbr == MPI_PROC_NULL || nbr == me) continue;
            MPI_Isend(sb.data(), n_send * 3, MPI_DOUBLE, nbr, f, dd->cart_comm,
                      &reqs[n_req++]);
            MPI_Irecv(rb.data(), n_recv * 3, MPI_DOUBLE, nbr, opposite_face(f),
                      dd->cart_comm, &reqs[n_req++]);
        }
        if (n_req > 0) MPI_Waitall(n_req, reqs, MPI_STATUSES_IGNORE);

        // Unpack in face order lo, hi — the slot layout the rebuild established.
        for (int s = 0; s < 2; ++s) {
            const int f = face[s];
            const int n_recv = dd->recv_counts[static_cast<std::size_t>(f)];
            const auto& rb = rbuf[static_cast<std::size_t>(s)];
            for (int k = 0; k < n_recv; ++k) {
                const std::size_t off = static_cast<std::size_t>(3 * k);
                sx[slot] = rb[off + 0];
                sy[slot] = rb[off + 1];
                sz[slot] = rb[off + 2];
                ++slot;
            }
        }
    }
}

void HaloRecv::build(const cpu::BuildContext& ctx) {
    auto* dd = ctx.runtime.md.dd;
    if (!dd || ctx.runtime.mpi_size <= 1)
        return;
    if (!bufs_)
        bufs_ = std::make_unique<detail::FaceBufs>();
    prealloc_bufs(*bufs_, dd->capacity);
}

void HaloRecv::launch(const cpu::LaunchContext& ctx) {
    if (ctx.runtime.mpi_size <= 1)
        return;
    auto* dd = ctx.runtime.md.dd;
    if (!dd || dd->cart_comm == MPI_COMM_NULL)
        return;
    if (dd->n_ghost == 0)
        return;

    Frame& frame = ctx.frame;
    FaceBufs& bufs = *bufs_;

    // Hoisted ONCE at the launch boundary; no `set_capacity`/`grow` happens
    // anywhere in this function (only HaloSend grows the frame), so these
    // stay valid across the whole reverse cascade below.
    AccReal* __restrict sp_fx = frame.span(kAtom_fx);
    AccReal* __restrict sp_fy = frame.span(kAtom_fy);
    AccReal* __restrict sp_fz = frame.span(kAtom_fz);

    // One-time flatten of bound participants' reverse-field declarations
    // (see HaloSend::launch for why this is deferred to first launch).
    if (!fields_cached_) {
        append_declared_fields(participants_, /*forward=*/false, rev_fields_);
        fields_cached_ = true;
    }
    const int per_ghost_doubles =
        kForceDoublesPerAtom + extra_width(rev_fields_);

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
            sbuf.resize(static_cast<std::size_t>(n)
                        * static_cast<std::size_t>(per_ghost_doubles));
            int slot = slot_start[face_s];
            for (int k = 0; k < n; ++k) {
                std::size_t off = static_cast<std::size_t>(k)
                                  * static_cast<std::size_t>(per_ghost_doubles);
                sbuf[off + 0] = sp_fx[slot];
                sbuf[off + 1] = sp_fy[slot];
                sbuf[off + 2] = sp_fz[slot];
                pack_extra_fields(rev_fields_, slot, sbuf,
                                  off + kForceDoublesPerAtom);
                ++slot;
            }
            bufs.recv_force[static_cast<std::size_t>(face_s)].assign(
                static_cast<std::size_t>(recv_force_counts[s])
                    * static_cast<std::size_t>(per_ghost_doubles),
                0.0);
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
            const auto& rbuf =
                bufs.recv_force[static_cast<std::size_t>(face_s)];
            const int n = static_cast<int>(list.size());
            for (int k = 0; k < n; ++k) {
                const int i = list[static_cast<std::size_t>(k)];
                std::size_t off = static_cast<std::size_t>(k)
                                  * static_cast<std::size_t>(per_ghost_doubles);
                sp_fx[i] += rbuf[off + 0];
                sp_fy[i] += rbuf[off + 1];
                sp_fz[i] += rbuf[off + 2];
                // Reverse extra fields: ACCUMULATE (exactly-once partial
                // sums), mirroring the intrinsic force semantics.
                off += 3;
                for (const auto& f : rev_fields_) {
                    const std::size_t base =
                        static_cast<std::size_t>(i)
                        * static_cast<std::size_t>(f.width);
                    for (int c = 0; c < f.width; ++c)
                        f.data[base + c] +=
                            static_cast<StorageReal>(rbuf[off++]);
                }
            }
        }
    }
}

#else // !ATV_HAS_MPI

void HaloSend::build(const cpu::BuildContext& /*ctx*/) {}
void HaloSend::launch(const cpu::LaunchContext& /*ctx*/) {}
void HaloRefresh::build(const cpu::BuildContext& /*ctx*/) {}
void HaloRefresh::launch(const cpu::LaunchContext& /*ctx*/) {}
void HaloRecv::build(const cpu::BuildContext& /*ctx*/) {}
void HaloRecv::launch(const cpu::LaunchContext& /*ctx*/) {}

#endif // ATV_HAS_MPI

} // namespace comm
} // namespace cpu
} // namespace ATV
