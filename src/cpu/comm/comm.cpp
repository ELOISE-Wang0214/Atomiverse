/// @file src/cpu/comm/comm.cpp
/// @brief cpu::comm::Comm -- impl of the LAMMPS-shaped four-method comm.

#include "atomiverse/cpu/comm/comm.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>

#ifdef ATV_HAS_MPI
#include <mpi.h>
#endif

#include "atomiverse/cpu/atom.h"
#include "atomiverse/cpu/box.h"
#include "atomiverse/cpu/comm/atom_pack.h"
#include "atomiverse/cpu/domain.h"

namespace ATV {
namespace cpu {
namespace comm {

// ── PIMPL state ────────────────────────────────────────────────────────────

struct CommState {
    /// Per-swap sendlist of local row indices that this rank ships across
    /// each of the 6 faces. Order: -x, +x, -y, +y, -z, +z. Refilled in
    /// `borders`; reused unchanged by `forward_comm` and `reverse_comm`.
    std::array<std::vector<int>, 6> sendlist;

    /// Number of ghosts received from each face in the most recent
    /// `borders` call. Defines the per-face ghost slot ranges that
    /// `reverse_comm` will ship forces back from.
    std::array<int, 6> recvnum{};

    /// Per-face send/recv buffers. Sized once on `build`, may grow on
    /// demand inside `borders` or `exchange`. Forward border uses
    /// `kBorderDoublesPerAtom` per atom; comm refresh uses
    /// `kCommDoublesPerAtom`; reverse uses 3 (force triple).
    std::array<std::vector<double>, 6> send_pos;
    std::array<std::vector<double>, 6> recv_pos;
    std::array<std::vector<double>, 6> send_force;
    std::array<std::vector<double>, 6> recv_force;
};

// ── Helpers ────────────────────────────────────────────────────────────────

namespace {

inline int face_dim(int f) noexcept { return f >> 1; }
inline int face_sign(int f) noexcept { return (f & 1) ? +1 : -1; }
inline int opposite_face(int f) noexcept { return f ^ 1; }

/// @brief PBC shift to apply when packing atoms across the periodic edge.
/// If this rank sits at the +dim edge and we're sending out the +dim
/// face, the receiver's coordinate frame is at the other periodic side,
/// so we pre-subtract `L`. Symmetric for the -dim edge.
inline double pbc_shift_for_face(const Domain& dom, int f, double L) noexcept {
    const int dim = face_dim(f);
    const int sign = face_sign(f);
    if (sign > 0 && dom.coords[dim] == dom.dims[dim] - 1)
        return -L;
    if (sign < 0 && dom.coords[dim] == 0)
        return +L;
    return 0.0;
}

/// @brief Copy every per-atom field from row `src` to row `dst`.
/// LAMMPS `AtomVec::copy(dst, src, delflag)` equivalent, minus the bonus
/// arrays (we don't have any).
void copy_atom_row(Atom& a, int dst, int src) {
    if (dst == src)
        return;
    a.x[dst] = a.x[src];
    a.y[dst] = a.y[src];
    a.z[dst] = a.z[src];
    a.vx[dst] = a.vx[src];
    a.vy[dst] = a.vy[src];
    a.vz[dst] = a.vz[src];
    a.fx[dst] = a.fx[src];
    a.fy[dst] = a.fy[src];
    a.fz[dst] = a.fz[src];
    a.mass[dst] = a.mass[src];
    a.charge[dst] = a.charge[src];
    a.type[dst] = a.type[src];
    a.tag[dst] = a.tag[src];
    a.image(dst, 0) = a.image(src, 0);
    a.image(dst, 1) = a.image(src, 1);
    a.image(dst, 2) = a.image(src, 2);

    a.num_bond[dst] = a.num_bond[src];
    for (int b = 0; b < Atom::kMaxBondPerAtom; ++b) {
        a.bond_atom(dst, b) = a.bond_atom(src, b);
        a.bond_type(dst, b) = a.bond_type(src, b);
    }
    a.num_angle[dst] = a.num_angle[src];
    for (int x = 0; x < 2 * Atom::kMaxAnglePerAtom; ++x)
        a.angle_atom(dst, x) = a.angle_atom(src, x);
    for (int x = 0; x < Atom::kMaxAnglePerAtom; ++x)
        a.angle_type(dst, x) = a.angle_type(src, x);
    a.num_dihedral[dst] = a.num_dihedral[src];
    for (int x = 0; x < 3 * Atom::kMaxDihedralPerAtom; ++x)
        a.dihedral_atom(dst, x) = a.dihedral_atom(src, x);
    for (int x = 0; x < Atom::kMaxDihedralPerAtom; ++x)
        a.dihedral_type(dst, x) = a.dihedral_type(src, x);
    a.nspecial(dst, 0) = a.nspecial(src, 0);
    a.nspecial(dst, 1) = a.nspecial(src, 1);
    a.nspecial(dst, 2) = a.nspecial(src, 2);
    for (int x = 0; x < Atom::kMaxSpecial; ++x)
        a.special(dst, x) = a.special(src, x);
}

/// @brief Orthorhombic PBC wrap of position `p` into `[0, L)`.
inline double wrap_to_cell(double p, double L) {
    if (L <= 0.0)
        return p;
    double q = p - std::floor(p / L) * L;
    if (q < 0.0)
        q += L;
    if (q >= L)
        q -= L;
    return q;
}

/// @brief Index of the slab along `dim` that contains wrapped position `q`,
/// under uniform brick decomposition.
inline int slab_uniform(double q, double L, int dims_d) {
    if (dims_d <= 1)
        return 0;
    const int r = static_cast<int>(q * static_cast<double>(dims_d) / L);
    return std::max(0, std::min(dims_d - 1, r));
}

/// @brief Index of the slab along `dim` under non-uniform RCB cuts.
inline int slab_cuts(double frac, const std::vector<double>& cuts,
                     int dims_d) {
    if (dims_d <= 1)
        return 0;
    auto it = std::upper_bound(cuts.begin(), cuts.end(), frac);
    int k = static_cast<int>(std::distance(cuts.begin(), it)) - 1;
    if (k < 0)
        k = 0;
    if (k > dims_d - 1)
        k = dims_d - 1;
    return k;
}

/// @brief Destination rank that owns wrapped position `(x, y, z)` under
/// the current `Domain` decomposition.
int dest_rank_for(const Domain& dom, double x, double y, double z,
                  const Box& box) {
    auto axis_coord = [&](int d, double q) {
        const double L = box.length(d);
        const double wrapped = wrap_to_cell(q, L);
        if (static_cast<int>(dom.cuts[d].size()) == dom.dims[d] + 1) {
            const double frac = (L > 0.0) ? (wrapped / L) : 0.0;
            return slab_cuts(frac, dom.cuts[d], dom.dims[d]);
        }
        return slab_uniform(wrapped, L, dom.dims[d]);
    };
    const int cx = axis_coord(0, x);
    const int cy = axis_coord(1, y);
    const int cz = axis_coord(2, z);
    return cx * (dom.dims[1] * dom.dims[2]) + cy * dom.dims[2] + cz;
}

} // namespace

// ── Comm ───────────────────────────────────────────────────────────────────

Comm::Comm(double cutghost)
    : cutghost_(cutghost), state_(std::make_unique<CommState>()) {}

Comm::~Comm() = default;
Comm::Comm(Comm&&) noexcept = default;
Comm& Comm::operator=(Comm&&) noexcept = default;

std::size_t Comm::send_buf_capacity_for_test() const {
    return state_ ? state_->send_pos[0].capacity() : 0;
}

void Comm::build(Atom& atom, Domain& /*domain*/) {
    if (!state_)
        state_ = std::make_unique<CommState>();
    // Conservative initial sizing: nlocal-sized send buffers per face.
    // Will grow on demand in `borders`.
    for (int f = 0; f < 6; ++f) {
        state_->sendlist[f].reserve(static_cast<std::size_t>(atom.nlocal));
        state_->send_pos[f].reserve(static_cast<std::size_t>(atom.nlocal)
                                    * kBorderDoublesPerAtom);
        state_->recv_pos[f].reserve(static_cast<std::size_t>(atom.nlocal)
                                    * kBorderDoublesPerAtom);
        state_->send_force[f].reserve(static_cast<std::size_t>(atom.nlocal)
                                      * 3);
        state_->recv_force[f].reserve(static_cast<std::size_t>(atom.nlocal)
                                      * 3);
    }
}

#ifdef ATV_HAS_MPI

void Comm::exchange(Atom& atom, Domain& domain, Box& box) {
    // Single-rank: nothing to migrate.
    int world_size = 1;
    int my_rank = 0;
    MPI_Comm world = domain.cart_comm;
    if (world == MPI_COMM_NULL) {
        // Single-rank fall-through. Still clear the map so borders can
        // rebuild it cleanly.
        atom.map_clear();
        atom.nghost = 0;
        return;
    }
    MPI_Comm_size(world, &world_size);
    MPI_Comm_rank(world, &my_rank);

    // map_clear before exchange (LAMMPS pattern). map will be rebuilt at
    // end of borders.
    atom.map_clear();
    atom.nghost = 0;

    // Classify each owner.
    std::vector<int> migrate_indices;
    std::vector<int> migrate_dest;
    migrate_indices.reserve(static_cast<std::size_t>(atom.nlocal) / 16 + 8);
    migrate_dest.reserve(migrate_indices.capacity());

    for (int i = 0; i < atom.nlocal; ++i) {
        const int dest = dest_rank_for(domain, atom.x[i], atom.y[i], atom.z[i],
                                       box);
        if (dest != my_rank) {
            migrate_indices.push_back(i);
            migrate_dest.push_back(dest);
        }
    }

    // Early-exit synchronization: if no rank has anything to migrate, skip
    // the collectives.
    int any_local = migrate_indices.empty() ? 0 : 1;
    int any_global = 0;
    MPI_Allreduce(&any_local, &any_global, 1, MPI_INT, MPI_LOR, world);
    if (any_global == 0)
        return;

    // Pack migrants into a single contiguous buffer. The per-rank payload
    // size in doubles is variable (per-atom topology varies), so we pack
    // into per-rank buffers first to know each rank's exact double count,
    // then concatenate.
    std::vector<std::vector<double>> per_rank_send(
        static_cast<std::size_t>(world_size));
    for (auto& v : per_rank_send)
        v.reserve(64); // typical small migration burst

    std::vector<double> scratch(kMaxExchangeDoubles);
    for (std::size_t k = 0; k < migrate_indices.size(); ++k) {
        const int i = migrate_indices[k];
        const int r = migrate_dest[k];
        const int n = pack_exchange(atom, i, scratch.data());
        auto& dst = per_rank_send[static_cast<std::size_t>(r)];
        dst.insert(dst.end(), scratch.begin(), scratch.begin() + n);
    }

    // Build the contiguous Alltoallv send buffer + send_counts (doubles).
    std::vector<int> send_counts(static_cast<std::size_t>(world_size), 0);
    std::vector<int> send_displs(static_cast<std::size_t>(world_size), 0);
    int total_send = 0;
    for (int r = 0; r < world_size; ++r) {
        send_displs[static_cast<std::size_t>(r)] = total_send;
        send_counts[static_cast<std::size_t>(r)] = static_cast<int>(
            per_rank_send[static_cast<std::size_t>(r)].size());
        total_send += send_counts[static_cast<std::size_t>(r)];
    }
    std::vector<double> send_buf(static_cast<std::size_t>(total_send));
    for (int r = 0; r < world_size; ++r) {
        const auto& src = per_rank_send[static_cast<std::size_t>(r)];
        std::memcpy(send_buf.data() + send_displs[static_cast<std::size_t>(r)],
                    src.data(), src.size() * sizeof(double));
    }

    // Exchange counts (in doubles).
    std::vector<int> recv_counts(static_cast<std::size_t>(world_size), 0);
    std::vector<int> recv_displs(static_cast<std::size_t>(world_size), 0);
    MPI_Alltoall(send_counts.data(), 1, MPI_INT, recv_counts.data(), 1, MPI_INT,
                 world);
    int total_recv = 0;
    for (int r = 0; r < world_size; ++r) {
        recv_displs[static_cast<std::size_t>(r)] = total_recv;
        total_recv += recv_counts[static_cast<std::size_t>(r)];
    }
    std::vector<double> recv_buf(static_cast<std::size_t>(total_recv));

    // Payload exchange.
    MPI_Alltoallv(send_buf.data(), send_counts.data(), send_displs.data(),
                  MPI_DOUBLE, recv_buf.data(), recv_counts.data(),
                  recv_displs.data(), MPI_DOUBLE, world);

    // Compact survivors via gap-fill (LAMMPS pattern: copy last atom to
    // freed slot, decrement nlocal). Iterate the migrate list in
    // descending order so each `copy_atom_row` source (currently last
    // owner) is always a non-migrating atom.
    std::sort(migrate_indices.rbegin(), migrate_indices.rend());
    for (int i : migrate_indices) {
        const int last = atom.nlocal - 1;
        copy_atom_row(atom, /*dst=*/i, /*src=*/last);
        --atom.nlocal;
    }

    // Unpack received atoms into the freed/extended tail.
    // First grow nmax if the post-receive total exceeds current capacity.
    // Estimate receive count by scanning recv_buf for atom-boundary
    // markers — but pack_exchange doesn't write a per-atom marker, so we
    // instead unpack incrementally and grow on demand.
    std::size_t off = 0;
    while (off < recv_buf.size()) {
        if (atom.nlocal + 1 > atom.nmax)
            atom.grow(atom.nlocal + 1);
        const int n_read = unpack_exchange(atom, atom.nlocal, recv_buf.data()
                                                                  + off);
        off += static_cast<std::size_t>(n_read);
        ++atom.nlocal;
    }
}

void Comm::borders(Atom& atom, Domain& domain, Box& box) {
    int world_size = 1;
    int my_rank = 0;
    MPI_Comm world = domain.cart_comm;
    if (world == MPI_COMM_NULL) {
        // Single-rank: no ghosts; just rebuild the map for the owner range.
        atom.nghost = 0;
        atom.map_clear();
        atom.map_set(0, atom.nlocal);
        // Clear per-swap state too.
        for (int f = 0; f < 6; ++f) {
            state_->sendlist[f].clear();
            state_->recvnum[f] = 0;
        }
        return;
    }
    MPI_Comm_size(world, &world_size);
    MPI_Comm_rank(world, &my_rank);

    // Cascade: each axis scans [nfirst, nlast) of the atom array. nlast
    // grows after each axis completes so the next axis can forward
    // earlier-axis ghosts to its neighbors. The scan must include all
    // atoms (owners + already-received ghosts), so nfirst stays 0.
    int slot = atom.nlocal;
    int nfirst = 0;
    int nlast = atom.nlocal;
    for (int dim = 0; dim < 3; ++dim) {
        for (int direction = 0; direction < 2; ++direction) {
            const int iswap = 2 * dim + direction;
            const int neighbor = domain.face_neighbors[iswap];

            // Rebuild sendlist[iswap] by scanning [nfirst, nlast).
            auto& sendlist = state_->sendlist[iswap];
            sendlist.clear();
            const double lo = domain.sublo[dim];
            const double hi = domain.subhi[dim];
            const double* q = (dim == 0)   ? atom.x.data()
                              : (dim == 1) ? atom.y.data()
                                           : atom.z.data();
            if (direction == 0) {
                for (int i = nfirst; i < nlast; ++i)
                    if (q[i] - lo < cutghost_)
                        sendlist.push_back(i);
            } else {
                for (int i = nfirst; i < nlast; ++i)
                    if (hi - q[i] < cutghost_)
                        sendlist.push_back(i);
            }
            const int n_send = static_cast<int>(sendlist.size());

            // Pack positions + tag + type + image for each atom in sendlist.
            // PBC shift applied for the periodic-edge faces.
            const double L = box.length(dim);
            const double shift = pbc_shift_for_face(domain, iswap, L);
            const double sx = (dim == 0) ? shift : 0.0;
            const double sy = (dim == 1) ? shift : 0.0;
            const double sz = (dim == 2) ? shift : 0.0;
            auto& sbuf = state_->send_pos[iswap];
            sbuf.resize(static_cast<std::size_t>(n_send)
                        * kBorderDoublesPerAtom);
            pack_border(atom, sendlist.data(), n_send, sx, sy, sz, sbuf.data());

            // Exchange counts.
            int n_recv = 0;
            if (neighbor != MPI_PROC_NULL && neighbor != my_rank) {
                MPI_Sendrecv(&n_send, 1, MPI_INT, neighbor, iswap, &n_recv, 1,
                             MPI_INT, neighbor, opposite_face(iswap), world,
                             MPI_STATUS_IGNORE);
            } else if (neighbor == my_rank) {
                // Self-neighbor (single rank along this axis). LAMMPS-style:
                // ship to self.
                n_recv = n_send;
            }
            state_->recvnum[iswap] = n_recv;

            // Grow atom storage if needed (LAMMPS unpack_border calls grow
            // internally; we do it explicitly here, same effect).
            if (slot + n_recv > atom.nmax)
                atom.grow(slot + n_recv);

            // Exchange payloads.
            auto& rbuf = state_->recv_pos[iswap];
            rbuf.resize(static_cast<std::size_t>(n_recv)
                        * kBorderDoublesPerAtom);
            if (neighbor != MPI_PROC_NULL && neighbor != my_rank) {
                MPI_Sendrecv(sbuf.data(),
                             static_cast<int>(sbuf.size()), MPI_DOUBLE,
                             neighbor, 100 + iswap, rbuf.data(),
                             static_cast<int>(rbuf.size()), MPI_DOUBLE,
                             neighbor, 100 + opposite_face(iswap), world,
                             MPI_STATUS_IGNORE);
            } else if (neighbor == my_rank) {
                std::memcpy(rbuf.data(), sbuf.data(),
                            sbuf.size() * sizeof(double));
            }

            // Unpack into ghost slots [slot, slot + n_recv).
            unpack_border(atom, n_recv, slot, rbuf.data());
            slot += n_recv;
        }
        // End of axis: include all ghosts received so far in the next
        // axis's scan range.
        nfirst = 0;
        nlast = slot;
    }
    atom.nghost = slot - atom.nlocal;

    // Rebuild the tag→local map over owners + ghosts (LAMMPS Atom::map_set).
    atom.map_clear();
    atom.map_set(0, atom.nlocal + atom.nghost);
}

void Comm::forward_comm(Atom& atom, Domain& domain, Box& box) {
    int my_rank = 0;
    MPI_Comm world = domain.cart_comm;
    if (world == MPI_COMM_NULL || atom.nghost == 0)
        return;
    MPI_Comm_rank(world, &my_rank);

    // Reuse the existing sendlist (set by the last `borders`). Walk the
    // same 6 swaps; ghost slot ranges are stable so the receiver writes
    // back into the same rows.
    int slot = atom.nlocal;
    for (int dim = 0; dim < 3; ++dim) {
        for (int direction = 0; direction < 2; ++direction) {
            const int iswap = 2 * dim + direction;
            const auto& sendlist = state_->sendlist[iswap];
            const int n_send = static_cast<int>(sendlist.size());
            const int n_recv = state_->recvnum[iswap];
            const int neighbor = domain.face_neighbors[iswap];

            const double L = box.length(dim);
            const double shift = pbc_shift_for_face(domain, iswap, L);
            const double sx = (dim == 0) ? shift : 0.0;
            const double sy = (dim == 1) ? shift : 0.0;
            const double sz = (dim == 2) ? shift : 0.0;
            auto& sbuf = state_->send_pos[iswap];
            sbuf.resize(static_cast<std::size_t>(n_send) * kCommDoublesPerAtom);
            pack_comm(atom, sendlist.data(), n_send, sx, sy, sz, sbuf.data());

            auto& rbuf = state_->recv_pos[iswap];
            rbuf.resize(static_cast<std::size_t>(n_recv) * kCommDoublesPerAtom);
            if (neighbor != MPI_PROC_NULL && neighbor != my_rank) {
                MPI_Sendrecv(sbuf.data(),
                             static_cast<int>(sbuf.size()), MPI_DOUBLE,
                             neighbor, 200 + iswap, rbuf.data(),
                             static_cast<int>(rbuf.size()), MPI_DOUBLE,
                             neighbor, 200 + opposite_face(iswap), world,
                             MPI_STATUS_IGNORE);
            } else if (neighbor == my_rank) {
                std::memcpy(rbuf.data(), sbuf.data(),
                            sbuf.size() * sizeof(double));
            }

            unpack_comm(atom, n_recv, slot, rbuf.data());
            slot += n_recv;
        }
    }
}

void Comm::reverse_comm(Atom& atom, Domain& domain) {
    int my_rank = 0;
    MPI_Comm world = domain.cart_comm;
    if (world == MPI_COMM_NULL || atom.nghost == 0)
        return;
    MPI_Comm_rank(world, &my_rank);

    // Per-face ghost slot starts: slot_start[f] = nlocal + sum_{g<f} recvnum[g].
    int slot_start[6];
    {
        int s = atom.nlocal;
        for (int f = 0; f < 6; ++f) {
            slot_start[f] = s;
            s += state_->recvnum[f];
        }
    }

    // Reverse the swap order: dim 2 → 0. For each swap, pack ghost forces
    // out of this rank's [slot_start[iswap], slot_start[iswap] + recvnum)
    // range; send to the rank we received them from (= face_neighbors[iswap]
    // = same as borders); accumulate into f[sendlist[iswap]] on the
    // receiver side.
    for (int dim = 2; dim >= 0; --dim) {
        for (int direction = 0; direction < 2; ++direction) {
            const int iswap = 2 * dim + direction;
            const int neighbor = domain.face_neighbors[iswap];
            const int n_ghost_local = state_->recvnum[iswap];
            const auto& owner_indices = state_->sendlist[iswap];
            const int n_owner_local = static_cast<int>(owner_indices.size());

            auto& sbuf = state_->send_force[iswap];
            sbuf.resize(static_cast<std::size_t>(n_ghost_local) * 3);
            {
                int slot = slot_start[iswap];
                int m = 0;
                for (int k = 0; k < n_ghost_local; ++k) {
                    sbuf[m++] = atom.fx[slot];
                    sbuf[m++] = atom.fy[slot];
                    sbuf[m++] = atom.fz[slot];
                    ++slot;
                }
            }

            auto& rbuf = state_->recv_force[iswap];
            rbuf.assign(static_cast<std::size_t>(n_owner_local) * 3, 0.0);
            if (neighbor != MPI_PROC_NULL && neighbor != my_rank) {
                MPI_Sendrecv(sbuf.data(),
                             static_cast<int>(sbuf.size()), MPI_DOUBLE,
                             neighbor, 300 + iswap, rbuf.data(),
                             static_cast<int>(rbuf.size()), MPI_DOUBLE,
                             neighbor, 300 + opposite_face(iswap), world,
                             MPI_STATUS_IGNORE);
            } else if (neighbor == my_rank) {
                std::memcpy(rbuf.data(), sbuf.data(),
                            sbuf.size() * sizeof(double));
            }

            // Accumulate received forces into the owner rows.
            int m = 0;
            for (int k = 0; k < n_owner_local; ++k) {
                const int i = owner_indices[static_cast<std::size_t>(k)];
                atom.fx[i] += rbuf[m++];
                atom.fy[i] += rbuf[m++];
                atom.fz[i] += rbuf[m++];
            }
        }
    }
}

#else // !ATV_HAS_MPI

void Comm::exchange(Atom& atom, Domain& /*domain*/, Box& /*box*/) {
    atom.map_clear();
    atom.nghost = 0;
}

void Comm::borders(Atom& atom, Domain& /*domain*/, Box& /*box*/) {
    atom.nghost = 0;
    atom.map_clear();
    atom.map_set(0, atom.nlocal);
    for (int f = 0; f < 6; ++f) {
        state_->sendlist[f].clear();
        state_->recvnum[f] = 0;
    }
}

void Comm::forward_comm(Atom& /*atom*/, Domain& /*domain*/, Box& /*box*/) {}
void Comm::reverse_comm(Atom& /*atom*/, Domain& /*domain*/) {}

#endif // ATV_HAS_MPI

} // namespace comm
} // namespace cpu
} // namespace ATV
