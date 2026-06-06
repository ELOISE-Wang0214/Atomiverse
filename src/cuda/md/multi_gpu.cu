/// @file cuda/md/multi_gpu.cu
/// @brief MultiGpu orchestrator — partition / gather (Phase-4 foundation).

#include "atomiverse/cuda/md/multi_gpu.cuh"

#include "atomiverse/cpu/frame.h"
#include "atomiverse/cpu/transfer.h"
#include "atomiverse/cuda/box.cuh"
#include "atomiverse/cuda/transfer.cuh"
#include "atomiverse/cuda/types.cuh"

#include <cuda_runtime.h>

#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

namespace ATV {
namespace cuda {
namespace md {

#define CKM(x)                                                                  \
    do {                                                                        \
        cudaError_t e_ = (x);                                                   \
        if (e_ != cudaSuccess)                                                  \
            throw std::runtime_error(std::string("MultiGpu CUDA: ") + #x        \
                                     + ": " + cudaGetErrorString(e_));          \
    } while (0)

namespace {
// Stream-ordered POOLED allocation. The per-reneighbor halo/migrate scratch was
// the multi-GPU bottleneck: ~1000 cudaMalloc/cudaFree calls/run = ~600 ms of
// host time (nsys). cudaMallocAsync/cudaFreeAsync use the device's default
// memory pool, which (with the release threshold raised — see
// set_pool_retain) RETAINS freed blocks, so after the first reneighbor these
// allocations hit the pool and are nearly free. Templated wrapper so the
// existing `&typed_ptr` call sites need no void** casts.
template <class T> inline cudaError_t mallocp(T** p, size_t bytes) {
    return cudaMallocAsync(reinterpret_cast<void**>(p), bytes, cudaStream_t{0});
}
inline void freep(void* p) {
    if (p)
        cudaFreeAsync(p, cudaStream_t{0});
}
// Raise the default mem-pool release threshold on every visible device so freed
// async blocks stay in the pool (default threshold 0 returns them to the OS on
// the next sync, defeating reuse). Idempotent; called from setup_force_engine.
inline void set_pool_retain(const DeviceManager& dm) {
    for (int i = 0; i < dm.count(); ++i) {
        cudaMemPool_t pool;
        if (cudaDeviceGetDefaultMemPool(&pool, dm.device_id(i)) != cudaSuccess)
            continue;
        std::uint64_t thr = ~0ull;
        cudaMemPoolSetAttribute(pool, cudaMemPoolAttrReleaseThreshold, &thr);
    }
}
} // namespace

MultiGpu::MultiGpu(int n_gpus) : dm_(n_gpus) {}

MultiGpu::~MultiGpu() { destroy_frames(); }

void MultiGpu::destroy_frames() {
    destroy_force_engine(); // frees per-subdomain NeighborList buffers (needs frames_)
    destroy_halo_plan();
    for (int s = 0; s < static_cast<int>(frames_.size()); ++s) {
        dm_.set_device(s);
        cuda::destroy(frames_[s]);
    }
    frames_.clear();
    n_local_.clear();
}

int MultiGpu::n_local(int slot) const { return n_local_[slot]; }
int MultiGpu::n_ghost(int slot) const { return frames_[slot].n_ghost_; }

namespace {
// Select atoms within `cutoff` of a subdomain boundary along `axis` and
// compact their x/y/z/tag/type into the send buffers (atomic compaction). Also
// records each selected atom's array INDEX in `oidx` so a later forward_comm
// can re-gather its position. `mode`: 0 = near the low (sublo) boundary only,
// 1 = near the high (subhi) boundary only, 2 = near EITHER boundary, sent ONCE
// (used when both face-neighbors are the same GPU, i.e. dims==2 along this
// axis: with the min-image force kernel a single raw-position ghost serves all
// periodic images, so the atom must not be sent twice or its force doubles).
__global__ void pack_face(const real_t* __restrict__ x,
                          const real_t* __restrict__ y,
                          const real_t* __restrict__ z,
                          const std::int64_t* __restrict__ tag,
                          const int* __restrict__ type, int n_local, int axis,
                          double lo, double hi, double cutoff, int mode,
                          real_t* __restrict__ ox, real_t* __restrict__ oy,
                          real_t* __restrict__ oz,
                          std::int64_t* __restrict__ otag,
                          int* __restrict__ otype, int* __restrict__ oidx,
                          int* __restrict__ counter, int cap) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n_local)
        return;
    real_t q = (axis == 0) ? x[i] : (axis == 1) ? y[i] : z[i];
    bool near_lo = q < lo + cutoff, near_hi = q > hi - cutoff;
    bool sel = (mode == 0) ? near_lo : (mode == 1) ? near_hi : (near_lo || near_hi);
    if (!sel)
        return;
    int s = atomicAdd(counter, 1);
    if (s < cap) {
        ox[s] = x[i];
        oy[s] = y[i];
        oz[s] = z[i];
        otag[s] = tag[i];
        otype[s] = type[i];
        oidx[s] = i;
    }
}

// Gather current positions of `cnt` owned atoms named by `idx` into the
// contiguous send buffers (forward_comm position refresh).
__global__ void gather_positions(const real_t* __restrict__ x,
                                 const real_t* __restrict__ y,
                                 const real_t* __restrict__ z,
                                 const int* __restrict__ idx, int cnt,
                                 real_t* __restrict__ ox,
                                 real_t* __restrict__ oy,
                                 real_t* __restrict__ oz) {
    int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= cnt)
        return;
    int i = idx[k];
    ox[k] = x[i];
    oy[k] = y[i];
    oz[k] = z[i];
}

// Reverse halo (force) — gather a destination's ghost-row forces [base,base+cnt)
// into a contiguous buffer for P2P shipment back to the owner subdomain. Mirror
// of gather_positions but reading the (contiguous) ghost force rows by offset.
__global__ void gather_ghost_forces(const real_t* __restrict__ fx,
                                    const real_t* __restrict__ fy,
                                    const real_t* __restrict__ fz, int base,
                                    int cnt, real_t* __restrict__ ofx,
                                    real_t* __restrict__ ofy,
                                    real_t* __restrict__ ofz) {
    int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= cnt)
        return;
    ofx[k] = fx[base + k];
    ofy[k] = fy[base + k];
    ofz[k] = fz[base + k];
}

// Reverse halo (force) — scatter-add received ghost forces into the source
// subdomain's owned (or prior-axis ghost) rows named by `idx` (the same src
// indices recorded for the forward position gather). atomicAdd because several
// transfers may target the same owned row, and because a prior-axis ghost row
// can receive contributions from multiple later-axis plans before being
// reversed further down the cascade.
__global__ void scatter_add_forces(real_t* __restrict__ fx,
                                   real_t* __restrict__ fy,
                                   real_t* __restrict__ fz,
                                   const int* __restrict__ idx, int cnt,
                                   const real_t* __restrict__ ifx,
                                   const real_t* __restrict__ ify,
                                   const real_t* __restrict__ ifz) {
    int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= cnt)
        return;
    int i = idx[k];
    atomicAdd(&fx[i], ifx[k]);
    atomicAdd(&fy[i], ify[k]);
    atomicAdd(&fz[i], ifz[k]);
}
} // namespace

void MultiGpu::destroy_halo_plan() {
    for (auto& p : halo_plan_) {
        dm_.set_device(p.src);
        freep(p.d_src_idx);
        freep(p.d_sx);
        freep(p.d_sy);
        freep(p.d_sz);
        // Reverse-halo gather buffers were allocated on the dst device.
        dm_.set_device(p.dst);
        freep(p.d_fx_dst);
        freep(p.d_fy_dst);
        freep(p.d_fz_dst);
    }
    halo_plan_.clear();
    // The freep's above are stream-ordered (cudaFreeAsync). Drain them so a
    // following teardown / next MultiGpu (the per-device pool is shared across
    // instances) can't reuse a block whose async free is still in flight.
    for (int s = 0; s < count(); ++s) {
        dm_.set_device(s);
        cudaDeviceSynchronize();
    }
}

void MultiGpu::build_halo(double cutoff) {
    const int N = count();
    destroy_halo_plan(); // free any prior reneighbor-interval plan
    // Reset ghost state.
    std::vector<int> goff(N, 0);
    for (int s = 0; s < N; ++s) {
        frames_[s].n_ghost_ = 0;
        n_local_[s] = frames_[s].n_local_;
    }

    // 3-axis cascade: x, then y, then z. The y-pass forwards owned+x-ghosts and
    // the z-pass forwards owned+x+y-ghosts, so edge/corner atoms reach
    // diagonal-neighbor GPUs. The sendable set per subdomain is snapshotted at
    // each axis start (owned + ghosts from prior axes); ghosts received during
    // this axis are appended beyond that snapshot and so are not re-forwarded.
    for (int axis = 0; axis < 3; ++axis) {
        std::vector<int> n_send(N);
        for (int s = 0; s < N; ++s)
            n_send[s] = frames_[s].n_local_ + goff[s];

        for (int src = 0; src < N; ++src) {
            const int ns = n_send[src];
            const int dlo = dd_->face_neighbor(src, 2 * axis);
            const int dhi = dd_->face_neighbor(src, 2 * axis + 1);
            if (ns == 0 || (dlo == src && dhi == src))
                continue; // empty, or dims==1 along this axis (periodic self)
            dm_.set_device(src);
            // Per-source scratch sized to the current sendable set.
            real_t *sx, *sy, *sz;
            std::int64_t* stag;
            int *stype, *sidx, *scnt;
            CKM(mallocp(&sx, ns * sizeof(real_t)));
            CKM(mallocp(&sy, ns * sizeof(real_t)));
            CKM(mallocp(&sz, ns * sizeof(real_t)));
            CKM(mallocp(&stag, ns * sizeof(std::int64_t)));
            CKM(mallocp(&stype, ns * sizeof(int)));
            CKM(mallocp(&sidx, ns * sizeof(int)));
            CKM(mallocp(&scnt, sizeof(int)));
            auto& a = frames_[src].atoms;
            const auto lo = dd_->sublo(src);
            const auto hi = dd_->subhi(src);

            // Distinct (dst, mode) sends along this axis. When the low and high
            // face-neighbors are the SAME GPU (dims==2 along the axis, periodic
            // self-pair), an atom near either boundary must be sent to it just
            // ONCE (mode 2) — the min-image force kernel resolves all periodic
            // images from a single raw-position ghost, so a second copy would
            // double that atom's force. Distinct neighbors (dims>=3) get
            // separate low/high sends.
            struct Send {
                int dst, mode;
            };
            Send sends[2];
            int n_sends = 0;
            if (dlo != src && dlo == dhi) {
                sends[n_sends++] = {dlo, 2};
            } else {
                if (dlo != src)
                    sends[n_sends++] = {dlo, 0};
                if (dhi != src)
                    sends[n_sends++] = {dhi, 1};
            }

            for (int si = 0; si < n_sends; ++si) {
                const int dst = sends[si].dst;
                const int mode = sends[si].mode;

                CKM(cudaMemset(scnt, 0, sizeof(int)));
                int block = 256, grid = (ns + block - 1) / block;
                pack_face<<<grid, block>>>(a.x, a.y, a.z, a.tag, a.type, ns, axis,
                                           lo[axis], hi[axis], cutoff, mode, sx,
                                           sy, sz, stag, stype, sidx, scnt, ns);
                CKM(cudaGetLastError());
                int cnt = 0; // D2H of the counter also syncs the pack kernel
                CKM(cudaMemcpy(&cnt, scnt, sizeof(int), cudaMemcpyDeviceToHost));
                if (cnt <= 0)
                    continue;

                // Append into dst's ghost region [n_local + goff, ...).
                auto& da = frames_[dst].atoms;
                const int base = frames_[dst].n_local_ + goff[dst];
                if (base + cnt > frames_[dst].capacity_)
                    throw std::runtime_error(
                        "MultiGpu::build_halo: ghost capacity exceeded on slot "
                        + std::to_string(dst) + " (raise ghost_headroom)");
                const int sd = dm_.device_id(src), dd_id = dm_.device_id(dst);
                CKM(cudaMemcpyPeer(da.x + base, dd_id, sx, sd,
                                   cnt * sizeof(real_t)));
                CKM(cudaMemcpyPeer(da.y + base, dd_id, sy, sd,
                                   cnt * sizeof(real_t)));
                CKM(cudaMemcpyPeer(da.z + base, dd_id, sz, sd,
                                   cnt * sizeof(real_t)));
                CKM(cudaMemcpyPeer(da.tag + base, dd_id, stag, sd,
                                   cnt * sizeof(std::int64_t)));
                CKM(cudaMemcpyPeer(da.type + base, dd_id, stype, sd,
                                   cnt * sizeof(int)));

                // Persist this transfer (src device): the selected src-array
                // indices and a contiguous position send buffer, for
                // forward_comm replay. Indices may reference prior-axis ghost
                // slots (the cascade), which forward_comm refreshes in order.
                HaloPlan p;
                p.src = src;
                p.dst = dst;
                p.axis = axis;
                p.cnt = cnt;
                p.base = base;
                CKM(mallocp(&p.d_src_idx, cnt * sizeof(int)));
                CKM(mallocp(&p.d_sx, cnt * sizeof(real_t)));
                CKM(mallocp(&p.d_sy, cnt * sizeof(real_t)));
                CKM(mallocp(&p.d_sz, cnt * sizeof(real_t)));
                CKM(cudaMemcpy(p.d_src_idx, sidx, cnt * sizeof(int),
                               cudaMemcpyDeviceToDevice));
                CKM(cudaMemcpy(p.d_sx, sx, cnt * sizeof(real_t),
                               cudaMemcpyDeviceToDevice));
                CKM(cudaMemcpy(p.d_sy, sy, cnt * sizeof(real_t),
                               cudaMemcpyDeviceToDevice));
                CKM(cudaMemcpy(p.d_sz, sz, cnt * sizeof(real_t),
                               cudaMemcpyDeviceToDevice));
                // Reverse-halo gather buffers live on the DST device (where the
                // ghost forces are read); allocate them there, then restore the
                // src device for the rest of this source's loop.
                dm_.set_device(dst);
                CKM(mallocp(&p.d_fx_dst, cnt * sizeof(real_t)));
                CKM(mallocp(&p.d_fy_dst, cnt * sizeof(real_t)));
                CKM(mallocp(&p.d_fz_dst, cnt * sizeof(real_t)));
                dm_.set_device(src);
                halo_plan_.push_back(p);

                goff[dst] += cnt;
                frames_[dst].n_ghost_ = goff[dst];
            }
            cudaSetDevice(dm_.device_id(src));
            freep(sx);
            freep(sy);
            freep(sz);
            freep(stag);
            freep(stype);
            freep(sidx);
            freep(scnt);
        }
        // Barrier between axes: received ghosts must be visible (and their
        // positions settled) before the next axis packs from them.
        for (int s = 0; s < N; ++s) {
            dm_.set_device(s);
            cudaDeviceSynchronize();
        }
    }
}

void MultiGpu::forward_comm() {
    // Replay the cascade IN AXIS ORDER with a barrier between axes: a y-axis
    // plan may gather from x-ghost slots, whose positions the x-axis plans must
    // have refreshed first. Plans were recorded axis-major, but barrier on the
    // axis boundary explicitly rather than rely on push order.
    for (int axis = 0; axis < 3; ++axis) {
        bool any = false;
        for (auto& p : halo_plan_) {
            if (p.axis != axis || p.cnt <= 0)
                continue;
            any = true;
            dm_.set_device(p.src);
            auto& a = frames_[p.src].atoms;
            int block = 256, grid = (p.cnt + block - 1) / block;
            gather_positions<<<grid, block>>>(a.x, a.y, a.z, p.d_src_idx, p.cnt,
                                              p.d_sx, p.d_sy, p.d_sz);
            CKM(cudaGetLastError());
            auto& da = frames_[p.dst].atoms;
            const int sd = dm_.device_id(p.src), dd_id = dm_.device_id(p.dst);
            CKM(cudaMemcpyPeer(da.x + p.base, dd_id, p.d_sx, sd,
                               p.cnt * sizeof(real_t)));
            CKM(cudaMemcpyPeer(da.y + p.base, dd_id, p.d_sy, sd,
                               p.cnt * sizeof(real_t)));
            CKM(cudaMemcpyPeer(da.z + p.base, dd_id, p.d_sz, sd,
                               p.cnt * sizeof(real_t)));
        }
        if (any)
            for (int s = 0; s < count(); ++s) {
                dm_.set_device(s);
                cudaDeviceSynchronize();
            }
    }
}

void MultiGpu::reverse_comm() {
    // Adjoint of forward_comm: replay the cascade in REVERSE axis order
    // (z, then y, then x) so a ghost force lands on its forwarding row before
    // that row is itself reversed one axis down. For each transfer: gather the
    // dst's ghost-row forces [base, base+cnt), P2P them back to the src device
    // (reusing the src-side d_s{x,y,z} buffers as the receive buffer), then
    // scatter-add into the src rows named by d_src_idx.
    for (int axis = 2; axis >= 0; --axis) {
        bool any = false;
        for (auto& p : halo_plan_) {
            if (p.axis != axis || p.cnt <= 0)
                continue;
            any = true;
            const int sd = dm_.device_id(p.src), dd_id = dm_.device_id(p.dst);
            int block = 256, grid = (p.cnt + block - 1) / block;

            // Gather the dst's ghost-row forces into the dst-device buffer.
            dm_.set_device(p.dst);
            auto& da = frames_[p.dst].atoms;
            gather_ghost_forces<<<grid, block>>>(da.fx, da.fy, da.fz, p.base,
                                                 p.cnt, p.d_fx_dst, p.d_fy_dst,
                                                 p.d_fz_dst);
            CKM(cudaGetLastError());

            // Ship dst → src (into the reused src-side receive buffers).
            CKM(cudaMemcpyPeer(p.d_sx, sd, p.d_fx_dst, dd_id,
                               p.cnt * sizeof(real_t)));
            CKM(cudaMemcpyPeer(p.d_sy, sd, p.d_fy_dst, dd_id,
                               p.cnt * sizeof(real_t)));
            CKM(cudaMemcpyPeer(p.d_sz, sd, p.d_fz_dst, dd_id,
                               p.cnt * sizeof(real_t)));

            // Scatter-add into the src owned (or prior-axis ghost) rows.
            dm_.set_device(p.src);
            auto& a = frames_[p.src].atoms;
            scatter_add_forces<<<grid, block>>>(a.fx, a.fy, a.fz, p.d_src_idx,
                                                p.cnt, p.d_sx, p.d_sy, p.d_sz);
            CKM(cudaGetLastError());
        }
        // Barrier between axes: this axis's scatter-adds into prior-axis ghost
        // rows must be visible before the next (lower) axis gathers from them.
        if (any)
            for (int s = 0; s < count(); ++s) {
                dm_.set_device(s);
                cudaDeviceSynchronize();
            }
    }
}

void MultiGpu::partition(const cpu::Frame& global, double ghost_headroom) {
    destroy_frames();
    const int N = dm_.count();
    const int ng = global.atoms.size();
    n_global_ = ng;

    // Global orthorhombic box edge lengths.
    const double Lx = global.box(0, 0);
    const double Ly = global.box(1, 1);
    const double Lz = global.box(2, 2);
    box_L_[0] = Lx;
    box_L_[1] = Ly;
    box_L_[2] = Lz;
    dd_ = std::make_unique<DomainDecomposition>(N, Lx, Ly, Lz);

    // Bucket global atom indices by owner slab.
    std::vector<std::vector<int>> buckets(N);
    for (int i = 0; i < ng; ++i) {
        int owner = dd_->owner_of(global.atoms.x[i], global.atoms.y[i],
                                  global.atoms.z[i]);
        buckets[owner].push_back(i);
    }

    // Detect bonded topology once: only then do we carry per-atom bonded SoA to
    // each subdomain (so a pure-LJ system keeps the device SoA null — the
    // single-GPU/LJ path is byte-unchanged).
    bool has_bonded = false;
    for (int i = 0; i < ng; ++i)
        if (global.atoms.num_bond[i] || global.atoms.num_angle[i]
            || global.atoms.num_dihedral[i]) {
            has_bonded = true;
            break;
        }

    frames_.resize(N);
    n_local_.assign(N, 0);
    for (int s = 0; s < N; ++s) {
        const auto& idx = buckets[s];
        const int cnt = static_cast<int>(idx.size());
        n_local_[s] = cnt;

        // Build the subdomain CPU frame (owned atoms only). The box stays the
        // GLOBAL cell — positions are global coordinates and PBC is global, so
        // each subdomain's neighbor list min-images against the full box.
        cpu::Frame sub(cnt > 0 ? cnt : 1);
        sub.atoms.nlocal = cnt;
        sub.box = global.box;
        for (int k = 0; k < cnt; ++k) {
            const int g = idx[k];
            sub.atoms.x[k] = global.atoms.x[g];
            sub.atoms.y[k] = global.atoms.y[g];
            sub.atoms.z[k] = global.atoms.z[g];
            sub.atoms.vx[k] = global.atoms.vx[g];
            sub.atoms.vy[k] = global.atoms.vy[g];
            sub.atoms.vz[k] = global.atoms.vz[g];
            sub.atoms.mass[k] = global.atoms.mass[g];
            sub.atoms.charge[k] = global.atoms.charge[g];
            sub.atoms.type[k] = global.atoms.type[g];
            sub.atoms.tag[k] = global.atoms.tag[g];

            // Owner-filtered per-atom bonded topology, partners as GLOBAL TAGS
            // (resolved to local rows each reneighbor via the tag→local map).
            if (has_bonded) {
                sub.atoms.num_bond[k] = global.atoms.num_bond[g];
                for (int m = 0; m < cpu::Atom::kMaxBondPerAtom; ++m) {
                    sub.atoms.bond_atom(k, m) = global.atoms.bond_atom(g, m);
                    sub.atoms.bond_type(k, m) = global.atoms.bond_type(g, m);
                }
                sub.atoms.num_angle[k] = global.atoms.num_angle[g];
                for (int m = 0; m < cpu::Atom::kMaxAnglePerAtom; ++m) {
                    sub.atoms.angle_atom(k, 2 * m) = global.atoms.angle_atom(g, 2 * m);
                    sub.atoms.angle_atom(k, 2 * m + 1) = global.atoms.angle_atom(g, 2 * m + 1);
                    sub.atoms.angle_type(k, m) = global.atoms.angle_type(g, m);
                }
                sub.atoms.num_dihedral[k] = global.atoms.num_dihedral[g];
                for (int m = 0; m < cpu::Atom::kMaxDihedralPerAtom; ++m) {
                    sub.atoms.dihedral_atom(k, 3 * m) = global.atoms.dihedral_atom(g, 3 * m);
                    sub.atoms.dihedral_atom(k, 3 * m + 1) = global.atoms.dihedral_atom(g, 3 * m + 1);
                    sub.atoms.dihedral_atom(k, 3 * m + 2) = global.atoms.dihedral_atom(g, 3 * m + 2);
                    sub.atoms.dihedral_type(k, m) = global.atoms.dihedral_type(g, m);
                }
                for (int t = 0; t < 3; ++t)
                    sub.atoms.nspecial(k, t) = global.atoms.nspecial(g, t);
                for (int t = 0; t < cpu::Atom::kMaxSpecial; ++t)
                    sub.atoms.special(k, t) = global.atoms.special(g, t);
            }
        }

        // Materialize on this slot's GPU with ghost headroom so the halo can
        // append ghosts in [n_local, capacity) without reallocation. Floor of
        // +64 covers tiny subdomains where cnt*headroom rounds to ~0.
        const int extra = std::max(64, static_cast<int>(cnt * ghost_headroom));
        dm_.set_device(s);
        frames_[s] = cpu::to_cuda(sub, cnt + extra);
        // Carry owner-filtered per-atom bonded topology onto the device (global
        // tags); the single-GPU flat lists from to_cuda are unused on this path.
        if (has_bonded)
            cuda::materialize_per_atom_bonded(frames_[s], sub);
    }
}

namespace {
// Brute-force LJ 12-6 force on owned atoms [0,n_local) from owned+ghost
// [0,n_all), global-box min-image. Each owned atom's force is complete because
// the halo supplied all cross-boundary neighbors as ghosts.
__global__ void lj_bf_kernel(const real_t* __restrict__ x,
                             const real_t* __restrict__ y,
                             const real_t* __restrict__ z,
                             real_t* __restrict__ fx, real_t* __restrict__ fy,
                             real_t* __restrict__ fz, int n_local, int n_all,
                             cuda::Box box, work_t eps, work_t sig_sq,
                             work_t cut_sq, work_t min_sq) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n_local)
        return;
    work_t xi = x[i], yi = y[i], zi = z[i];
    work_t fxi = 0, fyi = 0, fzi = 0;
    for (int j = 0; j < n_all; ++j) {
        if (j == i)
            continue;
        work_t dx = work_t(x[j]) - xi, dy = work_t(y[j]) - yi,
               dz = work_t(z[j]) - zi;
        box.min_image(dx, dy, dz);
        work_t r2 = dx * dx + dy * dy + dz * dz;
        if (r2 >= cut_sq || r2 < min_sq)
            continue;
        work_t inv = work_t(1) / r2;
        work_t s2 = sig_sq * inv, s6 = s2 * s2 * s2, s12 = s6 * s6;
        work_t fr = work_t(24) * eps * (work_t(2) * s12 - s6) * inv;
        fxi += -fr * dx;
        fyi += -fr * dy;
        fzi += -fr * dz;
    }
    fx[i] = static_cast<real_t>(fxi);
    fy[i] = static_cast<real_t>(fyi);
    fz[i] = static_cast<real_t>(fzi);
}
} // namespace

void MultiGpu::compute_lj_bruteforce(double cutoff, double epsilon,
                                     double sigma) {
    cuda::Box box{};
    box.ortho = true;
    box.cell[0] = static_cast<real_t>(box_L_[0]);
    box.cell[4] = static_cast<real_t>(box_L_[1]);
    box.cell[8] = static_cast<real_t>(box_L_[2]);
    box.inv_l[0] = static_cast<real_t>(1.0 / box_L_[0]);
    box.inv_l[1] = static_cast<real_t>(1.0 / box_L_[1]);
    box.inv_l[2] = static_cast<real_t>(1.0 / box_L_[2]);
    const work_t cut_sq = work_t(cutoff) * work_t(cutoff);
    const work_t sig_sq = work_t(sigma) * work_t(sigma);
    const work_t min_sq = work_t(1e-6);
    for (int s = 0; s < count(); ++s) {
        const int nl = frames_[s].n_local_;
        if (nl == 0)
            continue;
        const int na = nl + frames_[s].n_ghost_;
        dm_.set_device(s);
        auto& a = frames_[s].atoms;
        int block = 256, grid = (nl + block - 1) / block;
        lj_bf_kernel<<<grid, block>>>(a.x, a.y, a.z, a.fx, a.fy, a.fz, nl, na,
                                      box, work_t(epsilon), sig_sq, cut_sq,
                                      min_sq);
        CKM(cudaGetLastError());
    }
    for (int s = 0; s < count(); ++s) {
        dm_.set_device(s);
        CKM(cudaDeviceSynchronize());
    }
}

namespace {
// Velocity-Verlet first half: v += 0.5*dt*(f/m); x += dt*v   (owned only).
__global__ void vv_first(real_t* __restrict__ x, real_t* __restrict__ y,
                         real_t* __restrict__ z, real_t* __restrict__ vx,
                         real_t* __restrict__ vy, real_t* __restrict__ vz,
                         const real_t* __restrict__ fx,
                         const real_t* __restrict__ fy,
                         const real_t* __restrict__ fz,
                         const real_t* __restrict__ mass, int n_local,
                         real_t dt) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n_local)
        return;
    real_t inv_m = real_t(1) / mass[i];
    real_t h = real_t(0.5) * dt * inv_m;
    vx[i] += h * fx[i];
    vy[i] += h * fy[i];
    vz[i] += h * fz[i];
    x[i] += dt * vx[i];
    y[i] += dt * vy[i];
    z[i] += dt * vz[i];
}
// Velocity-Verlet second half: v += 0.5*dt*(f/m)   (owned only).
__global__ void vv_second(real_t* __restrict__ vx, real_t* __restrict__ vy,
                          real_t* __restrict__ vz, const real_t* __restrict__ fx,
                          const real_t* __restrict__ fy,
                          const real_t* __restrict__ fz,
                          const real_t* __restrict__ mass, int n_local,
                          real_t dt) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n_local)
        return;
    real_t h = real_t(0.5) * dt / mass[i];
    vx[i] += h * fx[i];
    vy[i] += h * fy[i];
    vz[i] += h * fz[i];
}
} // namespace

void MultiGpu::run_nve(int steps, double dt, double cutoff, double epsilon,
                       double sigma) {
    exchange_halo(cutoff);
    compute_lj_bruteforce(cutoff, epsilon, sigma);
    const real_t dtr = static_cast<real_t>(dt);
    for (int step = 0; step < steps; ++step) {
        for (int s = 0; s < count(); ++s) {
            const int nl = frames_[s].n_local_;
            if (nl == 0)
                continue;
            dm_.set_device(s);
            auto& a = frames_[s].atoms;
            int block = 256, grid = (nl + block - 1) / block;
            vv_first<<<grid, block>>>(a.x, a.y, a.z, a.vx, a.vy, a.vz, a.fx,
                                      a.fy, a.fz, a.mass, nl, dtr);
            CKM(cudaGetLastError());
        }
        exchange_halo(cutoff);
        compute_lj_bruteforce(cutoff, epsilon, sigma);
        for (int s = 0; s < count(); ++s) {
            const int nl = frames_[s].n_local_;
            if (nl == 0)
                continue;
            dm_.set_device(s);
            auto& a = frames_[s].atoms;
            int block = 256, grid = (nl + block - 1) / block;
            vv_second<<<grid, block>>>(a.vx, a.vy, a.vz, a.fx, a.fy, a.fz,
                                       a.mass, nl, dtr);
            CKM(cudaGetLastError());
        }
    }
    for (int s = 0; s < count(); ++s) {
        dm_.set_device(s);
        CKM(cudaDeviceSynchronize());
    }
}

// ── Cluster neighbor-list force path ──────────────────────────────────────

namespace {
cuda::Box make_ortho_box(const double L[3]) {
    cuda::Box box{};
    box.ortho = true;
    for (int d = 0; d < 3; ++d)
        box.periodic[d] = true;
    box.cell[0] = static_cast<real_t>(L[0]);
    box.cell[4] = static_cast<real_t>(L[1]);
    box.cell[8] = static_cast<real_t>(L[2]);
    box.half_l[0] = static_cast<real_t>(0.5 * L[0]);
    box.half_l[1] = static_cast<real_t>(0.5 * L[1]);
    box.half_l[2] = static_cast<real_t>(0.5 * L[2]);
    box.inv_l[0] = static_cast<real_t>(1.0 / L[0]);
    box.inv_l[1] = static_cast<real_t>(1.0 / L[1]);
    box.inv_l[2] = static_cast<real_t>(1.0 / L[2]);
    return box;
}
} // namespace

void MultiGpu::destroy_force_engine() {
    for (int s = 0; s < static_cast<int>(nlist_.size()); ++s) {
        if (nlist_[s] && s < static_cast<int>(nlist_built_.size())
            && nlist_built_[s]) {
            dm_.set_device(s);
            cuda::DestroyContext dc{ctx_[s], frames_[s]};
            nlist_[s]->destroy(dc);
            if (pair_[s])
                pair_[s]->destroy(dc); // free PairLJCut's sorted-pos scratch
        }
    }
    nlist_.clear();
    pair_.clear();
    ctx_.clear();
    nlist_built_.clear();
    engine_ready_ = false;
}

void MultiGpu::setup_force_engine(double cutoff, double epsilon, double sigma,
                                  double skin) {
    destroy_force_engine();
    const int N = count();
    // Retain freed pool blocks so the per-reneighbor halo/migrate async
    // allocations are reused instead of round-tripping to the OS.
    set_pool_retain(dm_);
    force_cutoff_ = cutoff;
    force_skin_ = (skin < 0) ? 0.3 * cutoff : skin;
    const cuda::Box box = make_ortho_box(box_L_);
    ctx_.resize(N);
    nlist_.resize(N);
    pair_.resize(N);
    nlist_built_.assign(N, 0);
    for (int s = 0; s < N; ++s) {
        ctx_[s].stream = dm_.stream(s);
        ctx_[s].defer_stream = dm_.defer_stream(s);
        ctx_[s].host_box = box;
        ctx_[s].pair_list = nullptr;
        ctx_[s].flat_pairs_requested = false;
        nlist_[s] = std::make_unique<cuda::NeighborList>(cutoff, force_skin_);
        pair_[s] = std::make_unique<cuda::PairLJCut>(epsilon, sigma, cutoff);
        // This force engine drives only the NVE neighbor-list path
        // (run_nve_nlist), which never consumes the pair virial. Gate it off
        // (LAMMPS-style vflag): skips the 6 inner-loop virial products and
        // frees their registers under the kernel's launch-bounds budget.
        pair_[s]->set_compute_stress(false);
    }
    nlist_cap_.assign(N, 0); // built lazily in reneighbor() at ~1.3x n_all
    engine_ready_ = true;
}

void MultiGpu::reneighbor() {
    if (!engine_ready_)
        throw std::runtime_error(
            "MultiGpu::reneighbor: call setup_force_engine first");
    // Atoms that drifted out of their subdomain change owner GPU first, THEN
    // the halo is rebuilt around the updated owned sets.
    migrate();
    // Halo over cutoff+skin so ghosts cover the skin region — no cross-boundary
    // neighbor can appear before the next reneighbor.
    build_halo(force_cutoff_ + force_skin_);
    const int N = count();
    for (int s = 0; s < N; ++s) {
        const int na = frames_[s].n_local_ + frames_[s].n_ghost_;
        dm_.set_device(s);
        // (Re)build the neighbor list only when the owned+ghost count first
        // exceeds the buffers' size; otherwise just reset the active count (no
        // realloc — the dominant per-reneighbor cost). Build at ~1.3x n_all
        // (capped at the frame's allocated capacity) so ghost-count fluctuation
        // between intervals rarely forces a rebuild.
        if (!nlist_built_[s] || na > nlist_cap_[s]) {
            if (nlist_built_[s]) {
                cuda::DestroyContext dc{ctx_[s], frames_[s]};
                nlist_[s]->destroy(dc);
            }
            // Single GPU has no ghosts and a fixed owned count, so size the
            // neighbor-list buffers (incl. the per-atom Verlet list) at exactly
            // na — the 1.3x ghost-fluctuation headroom is pure VRAM waste there.
            // Multi-GPU keeps the headroom so ghost-count drift between
            // reneighbor intervals rarely forces a costly rebuild.
            int want = (count() == 1) ? na : (static_cast<int>(na * 1.3) + 64);
            if (want > frames_[s].capacity_)
                want = frames_[s].capacity_;
            frames_[s].n_atoms_ = want;
            cuda::BuildContext bc{ctx_[s], frames_[s]};
            nlist_[s]->build(bc);
            nlist_cap_[s] = want;
            nlist_built_[s] = 1;
        }
        frames_[s].n_atoms_ = na; // active owned+ghost; force writes owned only
        nlist_[s]->set_active_count(na);
    }
    for (int s = 0; s < N; ++s) {
        dm_.set_device(s);
        CKM(cudaDeviceSynchronize());
    }
}

void MultiGpu::compute_lj_nlist() {
    if (!engine_ready_)
        throw std::runtime_error(
            "MultiGpu::compute_lj_nlist: call setup_force_engine + reneighbor");
    const int N = count();
    for (int s = 0; s < N; ++s) {
        const int na = frames_[s].n_local_ + frames_[s].n_ghost_;
        if (na == 0)
            continue;
        frames_[s].n_atoms_ = na;
        dm_.set_device(s);
        auto& a = frames_[s].atoms;
        cudaStream_t st = static_cast<cudaStream_t>(ctx_[s].stream);
        // The cluster kernel accumulates (fx[ai] += ...) over owned AND ghost
        // slots, so zero the whole owned+ghost range first.
        CKM(cudaMemsetAsync(a.fx, 0, na * sizeof(real_t), st));
        CKM(cudaMemsetAsync(a.fy, 0, na * sizeof(real_t), st));
        CKM(cudaMemsetAsync(a.fz, 0, na * sizeof(real_t), st));
        cuda::LaunchContext lc{ctx_[s], frames_[s]};
        nlist_[s]->launch(lc); // skin-gated CSR rebuild + republish pair_list
        pair_[s]->launch(lc);  // atom-centric cluster LJ force (calls compute)
    }
    for (int s = 0; s < N; ++s) {
        dm_.set_device(s);
        CKM(cudaDeviceSynchronize());
    }
}

void MultiGpu::run_nve_nlist(int steps, double dt, double cutoff,
                             double epsilon, double sigma, double skin,
                             int reneigh_every) {
    setup_force_engine(cutoff, epsilon, sigma, skin);
    reneighbor();
    compute_lj_nlist();
    const real_t dtr = static_cast<real_t>(dt);
    if (reneigh_every < 1)
        reneigh_every = 1;
    for (int step = 0; step < steps; ++step) {
        for (int s = 0; s < count(); ++s) {
            const int nl = frames_[s].n_local_;
            if (nl == 0)
                continue;
            dm_.set_device(s);
            auto& a = frames_[s].atoms;
            int block = 256, grid = (nl + block - 1) / block;
            vv_first<<<grid, block>>>(a.x, a.y, a.z, a.vx, a.vy, a.vz, a.fx,
                                      a.fy, a.fz, a.mass, nl, dtr);
            CKM(cudaGetLastError());
        }
        if ((step + 1) % reneigh_every == 0)
            reneighbor();    // re-select ghosts + rebuild neighbor list
        else
            forward_comm();  // refresh ghost positions only
        compute_lj_nlist();
        for (int s = 0; s < count(); ++s) {
            const int nl = frames_[s].n_local_;
            if (nl == 0)
                continue;
            dm_.set_device(s);
            auto& a = frames_[s].atoms;
            int block = 256, grid = (nl + block - 1) / block;
            vv_second<<<grid, block>>>(a.vx, a.vy, a.vz, a.fx, a.fy, a.fz,
                                       a.mass, nl, dtr);
            CKM(cudaGetLastError());
        }
    }
    for (int s = 0; s < count(); ++s) {
        dm_.set_device(s);
        CKM(cudaDeviceSynchronize());
    }
}

// ── Migration: atoms crossing subdomain boundaries change owner GPU ────────

namespace {
// Owner slot of a position under the global 3D decomposition (device-side
// mirror of DomainDecomposition::owner_of: wrap into [0,L), bin per axis,
// x-fastest slot index).
__device__ inline int dev_owner(double x, double y, double z, int d0, int d1,
                                int d2, double Lx, double Ly, double Lz) {
    double c[3] = {x, y, z}, L[3] = {Lx, Ly, Lz};
    int dims[3] = {d0, d1, d2}, ci[3];
    for (int d = 0; d < 3; ++d) {
        double w = c[d] - floor(c[d] / L[d]) * L[d]; // wrap to [0,L)
        int k = static_cast<int>(w * dims[d] / L[d]);
        if (k >= dims[d])
            k = dims[d] - 1;
        if (k < 0)
            k = 0;
        ci[d] = k;
    }
    return (ci[2] * dims[1] + ci[1]) * dims[0] + ci[0];
}

__global__ void classify_owner(const real_t* x, const real_t* y,
                               const real_t* z, int n, int d0, int d1, int d2,
                               double Lx, double Ly, double Lz, int* owner,
                               int* cnt) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n)
        return;
    int o = dev_owner(x[i], y[i], z[i], d0, d1, d2, Lx, Ly, Lz);
    owner[i] = o;
    atomicAdd(&cnt[o], 1);
}

// Stable counting-sort scatter: group owned indices by destination owner.
__global__ void scatter_order(const int* owner, int n, const int* boff,
                              int* fill, int* order) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n)
        return;
    int o = owner[i];
    int pos = boff[o] + atomicAdd(&fill[o], 1);
    order[pos] = i;
}

// Gather full per-atom state of `cnt` atoms named by idx[] into contiguous
// send buffers (one P2P hop to the new owner).
__global__ void gather_state(
    const int* idx, int cnt, const real_t* x, const real_t* y, const real_t* z,
    const real_t* vx, const real_t* vy, const real_t* vz, const real_t* mass,
    const real_t* charge, const int* type, const std::int64_t* tag, real_t* ox,
    real_t* oy, real_t* oz, real_t* ovx, real_t* ovy, real_t* ovz,
    real_t* omass, real_t* ocharge, int* otype, std::int64_t* otag) {
    int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= cnt)
        return;
    int i = idx[k];
    ox[k] = x[i];
    oy[k] = y[i];
    oz[k] = z[i];
    ovx[k] = vx[i];
    ovy[k] = vy[i];
    ovz[k] = vz[i];
    omass[k] = mass[i];
    ocharge[k] = charge[i];
    otype[k] = type[i];
    otag[k] = tag[i];
}

// 10 per-atom fields moved by migration.
struct MigBuf {
    real_t *x, *y, *z, *vx, *vy, *vz, *mass, *charge;
    int* type;
    std::int64_t* tag;
};
void mig_alloc(MigBuf& b, int m) {
    CKM(mallocp(&b.x, m * sizeof(real_t)));
    CKM(mallocp(&b.y, m * sizeof(real_t)));
    CKM(mallocp(&b.z, m * sizeof(real_t)));
    CKM(mallocp(&b.vx, m * sizeof(real_t)));
    CKM(mallocp(&b.vy, m * sizeof(real_t)));
    CKM(mallocp(&b.vz, m * sizeof(real_t)));
    CKM(mallocp(&b.mass, m * sizeof(real_t)));
    CKM(mallocp(&b.charge, m * sizeof(real_t)));
    CKM(mallocp(&b.type, m * sizeof(int)));
    CKM(mallocp(&b.tag, m * sizeof(std::int64_t)));
}
void mig_free(MigBuf& b) {
    freep(b.x);
    freep(b.y);
    freep(b.z);
    freep(b.vx);
    freep(b.vy);
    freep(b.vz);
    freep(b.mass);
    freep(b.charge);
    freep(b.type);
    freep(b.tag);
}
} // namespace

void MultiGpu::migrate() {
    const int N = count();
    if (N < 2)
        return; // single GPU: every atom already owned here
    const auto dims = dd_->dims();

    // Phase 1a — per source: classify owners + per-destination COUNTS only
    // (cheap; the counting-sort `order` is deferred to phase 1b so the common
    // no-migration case skips it entirely). `owner` is retained for phase 1b.
    std::vector<int*> d_owner(N, nullptr);
    std::vector<std::vector<int>> cnt(N, std::vector<int>(N, 0)); // cnt[src][dst]
    std::vector<std::vector<int>> boff(N, std::vector<int>(N, 0));
    for (int s = 0; s < N; ++s) {
        const int nl = frames_[s].n_local_;
        dm_.set_device(s);
        int *owner, *dcnt;
        CKM(mallocp(&owner, std::max(1, nl) * sizeof(int)));
        CKM(mallocp(&dcnt, N * sizeof(int)));
        CKM(cudaMemset(dcnt, 0, N * sizeof(int)));
        if (nl > 0) {
            auto& a = frames_[s].atoms;
            int blk = 256, grd = (nl + blk - 1) / blk;
            classify_owner<<<grd, blk>>>(a.x, a.y, a.z, nl, dims[0], dims[1],
                                         dims[2], box_L_[0], box_L_[1],
                                         box_L_[2], owner, dcnt);
            CKM(cudaGetLastError());
        }
        std::vector<int> hc(N, 0);
        CKM(cudaMemcpy(hc.data(), dcnt, N * sizeof(int), cudaMemcpyDeviceToHost));
        cnt[s] = hc;
        std::vector<int> hb(N, 0);
        for (int d = 1; d < N; ++d)
            hb[d] = hb[d - 1] + hc[d - 1];
        boff[s] = hb;
        d_owner[s] = owner;
        freep(dcnt);
        CKM(cudaDeviceSynchronize()); // counts (and owner[]) ready
    }

    // Fast path: if every owned atom stayed on its current subdomain
    // (cnt[s][s] == n_local[s] for all s), there is nothing to move — skip the
    // counting sort, gather, per-(src,dst) P2P and publish entirely. In a stable
    // run atoms cross a boundary only rarely, so the common reneighbor does just
    // the cheap classify above. build_halo() rebuilds ghosts regardless.
    {
        bool any_move = false;
        for (int s = 0; s < N; ++s)
            if (cnt[s][s] != frames_[s].n_local_) {
                any_move = true;
                break;
            }
        if (!any_move) {
            for (int s = 0; s < N; ++s) {
                dm_.set_device(s);
                freep(d_owner[s]);
            }
            for (int s = 0; s < N; ++s) {
                dm_.set_device(s);
                CKM(cudaDeviceSynchronize());
            }
            return; // owned sets unchanged; ghosts rebuilt by build_halo()
        }
    }

    // Phase 1b — atoms DID move: build the counting-sort order per source.
    std::vector<int*> d_order(N, nullptr);
    for (int s = 0; s < N; ++s) {
        const int nl = frames_[s].n_local_;
        dm_.set_device(s);
        int *order, *dfill, *dboff;
        CKM(mallocp(&order, std::max(1, nl) * sizeof(int)));
        CKM(mallocp(&dfill, N * sizeof(int)));
        CKM(mallocp(&dboff, N * sizeof(int)));
        CKM(cudaMemcpy(dboff, boff[s].data(), N * sizeof(int),
                       cudaMemcpyHostToDevice));
        CKM(cudaMemset(dfill, 0, N * sizeof(int)));
        if (nl > 0) {
            int blk = 256, grd = (nl + blk - 1) / blk;
            scatter_order<<<grd, blk>>>(d_owner[s], nl, dboff, dfill, order);
            CKM(cudaGetLastError());
        }
        d_order[s] = order;
        freep(d_owner[s]);
        freep(dfill);
        freep(dboff);
        CKM(cudaDeviceSynchronize()); // order[] ready for cross-device gather
    }

    // Phase 2 — new owned counts + per-(src) receive offsets on each dst.
    std::vector<int> new_nl(N, 0);
    std::vector<std::vector<int>> recv_off(N, std::vector<int>(N, 0));
    for (int dst = 0; dst < N; ++dst) {
        int off = 0;
        for (int src = 0; src < N; ++src) {
            recv_off[src][dst] = off;
            off += cnt[src][dst];
        }
        new_nl[dst] = off;
        if (off > frames_[dst].capacity_)
            throw std::runtime_error("MultiGpu::migrate: owned overflow on slot "
                                     + std::to_string(dst)
                                     + " (raise ghost_headroom)");
    }

    // Phase 3 — staging owned arrays per dst (avoids read/write hazard on the
    // live frame arrays while peers write into them).
    std::vector<MigBuf> stg(N);
    for (int dst = 0; dst < N; ++dst) {
        dm_.set_device(dst);
        mig_alloc(stg[dst], std::max(1, new_nl[dst]));
    }

    // Phase 4 — per (src,dst): gather state on src, P2P copy to dst staging.
    for (int src = 0; src < N; ++src) {
        dm_.set_device(src);
        auto& a = frames_[src].atoms;
        const int sd = dm_.device_id(src);
        for (int dst = 0; dst < N; ++dst) {
            const int c = cnt[src][dst];
            if (c <= 0)
                continue;
            MigBuf sb;
            mig_alloc(sb, c);
            int blk = 256, grd = (c + blk - 1) / blk;
            const int* idx = d_order[src] + boff[src][dst];
            gather_state<<<grd, blk>>>(idx, c, a.x, a.y, a.z, a.vx, a.vy, a.vz,
                                       a.mass, a.charge, a.type, a.tag, sb.x,
                                       sb.y, sb.z, sb.vx, sb.vy, sb.vz, sb.mass,
                                       sb.charge, sb.type, sb.tag);
            CKM(cudaGetLastError());
            const int dd_id = dm_.device_id(dst);
            const int o = recv_off[src][dst];
            MigBuf& D = stg[dst];
            CKM(cudaMemcpyPeer(D.x + o, dd_id, sb.x, sd, c * sizeof(real_t)));
            CKM(cudaMemcpyPeer(D.y + o, dd_id, sb.y, sd, c * sizeof(real_t)));
            CKM(cudaMemcpyPeer(D.z + o, dd_id, sb.z, sd, c * sizeof(real_t)));
            CKM(cudaMemcpyPeer(D.vx + o, dd_id, sb.vx, sd, c * sizeof(real_t)));
            CKM(cudaMemcpyPeer(D.vy + o, dd_id, sb.vy, sd, c * sizeof(real_t)));
            CKM(cudaMemcpyPeer(D.vz + o, dd_id, sb.vz, sd, c * sizeof(real_t)));
            CKM(cudaMemcpyPeer(D.mass + o, dd_id, sb.mass, sd,
                               c * sizeof(real_t)));
            CKM(cudaMemcpyPeer(D.charge + o, dd_id, sb.charge, sd,
                               c * sizeof(real_t)));
            CKM(cudaMemcpyPeer(D.type + o, dd_id, sb.type, sd, c * sizeof(int)));
            CKM(cudaMemcpyPeer(D.tag + o, dd_id, sb.tag, sd,
                               c * sizeof(std::int64_t)));
            mig_free(sb); // cudaMemcpyPeer is synchronous; safe to free
        }
    }

    // Phase 5 — publish staging into the frame arrays, set new owned counts,
    // drop ghosts (rebuilt by the following build_halo).
    for (int dst = 0; dst < N; ++dst) {
        dm_.set_device(dst);
        auto& a = frames_[dst].atoms;
        const int m = new_nl[dst];
        if (m > 0) {
            const size_t r = m * sizeof(real_t);
            CKM(cudaMemcpy(a.x, stg[dst].x, r, cudaMemcpyDeviceToDevice));
            CKM(cudaMemcpy(a.y, stg[dst].y, r, cudaMemcpyDeviceToDevice));
            CKM(cudaMemcpy(a.z, stg[dst].z, r, cudaMemcpyDeviceToDevice));
            CKM(cudaMemcpy(a.vx, stg[dst].vx, r, cudaMemcpyDeviceToDevice));
            CKM(cudaMemcpy(a.vy, stg[dst].vy, r, cudaMemcpyDeviceToDevice));
            CKM(cudaMemcpy(a.vz, stg[dst].vz, r, cudaMemcpyDeviceToDevice));
            CKM(cudaMemcpy(a.mass, stg[dst].mass, r, cudaMemcpyDeviceToDevice));
            CKM(cudaMemcpy(a.charge, stg[dst].charge, r, cudaMemcpyDeviceToDevice));
            CKM(cudaMemcpy(a.type, stg[dst].type, m * sizeof(int),
                           cudaMemcpyDeviceToDevice));
            CKM(cudaMemcpy(a.tag, stg[dst].tag, m * sizeof(std::int64_t),
                           cudaMemcpyDeviceToDevice));
        }
        frames_[dst].n_local_ = m;
        n_local_[dst] = m;
        frames_[dst].n_ghost_ = 0;
        mig_free(stg[dst]);
        freep(d_order[dst]);
    }
    for (int s = 0; s < N; ++s) {
        dm_.set_device(s);
        CKM(cudaDeviceSynchronize());
    }
}

void MultiGpu::gather_forces(cpu::Frame& out) const {
    if (out.atoms.size() < n_global_)
        throw std::runtime_error("MultiGpu::gather_forces: out frame too small");
    for (int s = 0; s < count(); ++s) {
        const int cnt = n_local_[s];
        if (cnt == 0)
            continue;
        cpu::Frame tmp(cnt);
        dm_.set_device(s);
        cuda::to_cpu(frames_[s], tmp);
        // Resolve each owned atom's global slot by its persistent TAG (robust
        // after migration, since an atom's owner subdomain changes).
        const int nout = out.atoms.size();
        for (int k = 0; k < cnt; ++k) {
            const int g = static_cast<int>(tmp.atoms.tag[k]);
            if (g < 0 || g >= nout)
                continue;
            out.atoms.fx[g] = tmp.atoms.fx[k];
            out.atoms.fy[g] = tmp.atoms.fy[k];
            out.atoms.fz[g] = tmp.atoms.fz[k];
        }
    }
}

void MultiGpu::gather(cpu::Frame& out) const {
    if (out.atoms.size() < n_global_)
        throw std::runtime_error("MultiGpu::gather: out frame too small");
    for (int s = 0; s < count(); ++s) {
        const int cnt = n_local_[s];
        if (cnt == 0)
            continue;
        cpu::Frame tmp(cnt);
        dm_.set_device(s);
        cuda::to_cpu(frames_[s], tmp);
        // Scatter owned atoms back to their global slots by persistent TAG
        // (to_cpu round-trips tags D2H); tag resolution survives migration,
        // where an atom's owner subdomain changes.
        const int nout = out.atoms.size();
        for (int k = 0; k < cnt; ++k) {
            const int g = static_cast<int>(tmp.atoms.tag[k]);
            if (g < 0 || g >= nout)
                continue;
            out.atoms.x[g] = tmp.atoms.x[k];
            out.atoms.y[g] = tmp.atoms.y[k];
            out.atoms.z[g] = tmp.atoms.z[k];
            out.atoms.vx[g] = tmp.atoms.vx[k];
            out.atoms.vy[g] = tmp.atoms.vy[k];
            out.atoms.vz[g] = tmp.atoms.vz[k];
        }
    }
}

} // namespace md
} // namespace cuda
} // namespace ATV
