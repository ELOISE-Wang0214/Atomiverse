/// @file cuda/md/multi_gpu.cu
/// @brief MultiGpu orchestrator — partition / gather (Phase-4 foundation).

#include "atomiverse/cuda/md/multi_gpu.cuh"

#include "atomiverse/cpu/frame.h"
#include "atomiverse/cpu/transfer.h"
#include "atomiverse/cuda/box.cuh"
#include "atomiverse/cuda/md/nve_limit.cuh"
#include "atomiverse/cuda/transfer.cuh"
#include "atomiverse/cuda/types.cuh"
// Concrete force styles built by the default factories (the header now carries
// only the abstract bases — Pair/Bond/Angle/Dihedral/KSpace).
#include "atomiverse/cuda/potentials/angle_harmonic.cuh"
#include "atomiverse/cuda/potentials/bond_harmonic.cuh"
#include "atomiverse/cuda/potentials/dihedral_multi_harmonic.cuh"
#include "atomiverse/cuda/potentials/pair_lj_cut.cuh"
#include "neighborlist/atom_cell_sort.cuh"

#include <cub/cub.cuh>
#include <cuda_runtime.h>

#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <map>

namespace ATV {
namespace cuda {
namespace md {

// Env-gated (AV_MG_PROFILE) per-phase wall-clock profiler for run_nve_pe.
// Off by default (zero overhead — no device syncs added). When on, each timed
// phase syncs all devices before/after so the accumulated time is the true
// device-inclusive cost of that phase; totals print to stderr at program exit.
namespace {
struct MgPhaseProf {
    bool on;
    std::map<std::string, double> sec;
    std::map<std::string, long> calls;
    MgPhaseProf() : on(std::getenv("AV_MG_PROFILE") != nullptr) {}
    void add(const char* k, double s) {
        sec[k] += s;
        calls[k] += 1;
    }
    ~MgPhaseProf() {
        if (!on)
            return;
        double tot = 0;
        for (auto& kv : sec)
            tot += kv.second;
        std::fprintf(stderr, "\n=== MultiGpu run_nve_pe phase profile ===\n");
        for (auto& kv : sec)
            std::fprintf(stderr, "  %-22s %9.4f s  %5.1f%%  (%ld)\n",
                         kv.first.c_str(), kv.second,
                         tot > 0 ? 100.0 * kv.second / tot : 0.0,
                         calls[kv.first]);
        std::fprintf(stderr, "  %-22s %9.4f s\n", "TOTAL(timed)", tot);
    }
};
inline MgPhaseProf& mg_prof() {
    static MgPhaseProf p;
    return p;
}
} // namespace

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
//
// CAPPED (not ~0ull/unbounded): a per-reneighbor migrate/halo working set is a
// few hundred MB, so retaining ~1 GB keeps every common reallocation a pool hit
// (the malloc/free-avoidance win) while letting anything larger return to the
// OS instead of permanently inflating the footprint. Unbounded retention made
// nvidia-smi "used" grow ~2x at 4 GPUs (freed-but-retained reneighbor scratch),
// which capped the addressable atom count. Override via AV_MG_POOL_RETAIN_MB.
inline void set_pool_retain(const DeviceManager& dm) {
    std::uint64_t thr = 1024ull * 1024 * 1024; // 1 GB default
    if (const char* e = std::getenv("AV_MG_POOL_RETAIN_MB")) {
        const long mb = std::atol(e);
        if (mb >= 0)
            thr = static_cast<std::uint64_t>(mb) * 1024 * 1024;
    }
    for (int i = 0; i < dm.count(); ++i) {
        cudaMemPool_t pool;
        if (cudaDeviceGetDefaultMemPool(&pool, dm.device_id(i)) != cudaSuccess)
            continue;
        cudaMemPoolSetAttribute(pool, cudaMemPoolAttrReleaseThreshold, &thr);
    }
}
} // namespace

MultiGpu::MultiGpu(int n_gpus) : dm_(n_gpus) {}

MultiGpu::~MultiGpu() { destroy_frames(); }

void MultiGpu::destroy_frames() {
    destroy_force_engine(); // frees per-subdomain NeighborList buffers (needs frames_)
    destroy_halo_plan();
    destroy_tag_maps(); // bonded tag→local maps (sized to frame capacity)
    for (int s = 0; s < static_cast<int>(frames_.size()); ++s) {
        dm_.set_device(s);
        cuda::destroy(frames_[s]);
    }
    frames_.clear();
    n_local_.clear();
    sort_buf_.clear(); // free reused cell-sort scratch with the frames
    bonded_ = false;
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

// Packed forward halo (multigpu-bonded-05): gather the positions of `cnt` atoms
// named by idx into ONE contiguous [3*cnt] buffer (x|y|z blocks), so the
// per-step position refresh is a single cudaMemcpyPeer per transfer.
__global__ void gather_positions_packed(const real_t* __restrict__ x,
                                        const real_t* __restrict__ y,
                                        const real_t* __restrict__ z,
                                        const int* __restrict__ idx, int cnt,
                                        real_t* __restrict__ out3) {
    int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= cnt)
        return;
    int i = idx[k];
    out3[k] = x[i];
    out3[cnt + k] = y[i];
    out3[2 * cnt + k] = z[i];
}

// Packed forward halo: scatter a received [3*cnt] block into the dst ghost rows
// [base, base+cnt).
__global__ void unpack_positions_packed(real_t* __restrict__ x,
                                        real_t* __restrict__ y,
                                        real_t* __restrict__ z, int base,
                                        int cnt, const real_t* __restrict__ in3) {
    int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= cnt)
        return;
    x[base + k] = in3[k];
    y[base + k] = in3[cnt + k];
    z[base + k] = in3[2 * cnt + k];
}

// Packed reverse halo (multigpu-bonded-05): gather a destination's ghost-row
// forces into ONE contiguous [3*cnt] buffer (x-block | y-block | z-block) so the
// shipment back to the owner is a single cudaMemcpyPeer instead of three.
__global__ void gather_ghost_forces_packed(const real_t* __restrict__ fx,
                                           const real_t* __restrict__ fy,
                                           const real_t* __restrict__ fz,
                                           int base, int cnt,
                                           real_t* __restrict__ out3) {
    int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= cnt)
        return;
    out3[k] = fx[base + k];
    out3[cnt + k] = fy[base + k];
    out3[2 * cnt + k] = fz[base + k];
}

// Packed reverse halo: scatter-add the received [3*cnt] block into the source's
// rows named by idx (atomicAdd; same contract as scatter_add_forces).
__global__ void scatter_add_forces_packed(real_t* __restrict__ fx,
                                          real_t* __restrict__ fy,
                                          real_t* __restrict__ fz,
                                          const int* __restrict__ idx, int cnt,
                                          const real_t* __restrict__ in3) {
    int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= cnt)
        return;
    int i = idx[k];
    atomicAdd(&fx[i], in3[k]);
    atomicAdd(&fy[i], in3[cnt + k]);
    atomicAdd(&fz[i], in3[2 * cnt + k]);
}

// Fill idx[k] = k for k in [0, n) — the values-in array for the CUB tag sort.
__global__ void iota_k(int* __restrict__ idx, int n) {
    int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k < n)
        idx[k] = k;
}

// ── Cross-GPU bonded re-materialization (multigpu-bonded P2.4/P2.5) ──────────
// Resolve a global tag to a local row via binary search over the per-subdomain
// sorted (tag → local-index) map covering owned+ghost. -1 if absent.
__device__ inline int resolve_tag_bsearch(const std::int64_t* __restrict__ st,
                                          const int* __restrict__ si, int n,
                                          std::int64_t tag) {
    int lo = 0, hi = n - 1;
    while (lo <= hi) {
        int mid = (lo + hi) >> 1;
        std::int64_t v = st[mid];
        if (v == tag)
            return si[mid];
        if (v < tag)
            lo = mid + 1;
        else
            hi = mid - 1;
    }
    return -1;
}

// One thread per owned atom (the min-tag bond owner). Emit (i_local, j_local,
// type) into the flat BondList, resolving the partner global tag → local row
// (owned or ghost). Unresolvable partner (halo too thin) → drop + count.
__global__ __launch_bounds__(256, 4) void rematerialize_bonds_k(
    const int* __restrict__ num_bond, const std::int64_t* __restrict__ bond_atom,
    const int* __restrict__ bond_type, int MB, int n_local,
    const std::int64_t* __restrict__ st, const int* __restrict__ si, int nmap,
    int* __restrict__ oi, int* __restrict__ oj, int* __restrict__ ot,
    int* __restrict__ counter, int* __restrict__ dropped) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n_local)
        return;
    const int nb = num_bond[i];
    for (int m = 0; m < nb; ++m) {
        int j = resolve_tag_bsearch(st, si, nmap, bond_atom[i * MB + m]);
        if (j < 0) {
            atomicAdd(dropped, 1);
            continue;
        }
        int slot = atomicAdd(counter, 1);
        oi[slot] = i;
        oj[slot] = j;
        ot[slot] = bond_type[i * MB + m];
    }
}

// One thread per owned center atom (angle owner). angle_atom(j, 2m+0/1) hold the
// i- and k-partner tags; the center j is the local owned row.
__global__ __launch_bounds__(256, 4) void rematerialize_angles_k(
    const int* __restrict__ num_angle,
    const std::int64_t* __restrict__ angle_atom,
    const int* __restrict__ angle_type, int MA, int n_local,
    const std::int64_t* __restrict__ st, const int* __restrict__ si, int nmap,
    int* __restrict__ oi, int* __restrict__ oj, int* __restrict__ ok,
    int* __restrict__ ot, int* __restrict__ counter, int* __restrict__ dropped) {
    int j = blockIdx.x * blockDim.x + threadIdx.x;
    if (j >= n_local)
        return;
    const int na = num_angle[j];
    for (int m = 0; m < na; ++m) {
        int ai = resolve_tag_bsearch(st, si, nmap, angle_atom[j * 2 * MA + 2 * m]);
        int ak = resolve_tag_bsearch(st, si, nmap,
                                     angle_atom[j * 2 * MA + 2 * m + 1]);
        if (ai < 0 || ak < 0) {
            atomicAdd(dropped, 1);
            continue;
        }
        int slot = atomicAdd(counter, 1);
        oi[slot] = ai;
        oj[slot] = j;
        ok[slot] = ak;
        ot[slot] = angle_type[j * MA + m];
    }
}

// One thread per owned atom-2 (dihedral owner). dihedral_atom(j, 3m+0/1/2) hold
// the i-, k-, l-partner tags; the owner j is the local owned row (atom 2 of 4).
// NOTE: no __launch_bounds__ — this kernel's register footprint (3 binary
// searches + 5 stores) can exceed what (256,4) permits; an over-tight bound
// made the 256-thread launch fail with cudaErrorInvalidValue on sm_90/CUDA 13.
__global__ void rematerialize_dihedrals_k(
    const int* __restrict__ num_dihedral,
    const std::int64_t* __restrict__ dihedral_atom,
    const int* __restrict__ dihedral_type, int MD, int n_local,
    const std::int64_t* __restrict__ st, const int* __restrict__ si, int nmap,
    int* __restrict__ oi, int* __restrict__ oj, int* __restrict__ ok,
    int* __restrict__ ol, int* __restrict__ ot, int* __restrict__ counter,
    int* __restrict__ dropped) {
    int j = blockIdx.x * blockDim.x + threadIdx.x;
    if (j >= n_local)
        return;
    const int nd = num_dihedral[j];
    for (int m = 0; m < nd; ++m) {
        int di = resolve_tag_bsearch(st, si, nmap, dihedral_atom[j * 3 * MD + 3 * m]);
        int dk = resolve_tag_bsearch(st, si, nmap,
                                     dihedral_atom[j * 3 * MD + 3 * m + 1]);
        int dl = resolve_tag_bsearch(st, si, nmap,
                                     dihedral_atom[j * 3 * MD + 3 * m + 2]);
        if (di < 0 || dk < 0 || dl < 0) {
            atomicAdd(dropped, 1);
            continue;
        }
        int slot = atomicAdd(counter, 1);
        oi[slot] = di;
        oj[slot] = j;
        ok[slot] = dk;
        ol[slot] = dl;
        ot[slot] = dihedral_type[j * MD + m];
    }
}
} // namespace

void MultiGpu::destroy_halo_plan() {
    for (auto& p : halo_plan_) {
        dm_.set_device(p.src);
        freep(p.d_src_idx);
        freep(p.d_pack_src); // packed forward/reverse staging on src
        if (p.ev_src)
            cudaEventDestroy(p.ev_src);
        dm_.set_device(p.dst);
        freep(p.d_pack_dst); // packed forward/reverse staging on dst
        if (p.ev_dst)
            cudaEventDestroy(p.ev_dst);
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
                CKM(cudaMemcpy(p.d_src_idx, sidx, cnt * sizeof(int),
                               cudaMemcpyDeviceToDevice));
                // Packed [3*cnt] staging shared by forward (positions) and
                // reverse (forces): one buffer on the SRC device + one on the DST
                // device, so each per-step exchange is a single P2P (vs three).
                CKM(mallocp(&p.d_pack_src, 3 * cnt * sizeof(real_t)));
                // Forward copy-done event on the SRC device (recorded on src's
                // stream after the forward P2P; dst waits before unpacking).
                CKM(cudaEventCreateWithFlags(&p.ev_src,
                                             cudaEventDisableTiming));
                dm_.set_device(dst);
                CKM(mallocp(&p.d_pack_dst, 3 * cnt * sizeof(real_t)));
                // Reverse copy-done event on the DST device (recorded on dst's
                // stream after the reverse P2P; src waits before scatter-add).
                CKM(cudaEventCreateWithFlags(&p.ev_dst,
                                             cudaEventDisableTiming));
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
            const int sd = dm_.device_id(p.src), dd_id = dm_.device_id(p.dst);
            int block = 256, grid = (p.cnt + block - 1) / block;
            // Gather positions into the packed src buffer, ship in ONE async P2P,
            // unpack into the dst ghost rows. All on the per-device default
            // stream (0): gather→copy→record are ordered on the src stream;
            // dst's unpack waits on the src copy-done event. Async (vs the prior
            // host-blocking cudaMemcpyPeer) so the axis's independent plans —
            // distinct (src,dst) device pairs — overlap rather than serialize.
            dm_.set_device(p.src);
            auto& a = frames_[p.src].atoms;
            gather_positions_packed<<<grid, block>>>(a.x, a.y, a.z, p.d_src_idx,
                                                     p.cnt, p.d_pack_src);
            CKM(cudaGetLastError());
            CKM(cudaMemcpyPeerAsync(p.d_pack_dst, dd_id, p.d_pack_src, sd,
                                    3 * p.cnt * sizeof(real_t), 0));
            CKM(cudaEventRecord(p.ev_src, 0));
            dm_.set_device(p.dst);
            CKM(cudaStreamWaitEvent(0, p.ev_src, 0));
            auto& da = frames_[p.dst].atoms;
            unpack_positions_packed<<<grid, block>>>(da.x, da.y, da.z, p.base,
                                                     p.cnt, p.d_pack_dst);
            CKM(cudaGetLastError());
        }
        // Barrier between axes: this axis's unpacks into ghost rows must complete
        // before the next axis gathers from them (cascade). The async P2P above
        // already lets the axis's independent plans overlap; this only bounds the
        // axis, so the host blocks once per axis instead of once per plan.
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

            // Gather the dst's ghost-row forces into ONE packed [3*cnt] buffer,
            // ship dst→src in a single async P2P, scatter-add into the src rows.
            // gather→copy→record on the dst default stream; src's scatter-add
            // waits on the dst copy-done event. Async so the axis's independent
            // plans overlap instead of host-serializing (adjoint of forward_comm).
            dm_.set_device(p.dst);
            auto& da = frames_[p.dst].atoms;
            gather_ghost_forces_packed<<<grid, block>>>(da.fx, da.fy, da.fz,
                                                        p.base, p.cnt,
                                                        p.d_pack_dst);
            CKM(cudaGetLastError());
            CKM(cudaMemcpyPeerAsync(p.d_pack_src, sd, p.d_pack_dst, dd_id,
                                    3 * p.cnt * sizeof(real_t), 0));
            CKM(cudaEventRecord(p.ev_dst, 0));

            // Scatter-add into the src owned (or prior-axis ghost) rows.
            dm_.set_device(p.src);
            CKM(cudaStreamWaitEvent(0, p.ev_dst, 0));
            auto& a = frames_[p.src].atoms;
            scatter_add_forces_packed<<<grid, block>>>(
                a.fx, a.fy, a.fz, p.d_src_idx, p.cnt, p.d_pack_src);
            CKM(cudaGetLastError());
        }
        // Barrier between axes: this axis's scatter-adds into prior-axis ghost
        // rows must be visible before the next (lower) axis gathers from them.
        // The async P2P above overlaps the axis's independent plans; this bounds
        // the axis (one host block per axis instead of per plan).
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
        // The partition writes the wide bonded arrays directly (not via
        // add_bond), so explicitly allocate them on the fresh sub frame when the
        // system carries topology — otherwise the lazy atom_style-molecular
        // storage stays empty and the per-atom writes below run out of bounds.
        if (has_bonded)
            sub.atoms.ensure_bonded_storage();
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
    // Remember whether to drive the cross-GPU bonded re-materialization each
    // reneighbor (set after the per-subdomain SoA fill; destroy_frames resets it).
    bonded_ = has_bonded;
    bonded_static_built_ = false; // force a first build on the next reneighbor
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
    // Bonded components destroy unconditionally (built in setup_bonded, not gated
    // on the neighbor-list build cadence) — free each BondHarmonic's K/r0.
    for (int s = 0; s < static_cast<int>(bond_.size()); ++s) {
        if (s >= static_cast<int>(ctx_.size()))
            continue;
        dm_.set_device(s);
        cuda::DestroyContext dc{ctx_[s], frames_[s]};
        if (bond_[s])
            bond_[s]->destroy(dc);
        if (s < static_cast<int>(angle_.size()) && angle_[s])
            angle_[s]->destroy(dc);
        if (s < static_cast<int>(dihedral_.size()) && dihedral_[s])
            dihedral_[s]->destroy(dc);
    }
    // KSpace destroys independently of the bonded cadence (built in set_kspace,
    // which needs neither setup_bonded nor a neighbor-list build).
    for (int s = 0; s < static_cast<int>(kspace_.size()); ++s) {
        if (s >= static_cast<int>(ctx_.size()) || !kspace_[s])
            continue;
        dm_.set_device(s);
        cuda::DestroyContext dc{ctx_[s], frames_[s]};
        kspace_[s]->destroy(dc);
    }
    nlist_.clear();
    pair_.clear();
    kspace_.clear();
    bond_.clear();
    angle_.clear();
    dihedral_.clear();
    ctx_.clear();
    nlist_built_.clear();
    engine_ready_ = false;
}

void MultiGpu::setup_force_engine(double cutoff, double epsilon, double sigma,
                                  double skin, PairFactory pair_factory) {
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
    kspace_.resize(N); // optional reciprocal-space family (set_kspace fills it)
    nlist_built_.assign(N, 0);
    // Default pair factory: LJ with the per-step vflag gated off. This force
    // engine drives only the NVE neighbor-list path (run_nve_nlist), which never
    // consumes the pair virial; gating it skips the 6 inner-loop virial products
    // and frees their registers under the kernel's launch-bounds budget. A
    // caller-supplied pair_factory replaces this wholesale (ignores eps/sigma).
    auto default_pair = [epsilon, sigma, cutoff]() -> std::unique_ptr<cuda::Pair> {
        auto p = std::make_unique<cuda::PairLJCut>(epsilon, sigma, cutoff);
        p->set_compute_stress(false);
        return p;
    };
    for (int s = 0; s < N; ++s) {
        ctx_[s].stream = dm_.stream(s);
        ctx_[s].defer_stream = dm_.defer_stream(s);
        ctx_[s].host_box = box;
        ctx_[s].pair_list = nullptr;
        ctx_[s].flat_pairs_requested = false;
        nlist_[s] = std::make_unique<cuda::NeighborList>(cutoff, force_skin_);
        // This driver reneighbors on a fixed cadence (reneigh_every) and re-arms
        // the list via set_active_count, so disable the per-step skin auto-gate:
        // it is redundant here and false-fires on ghost rows, causing ~4x the
        // necessary neighbor rebuilds (the dominant multi-GPU cost).
        nlist_[s]->set_rebuild_on_request(true);
        pair_[s] = pair_factory ? pair_factory() : default_pair();
    }
    nlist_cap_.assign(N, 0); // built lazily in reneighbor() at ~1.3x n_all
    engine_ready_ = true;
}

void MultiGpu::set_kspace(KSpaceFactory kspace_factory) {
    const int N = count();
    if (!engine_ready_ || static_cast<int>(kspace_.size()) != N)
        throw std::runtime_error(
            "MultiGpu::set_kspace: call setup_force_engine first");
    if (!kspace_factory)
        return;
    for (int s = 0; s < N; ++s) {
        dm_.set_device(s);
        kspace_[s] = kspace_factory();
        cuda::BuildContext bc{ctx_[s], frames_[s]};
        kspace_[s]->build(bc);
    }
}

void MultiGpu::build_tag_maps() {
    const int N = count();
    if (static_cast<int>(tagmap_tag_.size()) != N) {
        tagmap_tag_.assign(N, nullptr);
        tagmap_idx_.assign(N, nullptr);
        drop_ctr_.assign(N, nullptr);
        tagmap_iota_.assign(N, nullptr);
        tagmap_tmp_.assign(N, nullptr);
        tagmap_tmp_sz_.assign(N, 0);
    }
    for (int s = 0; s < N; ++s) {
        const int na = frames_[s].n_local_ + frames_[s].n_ghost_;
        dm_.set_device(s);
        if (!tagmap_tag_[s]) {
            const int cap = frames_[s].capacity_;
            CKM(cudaMalloc(&tagmap_tag_[s], cap * sizeof(std::int64_t)));
            CKM(cudaMalloc(&tagmap_idx_[s], cap * sizeof(int)));
            CKM(cudaMalloc(&tagmap_iota_[s], cap * sizeof(int)));
            CKM(cudaMalloc(&drop_ctr_[s], sizeof(int))); // persistent, reused each rebuild
            // Size the CUB temp storage once for the full capacity (the query is
            // a no-op that only computes bytes); reuse the blob every reneighbor.
            std::size_t need = 0;
            cub::DeviceRadixSort::SortPairs(
                nullptr, need, tagmap_tag_[s], tagmap_tag_[s], tagmap_iota_[s],
                tagmap_idx_[s], cap);
            CKM(cudaMalloc(&tagmap_tmp_[s], need));
            tagmap_tmp_sz_[s] = need;
        }
        if (na <= 0)
            continue;
        // Build the sorted (tag → local-index) map fully on device: iota the
        // value array, then CUB radix-sort pairs (keys = owned+ghost tags,
        // values = local indices). No D2H/host-sort/H2D round-trip — the host
        // sort was the dominant multi-GPU per-reneighbor overhead. Runs on the
        // default stream (0), same as the rematerialize kernels that consume the
        // map in rebuild_bonded_tuples — so stream ordering carries the
        // dependency with no extra host sync.
        const cudaStream_t st = cudaStream_t(0);
        const int block = 256, grid = (na + block - 1) / block;
        iota_k<<<grid, block, 0, st>>>(tagmap_iota_[s], na);
        CKM(cudaGetLastError());
        std::size_t avail = tagmap_tmp_sz_[s];
        cub::DeviceRadixSort::SortPairs(
            tagmap_tmp_[s], avail, frames_[s].atoms.tag, tagmap_tag_[s],
            tagmap_iota_[s], tagmap_idx_[s], na, 0, sizeof(std::int64_t) * 8, st);
        CKM(cudaGetLastError());
    }
}

void MultiGpu::rebuild_bonded_tuples() {
    if (!bonded_)
        return;
    // Single-GPU: tags, owned set, and ghost set (none) never change, so the
    // tag→local map, bonded tuples, and special-bonds CSR are identical every
    // reneighbor. Build once; skip the per-reneighbor host-side tag sort +
    // rematerialize + CSR rebuild (the dominant single-GPU host cost).
    if (count() == 1 && bonded_static_built_)
        return;
    build_tag_maps();
    const int N = count();
    constexpr int MB = cpu::Atom::kMaxBondPerAtom;
    constexpr int MA = cpu::Atom::kMaxAnglePerAtom;
    constexpr int MD = cpu::Atom::kMaxDihedralPerAtom;
    for (int s = 0; s < N; ++s) {
        const int nl = frames_[s].n_local_;
        const int na = nl + frames_[s].n_ghost_;
        dm_.set_device(s);
        auto& a = frames_[s].atoms;
        auto& F = frames_[s];
        int* d_drop = drop_ctr_[s]; // persistent (allocated in build_tag_maps)
        CKM(cudaMemset(d_drop, 0, sizeof(int)));
        const int block = 256, grid = (nl + block - 1) / block;
        // Each tuple type's flat arrays + count cell exist only if that subdomain
        // owns ≥1 such tuple (to_cuda allocates them from the owned count). Guard
        // on the count pointer; the device count cell doubles as the append
        // counter (memset 0 → atomicAdd in the kernel → host-cached count).
        int nb = 0, nan = 0, nd = 0;
        if (nl > 0 && F.bonds.n_bonds) {
            CKM(cudaMemset(F.bonds.n_bonds, 0, sizeof(int)));
            rematerialize_bonds_k<<<grid, block>>>(
                a.num_bond, a.bond_atom, a.bond_type, MB, nl, tagmap_tag_[s],
                tagmap_idx_[s], na, F.bonds.atomi, F.bonds.atomj, F.bonds.type,
                F.bonds.n_bonds, d_drop);
            CKM(cudaGetLastError());
            CKM(cudaMemcpy(&nb, F.bonds.n_bonds, sizeof(int),
                           cudaMemcpyDeviceToHost));
        }
        if (nl > 0 && F.angles.n_angles) {
            CKM(cudaMemset(F.angles.n_angles, 0, sizeof(int)));
            rematerialize_angles_k<<<grid, block>>>(
                a.num_angle, a.angle_atom, a.angle_type, MA, nl, tagmap_tag_[s],
                tagmap_idx_[s], na, F.angles.atomi, F.angles.atomj,
                F.angles.atomk, F.angles.type, F.angles.n_angles, d_drop);
            CKM(cudaGetLastError());
            CKM(cudaMemcpy(&nan, F.angles.n_angles, sizeof(int),
                           cudaMemcpyDeviceToHost));
        }
        if (nl > 0 && F.dihedrals.n_dihedrals) {
            CKM(cudaMemset(F.dihedrals.n_dihedrals, 0, sizeof(int)));
            rematerialize_dihedrals_k<<<grid, block>>>(
                a.num_dihedral, a.dihedral_atom, a.dihedral_type, MD, nl,
                tagmap_tag_[s], tagmap_idx_[s], na, F.dihedrals.atomi,
                F.dihedrals.atomj, F.dihedrals.atomk, F.dihedrals.atoml,
                F.dihedrals.type, F.dihedrals.n_dihedrals, d_drop);
            CKM(cudaGetLastError());
            CKM(cudaMemcpy(&nd, F.dihedrals.n_dihedrals, sizeof(int),
                           cudaMemcpyDeviceToHost));
        }
        F.n_bonds_ = nb;
        F.n_angles_ = nan;
        F.n_dihedrals_ = nd;
        (void)d_drop; // persistent counter; not freed here
        // Rebuild the LJ special-bonds exclusion CSR from the fresh flat lists
        // (local/ghost indices), so 1-2/1-3/1-4 partners — including cross-GPU
        // ghosts — are excluded from the pair force this interval.
        if (special_bonds_ && s < static_cast<int>(pair_.size()) && pair_[s]) {
            F.n_atoms_ = na;
            cuda::BuildContext bc{ctx_[s], F};
            pair_[s]->build(bc);
        }
    }
    bonded_static_built_ = true; // single-GPU: skip subsequent rebuilds
}

void MultiGpu::destroy_tag_maps() {
    for (std::size_t s = 0; s < tagmap_tag_.size(); ++s) {
        dm_.set_device(static_cast<int>(s));
        cudaFree(tagmap_tag_[s]);
        cudaFree(tagmap_idx_[s]);
        if (s < tagmap_iota_.size())
            cudaFree(tagmap_iota_[s]);
        if (s < tagmap_tmp_.size())
            cudaFree(tagmap_tmp_[s]);
        if (s < drop_ctr_.size())
            cudaFree(drop_ctr_[s]);
    }
    tagmap_tag_.clear();
    tagmap_idx_.clear();
    tagmap_iota_.clear();
    tagmap_tmp_.clear();
    tagmap_tmp_sz_.clear();
    drop_ctr_.clear();
}

void MultiGpu::setup_bonded(std::vector<double> bond_K,
                            std::vector<double> bond_r0,
                            std::vector<double> angle_K,
                            std::vector<double> angle_theta0,
                            std::vector<std::array<double, 5>> dihedral_C,
                            BondFactory bond_factory, AngleFactory angle_factory,
                            DihedralFactory dihedral_factory) {
    if (!engine_ready_)
        throw std::runtime_error(
            "MultiGpu::setup_bonded: call setup_force_engine first");
    const int N = count();
    bond_.resize(N);
    angle_.resize(N);
    dihedral_.resize(N);
    const bool do_angle = !angle_K.empty();
    const bool do_dih = !dihedral_C.empty();
    // Default bonded factories: harmonic bond/angle + multi/harmonic dihedral
    // from the coefficient lists. A caller factory overrides the concrete style.
    auto default_bond = [&bond_K, &bond_r0]() -> std::unique_ptr<cuda::Bond> {
        return std::make_unique<cuda::BondHarmonic>(bond_K, bond_r0);
    };
    auto default_angle = [&angle_K, &angle_theta0]() -> std::unique_ptr<cuda::Angle> {
        return std::make_unique<cuda::AngleHarmonic>(angle_K, angle_theta0);
    };
    auto default_dih = [&dihedral_C]() -> std::unique_ptr<cuda::Dihedral> {
        return std::make_unique<cuda::DihedralMultiHarmonic>(dihedral_C);
    };
    for (int s = 0; s < N; ++s) {
        dm_.set_device(s);
        cuda::BuildContext bc{ctx_[s], frames_[s]};
        bond_[s] = bond_factory ? bond_factory() : default_bond();
        bond_[s]->build(bc);
        if (do_angle) {
            angle_[s] = angle_factory ? angle_factory() : default_angle();
            angle_[s]->build(bc);
        }
        if (do_dih) {
            dihedral_[s] = dihedral_factory ? dihedral_factory() : default_dih();
            dihedral_[s]->build(bc);
        }
    }
}

void MultiGpu::compute_bonded_forces() {
    const int N = count();
    if (static_cast<int>(bond_.size()) != N)
        throw std::runtime_error(
            "MultiGpu::compute_bonded_forces: call setup_bonded first");
    // Zero owned+ghost forces and the scalar accumulators (bond-only here; the
    // combined LJ+bonded path zeros once and runs both before reverse_comm).
    for (int s = 0; s < N; ++s) {
        const int na = frames_[s].n_local_ + frames_[s].n_ghost_;
        dm_.set_device(s);
        auto& a = frames_[s].atoms;
        if (na > 0) {
            CKM(cudaMemset(a.fx, 0, na * sizeof(real_t)));
            CKM(cudaMemset(a.fy, 0, na * sizeof(real_t)));
            CKM(cudaMemset(a.fz, 0, na * sizeof(real_t)));
        }
        CKM(cudaMemset(frames_[s].scalars, 0, sizeof(SystemScalars)));
    }
    // BondHarmonic over each subdomain's re-materialized tuples: deposits on
    // owned AND ghost partner rows (atomicAdd); energy counted once per bond on
    // its owner GPU.
    for (int s = 0; s < N; ++s) {
        const int na = frames_[s].n_local_ + frames_[s].n_ghost_;
        if (na == 0)
            continue;
        dm_.set_device(s);
        frames_[s].n_atoms_ = na;
        cuda::LaunchContext lc{ctx_[s], frames_[s]};
        bond_[s]->launch(lc);
        if (angle_[s])
            angle_[s]->launch(lc);
        if (dihedral_[s])
            dihedral_[s]->launch(lc);
    }
    for (int s = 0; s < N; ++s) {
        dm_.set_device(s);
        CKM(cudaDeviceSynchronize());
    }
    // Accumulate ghost-row bond forces onto their owners.
    reverse_comm();
    // Sum the per-subdomain bond energy (each bond owned by one GPU → no double
    // count) into the host-side total.
    bonded_energy_ = 0.0;
    for (int s = 0; s < N; ++s) {
        dm_.set_device(s);
        SystemScalars sc{};
        CKM(cudaMemcpy(&sc, frames_[s].scalars, sizeof(SystemScalars),
                       cudaMemcpyDeviceToHost));
        bonded_energy_ += static_cast<double>(sc.energy);
    }
}

void MultiGpu::enable_special_bonds() {
    const int N = count();
    if (static_cast<int>(pair_.size()) != N)
        throw std::runtime_error(
            "MultiGpu::enable_special_bonds: call setup_force_engine first");
    for (int s = 0; s < N; ++s)
        pair_[s]->enable_full_exclusions(); // full 1-2/1-3/1-4 exclusion (Pair hook)
    special_bonds_ = true; // rebuild_bonded_tuples will (re)build the CSR
}

void MultiGpu::compute_forces(bool read_energy) {
    const int N = count();
    if (!engine_ready_ || static_cast<int>(bond_.size()) != N)
        throw std::runtime_error(
            "MultiGpu::compute_forces: call setup_force_engine + setup_bonded");
    // Zero owned+ghost forces + scalars once; LJ writes owned (full), bonded
    // deposits owned+ghost, then one reverse_comm brings ghost contributions
    // home. (LJ leaves ghost forces zero, so only bonded ghost forces travel.)
    // ASYNC memset on each subdomain's stream — synchronous cudaMemset here
    // stalled the host 4×/step (the dominant per-step single-GPU overhead).
    for (int s = 0; s < N; ++s) {
        const int na = frames_[s].n_local_ + frames_[s].n_ghost_;
        dm_.set_device(s);
        auto& a = frames_[s].atoms;
        auto st = static_cast<cudaStream_t>(ctx_[s].stream);
        if (na > 0) {
            CKM(cudaMemsetAsync(a.fx, 0, na * sizeof(real_t), st));
            CKM(cudaMemsetAsync(a.fy, 0, na * sizeof(real_t), st));
            CKM(cudaMemsetAsync(a.fz, 0, na * sizeof(real_t), st));
        }
        CKM(cudaMemsetAsync(frames_[s].scalars, 0, sizeof(SystemScalars), st));
    }
    for (int s = 0; s < N; ++s) {
        const int na = frames_[s].n_local_ + frames_[s].n_ghost_;
        if (na == 0)
            continue;
        dm_.set_device(s);
        frames_[s].n_atoms_ = na;
        cuda::LaunchContext lc{ctx_[s], frames_[s]};
        // Multi-GPU: the neighbor list was rebuilt CONCURRENTLY in reneighbor()
        // (rebuild_issue/finalize), so first_launch is cleared and this launch
        // only republishes the pair list — no per-step rebuild here (which, in
        // this per-device loop, would serialize the GPUs on its overflow syncs).
        nlist_[s]->launch(lc); // neighbor refresh (driver-cadence rebuild)
        pair_[s]->launch(lc);  // short-range (owned, full, special-bonds-excluded)
        bond_[s]->launch(lc);
        if (angle_[s])
            angle_[s]->launch(lc);
        if (dihedral_[s])
            dihedral_[s]->launch(lc);
        // Reciprocal-space (LAMMPS order: pair → bonded → kspace). Optional;
        // single-GPU only until force-infra-04 adds the cross-GPU S(k) reduce.
        if (kspace_[s])
            kspace_[s]->launch(lc);
    }
    // Multi-GPU: sync all subdomains, then bring ghost (bonded) forces home.
    // Single-GPU has no ghosts (reverse_comm is a no-op) and runs the integrator
    // on the same stream as the force kernels, so no per-step host sync is
    // needed — stream ordering carries forces into the next vv_second. Dropping
    // this per-step cudaDeviceSynchronize is the last big single-GPU host win.
    if (N > 1) {
        auto tk0 = std::chrono::steady_clock::now();
        for (int s = 0; s < N; ++s) {
            dm_.set_device(s);
            CKM(cudaDeviceSynchronize());
        }
        if (mg_prof().on)
            mg_prof().add("C_force_kernels",
                          std::chrono::duration<double>(
                              std::chrono::steady_clock::now() - tk0)
                              .count());
        auto tr0 = std::chrono::steady_clock::now();
        reverse_comm(); // internally syncs per axis
        if (mg_prof().on)
            mg_prof().add("C_reverse_comm",
                          std::chrono::duration<double>(
                              std::chrono::steady_clock::now() - tr0)
                              .count());
    }
    if (read_energy) {
        // Energy readback syncs implicitly (blocking D2H); ensure single-GPU
        // force kernels are done first.
        if (N == 1)
            CKM(cudaDeviceSynchronize());
        bonded_energy_ = 0.0; // total PE (LJ + bonded), summed across subdomains
        for (int s = 0; s < N; ++s) {
            dm_.set_device(s);
            SystemScalars sc{};
            CKM(cudaMemcpy(&sc, frames_[s].scalars, sizeof(SystemScalars),
                           cudaMemcpyDeviceToHost));
            bonded_energy_ += static_cast<double>(sc.energy);
        }
    }
}

void MultiGpu::reneighbor() {
    if (!engine_ready_)
        throw std::runtime_error(
            "MultiGpu::reneighbor: call setup_force_engine first");
    auto timed = [&](const char* name, auto&& fn) {
        auto& P = mg_prof();
        if (!P.on) {
            fn();
            return;
        }
        for (int s = 0; s < count(); ++s) {
            dm_.set_device(s);
            cudaDeviceSynchronize();
        }
        auto t0 = std::chrono::steady_clock::now();
        fn();
        for (int s = 0; s < count(); ++s) {
            dm_.set_device(s);
            cudaDeviceSynchronize();
        }
        P.add(name, std::chrono::duration<double>(
                        std::chrono::steady_clock::now() - t0)
                        .count());
    };
    // Atoms that drifted out of their subdomain change owner GPU first, THEN
    // the halo is rebuilt around the updated owned sets.
    timed("R_migrate", [&] { migrate(); });
    // Cell-sort each subdomain's OWNED atoms by neighbor-grid cell index so a
    // warp's atoms share cells — coalesces the Verlet build's bin_pos scan and
    // the force gather (the multi-GPU even-count memory cliff). Runs BEFORE
    // build_halo so HaloPlan (built against owned), build_tag_maps, and
    // rebuild_bonded_tuples all see the sorted owned set. Grid dims mirror
    // nlist setup_cell_body: binsize = 0.5*(cutoff+skin), dims=max(3,L/binsize).
    // Beneficial-or-skip gate: single GPU has no halo/ghost gather and never
    // hits the multi-GPU build cliff, so the sort is pure overhead — skip it.
    // (Tiny multi-GPU subdomains below the cache cliff are handled by the
    // per-GPU size guard below.) Reused per-slot buffers => no per-reneighbor
    // cudaMalloc/Free in steady state.
    if (count() > 1) {
        timed("R_cell_sort", [&] {
            const double rs = force_cutoff_ + force_skin_;
            const double binsize = 0.5 * rs;
            if (static_cast<int>(sort_buf_.size()) != count()) {
                sort_buf_.clear();
                sort_buf_.resize(count());
                for (auto& p : sort_buf_)
                    p = std::make_unique<nlist::AtomCellSortBuffers>();
            }
            // Skip when the largest subdomain is below the cliff threshold:
            // the uncoalesced-build penalty only bites large per-GPU counts;
            // small systems (e.g. the parity tests) don't benefit and the sort
            // would just cost. ~2M owned atoms/GPU is well below the ~8M+ where
            // the cliff appears.
            int max_nl = 0;
            for (int s = 0; s < count(); ++s)
                max_nl = std::max(max_nl, frames_[s].n_local_);
            constexpr int kCellSortMinOwned = 2'000'000;
            // AV_CELL_SORT_FORCE bypasses the size gate so the parity suite can
            // exercise the (otherwise large-system-only) sort path on small
            // test frames.
            if (std::getenv("AV_CELL_SORT_OFF") != nullptr)
                return; // A/B kill-switch (regression measurement)
            if (max_nl < kCellSortMinOwned && std::getenv("AV_CELL_SORT_FORCE") == nullptr)
                return;
            for (int s = 0; s < count(); ++s) {
                if (frames_[s].n_local_ <= 0)
                    continue;
                dm_.set_device(s);
                const int dx = std::max(3, static_cast<int>(box_L_[0] / binsize));
                const int dy = std::max(3, static_cast<int>(box_L_[1] / binsize));
                const int dz = std::max(3, static_cast<int>(box_L_[2] / binsize));
                nlist::cell_sort_owned(
                    frames_[s], static_cast<real_t>(box_L_[0]),
                    static_cast<real_t>(box_L_[1]),
                    static_cast<real_t>(box_L_[2]), dx, dy, dz, /*periodic=*/true,
                    *sort_buf_[s], static_cast<cudaStream_t>(ctx_[s].stream));
            }
            for (int s = 0; s < count(); ++s) {
                dm_.set_device(s);
                CKM(cudaDeviceSynchronize());
            }
        });
    }
    // Halo over cutoff+skin so ghosts cover the skin region — no cross-boundary
    // neighbor can appear before the next reneighbor.
    timed("R_build_halo", [&] { build_halo(force_cutoff_ + force_skin_); });
    // Resolve cross-GPU bonded tuples against the fresh owned+ghost set.
    timed("R_rebuild_bonded", [&] { rebuild_bonded_tuples(); });
    auto t_nlist0 = std::chrono::steady_clock::now();
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
            // na. Multi-GPU keeps a small headroom so ghost-count drift between
            // reneighbor intervals rarely forces a costly rebuild — but the
            // ghost set is only ~2% of owned and fluctuates <1%, so the old
            // 1.3x (30% over-allocation of the single largest buffer, the
            // per-atom Verlet list) was pure waste that inflated the multi-GPU
            // footprint to ~3.9 GB/Matom and OOM'd large per-GPU systems. 1.08x
            // covers the observed fluctuation with ~25% less neighbor-list VRAM.
            int want =
                (count() == 1) ? na : (static_cast<int>(na * 1.08) + 64);
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
    // Multi-GPU: rebuild every subdomain's neighbor list CONCURRENTLY. The
    // per-device launch() rebuild has two internal cudaStreamSynchronize
    // overflow checks; letting compute_forces' per-device launch loop trigger
    // the rebuild (via set_active_count → first_launch) serialized the GPUs —
    // each device's ~75 ms build ran one-after-another (the dominant per-
    // interval multi-GPU cost, ~N× the single-GPU build). Instead: ISSUE the
    // rebuild on every GPU (no host sync), block on ONE cross-device barrier,
    // then FINALIZE (handle the rare overflow). The builds now overlap, so the
    // wall cost is ~one build regardless of N. compute_forces' following
    // launch() sees first_launch cleared and republishes only. (Single-GPU
    // keeps the launch()-driven rebuild: first_launch stays armed here.)
    if (N > 1) {
        for (int s = 0; s < N; ++s) {
            const int na = frames_[s].n_local_ + frames_[s].n_ghost_;
            if (na == 0)
                continue;
            dm_.set_device(s);
            frames_[s].n_atoms_ = na;
            cuda::LaunchContext lc{ctx_[s], frames_[s], 0, 0};
            nlist_[s]->rebuild_issue(lc);
        }
        for (int s = 0; s < N; ++s) {
            dm_.set_device(s);
            CKM(cudaDeviceSynchronize());
        }
        for (int s = 0; s < N; ++s) {
            const int na = frames_[s].n_local_ + frames_[s].n_ghost_;
            if (na == 0)
                continue;
            dm_.set_device(s);
            cuda::LaunchContext lc{ctx_[s], frames_[s], 0, 0};
            nlist_[s]->finalize_rebuild(lc);
        }
    }
    // Cap the retained-pool ratchet. migrate()/build_halo() allocate large
    // transient scratch (incl. packing per-atom bonded arrays) via
    // cudaMallocAsync; with the release threshold pinned at ~0ull
    // (set_pool_retain) the pool's reserved size climbs to the PEAK transient —
    // where old + new buffers briefly coexist — and never shrinks, so VRAM
    // ratcheted to ~2x the steady footprint and OOM'd large per-GPU systems
    // (e.g. 31->61 GB at 12M atoms/GPU, constant atom counts). The frees above
    // are complete after the sync, so release that unused excess back to the OS
    // now, keeping a small reserve for intra-interval reuse. Per-step (non-
    // reneighbor) allocations still hit the warm pool; only the once-per-
    // interval transient peak is reclaimed.
    constexpr std::size_t kPoolKeepBytes = std::size_t(256) << 20; // 256 MiB
    for (int s = 0; s < N; ++s) {
        cudaMemPool_t pool;
        if (cudaDeviceGetDefaultMemPool(&pool, dm_.device_id(s)) == cudaSuccess)
            cudaMemPoolTrimTo(pool, kPoolKeepBytes);
    }
    if (mg_prof().on) {
        mg_prof().add("R_nlist", std::chrono::duration<double>(
                                     std::chrono::steady_clock::now() - t_nlist0)
                                     .count());
        for (int s = 0; s < N; ++s) {
            dm_.set_device(s);
            size_t fr = 0, tt = 0;
            cudaMemGetInfo(&fr, &tt);
            std::fprintf(stderr,
                         "[MG_VRAM] reneigh dev%d used=%.1f/%.1f GB "
                         "n_local=%d n_ghost=%d nlist_cap=%d\n",
                         s, (tt - fr) / 1e9, tt / 1e9, frames_[s].n_local_,
                         frames_[s].n_ghost_, nlist_cap_[s]);
        }
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

void MultiGpu::run_nve_pe(int steps, double dt, int reneigh_every) {
    if (!engine_ready_ || bond_.empty())
        throw std::runtime_error(
            "MultiGpu::run_nve_pe: call setup_force_engine + setup_bonded first");
    reneighbor();      // migrate + halo + re-materialize topology + nlist
    compute_forces(/*read_energy=*/false); // initial LJ + bonded force
    const work_t dtw = static_cast<work_t>(dt);
    // The MSU PE stage-1 deck integrates the UNRELAXED amorphous start-state
    // (E_pair ~3.4e8) with LAMMPS `fix nve/limit 0.05` — a per-step
    // displacement cap that tames the severe initial atom overlap. Plain
    // velocity-Verlet diverges within a few steps here: atoms accelerate under
    // the overlap force spikes, the system clumps, and the neighbor count grows
    // without bound until it overflows int32 at ~1M atoms. Match LAMMPS:
    // displacement-capped integration (kFtm2v units + orthorhombic PBC wrap).
    constexpr work_t kNveLimitXmax = work_t(0.05); // Å, == `fix nve/limit 0.05`
    if (reneigh_every < 1)
        reneigh_every = 1;
    // Phase timer (no-op unless AV_MG_PROFILE is set): syncs all devices around
    // the phase so the recorded time is device-inclusive.
    auto timed = [&](const char* name, auto&& fn) {
        auto& P = mg_prof();
        if (!P.on) {
            fn();
            return;
        }
        for (int s = 0; s < count(); ++s) {
            dm_.set_device(s);
            cudaDeviceSynchronize();
        }
        auto t0 = std::chrono::steady_clock::now();
        fn();
        for (int s = 0; s < count(); ++s) {
            dm_.set_device(s);
            cudaDeviceSynchronize();
        }
        P.add(name, std::chrono::duration<double>(
                        std::chrono::steady_clock::now() - t0)
                        .count());
    };
    for (int step = 0; step < steps; ++step) {
      timed("1_integrate_first", [&] {
        for (int s = 0; s < count(); ++s) {
            const int nl = frames_[s].n_local_;
            if (nl == 0)
                continue;
            dm_.set_device(s);
            auto& a = frames_[s].atoms;
            // Stream choice: single-GPU runs the integrator on the SAME stream
            // as the force kernels (ctx_[s].stream) so compute_forces needs no
            // per-step host sync — stream ordering carries forces → vv_second.
            // Multi-GPU must use the per-device default stream (0): migrate() /
            // forward_comm() read these positions on the default stream, so the
            // integrator shares it (else the position write races the migrate
            // read — a Heisenbug that vanishes under compute-sanitizer).
            auto st = (count() == 1) ? static_cast<cudaStream_t>(ctx_[s].stream)
                                     : cudaStream_t(0);
            launch_nve_limit_initial(a.x, a.y, a.z, a.vx, a.vy, a.vz, a.fx, a.fy,
                                     a.fz, a.mass, ctx_[s].host_box, dtw,
                                     kNveLimitXmax, nl, st);
            CKM(cudaGetLastError());
        }
      });
        if ((step + 1) % reneigh_every == 0)
            timed("2_reneighbor", [&] { reneighbor(); });
        else
            timed("2_forward_comm", [&] { forward_comm(); });
        timed("3_compute_forces",
              [&] { compute_forces(/*read_energy=*/false); });
      timed("4_integrate_final", [&] {
        for (int s = 0; s < count(); ++s) {
            const int nl = frames_[s].n_local_;
            if (nl == 0)
                continue;
            dm_.set_device(s);
            auto& a = frames_[s].atoms;
            auto st = (count() == 1) ? static_cast<cudaStream_t>(ctx_[s].stream)
                                     : cudaStream_t(0);
            launch_nve_limit_final(a.vx, a.vy, a.vz, a.fx, a.fy, a.fz, a.mass,
                                   dtw, kNveLimitXmax, nl, st);
            CKM(cudaGetLastError());
        }
      });
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

// Per-block shared-memory partial counts: every atom atomicAdds into a block-
// local s_cnt[owner], then each block contributes ONE global atomicAdd per
// destination slot. Cuts global-atomic traffic from O(n) (all n atoms hammering
// the few cnt[] counters — the dominant migrate cost at ~16M atoms/GPU) to
// O(blocks * n_slots). n_slots is tiny (= #subdomains), so the shared array fits.
__global__ void classify_owner(const real_t* x, const real_t* y,
                               const real_t* z, int n, int d0, int d1, int d2,
                               double Lx, double Ly, double Lz, int* owner,
                               int* cnt, int n_slots) {
    extern __shared__ int s_cnt[];
    for (int t = threadIdx.x; t < n_slots; t += blockDim.x)
        s_cnt[t] = 0;
    __syncthreads();
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) {
        int o = dev_owner(x[i], y[i], z[i], d0, d1, d2, Lx, Ly, Lz);
        owner[i] = o;
        atomicAdd(&s_cnt[o], 1);
    }
    __syncthreads();
    for (int t = threadIdx.x; t < n_slots; t += blockDim.x)
        if (s_cnt[t] != 0)
            atomicAdd(&cnt[t], s_cnt[t]);
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

// Per-atom bonded topology carried alongside the core MigBuf when the system
// has bonds (multigpu-bonded P2.3). Row widths mirror cpu::Atom; partners are
// global tags (resolved to local rows on the new owner each reneighbor).
struct MigBondedBuf {
    int *num_bond, *bond_type, *num_angle, *angle_type, *num_dihedral,
        *dihedral_type, *nspecial;
    std::int64_t *bond_atom, *angle_atom, *dihedral_atom, *special;
};
void migb_alloc(MigBondedBuf& b, int m) {
    const int MB = cpu::Atom::kMaxBondPerAtom, MA = cpu::Atom::kMaxAnglePerAtom,
              MD = cpu::Atom::kMaxDihedralPerAtom, MS = cpu::Atom::kMaxSpecial;
    CKM(mallocp(&b.num_bond, m * sizeof(int)));
    CKM(mallocp(&b.bond_atom, m * MB * sizeof(std::int64_t)));
    CKM(mallocp(&b.bond_type, m * MB * sizeof(int)));
    CKM(mallocp(&b.num_angle, m * sizeof(int)));
    CKM(mallocp(&b.angle_atom, m * 2 * MA * sizeof(std::int64_t)));
    CKM(mallocp(&b.angle_type, m * MA * sizeof(int)));
    CKM(mallocp(&b.num_dihedral, m * sizeof(int)));
    CKM(mallocp(&b.dihedral_atom, m * 3 * MD * sizeof(std::int64_t)));
    CKM(mallocp(&b.dihedral_type, m * MD * sizeof(int)));
    CKM(mallocp(&b.nspecial, m * 3 * sizeof(int)));
    CKM(mallocp(&b.special, m * MS * sizeof(std::int64_t)));
}
void migb_free(MigBondedBuf& b) {
    freep(b.num_bond);
    freep(b.bond_atom);
    freep(b.bond_type);
    freep(b.num_angle);
    freep(b.angle_atom);
    freep(b.angle_type);
    freep(b.num_dihedral);
    freep(b.dihedral_atom);
    freep(b.dihedral_type);
    freep(b.nspecial);
    freep(b.special);
}

// Gather the per-atom bonded SoA rows of `cnt` atoms named by idx[] into
// contiguous send buffers (one row per atom, fixed widths). Mirrors gather_state.
__global__ void gather_bonded_state(
    const int* __restrict__ idx, int cnt, int MB, int MA, int MD, int MS,
    const int* num_bond, const std::int64_t* bond_atom, const int* bond_type,
    const int* num_angle, const std::int64_t* angle_atom, const int* angle_type,
    const int* num_dihedral, const std::int64_t* dihedral_atom,
    const int* dihedral_type, const int* nspecial, const std::int64_t* special,
    int* onb, std::int64_t* oba, int* obt, int* ona, std::int64_t* oaa,
    int* oat, int* ond, std::int64_t* oda, int* odt, int* ons,
    std::int64_t* osp) {
    int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= cnt)
        return;
    int i = idx[k];
    onb[k] = num_bond[i];
    for (int w = 0; w < MB; ++w) {
        oba[k * MB + w] = bond_atom[i * MB + w];
        obt[k * MB + w] = bond_type[i * MB + w];
    }
    ona[k] = num_angle[i];
    for (int w = 0; w < 2 * MA; ++w)
        oaa[k * 2 * MA + w] = angle_atom[i * 2 * MA + w];
    for (int w = 0; w < MA; ++w)
        oat[k * MA + w] = angle_type[i * MA + w];
    ond[k] = num_dihedral[i];
    for (int w = 0; w < 3 * MD; ++w)
        oda[k * 3 * MD + w] = dihedral_atom[i * 3 * MD + w];
    for (int w = 0; w < MD; ++w)
        odt[k * MD + w] = dihedral_type[i * MD + w];
    for (int w = 0; w < 3; ++w)
        ons[k * 3 + w] = nspecial[i * 3 + w];
    for (int w = 0; w < MS; ++w)
        osp[k * MS + w] = special[i * MS + w];
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
    // ISSUE the classify on every GPU first (no host sync), then COLLECT the
    // per-destination counts. The previous per-device "classify → blocking D2H
    // → sync" loop serialized the GPUs (each D2H stalled the host until that
    // device's classify finished); on the common no-migration reneighbor this
    // was pure serial overhead × N. Splitting the issue from the collect lets
    // the N classifies overlap; the collect-loop D2H then waits on
    // already-finished kernels.
    std::vector<int*> d_dcnt(N, nullptr);
    for (int s = 0; s < N; ++s) {
        const int nl = frames_[s].n_local_;
        dm_.set_device(s);
        int *owner, *dcnt;
        CKM(mallocp(&owner, std::max(1, nl) * sizeof(int)));
        CKM(mallocp(&dcnt, N * sizeof(int)));
        CKM(cudaMemsetAsync(dcnt, 0, N * sizeof(int),
                            static_cast<cudaStream_t>(0)));
        if (nl > 0) {
            auto& a = frames_[s].atoms;
            int blk = 256, grd = (nl + blk - 1) / blk;
            const int nslots = dims[0] * dims[1] * dims[2]; // == N
            classify_owner<<<grd, blk, nslots * sizeof(int)>>>(
                a.x, a.y, a.z, nl, dims[0], dims[1], dims[2], box_L_[0],
                box_L_[1], box_L_[2], owner, dcnt, nslots);
            CKM(cudaGetLastError());
        }
        d_owner[s] = owner;
        d_dcnt[s] = dcnt;
    }
    for (int s = 0; s < N; ++s) {
        dm_.set_device(s);
        std::vector<int> hc(N, 0);
        CKM(cudaMemcpy(hc.data(), d_dcnt[s], N * sizeof(int),
                       cudaMemcpyDeviceToHost)); // blocks on s's classify only
        cnt[s] = hc;
        std::vector<int> hb(N, 0);
        for (int d = 1; d < N; ++d)
            hb[d] = hb[d - 1] + hc[d - 1];
        boff[s] = hb;
        freep(d_dcnt[s]);
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
    std::vector<MigBondedBuf> stgb(N); // per-atom bonded SoA staging (if bonded_)
    for (int dst = 0; dst < N; ++dst) {
        dm_.set_device(dst);
        mig_alloc(stg[dst], std::max(1, new_nl[dst]));
        if (bonded_)
            migb_alloc(stgb[dst], std::max(1, new_nl[dst]));
    }
    // These staging buffers come from cudaMallocAsync on each dst's stream 0
    // (stream-ordered). Phase 4 writes into them from a DIFFERENT device via
    // cudaMemcpyPeer, which carries no ordering against dst's pending async
    // allocation — using the memory before the alloc completes on dst's stream
    // is undefined (manifests as a Heisenbug "invalid argument" that vanishes
    // under compute-sanitizer's serialization). Drain every dst's stream so the
    // allocations are real before any cross-device peer write touches them.
    for (int dst = 0; dst < N; ++dst) {
        dm_.set_device(dst);
        CKM(cudaDeviceSynchronize());
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
            if (bonded_) {
                // Carry the per-atom bonded SoA rows for these atoms too.
                const int MB = cpu::Atom::kMaxBondPerAtom,
                          MA = cpu::Atom::kMaxAnglePerAtom,
                          MD = cpu::Atom::kMaxDihedralPerAtom,
                          MS = cpu::Atom::kMaxSpecial;
                MigBondedBuf sbb;
                migb_alloc(sbb, c);
                gather_bonded_state<<<grd, blk>>>(
                    idx, c, MB, MA, MD, MS, a.num_bond, a.bond_atom, a.bond_type,
                    a.num_angle, a.angle_atom, a.angle_type, a.num_dihedral,
                    a.dihedral_atom, a.dihedral_type, a.nspecial, a.special,
                    sbb.num_bond, sbb.bond_atom, sbb.bond_type, sbb.num_angle,
                    sbb.angle_atom, sbb.angle_type, sbb.num_dihedral,
                    sbb.dihedral_atom, sbb.dihedral_type, sbb.nspecial,
                    sbb.special);
                CKM(cudaGetLastError());
                MigBondedBuf& DB = stgb[dst];
                auto peer_i = [&](int* dstp, int* srcp, int width) {
                    CKM(cudaMemcpyPeer(dstp + static_cast<std::size_t>(o) * width,
                                       dd_id, srcp, sd,
                                       static_cast<std::size_t>(c) * width
                                           * sizeof(int)));
                };
                auto peer_t = [&](std::int64_t* dstp, std::int64_t* srcp,
                                  int width) {
                    CKM(cudaMemcpyPeer(dstp + static_cast<std::size_t>(o) * width,
                                       dd_id, srcp, sd,
                                       static_cast<std::size_t>(c) * width
                                           * sizeof(std::int64_t)));
                };
                peer_i(DB.num_bond, sbb.num_bond, 1);
                peer_t(DB.bond_atom, sbb.bond_atom, MB);
                peer_i(DB.bond_type, sbb.bond_type, MB);
                peer_i(DB.num_angle, sbb.num_angle, 1);
                peer_t(DB.angle_atom, sbb.angle_atom, 2 * MA);
                peer_i(DB.angle_type, sbb.angle_type, MA);
                peer_i(DB.num_dihedral, sbb.num_dihedral, 1);
                peer_t(DB.dihedral_atom, sbb.dihedral_atom, 3 * MD);
                peer_i(DB.dihedral_type, sbb.dihedral_type, MD);
                peer_i(DB.nspecial, sbb.nspecial, 3);
                peer_t(DB.special, sbb.special, MS);
                migb_free(sbb);
            }
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
            if (bonded_) {
                const int MB = cpu::Atom::kMaxBondPerAtom,
                          MA = cpu::Atom::kMaxAnglePerAtom,
                          MD = cpu::Atom::kMaxDihedralPerAtom,
                          MS = cpu::Atom::kMaxSpecial;
                MigBondedBuf& B = stgb[dst];
                auto d2d_i = [&](int* dstp, int* srcp, int width) {
                    CKM(cudaMemcpy(dstp, srcp,
                                   static_cast<std::size_t>(m) * width
                                       * sizeof(int),
                                   cudaMemcpyDeviceToDevice));
                };
                auto d2d_t = [&](std::int64_t* dstp, std::int64_t* srcp,
                                 int width) {
                    CKM(cudaMemcpy(dstp, srcp,
                                   static_cast<std::size_t>(m) * width
                                       * sizeof(std::int64_t),
                                   cudaMemcpyDeviceToDevice));
                };
                d2d_i(a.num_bond, B.num_bond, 1);
                d2d_t(a.bond_atom, B.bond_atom, MB);
                d2d_i(a.bond_type, B.bond_type, MB);
                d2d_i(a.num_angle, B.num_angle, 1);
                d2d_t(a.angle_atom, B.angle_atom, 2 * MA);
                d2d_i(a.angle_type, B.angle_type, MA);
                d2d_i(a.num_dihedral, B.num_dihedral, 1);
                d2d_t(a.dihedral_atom, B.dihedral_atom, 3 * MD);
                d2d_i(a.dihedral_type, B.dihedral_type, MD);
                d2d_i(a.nspecial, B.nspecial, 3);
                d2d_t(a.special, B.special, MS);
            }
        }
        frames_[dst].n_local_ = m;
        n_local_[dst] = m;
        frames_[dst].n_ghost_ = 0;
        mig_free(stg[dst]);
        if (bonded_)
            migb_free(stgb[dst]);
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
