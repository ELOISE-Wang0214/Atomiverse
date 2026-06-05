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
} // namespace

void MultiGpu::destroy_halo_plan() {
    for (auto& p : halo_plan_) {
        dm_.set_device(p.src);
        cudaFree(p.d_src_idx);
        cudaFree(p.d_sx);
        cudaFree(p.d_sy);
        cudaFree(p.d_sz);
    }
    halo_plan_.clear();
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
            cudaMalloc(&sx, ns * sizeof(real_t));
            cudaMalloc(&sy, ns * sizeof(real_t));
            cudaMalloc(&sz, ns * sizeof(real_t));
            cudaMalloc(&stag, ns * sizeof(std::int64_t));
            cudaMalloc(&stype, ns * sizeof(int));
            cudaMalloc(&sidx, ns * sizeof(int));
            cudaMalloc(&scnt, sizeof(int));
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
                cudaMalloc(&p.d_src_idx, cnt * sizeof(int));
                cudaMalloc(&p.d_sx, cnt * sizeof(real_t));
                cudaMalloc(&p.d_sy, cnt * sizeof(real_t));
                cudaMalloc(&p.d_sz, cnt * sizeof(real_t));
                CKM(cudaMemcpy(p.d_src_idx, sidx, cnt * sizeof(int),
                               cudaMemcpyDeviceToDevice));
                CKM(cudaMemcpy(p.d_sx, sx, cnt * sizeof(real_t),
                               cudaMemcpyDeviceToDevice));
                CKM(cudaMemcpy(p.d_sy, sy, cnt * sizeof(real_t),
                               cudaMemcpyDeviceToDevice));
                CKM(cudaMemcpy(p.d_sz, sz, cnt * sizeof(real_t),
                               cudaMemcpyDeviceToDevice));
                halo_plan_.push_back(p);

                goff[dst] += cnt;
                frames_[dst].n_ghost_ = goff[dst];
            }
            cudaSetDevice(dm_.device_id(src));
            cudaFree(sx);
            cudaFree(sy);
            cudaFree(sz);
            cudaFree(stag);
            cudaFree(stype);
            cudaFree(sidx);
            cudaFree(scnt);
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

    frames_.resize(N);
    n_local_.assign(N, 0);
    owned_idx_ = buckets; // remember each owned atom's global index (gather map)
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
        }

        // Materialize on this slot's GPU with ghost headroom so the halo can
        // append ghosts in [n_local, capacity) without reallocation. Floor of
        // +64 covers tiny subdomains where cnt*headroom rounds to ~0.
        const int extra = std::max(64, static_cast<int>(cnt * ghost_headroom));
        dm_.set_device(s);
        frames_[s] = cpu::to_cuda(sub, cnt + extra);
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
    }
    engine_ready_ = true;
}

void MultiGpu::reneighbor() {
    if (!engine_ready_)
        throw std::runtime_error(
            "MultiGpu::reneighbor: call setup_force_engine first");
    // Halo over cutoff+skin so ghosts cover the skin region — no cross-boundary
    // neighbor can appear before the next reneighbor.
    build_halo(force_cutoff_ + force_skin_);
    const int N = count();
    for (int s = 0; s < N; ++s) {
        const int na = frames_[s].n_local_ + frames_[s].n_ghost_;
        frames_[s].n_atoms_ = na; // cluster owned+ghost; force writes owned only
        if (na == 0)
            continue;
        dm_.set_device(s);
        // n_all changed vs the prior interval, so the neighbor list (sized at
        // build) must be torn down and reallocated. Reneighbor is infrequent.
        if (nlist_built_[s]) {
            cuda::DestroyContext dc{ctx_[s], frames_[s]};
            nlist_[s]->destroy(dc);
        }
        cuda::BuildContext bc{ctx_[s], frames_[s]};
        nlist_[s]->build(bc); // publishes pair_list to ctx_[s]
        nlist_built_[s] = 1;
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
        const auto& gidx = owned_idx_[s];
        for (int k = 0; k < cnt; ++k) {
            const int g = gidx[k];
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
        // Scatter owned atoms back to their global slots. Use the global
        // indices captured at partition (atom k of subdomain s came from
        // global atom owned_idx_[s][k]) rather than the device tags, which is
        // robust and independent of tag round-trip.
        const auto& gidx = owned_idx_[s];
        for (int k = 0; k < cnt; ++k) {
            const int g = gidx[k];
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
