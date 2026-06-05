/// @file cuda/md/multi_gpu.cu
/// @brief MultiGpu orchestrator — partition / gather (Phase-4 foundation).

#include "atomiverse/cuda/md/multi_gpu.cuh"

#include "atomiverse/cpu/frame.h"
#include "atomiverse/cpu/transfer.h"
#include "atomiverse/cuda/transfer.cuh"

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
// Select owned atoms within `cutoff` of a face plane and compact their
// x/y/z/tag/type into the send buffers (atomic compaction). `dir > 0`: atoms
// near the high (subhi) boundary; `dir < 0`: near the low (sublo) boundary.
__global__ void pack_face(const real_t* __restrict__ x,
                          const real_t* __restrict__ y,
                          const real_t* __restrict__ z,
                          const std::int64_t* __restrict__ tag,
                          const int* __restrict__ type, int n_local, int axis,
                          int dir, double plane, double cutoff,
                          real_t* __restrict__ ox, real_t* __restrict__ oy,
                          real_t* __restrict__ oz,
                          std::int64_t* __restrict__ otag,
                          int* __restrict__ otype, int* __restrict__ counter,
                          int cap) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n_local)
        return;
    real_t q = (axis == 0) ? x[i] : (axis == 1) ? y[i] : z[i];
    bool sel = (dir > 0) ? (q > plane - cutoff) : (q < plane + cutoff);
    if (!sel)
        return;
    int s = atomicAdd(counter, 1);
    if (s < cap) {
        ox[s] = x[i];
        oy[s] = y[i];
        oz[s] = z[i];
        otag[s] = tag[i];
        otype[s] = type[i];
    }
}
} // namespace

void MultiGpu::exchange_halo(double cutoff) {
    const int N = count();
    // Reset ghost state.
    std::vector<int> goff(N, 0);
    for (int s = 0; s < N; ++s) {
        frames_[s].n_ghost_ = 0;
        n_local_[s] = frames_[s].n_local_;
    }

    for (int src = 0; src < N; ++src) {
        const int nl = frames_[src].n_local_;
        if (nl == 0)
            continue;
        dm_.set_device(src);
        // Per-source send buffers (worst case = all owned atoms).
        real_t *sx, *sy, *sz;
        std::int64_t* stag;
        int *stype, *scnt;
        cudaMalloc(&sx, nl * sizeof(real_t));
        cudaMalloc(&sy, nl * sizeof(real_t));
        cudaMalloc(&sz, nl * sizeof(real_t));
        cudaMalloc(&stag, nl * sizeof(std::int64_t));
        cudaMalloc(&stype, nl * sizeof(int));
        cudaMalloc(&scnt, sizeof(int));
        auto& a = frames_[src].atoms;
        const auto lo = dd_->sublo(src);
        const auto hi = dd_->subhi(src);

        for (int face = 0; face < 6; ++face) {
            const int dst = dd_->face_neighbor(src, face);
            if (dst == src)
                continue; // dims==1 along this axis: periodic self, no ghost
            const int axis = face / 2;
            const int dir = (face & 1) ? +1 : -1;
            const double plane = (dir > 0) ? hi[axis] : lo[axis];

            CKM(cudaMemset(scnt, 0, sizeof(int)));
            int block = 256, grid = (nl + block - 1) / block;
            pack_face<<<grid, block>>>(a.x, a.y, a.z, a.tag, a.type, nl, axis,
                                       dir, plane, cutoff, sx, sy, sz, stag,
                                       stype, scnt, nl);
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
                    "MultiGpu::exchange_halo: ghost capacity exceeded on slot "
                    + std::to_string(dst) + " (raise ghost_headroom)");
            const int sd = dm_.device_id(src), dd_id = dm_.device_id(dst);
            CKM(cudaMemcpyPeer(da.x + base, dd_id, sx, sd, cnt * sizeof(real_t)));
            CKM(cudaMemcpyPeer(da.y + base, dd_id, sy, sd, cnt * sizeof(real_t)));
            CKM(cudaMemcpyPeer(da.z + base, dd_id, sz, sd, cnt * sizeof(real_t)));
            CKM(cudaMemcpyPeer(da.tag + base, dd_id, stag, sd,
                               cnt * sizeof(std::int64_t)));
            CKM(cudaMemcpyPeer(da.type + base, dd_id, stype, sd,
                               cnt * sizeof(int)));
            goff[dst] += cnt;
            frames_[dst].n_ghost_ = goff[dst];
        }
        cudaSetDevice(dm_.device_id(src));
        cudaFree(sx);
        cudaFree(sy);
        cudaFree(sz);
        cudaFree(stag);
        cudaFree(stype);
        cudaFree(scnt);
    }
    // Ensure all peer copies completed before callers read ghosts.
    for (int s = 0; s < N; ++s) {
        dm_.set_device(s);
        cudaDeviceSynchronize();
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
