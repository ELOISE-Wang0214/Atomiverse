/// @file md/neighborlist/nlist_direct.cu
/// @brief Direct cell-bin FULL per-atom Verlet-list build (LAMMPS-style).
///
/// Replaces the cluster-CSR → build_verlet_list pipeline (which spent ~65% of
/// GPU time building a redundant second neighbor structure). One pass:
///   1. bin_atoms_k       — bin each atom into its cell (atomic slot claim).
///   2. build_verlet_direct_k — one thread per atom: scan the precomputed
///      stencil cells, append every j within (cutoff+skin) to the FULL
///      per-atom list (column-major vl_nbr[k*stride+i], coalesced across atoms).
/// Both auto-grow (bins / verlet columns) on overflow, mirroring LAMMPS'
/// resize-and-rebuild. Orthorhombic fast path; triclinic via Box::min_image.

#include "nlist_impl.cuh"

#include "atomiverse/cuda/box.cuh"
#include "atomiverse/cuda/launch.cuh"
#include "atomiverse/cuda/types.cuh"
#include <cuda_runtime.h>

namespace ATV::cuda::md::nlist {

namespace {

// Uniform orthorhombic cell index, matching setup_cell's convention:
//   cell = cx*(dy*dz) + cy*dz + cz,  cx the x-index (slowest), cz fastest.
// Positions are wrapped to [0,L) so the binning agrees with the stencil's
// periodic wrap. `cx/cy/cz` are returned for stencil neighbor enumeration.
__device__ inline int cell_of(real_t x, real_t y, real_t z, real_t Lx,
                              real_t Ly, real_t Lz, int dx, int dy, int dz,
                              bool periodic, int& cx, int& cy, int& cz) {
    real_t fx = x / Lx, fy = y / Ly, fz = z / Lz;
    if (periodic) {
        fx -= floorf(fx);
        fy -= floorf(fy);
        fz -= floorf(fz);
    }
    cx = static_cast<int>(fx * dx);
    cy = static_cast<int>(fy * dy);
    cz = static_cast<int>(fz * dz);
    cx = cx < 0 ? 0 : (cx >= dx ? dx - 1 : cx);
    cy = cy < 0 ? 0 : (cy >= dy ? dy - 1 : cy);
    cz = cz < 0 ? 0 : (cz >= dz ? dz - 1 : cz);
    return (cx * dy + cy) * dz + cz;
}

__global__ void bin_atoms_k(const real_t* __restrict__ x,
                            const real_t* __restrict__ y,
                            const real_t* __restrict__ z, int n, cuda::Box box,
                            int dx, int dy, int dz, bool periodic,
                            int* __restrict__ bincount, int* __restrict__ bins,
                            int bin_cap, int* __restrict__ overflow) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n)
        return;
    int cx, cy, cz;
    int cell = cell_of(x[i], y[i], z[i], real_t(box.cell[0]), real_t(box.cell[4]),
                       real_t(box.cell[8]), dx, dy, dz, periodic, cx, cy, cz);
    int slot = atomicAdd(&bincount[cell], 1);
    if (slot < bin_cap)
        bins[static_cast<std::size_t>(cell) * bin_cap + slot] = i;
    else
        atomicMax(overflow, slot + 1);
}

template <bool Ortho, typename work_t>
__global__ void build_verlet_direct_k(
    const real_t* __restrict__ x, const real_t* __restrict__ y,
    const real_t* __restrict__ z, int n, cuda::Box box, int dx, int dy, int dz,
    bool periodic, const int* __restrict__ sdx, const int* __restrict__ sdy,
    const int* __restrict__ sdz, const int* __restrict__ n_stencil_ptr,
    const int* __restrict__ bincount, const int* __restrict__ bins, int bin_cap,
    work_t rsq_search, const int* __restrict__ excl_offsets,
    const int* __restrict__ excl_partners, int* __restrict__ vl_nbr,
    int* __restrict__ vl_cnt, int vl_stride, int vl_max,
    int* __restrict__ overflow) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n)
        return;
    const real_t xi = x[i], yi = y[i], zi = z[i];
    int cx, cy, cz;
    cell_of(xi, yi, zi, real_t(box.cell[0]), real_t(box.cell[4]),
            real_t(box.cell[8]), dx, dy, dz, periodic, cx, cy, cz);
    const int nst = *n_stencil_ptr;
    const int ncells_yz = dy * dz;
    const work_t rsq = rsq_search;
    // LAMMPS-style special-bonds baking: drop 1-2/1-3/1-4 partners of i here
    // (build time) so the force kernel needs no per-pair exclusion scan. Only
    // reached when the caller enables it (excl_offsets != nullptr, s_excl == 0).
    const bool has_excl = (excl_offsets != nullptr);
    const int ex_beg = has_excl ? excl_offsets[i] : 0;
    const int ex_end = has_excl ? excl_offsets[i + 1] : 0;
    int cnt = 0;
    for (int s = 0; s < nst; ++s) {
        int ncx = cx + sdx[s], ncy = cy + sdy[s], ncz = cz + sdz[s];
        if (periodic) {
            ncx = ((ncx % dx) + dx) % dx;
            ncy = ((ncy % dy) + dy) % dy;
            ncz = ((ncz % dz) + dz) % dz;
        } else if (ncx < 0 || ncx >= dx || ncy < 0 || ncy >= dy || ncz < 0
                   || ncz >= dz) {
            continue;
        }
        const int ncell = ncx * ncells_yz + ncy * dz + ncz;
        int bc = bincount[ncell];
        if (bc > bin_cap)
            bc = bin_cap;
        for (int t = 0; t < bc; ++t) {
            int j = bins[static_cast<std::size_t>(ncell) * bin_cap + t];
            if (j == i)
                continue;
            work_t ddx = work_t(x[j]) - work_t(xi);
            work_t ddy = work_t(y[j]) - work_t(yi);
            work_t ddz = work_t(z[j]) - work_t(zi);
            if (periodic) {
                if constexpr (Ortho) {
                    ddx -= work_t(box.cell[0])
                           * work_t(nearbyint(ddx * work_t(box.inv_l[0])));
                    ddy -= work_t(box.cell[4])
                           * work_t(nearbyint(ddy * work_t(box.inv_l[1])));
                    ddz -= work_t(box.cell[8])
                           * work_t(nearbyint(ddz * work_t(box.inv_l[2])));
                } else {
                    box.min_image(ddx, ddy, ddz);
                }
            }
            work_t r2 = ddx * ddx + ddy * ddy + ddz * ddz;
            if (r2 < rsq) {
                // Skip excluded (special-bonds) partners: short linear scan of
                // i's exclusion list (~6 entries for a chain) — done once per
                // rebuild here instead of every step in the force kernel.
                if (has_excl) {
                    bool excluded = false;
                    for (int e = ex_beg; e < ex_end; ++e)
                        if (excl_partners[e] == j) {
                            excluded = true;
                            break;
                        }
                    if (excluded)
                        continue;
                }
                if (cnt < vl_max)
                    vl_nbr[static_cast<std::size_t>(cnt) * vl_stride + i] = j;
                ++cnt;
            }
        }
    }
    vl_cnt[i] = cnt;
    atomicMax(overflow, cnt);
}

} // namespace

void launch_bin_atoms(const real_t* x, const real_t* y, const real_t* z, int n,
                      ATV::cuda::Box box, int dx, int dy, int dz, bool periodic,
                      int* bincount, int* bins, int bin_cap, int* overflow,
                      cudaStream_t stream) {
    if (n <= 0)
        return;
    int block = ATV::cuda::kDefaultBlock;
    int grid = (n + block - 1) / block;
    bin_atoms_k<<<grid, block, 0, stream>>>(x, y, z, n, box, dx, dy, dz,
                                            periodic, bincount, bins, bin_cap,
                                            overflow);
    ATV_CUDA_CHECK(cudaGetLastError());
}

template <typename work_t>
void launch_build_verlet_direct(
    const real_t* x, const real_t* y, const real_t* z, int n, ATV::cuda::Box box,
    int dx, int dy, int dz, bool periodic, const int* sdx, const int* sdy,
    const int* sdz, const int* n_stencil_ptr, const int* bincount,
    const int* bins, int bin_cap, work_t rsq_search, const int* excl_offsets,
    const int* excl_partners, int* vl_nbr, int* vl_cnt, int vl_stride,
    int vl_max, int* overflow, cudaStream_t stream) {
    if (n <= 0)
        return;
    int block = ATV::cuda::kDefaultBlock;
    int grid = (n + block - 1) / block;
    if (box.ortho) {
        build_verlet_direct_k<true, work_t><<<grid, block, 0, stream>>>(
            x, y, z, n, box, dx, dy, dz, periodic, sdx, sdy, sdz, n_stencil_ptr,
            bincount, bins, bin_cap, rsq_search, excl_offsets, excl_partners,
            vl_nbr, vl_cnt, vl_stride, vl_max, overflow);
    } else {
        build_verlet_direct_k<false, work_t><<<grid, block, 0, stream>>>(
            x, y, z, n, box, dx, dy, dz, periodic, sdx, sdy, sdz, n_stencil_ptr,
            bincount, bins, bin_cap, rsq_search, excl_offsets, excl_partners,
            vl_nbr, vl_cnt, vl_stride, vl_max, overflow);
    }
    ATV_CUDA_CHECK(cudaGetLastError());
}

namespace {
[[maybe_unused]] const void* _inst_launch_build_verlet_direct =
    reinterpret_cast<const void*>(&launch_build_verlet_direct<work_t>);
} // namespace

} // namespace ATV::cuda::md::nlist
