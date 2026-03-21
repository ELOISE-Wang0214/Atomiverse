/// @file neighborlist/setup.cu
/// @brief GPU-side cell setup: inverse cell, grid dimensions, stencil computation.
///
/// All computation runs on the GPU. No host-side math, no H2D transfers.

#include "nlist_impl.cuh"
#include "triclinic.cuh"

#include <cuda_runtime.h>
#include <cmath>

namespace ATV::md::nlist {

// ---------------------------------------------------------------------------
// Shared device logic (called by both kernel entry points)
// ---------------------------------------------------------------------------

/// @brief Compute inverse, dims, stencil from d_cell already populated.
__device__ void setup_cell_body(
    Real* __restrict__ d_cell,
    Real* __restrict__ d_inv_cell,
    int* __restrict__ d_dims,
    int* __restrict__ d_stencil_dx,
    int* __restrict__ d_stencil_dy,
    int* __restrict__ d_stencil_dz,
    int* __restrict__ d_n_stencil,
    int* __restrict__ d_ncells,
    Real cutoff_plus_skin,
    bool periodic,
    int max_stencil) {

    invert_3x3(d_cell, d_inv_cell);

    // Perpendicular widths via cross products
    Real widths[3];
    {
        Real inv[9];
        Real det = invert_3x3(d_cell, inv);
        Real volume = (det < Real(0)) ? -det : det;
        auto cross_norm = [&](int i, int j) -> Real {
            Real cx = d_cell[i*3+1]*d_cell[j*3+2] - d_cell[i*3+2]*d_cell[j*3+1];
            Real cy = d_cell[i*3+2]*d_cell[j*3+0] - d_cell[i*3+0]*d_cell[j*3+2];
            Real cz = d_cell[i*3+0]*d_cell[j*3+1] - d_cell[i*3+1]*d_cell[j*3+0];
            return sqrt(cx*cx + cy*cy + cz*cz);
        };
        widths[0] = volume / cross_norm(1, 2);
        widths[1] = volume / cross_norm(0, 2);
        widths[2] = volume / cross_norm(0, 1);
    }

    // Grid dimensions
    Real rs = cutoff_plus_skin;
    int dx = max(3, int(widths[0] / rs));
    int dy = max(3, int(widths[1] / rs));
    int dz = max(3, int(widths[2] / rs));
    d_dims[0] = dx; d_dims[1] = dy; d_dims[2] = dz;
    *d_ncells = dx * dy * dz;

    // Stencil offsets with periodic deduplication
    int range[3];
    int dim_arr[3] = {dx, dy, dz};
    for (int d = 0; d < 3; ++d) {
        Real row_norm = sqrt(d_inv_cell[d*3+0]*d_inv_cell[d*3+0] +
                             d_inv_cell[d*3+1]*d_inv_cell[d*3+1] +
                             d_inv_cell[d*3+2]*d_inv_cell[d*3+2]);
        range[d] = max(1, int(rs * row_norm * Real(dim_arr[d])) + 1);
    }
    Real fx = Real(1)/Real(dx), fy = Real(1)/Real(dy), fz = Real(1)/Real(dz);
    Real cut_sq = rs * rs;
    int count = 0;

    for (int di = -range[0]; di <= range[0]; ++di)
      for (int dj = -range[1]; dj <= range[1]; ++dj)
        for (int dk = -range[2]; dk <= range[2]; ++dk) {
            Real cfx = Real(di)*fx, cfy = Real(dj)*fy, cfz = Real(dk)*fz;
            Real cx = d_cell[0]*cfx + d_cell[1]*cfy + d_cell[2]*cfz;
            Real cy = d_cell[3]*cfx + d_cell[4]*cfy + d_cell[5]*cfz;
            Real cz = d_cell[6]*cfx + d_cell[7]*cfy + d_cell[8]*cfz;
            Real hwx = Real(0.5)*(fabs(d_cell[0]*fx) + fabs(d_cell[1]*fy) + fabs(d_cell[2]*fz));
            Real hwy = Real(0.5)*(fabs(d_cell[3]*fx) + fabs(d_cell[4]*fy) + fabs(d_cell[5]*fz));
            Real hwz = Real(0.5)*(fabs(d_cell[6]*fx) + fabs(d_cell[7]*fy) + fabs(d_cell[8]*fz));
            Real gx = fmax(Real(0), fabs(cx) - Real(2)*hwx);
            Real gy = fmax(Real(0), fabs(cy) - Real(2)*hwy);
            Real gz = fmax(Real(0), fabs(cz) - Real(2)*hwz);
            if (gx*gx + gy*gy + gz*gz <= cut_sq) {
                bool dup = false;
                if (periodic) {
                    int wx = ((di%dx)+dx)%dx, wy = ((dj%dy)+dy)%dy, wz = ((dk%dz)+dz)%dz;
                    for (int e = 0; e < count; ++e) {
                        if (((d_stencil_dx[e]%dx)+dx)%dx == wx &&
                            ((d_stencil_dy[e]%dy)+dy)%dy == wy &&
                            ((d_stencil_dz[e]%dz)+dz)%dz == wz) { dup = true; break; }
                    }
                }
                if (!dup && count < max_stencil) {
                    d_stencil_dx[count] = di;
                    d_stencil_dy[count] = dj;
                    d_stencil_dz[count] = dk;
                    ++count;
                }
            }
        }
    *d_n_stencil = count;
}

// ---------------------------------------------------------------------------
// Kernel entry points
// ---------------------------------------------------------------------------

/// @brief Setup from host-passed box scalars (used by tests/host-side setup).
__global__ void setup_cell_from_scalars(
    Real b00, Real b01, Real b02,
    Real b10, Real b11, Real b12,
    Real b20, Real b21, Real b22,
    Real* d_cell, Real* d_inv_cell,
    int* d_dims, int* d_stencil_dx, int* d_stencil_dy, int* d_stencil_dz,
    int* d_n_stencil, int* d_ncells,
    Real cutoff_plus_skin, bool periodic, int max_stencil) {

    d_cell[0] = b00; d_cell[1] = b01; d_cell[2] = b02;
    d_cell[3] = b10; d_cell[4] = b11; d_cell[5] = b12;
    d_cell[6] = b20; d_cell[7] = b21; d_cell[8] = b22;

    setup_cell_body(d_cell, d_inv_cell, d_dims,
                    d_stencil_dx, d_stencil_dy, d_stencil_dz,
                    d_n_stencil, d_ncells, cutoff_plus_skin, periodic, max_stencil);
}

/// @brief Setup from device box pointer (production path).
__global__ void setup_cell_from_device(
    const Real* __restrict__ d_box,
    Real* __restrict__ d_cell, Real* __restrict__ d_inv_cell,
    int* __restrict__ d_dims,
    int* __restrict__ d_stencil_dx, int* __restrict__ d_stencil_dy,
    int* __restrict__ d_stencil_dz,
    int* __restrict__ d_n_stencil, int* __restrict__ d_ncells,
    Real cutoff_plus_skin, bool periodic, int max_stencil) {

    for (int i = 0; i < 9; ++i) d_cell[i] = d_box[i];

    setup_cell_body(d_cell, d_inv_cell, d_dims,
                    d_stencil_dx, d_stencil_dy, d_stencil_dz,
                    d_n_stencil, d_ncells, cutoff_plus_skin, periodic, max_stencil);
}

/// @brief Single-thread kernel: write sentinel value at arr[idx].
__global__ void write_int_value(int* __restrict__ arr, int idx, int val) {
    if (threadIdx.x == 0) arr[idx] = val;
}

// ---------------------------------------------------------------------------
// Launch wrappers
// ---------------------------------------------------------------------------

void launch_setup_cell(
    const Real* box_data,
    Real* d_cell, Real* d_inv_cell,
    int* d_dims, int* d_stencil_dx, int* d_stencil_dy, int* d_stencil_dz,
    int* d_n_stencil, int* d_ncells,
    Real cutoff_plus_skin, bool periodic, int max_stencil,
    cudaStream_t stream) {
    setup_cell_from_scalars<<<1, 1, 0, stream>>>(
        box_data[0], box_data[1], box_data[2],
        box_data[3], box_data[4], box_data[5],
        box_data[6], box_data[7], box_data[8],
        d_cell, d_inv_cell, d_dims,
        d_stencil_dx, d_stencil_dy, d_stencil_dz,
        d_n_stencil, d_ncells,
        cutoff_plus_skin, periodic, max_stencil);
    ATV_CUDA_CHECK(cudaGetLastError());
}

void launch_setup_cell_device(
    const Real* d_box,
    Real* d_cell, Real* d_inv_cell,
    int* d_dims, int* d_stencil_dx, int* d_stencil_dy, int* d_stencil_dz,
    int* d_n_stencil, int* d_ncells,
    Real cutoff_plus_skin, bool periodic, int max_stencil,
    cudaStream_t stream) {
    setup_cell_from_device<<<1, 1, 0, stream>>>(
        d_box, d_cell, d_inv_cell, d_dims,
        d_stencil_dx, d_stencil_dy, d_stencil_dz,
        d_n_stencil, d_ncells,
        cutoff_plus_skin, periodic, max_stencil);
    ATV_CUDA_CHECK(cudaGetLastError());
}

void launch_write_int(int* arr, int idx, int val, cudaStream_t stream) {
    write_int_value<<<1, 1, 0, stream>>>(arr, idx, val);
    ATV_CUDA_CHECK(cudaGetLastError());
}

} // namespace ATV::md::nlist
