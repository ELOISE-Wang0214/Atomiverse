/// @file cuda/potentials/dihedral_multi_harmonic.cu
/// @brief CUDA DihedralMultiHarmonic — V = Σₙ Cₙ cosⁿ(φ).

#include "atomiverse/cuda/launch.cuh"
#include "atomiverse/cuda/potentials/dihedral_multi_harmonic.cuh"
#include "atomiverse/cuda/types.cuh"
#include <cmath>
#include <cuda_runtime.h>

namespace ATV {
namespace cuda {

struct DihedralMultiHarmonic::State {
    work_t* d_C = nullptr; // n_types * 5 coefficients, flat
    int n_types = 0;
};

DihedralMultiHarmonic::DihedralMultiHarmonic(
    std::vector<std::array<double, 5>> C)
    : C_(std::move(C)), state_(std::make_unique<State>()) {}

DihedralMultiHarmonic::~DihedralMultiHarmonic() = default;
DihedralMultiHarmonic::DihedralMultiHarmonic(DihedralMultiHarmonic&&) noexcept =
    default;
DihedralMultiHarmonic&
DihedralMultiHarmonic::operator=(DihedralMultiHarmonic&&) noexcept = default;

namespace {

__global__ void multi_harmonic_dihedral_k(
    const real_t* __restrict__ x, const real_t* __restrict__ y,
    const real_t* __restrict__ z, real_t* __restrict__ fx,
    real_t* __restrict__ fy, real_t* __restrict__ fz,
    prec_t* __restrict__ energy, const int* __restrict__ di,
    const int* __restrict__ dj, const int* __restrict__ dk,
    const int* __restrict__ dl, const int* __restrict__ dt,
    const work_t* __restrict__ C, // [n_types * 5]
    const prec_t* __restrict__ box, int n_dihed) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_dihed)
        return;

    int i = di[idx], j = dj[idx], k = dk[idx], l = dl[idx];
    int type = dt[idx];
    const work_t* Cn = &C[type * 5];

    // Vectors along the backbone
    work_t b1x = x[j] - x[i], b1y = y[j] - y[i], b1z = z[j] - z[i];
    work_t b2x = x[k] - x[j], b2y = y[k] - y[j], b2z = z[k] - z[j];
    work_t b3x = x[l] - x[k], b3y = y[l] - y[k], b3z = z[l] - z[k];

    // Minimum image
    work_t bx = box[0], by = box[4], bz = box[8];
    b1x -= bx * round(b1x / bx);
    b1y -= by * round(b1y / by);
    b1z -= bz * round(b1z / bz);
    b2x -= bx * round(b2x / bx);
    b2y -= by * round(b2y / by);
    b2z -= bz * round(b2z / bz);
    b3x -= bx * round(b3x / bx);
    b3y -= by * round(b3y / by);
    b3z -= bz * round(b3z / bz);

    // Normal vectors: n1 = b1 × b2, n2 = b2 × b3
    work_t n1x = b1y * b2z - b1z * b2y;
    work_t n1y = b1z * b2x - b1x * b2z;
    work_t n1z = b1x * b2y - b1y * b2x;

    work_t n2x = b2y * b3z - b2z * b3y;
    work_t n2y = b2z * b3x - b2x * b3z;
    work_t n2z = b2x * b3y - b2y * b3x;

    work_t n1_len = sqrt(n1x * n1x + n1y * n1y + n1z * n1z);
    work_t n2_len = sqrt(n2x * n2x + n2y * n2y + n2z * n2z);
    constexpr work_t kVecEpsilon = work_t(1e-10);
    if (n1_len < kVecEpsilon || n2_len < kVecEpsilon)
        return;

    work_t n1_inv = work_t(1) / n1_len;
    work_t n2_inv = work_t(1) / n2_len;
    n1x *= n1_inv;
    n1y *= n1_inv;
    n1z *= n1_inv;
    n2x *= n2_inv;
    n2y *= n2_inv;
    n2z *= n2_inv;

    work_t cos_phi = n1x * n2x + n1y * n2y + n1z * n2z;
    cos_phi = fmin(fmax(cos_phi, work_t(-1)), work_t(1));

    // V = Σ C_n cos^n(φ)
    work_t cp = cos_phi;
    work_t v = Cn[0] + Cn[1] * cp;
    work_t cp2 = cp * cp;
    v += Cn[2] * cp2;
    v += Cn[3] * cp2 * cp;
    v += Cn[4] * cp2 * cp2;

    // dV/d(cos φ) = C1 + 2*C2*cos + 3*C3*cos^2 + 4*C4*cos^3
    work_t dv_dcos = Cn[1] + work_t(2) * Cn[2] * cp + work_t(3) * Cn[3] * cp2
                     + work_t(4) * Cn[4] * cp2 * cp;

    // Forces via the standard dihedral force decomposition
    // F_i = -dV/d(cos φ) * d(cos φ)/dr_i
    work_t b2_len = sqrt(b2x * b2x + b2y * b2y + b2z * b2z);
    if (b2_len < kVecEpsilon)
        return;

    // Force on atom i: F_i = -(dV/dcos) * |b2| / |n1|^2 * n1
    work_t fi_coeff = -dv_dcos * b2_len * n1_inv * n1_inv * n1_len;
    // This simplifies. Using standard LAMMPS-style force decomposition:
    // F_i = (dV/dcos) * |b2| * n1_hat / |n1|
    // But we already normalized n1. Let me use the standard formula.

    work_t f1 = -dv_dcos * b2_len / (n1_len * n1_len) * n1_inv;
    // Actually, let me use a simpler correct formula:
    // See Allen & Tildesley, or LAMMPS dihedral_multi_harmonic.cpp

    work_t a11 = -dv_dcos * b2_len / (n1_len * n1_len);
    work_t a22 = dv_dcos * b2_len / (n2_len * n2_len);

    // b1 dot b2 and b3 dot b2
    work_t b1b2 = b1x * b2x + b1y * b2y + b1z * b2z;
    work_t b3b2 = b3x * b2x + b3y * b2y + b3z * b2z;
    work_t b2_len2 = b2_len * b2_len;

    work_t a12 = a11 * b1b2 / b2_len2;
    work_t a23 = a22 * b3b2 / b2_len2;

    // Force on i
    work_t fix = a11 * n1x;
    work_t fiy = a11 * n1y;
    work_t fiz = a11 * n1z;

    // Force on l
    work_t flx = a22 * n2x;
    work_t fly = a22 * n2y;
    work_t flz = a22 * n2z;

    // Force on j = -F_i + a12*n1 - a23*n2  (from chain rule)
    work_t fjx = -fix + a12 * n1x - a23 * n2x;
    work_t fjy = -fiy + a12 * n1y - a23 * n2y;
    work_t fjz = -fiz + a12 * n1z - a23 * n2z;

    // Force on k = -F_l - a12*n1 + a23*n2
    work_t fkx = -flx - a12 * n1x + a23 * n2x;
    work_t fky = -fly - a12 * n1y + a23 * n2y;
    work_t fkz = -flz - a12 * n1z + a23 * n2z;

    atomicAdd(&fx[i], static_cast<real_t>(fix));
    atomicAdd(&fy[i], static_cast<real_t>(fiy));
    atomicAdd(&fz[i], static_cast<real_t>(fiz));
    atomicAdd(&fx[j], static_cast<real_t>(fjx));
    atomicAdd(&fy[j], static_cast<real_t>(fjy));
    atomicAdd(&fz[j], static_cast<real_t>(fjz));
    atomicAdd(&fx[k], static_cast<real_t>(fkx));
    atomicAdd(&fy[k], static_cast<real_t>(fky));
    atomicAdd(&fz[k], static_cast<real_t>(fkz));
    atomicAdd(&fx[l], static_cast<real_t>(flx));
    atomicAdd(&fy[l], static_cast<real_t>(fly));
    atomicAdd(&fz[l], static_cast<real_t>(flz));

    // Energy reduction
    extern __shared__ work_t sdata[];
    sdata[threadIdx.x] = v;
    __syncthreads();
    for (int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (threadIdx.x < s)
            sdata[threadIdx.x] += sdata[threadIdx.x + s];
        __syncthreads();
    }
    if (threadIdx.x == 0)
        atomicAdd(energy, static_cast<prec_t>(sdata[0]));
}

} // anonymous namespace

void DihedralMultiHarmonic::build(const BuildContext& /*ctx*/) {
    int nt = static_cast<int>(C_.size());
    state_->n_types = nt;
    // Flatten coefficients: [C0_0, C1_0, ..., C4_0, C0_1, ...]
    std::vector<work_t> flat(nt * 5);
    for (int t = 0; t < nt; ++t)
        for (int c = 0; c < 5; ++c)
            flat[t * 5 + c] = C_[t][c];
    cudaMalloc(&state_->d_C, nt * 5 * sizeof(work_t));
    cudaMemcpy(state_->d_C, flat.data(), nt * 5 * sizeof(work_t),
               cudaMemcpyHostToDevice);
}

void DihedralMultiHarmonic::compute(const LaunchContext& ctx) {
    auto stream = static_cast<cudaStream_t>(ctx.runtime.stream);
    auto& f = ctx.frame;
    int nd = f.n_dihedrals();
    if (nd == 0)
        return;

    int block = kDefaultBlock;
    int grid = grid_size(static_cast<int>(nd));
    multi_harmonic_dihedral_k<<<grid, block, block * sizeof(work_t), stream>>>(
        f.atoms.x, f.atoms.y, f.atoms.z, f.atoms.fx, f.atoms.fy, f.atoms.fz,
        f.energy, f.dihedrals.atomi, f.dihedrals.atomj, f.dihedrals.atomk,
        f.dihedrals.atoml, f.dihedrals.type, state_->d_C, f.box, nd);
}

void DihedralMultiHarmonic::destroy(const DestroyContext& /*ctx*/) {
    cudaFree(state_->d_C);
    *state_ = {};
}

} // namespace cuda
} // namespace ATV
