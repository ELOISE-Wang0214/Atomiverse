/// @file bonded/multi_harmonic_dihedral.cu
/// @brief MultiHarmonicDihedral -- V(φ) = Σ C_n cos^n(φ).

#include "atomiverse/cuda/potentials/multi_harmonic_dihedral.cuh"
#include "atomiverse/cuda/launch.cuh"
#include <cuda_runtime.h>
#include <cmath>

namespace ATV {

struct MultiHarmonicDihedral::State {
    Real* d_C = nullptr;  // n_types * 5 coefficients, flat
    int n_types = 0;
};

MultiHarmonicDihedral::MultiHarmonicDihedral(std::vector<std::array<Real, 5>> C)
    : C_(std::move(C)), state_(std::make_unique<State>()) {}

MultiHarmonicDihedral::~MultiHarmonicDihedral() = default;
MultiHarmonicDihedral::MultiHarmonicDihedral(MultiHarmonicDihedral&&) noexcept = default;
MultiHarmonicDihedral& MultiHarmonicDihedral::operator=(MultiHarmonicDihedral&&) noexcept = default;

namespace {

__global__ void multi_harmonic_dihedral_k(
    const Real* __restrict__ x, const Real* __restrict__ y, const Real* __restrict__ z,
    Real* __restrict__ fx, Real* __restrict__ fy, Real* __restrict__ fz,
    Real* __restrict__ energy,
    const Int* __restrict__ di, const Int* __restrict__ dj,
    const Int* __restrict__ dk, const Int* __restrict__ dl,
    const Int* __restrict__ dt,
    const Real* __restrict__ C,  // [n_types * 5]
    const Real* __restrict__ box, Int n_dihed)
{
    Int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_dihed) return;

    Int i = di[idx], j = dj[idx], k = dk[idx], l = dl[idx];
    Int type = dt[idx];
    const Real* Cn = &C[type * 5];

    // Vectors along the backbone
    Real b1x = x[j] - x[i], b1y = y[j] - y[i], b1z = z[j] - z[i];
    Real b2x = x[k] - x[j], b2y = y[k] - y[j], b2z = z[k] - z[j];
    Real b3x = x[l] - x[k], b3y = y[l] - y[k], b3z = z[l] - z[k];

    // Minimum image
    Real bx = box[0], by = box[4], bz = box[8];
    b1x -= bx * round(b1x / bx); b1y -= by * round(b1y / by); b1z -= bz * round(b1z / bz);
    b2x -= bx * round(b2x / bx); b2y -= by * round(b2y / by); b2z -= bz * round(b2z / bz);
    b3x -= bx * round(b3x / bx); b3y -= by * round(b3y / by); b3z -= bz * round(b3z / bz);

    // Normal vectors: n1 = b1 × b2, n2 = b2 × b3
    Real n1x = b1y * b2z - b1z * b2y;
    Real n1y = b1z * b2x - b1x * b2z;
    Real n1z = b1x * b2y - b1y * b2x;

    Real n2x = b2y * b3z - b2z * b3y;
    Real n2y = b2z * b3x - b2x * b3z;
    Real n2z = b2x * b3y - b2y * b3x;

    Real n1_len = sqrt(n1x * n1x + n1y * n1y + n1z * n1z);
    Real n2_len = sqrt(n2x * n2x + n2y * n2y + n2z * n2z);
    constexpr Real kVecEpsilon = Real(1e-10);
    if (n1_len < kVecEpsilon || n2_len < kVecEpsilon) return;

    Real n1_inv = Real(1) / n1_len;
    Real n2_inv = Real(1) / n2_len;
    n1x *= n1_inv; n1y *= n1_inv; n1z *= n1_inv;
    n2x *= n2_inv; n2y *= n2_inv; n2z *= n2_inv;

    Real cos_phi = n1x * n2x + n1y * n2y + n1z * n2z;
    cos_phi = fmin(fmax(cos_phi, Real(-1)), Real(1));

    // V = Σ C_n cos^n(φ)
    Real cp = cos_phi;
    Real v = Cn[0] + Cn[1] * cp;
    Real cp2 = cp * cp;
    v += Cn[2] * cp2;
    v += Cn[3] * cp2 * cp;
    v += Cn[4] * cp2 * cp2;

    // dV/d(cos φ) = C1 + 2*C2*cos + 3*C3*cos^2 + 4*C4*cos^3
    Real dv_dcos = Cn[1] + Real(2) * Cn[2] * cp + Real(3) * Cn[3] * cp2 +
                   Real(4) * Cn[4] * cp2 * cp;

    // Forces via the standard dihedral force decomposition
    // F_i = -dV/d(cos φ) * d(cos φ)/dr_i
    Real b2_len = sqrt(b2x * b2x + b2y * b2y + b2z * b2z);
    if (b2_len < kVecEpsilon) return;

    // Force on atom i: F_i = -(dV/dcos) * |b2| / |n1|^2 * n1
    Real fi_coeff = -dv_dcos * b2_len * n1_inv * n1_inv * n1_len;
    // This simplifies. Using standard LAMMPS-style force decomposition:
    // F_i = (dV/dcos) * |b2| * n1_hat / |n1|
    // But we already normalized n1. Let me use the standard formula.

    Real f1 = -dv_dcos * b2_len / (n1_len * n1_len) * n1_inv;
    // Actually, let me use a simpler correct formula:
    // See Allen & Tildesley, or LAMMPS dihedral_multi_harmonic.cpp

    Real a11 = -dv_dcos * b2_len / (n1_len * n1_len);
    Real a22 = dv_dcos * b2_len / (n2_len * n2_len);

    // b1 dot b2 and b3 dot b2
    Real b1b2 = b1x * b2x + b1y * b2y + b1z * b2z;
    Real b3b2 = b3x * b2x + b3y * b2y + b3z * b2z;
    Real b2_len2 = b2_len * b2_len;

    Real a12 = a11 * b1b2 / b2_len2;
    Real a23 = a22 * b3b2 / b2_len2;

    // Force on i
    Real fix = a11 * n1x;
    Real fiy = a11 * n1y;
    Real fiz = a11 * n1z;

    // Force on l
    Real flx = a22 * n2x;
    Real fly = a22 * n2y;
    Real flz = a22 * n2z;

    // Force on j = -F_i + a12*n1 - a23*n2  (from chain rule)
    Real fjx = -fix + a12 * n1x - a23 * n2x;
    Real fjy = -fiy + a12 * n1y - a23 * n2y;
    Real fjz = -fiz + a12 * n1z - a23 * n2z;

    // Force on k = -F_l - a12*n1 + a23*n2
    Real fkx = -flx - a12 * n1x + a23 * n2x;
    Real fky = -fly - a12 * n1y + a23 * n2y;
    Real fkz = -flz - a12 * n1z + a23 * n2z;

    atomicAdd(&fx[i], fix); atomicAdd(&fy[i], fiy); atomicAdd(&fz[i], fiz);
    atomicAdd(&fx[j], fjx); atomicAdd(&fy[j], fjy); atomicAdd(&fz[j], fjz);
    atomicAdd(&fx[k], fkx); atomicAdd(&fy[k], fky); atomicAdd(&fz[k], fkz);
    atomicAdd(&fx[l], flx); atomicAdd(&fy[l], fly); atomicAdd(&fz[l], flz);

    // Energy reduction
    extern __shared__ Real sdata[];
    sdata[threadIdx.x] = v;
    __syncthreads();
    for (int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (threadIdx.x < s) sdata[threadIdx.x] += sdata[threadIdx.x + s];
        __syncthreads();
    }
    if (threadIdx.x == 0) atomicAdd(energy, sdata[0]);
}

} // anonymous namespace

void MultiHarmonicDihedral::build(const BuildContext& ctx) {
    int nt = static_cast<int>(C_.size());
    state_->n_types = nt;
    // Flatten coefficients: [C0_0, C1_0, ..., C4_0, C0_1, ...]
    std::vector<Real> flat(nt * 5);
    for (int t = 0; t < nt; ++t)
        for (int c = 0; c < 5; ++c)
            flat[t * 5 + c] = C_[t][c];
    cudaMalloc(&state_->d_C, nt * 5 * sizeof(Real));
    cudaMemcpy(state_->d_C, flat.data(), nt * 5 * sizeof(Real), cudaMemcpyHostToDevice);
}

void MultiHarmonicDihedral::launch(const LaunchContext& ctx) {
    auto stream = static_cast<cudaStream_t>(ctx.runtime.stream);
    auto& f = ctx.frame;
    Int nd = f.n_dihedrals();
    if (nd == 0) return;

    int block = cuda::kDefaultBlock;
    int grid = cuda::grid_size(static_cast<int>(nd));
    multi_harmonic_dihedral_k<<<grid, block, block * sizeof(Real), stream>>>(
        f.x, f.y, f.z, f.fx, f.fy, f.fz, f.energy,
        f.dihedrals.i, f.dihedrals.j, f.dihedrals.k, f.dihedrals.l,
        f.dihedrals.type, state_->d_C, f.box, nd);
}

void MultiHarmonicDihedral::destroy(const DestroyContext& ctx) {
    cudaFree(state_->d_C);
    *state_ = {};
}

} // namespace ATV
