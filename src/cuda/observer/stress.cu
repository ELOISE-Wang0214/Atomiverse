/// @file observer/stress.cu
/// @brief VirialStress -- kinetic + pair virial stress tensor.

#include "atomiverse/core/units.h"
#include "atomiverse/cuda/launch.cuh"
#include "atomiverse/cuda/observer/stress.cuh"
#include "atomiverse/cuda/types.cuh"
#include <cuda_runtime.h>

namespace ATV {
namespace cuda {

using cuda::prec_t;
using cuda::real_t;
using cuda::work_t;

struct VirialStress::State {
    prec_t* d_stress_tmp = nullptr; // [6] accumulation buffer
    int n_atoms = 0;
};

VirialStress::VirialStress() : state_(std::make_unique<State>()) {}
VirialStress::~VirialStress() = default;
VirialStress::VirialStress(VirialStress&&) noexcept = default;
VirialStress& VirialStress::operator=(VirialStress&&) noexcept = default;

namespace {

// Kinetic stress tensor Σ m·v·v in eV units (kMvv2e folded in).
// Convention matches `cpu::md::VirialStress`: writes RAW Σ m·v·v · kMvv2e
// (NOT divided by V) to `f.scalars->stress_kinetic[6]`. Volume division
// happens downstream in pressure-computing components (NoseHooverNPT,
// thermo recorders).
__global__ void kinetic_stress_k(const real_t* __restrict__ vx,
                                 const real_t* __restrict__ vy,
                                 const real_t* __restrict__ vz,
                                 const real_t* __restrict__ mass,
                                 prec_t* __restrict__ stress_kinetic,
                                 int n) {
    extern __shared__ work_t sdata[]; // 6 * blockDim.x
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    work_t s[6] = {};
    if (idx < n) {
        work_t m = static_cast<work_t>(mass[idx]);
        work_t vxi = static_cast<work_t>(vx[idx]),
               vyi = static_cast<work_t>(vy[idx]),
               vzi = static_cast<work_t>(vz[idx]);
        s[0] = m * vxi * vxi; // xx
        s[1] = m * vyi * vyi; // yy
        s[2] = m * vzi * vzi; // zz
        s[3] = m * vxi * vyi; // xy
        s[4] = m * vxi * vzi; // xz
        s[5] = m * vyi * vzi; // yz
    }

    for (int c = 0; c < 6; ++c)
        sdata[c * blockDim.x + threadIdx.x] = s[c];
    __syncthreads();
    for (int stride = blockDim.x / 2; stride > 0; stride >>= 1) {
        if (threadIdx.x < stride)
            for (int c = 0; c < 6; ++c)
                sdata[c * blockDim.x + threadIdx.x] +=
                    sdata[c * blockDim.x + threadIdx.x + stride];
        __syncthreads();
    }
    if (threadIdx.x == 0) {
        // Apply mvv2e conversion (amu·(Å/fs)² → eV). No volume division —
        // that's the consumer's job. Result is Σ m·v_a·v_b in eV.
        const work_t mvv2e = static_cast<work_t>(kMvv2e);
        for (int c = 0; c < 6; ++c)
            atomicAdd(&stress_kinetic[c],
                      static_cast<prec_t>(sdata[c * blockDim.x] * mvv2e));
    }
}

} // anonymous namespace

void VirialStress::build(const BuildContext& ctx) {
    state_->n_atoms = ctx.frame.n_atoms();
    cudaMalloc(&state_->d_stress_tmp, 6 * sizeof(prec_t));
}

void VirialStress::launch(const LaunchContext& ctx) {
    auto stream = static_cast<cudaStream_t>(ctx.runtime.stream);
    auto& f = ctx.frame;
    int n = f.n_atoms();
    int block = cuda::kDefaultBlock;
    int grid = cuda::grid_size(static_cast<int>(n));

    // Overwrite-only: this component is the sole writer of
    // `stress_kinetic[]`; zero it first and then add the per-atom
    // contributions. The pair / bonded virial (`stress[]`) is owned by
    // the force-phase evaluators (`pair_lj_cut`, `bond_harmonic`,
    // `angle_harmonic`, `dihedral_multi_harmonic`) and is NOT touched
    // here — clearing `stress` here previously wiped the virial just
    // accumulated by the force kernels, leaving NoseHooverNPT with a
    // kinetic-only pressure that was always positive (no LJ cohesion
    // to balance it), so the barostat scaled velocities monotonically
    // down to zero across stage 2's 50000 NPT steps.
    cudaMemsetAsync(f.scalars->stress_kinetic, 0, 6 * sizeof(prec_t), stream);

    kinetic_stress_k<<<grid, block, 6 * block * sizeof(work_t), stream>>>(
        f.atoms.vx, f.atoms.vy, f.atoms.vz, f.atoms.mass,
        f.scalars->stress_kinetic, n);
}

void VirialStress::destroy(const DestroyContext& ctx) {
    cudaFree(state_->d_stress_tmp);
    *state_ = {};
}

} // namespace cuda
} // namespace ATV
