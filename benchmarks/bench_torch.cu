/// @file bench_torch.cu
/// @brief LibTorch backend benchmarks: LJ forward/backward, transfer round-trip.
///
/// Compares the torch backend against the native CUDA backend for:
///   1. LJ forward: cuda::LJ (via MD Driver) vs ATV::torch::lj
///   2. LJ backward: torch::lj forward + backward cost
///   3. Transfer round-trip: from_cpu + to_cpu for varying atom counts
///
/// All benchmarks are gated behind ATV_HAS_TORCH. When Torch is not available,
/// the file compiles to an empty main().

#ifdef ATV_HAS_TORCH

#include <benchmark/benchmark.h>
#include <cmath>
#include <cuda_runtime.h>

#include "atomiverse/torch/frame.cuh"
#include "atomiverse/torch/transfer.cuh"
#include "atomiverse/torch/potentials/lj.cuh"
#include "atomiverse/torch/neighbor_list.cuh"

#include "atomiverse/cpu/transfer.h"
#include "atomiverse/cuda/transfer.cuh"
#include "atomiverse/cuda/potentials/lj.cuh"
#include "atomiverse/cuda/md/md.cuh"
#include "atomiverse/cuda/neighbor_list.cuh"
#include "atomiverse/cuda/md/velocity_verlet.cuh"

using ATV::Real;
using ATV::Int;
using ATV::Mat3;
namespace atorch = ATV::torch;

// ── Helpers ──────────────────────────────────────────────────────────────

/// @brief Build an FCC lattice of n_cells^3 * 4 atoms in a cubic box.
static ATV::cpu::Frame make_fcc(int n_cells, Real lattice_const) {
    int n = n_cells * n_cells * n_cells * 4;
    Real box_len = lattice_const * Real(n_cells);

    ATV::cpu::Frame frame(n);

    const Real basis[4][3] = {
        {Real(0),   Real(0),   Real(0)},
        {Real(0.5), Real(0.5), Real(0)},
        {Real(0.5), Real(0),   Real(0.5)},
        {Real(0),   Real(0.5), Real(0.5)},
    };

    int idx = 0;
    for (int ix = 0; ix < n_cells; ++ix)
        for (int iy = 0; iy < n_cells; ++iy)
            for (int iz = 0; iz < n_cells; ++iz)
                for (int b = 0; b < 4; ++b) {
                    frame.x[idx] = (Real(ix) + basis[b][0]) * lattice_const;
                    frame.y[idx] = (Real(iy) + basis[b][1]) * lattice_const;
                    frame.z[idx] = (Real(iz) + basis[b][2]) * lattice_const;
                    frame.mass[idx] = Real(39.948); // Argon [amu]
                    frame.type_id[idx] = 0;
                    idx++;
                }

    frame.box = Mat3::cubic(static_cast<double>(box_len));
    return frame;
}

/// @brief Map atom count to FCC n_cells (rounded up so n_cells^3*4 >= target).
static int ncells_for(int target_atoms) {
    int nc = static_cast<int>(std::ceil(std::cbrt(target_atoms / 4.0)));
    return std::max(nc, 1);
}

// ── LJ parameters ────────────────────────────────────────────────────────

static constexpr Real kEpsilon = Real(0.0104);   // LJ well depth [eV]
static constexpr Real kSigma   = Real(3.4);      // LJ sigma [Angstrom]
static constexpr Real kCutoff  = Real(6.8);      // 2.0 * sigma [Angstrom]
static constexpr Real kLatConst = Real(5.26);    // FCC lattice constant for Ar [Angstrom]

// ═════════════════════════════════════════════════════════════════════════
// 1. LJ Forward: CUDA backend (via MD Driver) vs torch backend
// ═════════════════════════════════════════════════════════════════════════

/// @brief Torch LJ forward only (neighbor list pre-built, measures kernel cost).
static void BM_LJ_KernelOnly_Torch(benchmark::State& state) {
    int n_cells = ncells_for(state.range(0));
    auto cpu_frame = make_fcc(n_cells, kLatConst);
    int n = cpu_frame.n_atoms();

    auto tframe = atorch::from_cpu(cpu_frame);
    auto real_opts = ::torch::TensorOptions().dtype(atorch::kRealDtype).device(::torch::kCUDA);
    auto t_eps = ::torch::tensor({kEpsilon}, real_opts);
    auto t_sig = ::torch::tensor({kSigma}, real_opts);

    // Build neighbor list once (reuse across iterations)
    atorch::neighbor_list(tframe, kCutoff);
    tframe.energy = ::torch::zeros({1}, real_opts);
    tframe.fx = ::torch::zeros({n}, real_opts);
    tframe.fy = ::torch::zeros({n}, real_opts);
    tframe.fz = ::torch::zeros({n}, real_opts);
    atorch::lj(tframe, t_eps, t_sig, kCutoff);
    cudaDeviceSynchronize();

    for (auto _ : state) {
        tframe.energy.zero_();
        tframe.fx.zero_();
        tframe.fy.zero_();
        tframe.fz.zero_();
        atorch::lj(tframe, t_eps, t_sig, kCutoff);
        cudaDeviceSynchronize();
    }

    state.SetItemsProcessed(state.iterations());
    state.counters["atoms"] = n;
    state.counters["pairs"] = tframe.n_pairs();

    atorch::destroy(tframe);
}

BENCHMARK(BM_LJ_KernelOnly_Torch)
    ->Arg(100)->Arg(1000)->Arg(10000)
    ->Unit(benchmark::kMicrosecond);

// ═════════════════════════════════════════════════════════════════════════
// 2. LJ Backward: forward + backward through autograd
// ═════════════════════════════════════════════════════════════════════════

/// @brief Torch LJ forward + backward (positions require grad).
static void BM_LJ_Backward_Positions(benchmark::State& state) {
    int n_cells = ncells_for(state.range(0));
    auto cpu_frame = make_fcc(n_cells, kLatConst);
    int n = cpu_frame.n_atoms();

    auto tframe = atorch::from_cpu(cpu_frame);
    auto real_opts = ::torch::TensorOptions().dtype(atorch::kRealDtype).device(::torch::kCUDA);

    // Enable grad on positions
    tframe.x = tframe.x.detach().requires_grad_(true);
    tframe.y = tframe.y.detach().requires_grad_(true);
    tframe.z = tframe.z.detach().requires_grad_(true);

    auto t_eps = ::torch::tensor({kEpsilon}, real_opts);
    auto t_sig = ::torch::tensor({kSigma}, real_opts);

    // Build neighbor list (detached -- non-differentiable)
    atorch::neighbor_list(tframe, kCutoff);
    tframe.energy = ::torch::zeros({1}, real_opts);
    tframe.fx = ::torch::zeros({n}, real_opts);
    tframe.fy = ::torch::zeros({n}, real_opts);
    tframe.fz = ::torch::zeros({n}, real_opts);

    // Warm up
    {
        auto energy = atorch::lj(tframe, t_eps, t_sig, kCutoff);
        energy.backward();
        cudaDeviceSynchronize();
    }

    for (auto _ : state) {
        // Zero accumulators and gradients
        tframe.energy.zero_();
        tframe.fx.zero_();
        tframe.fy.zero_();
        tframe.fz.zero_();
        if (tframe.x.grad().defined()) tframe.x.grad().zero_();
        if (tframe.y.grad().defined()) tframe.y.grad().zero_();
        if (tframe.z.grad().defined()) tframe.z.grad().zero_();

        auto energy = atorch::lj(tframe, t_eps, t_sig, kCutoff);
        energy.backward();
        cudaDeviceSynchronize();
    }

    state.SetItemsProcessed(state.iterations());
    state.counters["atoms"] = n;
    state.counters["pairs"] = tframe.n_pairs();

    atorch::destroy(tframe);
}

BENCHMARK(BM_LJ_Backward_Positions)
    ->Arg(100)->Arg(1000)->Arg(10000)
    ->Unit(benchmark::kMicrosecond);

/// @brief Torch LJ forward + backward (epsilon and sigma require grad).
static void BM_LJ_Backward_Params(benchmark::State& state) {
    int n_cells = ncells_for(state.range(0));
    auto cpu_frame = make_fcc(n_cells, kLatConst);
    int n = cpu_frame.n_atoms();

    auto tframe = atorch::from_cpu(cpu_frame);
    auto real_opts = ::torch::TensorOptions().dtype(atorch::kRealDtype).device(::torch::kCUDA);

    auto t_eps = ::torch::tensor({kEpsilon}, real_opts).requires_grad_(true);
    auto t_sig = ::torch::tensor({kSigma}, real_opts).requires_grad_(true);

    // Build neighbor list
    atorch::neighbor_list(tframe, kCutoff);
    tframe.energy = ::torch::zeros({1}, real_opts);
    tframe.fx = ::torch::zeros({n}, real_opts);
    tframe.fy = ::torch::zeros({n}, real_opts);
    tframe.fz = ::torch::zeros({n}, real_opts);

    // Warm up
    {
        auto energy = atorch::lj(tframe, t_eps, t_sig, kCutoff);
        energy.backward();
        cudaDeviceSynchronize();
    }

    for (auto _ : state) {
        tframe.energy.zero_();
        tframe.fx.zero_();
        tframe.fy.zero_();
        tframe.fz.zero_();
        if (t_eps.grad().defined()) t_eps.grad().zero_();
        if (t_sig.grad().defined()) t_sig.grad().zero_();

        auto energy = atorch::lj(tframe, t_eps, t_sig, kCutoff);
        energy.backward();
        cudaDeviceSynchronize();
    }

    state.SetItemsProcessed(state.iterations());
    state.counters["atoms"] = n;
    state.counters["pairs"] = tframe.n_pairs();

    atorch::destroy(tframe);
}

BENCHMARK(BM_LJ_Backward_Params)
    ->Arg(100)->Arg(1000)->Arg(10000)
    ->Unit(benchmark::kMicrosecond);

// ═════════════════════════════════════════════════════════════════════════
// 3. Transfer round-trip: from_cpu -> to_cpu
// ═════════════════════════════════════════════════════════════════════════

/// @brief H2D transfer: cpu::Frame -> torch::Frame via from_cpu.
static void BM_Transfer_FromCpu(benchmark::State& state) {
    int n_cells = ncells_for(state.range(0));
    auto cpu_frame = make_fcc(n_cells, kLatConst);
    int n = cpu_frame.n_atoms();

    // Warm up
    {
        auto tframe = atorch::from_cpu(cpu_frame);
        cudaDeviceSynchronize();
        atorch::destroy(tframe);
    }

    for (auto _ : state) {
        auto tframe = atorch::from_cpu(cpu_frame);
        cudaDeviceSynchronize();
        benchmark::DoNotOptimize(tframe.x.data_ptr());
        atorch::destroy(tframe);
    }

    state.SetItemsProcessed(state.iterations());
    state.counters["atoms"] = n;
    // Approximate bytes transferred: 3 positions + 3 velocities + 3 forces + mass + type_id
    // = 10 Real arrays + 1 Int array ~ (10 * sizeof(Real) + sizeof(Int)) * n
    state.SetBytesProcessed(
        state.iterations() * n *
        static_cast<int64_t>(10 * sizeof(Real) + sizeof(Int)));
}

BENCHMARK(BM_Transfer_FromCpu)
    ->Arg(100)->Arg(1000)->Arg(10000)
    ->Unit(benchmark::kMicrosecond);

/// @brief D2H transfer: torch::Frame -> cpu::Frame via to_cpu.
static void BM_Transfer_ToCpu(benchmark::State& state) {
    int n_cells = ncells_for(state.range(0));
    auto cpu_frame = make_fcc(n_cells, kLatConst);
    int n = cpu_frame.n_atoms();

    auto tframe = atorch::from_cpu(cpu_frame);
    cudaDeviceSynchronize();

    ATV::cpu::Frame dst(n);

    // Warm up
    {
        atorch::to_cpu(tframe, dst);
        cudaDeviceSynchronize();
    }

    for (auto _ : state) {
        atorch::to_cpu(tframe, dst);
        cudaDeviceSynchronize();
        benchmark::DoNotOptimize(dst.x.data());
    }

    state.SetItemsProcessed(state.iterations());
    state.counters["atoms"] = n;
    state.SetBytesProcessed(
        state.iterations() * n *
        static_cast<int64_t>(10 * sizeof(Real) + sizeof(Int)));

    atorch::destroy(tframe);
}

BENCHMARK(BM_Transfer_ToCpu)
    ->Arg(100)->Arg(1000)->Arg(10000)
    ->Unit(benchmark::kMicrosecond);

/// @brief Full round-trip: from_cpu -> to_cpu.
static void BM_Transfer_RoundTrip(benchmark::State& state) {
    int n_cells = ncells_for(state.range(0));
    auto cpu_frame = make_fcc(n_cells, kLatConst);
    int n = cpu_frame.n_atoms();

    ATV::cpu::Frame dst(n);

    // Warm up
    {
        auto tframe = atorch::from_cpu(cpu_frame);
        cudaDeviceSynchronize();
        atorch::to_cpu(tframe, dst);
        cudaDeviceSynchronize();
        atorch::destroy(tframe);
    }

    for (auto _ : state) {
        auto tframe = atorch::from_cpu(cpu_frame);
        atorch::to_cpu(tframe, dst);
        cudaDeviceSynchronize();
        benchmark::DoNotOptimize(dst.x.data());
        atorch::destroy(tframe);
    }

    state.SetItemsProcessed(state.iterations());
    state.counters["atoms"] = n;
    // Round-trip: both directions
    state.SetBytesProcessed(
        state.iterations() * n * 2 *
        static_cast<int64_t>(10 * sizeof(Real) + sizeof(Int)));
}

BENCHMARK(BM_Transfer_RoundTrip)
    ->Arg(100)->Arg(1000)->Arg(10000)
    ->Unit(benchmark::kMicrosecond);

/// @brief Zero-copy wrap: cuda::Frame -> torch::Frame via to_torch.
static void BM_Transfer_ToTorch_ZeroCopy(benchmark::State& state) {
    int n_cells = ncells_for(state.range(0));
    auto cpu_frame = make_fcc(n_cells, kLatConst);
    int n = cpu_frame.n_atoms();

    auto cuda_frame = ATV::cpu::to_cuda(cpu_frame);
    cudaDeviceSynchronize();

    // Warm up
    {
        auto tframe = atorch::to_torch(cuda_frame);
        benchmark::DoNotOptimize(tframe.x.data_ptr());
    }

    for (auto _ : state) {
        auto tframe = atorch::to_torch(cuda_frame);
        benchmark::DoNotOptimize(tframe.x.data_ptr());
    }

    state.SetItemsProcessed(state.iterations());
    state.counters["atoms"] = n;

    ATV::cuda::destroy(cuda_frame);
}

BENCHMARK(BM_Transfer_ToTorch_ZeroCopy)
    ->Arg(100)->Arg(1000)->Arg(10000)
    ->Unit(benchmark::kNanosecond);

BENCHMARK_MAIN();

#else
int main() { return 0; }
#endif // ATV_HAS_TORCH
