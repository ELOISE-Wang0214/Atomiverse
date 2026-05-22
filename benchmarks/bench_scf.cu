/// @file bench_scf.cu
/// @brief SCF benchmarks: CPU (xtensor-blas) vs CUDA (cuSOLVER/cuBLAS).
///
/// Measures end-to-end SCF solve time for different matrix sizes.
/// CPU uses LAPACK via xtensor-blas; CUDA uses cuSOLVER + cuBLAS via the
/// Driver/Component architecture (includes transfer overhead).

#include <benchmark/benchmark.h>

#include "atomiverse/cpu/scf/solve.h"
#include "atomiverse/cuda/evaluate.cuh"
#include "atomiverse/cuda/scf/scf_evaluator.cuh"

using namespace ATV;

static cpu::Frame make_frame(int n) {
    cpu::Frame frame(n);
    for (int i = 0; i < n; ++i) {
        frame.atoms.x[i] = double(i);
        frame.atoms.mass[i] = double(1);
    }
    frame.box = Mat3::cubic(100.0);
    return frame;
}

// ── CPU solver (xtensor-blas / LAPACK) ────────────────────────────────

static void BM_SCF_CPU(benchmark::State& state) {
    int n = state.range(0);
    int n_occ = n / 4;

    for (auto _ : state) {
        auto result = cpu::scf::solve(n, n_occ, double(0.1), double(1e-5), 50);
        benchmark::DoNotOptimize(result.energy);
    }

    auto result = cpu::scf::solve(n, n_occ, double(0.1), double(1e-5), 50);
    state.counters["iters"] = result.iterations;
}

BENCHMARK(BM_SCF_CPU)
    ->Arg(16)
    ->Arg(32)
    ->Arg(64)
    ->Arg(128)
    ->Unit(benchmark::kMillisecond);

// ── CUDA solver (end-to-end: alloc + transfer + solve + retrieve) ─────

static void BM_SCF_CUDA(benchmark::State& state) {
    int n = state.range(0);
    int n_occ = n / 4;

    // Warm up GPU
    {
        auto frame = make_frame(n);
        cuda::Evaluate eval;
        eval.add(cuda::scf::SCFEvaluator{n_occ, double(0.1), double(1e-5), 50});
        eval.build();
        eval.run(frame);
    }

    for (auto _ : state) {
        auto frame = make_frame(n);
        cuda::Evaluate eval;
        eval.add(cuda::scf::SCFEvaluator{n_occ, double(0.1), double(1e-5), 50});
        eval.build();
        eval.run(frame);
        benchmark::DoNotOptimize(frame.energy);
    }
}

BENCHMARK(BM_SCF_CUDA)
    ->Arg(16)
    ->Arg(32)
    ->Arg(64)
    ->Arg(128)
    ->Unit(benchmark::kMillisecond);

BENCHMARK_MAIN();
