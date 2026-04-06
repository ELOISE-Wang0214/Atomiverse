#include <benchmark/benchmark.h>
#include <cmath>
#include <cstdlib>
#include <vector>

#include "atomiverse/cuda/potentials/lj.cuh"
#include "atomiverse/cuda/md/md.cuh"
#include "atomiverse/cuda/neighbor_list.cuh"
#include "atomiverse/cuda/md/velocity_verlet.cuh"

using namespace ATV;

static void BM_MDStep(benchmark::State& state) {
    int n = state.range(0);
    Real cutoff = Real(3.0);
    Real box_len = Real(std::cbrt(n / 0.8));

    cpu::Frame frame(n);
    std::srand(42);
    for (int i = 0; i < n; ++i) {
        frame.atoms.x[i] = box_len * Real(std::rand()) / Real(RAND_MAX);
        frame.atoms.y[i] = box_len * Real(std::rand()) / Real(RAND_MAX);
        frame.atoms.z[i] = box_len * Real(std::rand()) / Real(RAND_MAX);
        frame.atoms.mass[i] = Real(1.0);
    }
    frame.box = Mat3::cubic(static_cast<double>(box_len));

    md::MD sim;
    sim.add(md::VelocityVerlet{Real(0.001)});
    sim.on(md::Phase::force)
        .add(NeighborList{cutoff}, -10)
        .add(LJ{Real(1.0), Real(1.0), cutoff});

    sim.build();
    sim.run(frame, 10); // warm up

    for (auto _ : state) {
        sim.run(frame, 1000);
        cudaDeviceSynchronize();
    }

    state.SetItemsProcessed(state.iterations() * 1000);
}

BENCHMARK(BM_MDStep)->Arg(100)->Arg(1000)->Arg(10000)->Unit(benchmark::kMillisecond);

BENCHMARK_MAIN();
