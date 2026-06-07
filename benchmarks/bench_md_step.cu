#include <benchmark/benchmark.h>
#include <cmath>
#include <cstdlib>
#include <vector>

#include "atomiverse/cuda/md/md.cuh"
#include "atomiverse/cuda/md/velocity_verlet.cuh"
#include "atomiverse/cuda/neighbor_list.cuh"
#include "atomiverse/cuda/potentials/lj.cuh"

using namespace ATV;
// `cuda::md::MD`, `cuda::NeighborList`, `cuda::LJ` are fully qualified
// below to avoid an `md` ambiguity between `ATV::md` (phase enum) and
// `ATV::cuda::md` (driver namespace).

static void BM_MDStep(benchmark::State& state) {
    int n = state.range(0);
    double cutoff = double(3.0);
    double box_len = double(std::cbrt(n / 0.8));

    cpu::Frame frame(n);
    std::srand(42);
    for (int i = 0; i < n; ++i) {
        frame.atoms.x[i] = box_len * double(std::rand()) / double(RAND_MAX);
        frame.atoms.y[i] = box_len * double(std::rand()) / double(RAND_MAX);
        frame.atoms.z[i] = box_len * double(std::rand()) / double(RAND_MAX);
        frame.atoms.mass[i] = double(1.0);
    }
    frame.box = Mat3::cubic(static_cast<double>(box_len));

    cuda::md::MD sim;
    sim.add(cuda::md::VelocityVerlet{double(0.001)});
    sim.on(ATV::md::Phase::force)
        .add(cuda::NeighborList{cutoff}, -10)
        .add(cuda::LJ{double(1.0), double(1.0), cutoff});

    sim.build();
    sim.run(frame, 10); // warm up

    for (auto _ : state) {
        sim.run(frame, 1000);
        cudaDeviceSynchronize();
    }

    state.SetItemsProcessed(state.iterations() * 1000);
}

BENCHMARK(BM_MDStep)->Arg(100)->Arg(1000)->Arg(10000)->Unit(
    benchmark::kMillisecond);

// Large FCC-lattice LJ system at rho*=0.8 — stable, realistic neighbor
// density (~55-80 neighbors at rc=3 sigma), and large enough that the atom
// arrays exceed L2, so the uncoalesced neighbor gather is actually exercised
// (the 10000-atom case fits in L2 and hides it). target_n is rounded to the
// nearest 4*nc^3 FCC atom count.
static void BM_MDStepFCC(benchmark::State& state) {
    const long target_n = state.range(0);
    const double rho = 0.8;
    const double cutoff = double(3.0);
    const double a = std::cbrt(4.0 / rho); // FCC lattice constant (sigma)
    int nc = std::max(1, (int)std::lround(std::cbrt(target_n / 4.0)));
    const int n = 4 * nc * nc * nc;
    const double box_len = nc * a;

    // 4-atom FCC conventional-cell basis (fractional)
    const double bas[4][3] = {
        {0, 0, 0}, {0.5, 0.5, 0}, {0.5, 0, 0.5}, {0, 0.5, 0.5}};

    cpu::Frame frame(n);
    int idx = 0;
    for (int cx = 0; cx < nc; ++cx)
        for (int cy = 0; cy < nc; ++cy)
            for (int cz = 0; cz < nc; ++cz)
                for (int b = 0; b < 4; ++b) {
                    frame.atoms.x[idx] = (cx + bas[b][0]) * a;
                    frame.atoms.y[idx] = (cy + bas[b][1]) * a;
                    frame.atoms.z[idx] = (cz + bas[b][2]) * a;
                    frame.atoms.mass[idx] = double(1.0);
                    ++idx;
                }
    frame.box = Mat3::cubic(static_cast<double>(box_len));

    cuda::md::MD sim;
    sim.add(cuda::md::VelocityVerlet{double(0.001)});
    sim.on(ATV::md::Phase::force)
        .add(cuda::NeighborList{cutoff}, -10)
        .add(cuda::LJ{double(1.0), double(1.0), cutoff});
    // No barostat / stress observer → the driver auto-gates the pair virial off
    // (NVE fast path). This is the steady-state kernel that matters.

    sim.build();
    sim.run(frame, 10); // warm up

    for (auto _ : state) {
        sim.run(frame, 100);
        cudaDeviceSynchronize();
    }
    state.SetItemsProcessed(state.iterations() * 100 * (long)n);
    state.counters["atoms"] = n;
    state.counters["Matom-step/s"] = benchmark::Counter(
        (double)n * 100, benchmark::Counter::kIsIterationInvariantRate,
        benchmark::Counter::kIs1000);
}

BENCHMARK(BM_MDStepFCC)
    ->Arg(1000000)
    ->Arg(4000000)
    ->Unit(benchmark::kMillisecond)
    ->UseRealTime()
    ->MinTime(2.0);

BENCHMARK_MAIN();
