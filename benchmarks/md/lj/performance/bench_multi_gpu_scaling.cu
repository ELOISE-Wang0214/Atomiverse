/// @file benchmarks/bench_multi_gpu_scaling.cu
/// @brief Single-process multi-GPU LJ MD strong-scaling benchmark.
///
/// Builds a cubic LJ system, partitions it across N GPUs (CUDA P2P over
/// NVLink — no MPI), and times a velocity-Verlet NVE run using the cluster
/// neighbor-list force path (cuda::md::MultiGpu::run_nve_nlist). Emits a CSV
/// line so an external driver can collect 1/2/4-GPU points and report speedup.
///
/// Usage:  bench_multi_gpu_scaling [n_atoms_target] [n_gpus] [steps] [reneigh]
///   n_atoms_target  approx atom count (rounded to a cubic lattice); default 262144
///   n_gpus          GPUs to use (<=0 => all visible); default -1
///   steps           timed NVE steps; default 100
///   reneigh         steps between halo+neighbor rebuilds; default 10
///
/// Timing includes the one-time force-engine setup + first reneighbor inside
/// run_nve_nlist; a warmup run precedes the timed run so allocation/JIT costs
/// do not pollute the measurement. Over the default step count the residual
/// setup share is small and identical in spirit across GPU counts, so the
/// 1->2->4 GPU speedup ratio is representative.

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "atomiverse/core/mat3.h"
#include "atomiverse/cpu/frame.h"
#include "atomiverse/cuda/comm/multi_gpu.cuh"

using namespace ATV;

namespace {
cpu::Frame make_cubic_lj(int cells, double spacing) {
    const int n = cells * cells * cells;
    cpu::Frame f(n);
    f.box = Mat3::cubic(cells * spacing);
    int idx = 0;
    for (int ix = 0; ix < cells; ++ix)
        for (int iy = 0; iy < cells; ++iy)
            for (int iz = 0; iz < cells; ++iz) {
                // Mild deterministic jitter so it isn't a frozen lattice.
                f.atoms.x[idx] = (ix + 0.5) * spacing + 0.1 * std::sin(0.7 * idx);
                f.atoms.y[idx] = (iy + 0.5) * spacing + 0.1 * std::sin(1.3 * idx);
                f.atoms.z[idx] = (iz + 0.5) * spacing + 0.1 * std::sin(1.9 * idx);
                f.atoms.vx[idx] = 0.01 * std::sin(0.5 * idx);
                f.atoms.vy[idx] = 0.01 * std::sin(0.9 * idx);
                f.atoms.vz[idx] = 0.01 * std::sin(1.1 * idx);
                f.atoms.mass[idx] = 1.0;
                f.atoms.type[idx] = 0;
                f.atoms.tag[idx] = idx;
                ++idx;
            }
    return f;
}
} // namespace

int main(int argc, char** argv) {
    int n_target = (argc > 1) ? std::atoi(argv[1]) : 262144;
    int n_gpus = (argc > 2) ? std::atoi(argv[2]) : -1;
    int steps = (argc > 3) ? std::atoi(argv[3]) : 100;
    int reneigh = (argc > 4) ? std::atoi(argv[4]) : 10;

    const double sigma = 2.5, eps = 1.0, cut = 5.0;
    const double spacing = 1.12 * sigma; // ~LJ liquid spacing
    const double dt = 0.001;
    int cells = static_cast<int>(std::round(std::cbrt((double)n_target)));
    if (cells < 2)
        cells = 2;
    const int n = cells * cells * cells;

    cpu::Frame g = make_cubic_lj(cells, spacing);

    cuda::md::MultiGpu mg(n_gpus);
    const int ng = mg.count();
    mg.partition(g, 1.0);

    // Warmup (untimed): pays allocation / first-build / JIT.
    mg.run_nve_nlist(5, dt, cut, eps, sigma, -1.0, reneigh);

    auto t0 = std::chrono::high_resolution_clock::now();
    mg.run_nve_nlist(steps, dt, cut, eps, sigma, -1.0, reneigh);
    auto t1 = std::chrono::high_resolution_clock::now();
    double wall = std::chrono::duration<double>(t1 - t0).count();

    // Peak VRAM in use per device (total - free) after the run. peak_mb = the
    // busiest single GPU (the meaningful per-card figure); total_mb = summed
    // across the GPUs used.
    double peak_mb = 0.0, total_mb = 0.0;
    for (int d = 0; d < ng; ++d) {
        cudaSetDevice(mg.devices().device_id(d));
        size_t freeb = 0, totb = 0;
        cudaMemGetInfo(&freeb, &totb);
        double used = double(totb - freeb) / 1.0e6;
        total_mb += used;
        if (used > peak_mb)
            peak_mb = used;
    }

    double atom_steps = double(n) * double(steps);
    double aps = atom_steps / wall;
    double steps_per_s = double(steps) / wall;
    // CSV: tag,atoms,gpus,steps,wall_s,steps_per_s,Matom_steps_per_s,peak_vram_mb,total_vram_mb
    std::printf("CSV,%d,%d,%d,%.4f,%.2f,%.2f,%.1f,%.1f\n", n, ng, steps, wall,
                steps_per_s, aps / 1e6, peak_mb, total_mb);
    std::printf("# %d atoms on %d GPU(s): %.4f s, %d steps -> %.2f Matom-step/s,"
                " peak %.0f MB/GPU, total %.0f MB\n",
                n, ng, wall, steps, aps / 1e6, peak_mb, total_mb);
    return 0;
}
