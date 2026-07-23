/// @file tools/cuda_nonbonded_style_bench.cu
/// @brief Single-GPU LJ NVE wall timing for the OpenMM cluster MD path.
///
/// Usage:
///   ./cuda_nonbonded_style_bench [n_side] [n_steps] [warmup]
///
/// Prints ms/step after warmup. Not a substitute for ncu; relative only.
/// Dual-style A/B (openmm|lammps) was removed with cuda-openmm-md-01-contract.
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cuda/dynamics/md.cuh"
#include "atomiverse/cuda/dynamics/velocity_verlet.cuh"
#include "atomiverse/cuda/neighbor/neighbor_list.cuh"
#include "atomiverse/cuda/potentials/pair_lj_cut.cuh"

using namespace ATV;

static cpu::Frame make_lj(int n_side, double density) {
    int n = n_side * n_side * n_side;
    double box_len = std::cbrt(double(n) / density);
    double spacing = box_len / double(n_side);
    cpu::Frame frame(n);
    int idx = 0;
    for (int ix = 0; ix < n_side; ++ix)
        for (int iy = 0; iy < n_side; ++iy)
            for (int iz = 0; iz < n_side; ++iz) {
                frame.span(kAtom_x)[idx] = (ix + 0.5) * spacing;
                frame.span(kAtom_y)[idx] = (iy + 0.5) * spacing;
                frame.span(kAtom_z)[idx] = (iz + 0.5) * spacing;
                frame.span(kAtom_mass)[idx] = 1.0;
                frame.span(kAtom_type)[idx] = 18;
                frame.span(kAtom_vx)[idx] = 0.01 * ((ix + iy + iz) % 3 - 1);
                frame.span(kAtom_vy)[idx] = 0.01 * ((ix + 2 * iy) % 3 - 1);
                frame.span(kAtom_vz)[idx] = 0.01 * ((2 * ix + iz) % 3 - 1);
                ++idx;
            }
    frame.simbox = cpu::Box::cubic(box_len);
    return frame;
}

int main(int argc, char** argv) {
    int n_side = (argc > 1) ? std::atoi(argv[1]) : 16;
    int steps = (argc > 2) ? std::atoi(argv[2]) : 50;
    int warm = (argc > 3) ? std::atoi(argv[3]) : 10;
    if (n_side < 4)
        n_side = 4;
    if (steps < 5)
        steps = 5;

    auto frame = make_lj(n_side, 0.8);
    const double cutoff = 2.5;
    const double skin = 0.3; // production-like skin (not hostile rebuild spam)

    cuda::md::MD sim;
    sim.set_neighbor(cuda::NeighborList{cutoff, skin});
    sim.add_potential(cuda::PairLJCut{1.0, 1.0, cutoff});
    sim.set_integrator(cuda::md::VelocityVerlet{0.005});
    sim.build();
    sim.run(frame, warm);
    auto t0 = std::chrono::steady_clock::now();
    sim.run(frame, steps);
    auto t1 = std::chrono::steady_clock::now();
    double ms =
        std::chrono::duration<double, std::milli>(t1 - t0).count() / steps;
    std::printf(
        "path=openmm_cluster n_side=%d N=%d steps=%d  ms/step=%.4f  energy=%.6f\n",
        n_side, n_side * n_side * n_side, steps, ms,
        static_cast<double>(frame.meta.get(kEnergy)));
    return 0;
}
