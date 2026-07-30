/// @file tools/bench_cuda_rhf.cu
/// @brief Micro-bench for cuda::scf::RHF: wall time of energy / force paths.
///
/// Build (on GH200 node, after configuring ATV with CUDA+HF):
///   g++/nvcc not needed standalone — link against atomiverse like regressions.
/// Prefer: add as a one-off binary or run via timed gtest.
///
/// Standalone usage when linked:
///   ./bench_cuda_rhf [n_repeat=20]
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>

#include <cuda_runtime.h>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/core/units.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/scf/scf_enums.h"
#include "atomiverse/cuda/scf/rhf.cuh"

using namespace ATV;
using clock_t_ = std::chrono::steady_clock;

static cpu::Frame make_h2() {
    cpu::Frame f(2);
    f.span(kAtom_x)[0] = 0.0;
    f.span(kAtom_y)[0] = 0.0;
    f.span(kAtom_z)[0] = 0.0;
    f.span(kAtom_x)[1] = 0.74;
    f.span(kAtom_y)[1] = 0.0;
    f.span(kAtom_z)[1] = 0.0;
    f.span(kAtom_type)[0] = 1;
    f.span(kAtom_type)[1] = 1;
    f.span(kAtom_mass)[0] = 1.008;
    f.span(kAtom_mass)[1] = 1.008;
    f.box = cpu::Box::cubic(100.0);
    return f;
}

static cpu::Frame make_h2o() {
    // Standard water geometry (Å), STO-3G n_AO=7.
    cpu::Frame f(3);
    f.span(kAtom_x)[0] = 0.0;
    f.span(kAtom_y)[0] = 0.0;
    f.span(kAtom_z)[0] = 0.1173;
    f.span(kAtom_x)[1] = 0.0;
    f.span(kAtom_y)[1] = 0.7572;
    f.span(kAtom_z)[1] = -0.4692;
    f.span(kAtom_x)[2] = 0.0;
    f.span(kAtom_y)[2] = -0.7572;
    f.span(kAtom_z)[2] = -0.4692;
    f.span(kAtom_type)[0] = 8;
    f.span(kAtom_type)[1] = 1;
    f.span(kAtom_type)[2] = 1;
    f.span(kAtom_mass)[0] = 15.999;
    f.span(kAtom_mass)[1] = 1.008;
    f.span(kAtom_mass)[2] = 1.008;
    f.box = cpu::Box::cubic(100.0);
    return f;
}

static double ms_since(clock_t_::time_point t0) {
    return std::chrono::duration<double, std::milli>(clock_t_::now() - t0)
        .count();
}

template <typename MakeFrame>
static void bench(const char* name, MakeFrame make, int n_rep, bool forces) {
    // Warmup
    {
        auto f = make();
        cuda::scf::RHF rhf{cpu::scf::BasisSet::STO3G, 0, 1e-10, 100};
        if (!forces)
            rhf.with_gradient(false);
        rhf.run(f);
        cudaDeviceSynchronize();
    }

    double t_sum = 0, t_min = 1e300, t_max = 0;
    double e_ha = 0;
    int iters = 0;
    for (int i = 0; i < n_rep; ++i) {
        auto f = make();
        cuda::scf::RHF rhf{cpu::scf::BasisSet::STO3G, 0, 1e-10, 100};
        if (!forces)
            rhf.with_gradient(false);
        cudaDeviceSynchronize();
        auto t0 = clock_t_::now();
        rhf.run(f);
        cudaDeviceSynchronize();
        double dt = ms_since(t0);
        t_sum += dt;
        t_min = std::min(t_min, dt);
        t_max = std::max(t_max, dt);
        e_ha = rhf.energy();
        iters = rhf.iterations();
    }
    std::printf("%-28s  n=%2d  mean=%8.3f ms  min=%8.3f  max=%8.3f  "
                "iters=%d  E=%.6f Ha\n",
                name, n_rep, t_sum / n_rep, t_min, t_max, iters, e_ha);
}

int main(int argc, char** argv) {
    int n_rep = 20;
    if (argc > 1)
        n_rep = std::atoi(argv[1]);
    if (n_rep < 1)
        n_rep = 1;

    int dev = 0;
    cudaGetDevice(&dev);
    cudaDeviceProp prop{};
    cudaGetDeviceProperties(&prop, dev);
    std::printf("GPU: %s  CC %d.%d  n_rep=%d\n", prop.name, prop.major,
                prop.minor, n_rep);

    bench("H2 energy-only", make_h2, n_rep, false);
    bench("H2 energy+forces", make_h2, n_rep, true);
    bench("H2O energy-only", make_h2o, n_rep, false);
    // H2O forces still host path (p shells) — still time it.
    bench("H2O energy+forces", make_h2o, std::max(3, n_rep / 4), true);
    return 0;
}
