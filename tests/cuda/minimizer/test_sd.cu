/// @file test_minimizer.cu
/// @brief Integration test: steepest descent minimization of LJ system.

#include <gtest/gtest.h>
#include <cmath>

#include "atomiverse/cuda/potentials/lj.cuh"
#include "atomiverse/cuda/neighbor_list.cuh"
#include "atomiverse/cuda/minimizer/sd.cuh"
#include "atomiverse/cuda/recorder/thermo.cuh"

using namespace ATV;

namespace {

/// Create a simple cubic lattice slightly perturbed from equilibrium.
cpu::Frame make_displaced_system(int n_side, double density, double displacement) {
    int n = n_side * n_side * n_side;
    double box_len = std::cbrt(double(n) / density);
    double spacing = box_len / double(n_side);

    cpu::Frame frame(n);
    std::srand(123);
    int idx = 0;
    for (int ix = 0; ix < n_side; ++ix)
        for (int iy = 0; iy < n_side; ++iy)
            for (int iz = 0; iz < n_side; ++iz) {
                frame.atoms.x[idx] = (double(ix) + double(0.5)) * spacing +
                    displacement * (double(std::rand()) / double(RAND_MAX) - double(0.5));
                frame.atoms.y[idx] = (double(iy) + double(0.5)) * spacing +
                    displacement * (double(std::rand()) / double(RAND_MAX) - double(0.5));
                frame.atoms.z[idx] = (double(iz) + double(0.5)) * spacing +
                    displacement * (double(std::rand()) / double(RAND_MAX) - double(0.5));
                frame.atoms.mass[idx] = double(1.0);
                frame.atoms.type[idx] = 0;
                idx++;
            }
    frame.box = Mat3::cubic(static_cast<double>(box_len));
    return frame;
}

} // anonymous namespace

TEST(Minimizer, SDEnergyDecreases) {
    auto frame = make_displaced_system(4, double(0.8), double(0.3));
    double cutoff = double(3.0);

    // Get initial energy with a single force evaluation
    cpu::Frame frame_init = frame;
    {
        minimizer::SD min;
        min.on(minimizer::Phase::force)
            .add(NeighborList{cutoff}, -10)
            .add(LJ{double(1.0), double(1.0), cutoff});
        min.build();
        min.run(frame_init, 1);  // 1 step: force eval only (step 0 has zero forces → step update is no-op)
    }
    double e_after_1 = frame_init.energy;

    // Run more steps of minimization
    cpu::Frame frame_min = frame;
    {
        minimizer::SD min;
        min.add(minimizer::SteepestDescentStep{double(0.01)});
        min.on(minimizer::Phase::force)
            .add(NeighborList{cutoff}, -10)
            .add(LJ{double(1.0), double(1.0), cutoff});
        min.build();
        min.run(frame_min, 50);
    }
    double e_final = frame_min.energy;

    // Energy should have decreased (system relaxed)
    EXPECT_LT(static_cast<double>(e_final), static_cast<double>(e_after_1));
}

TEST(Minimizer, SDWithThermo) {
    auto frame = make_displaced_system(3, double(0.8), double(0.2));
    double cutoff = double(3.0);

    minimizer::SD min;
    min.add(minimizer::SteepestDescentStep{double(0.01)});
    min.on(minimizer::Phase::force)
        .add(NeighborList{cutoff}, -10)
        .add(LJ{double(1.0), double(1.0), cutoff});
    min.on(minimizer::Phase::end_of_step)
        .add(Thermo{10});
    min.build();

    // Should print thermo without crashing
    testing::internal::CaptureStdout();
    min.run(frame, 30);
    std::string output = testing::internal::GetCapturedStdout();

    EXPECT_NE(output.find("Step"), std::string::npos);
    EXPECT_NE(output.find("PotEng"), std::string::npos);

    // Energy should be finite
    EXPECT_TRUE(std::isfinite(static_cast<double>(frame.energy)));
}
