/// @file test_md.cu
/// @brief Integration test: NVE MD with LJ + output components.

#include <gtest/gtest.h>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <string>

#include "atomiverse/cuda/md/init.cuh"
#include "atomiverse/cuda/potentials/lj.cuh"
#include "atomiverse/cuda/md/md.cuh"
#include "atomiverse/cuda/neighbor_list.cuh"
#include "atomiverse/cuda/md/velocity_verlet.cuh"
#include "atomiverse/cuda/observer/rdf.cuh"
#include "atomiverse/cuda/recorder/recorder.cuh"
#include "atomiverse/cuda/recorder/thermo.cuh"
#include "atomiverse/core/units.h"

using namespace ATV;

namespace {

/// Create a simple cubic lattice of LJ particles.
cpu::Frame make_lj_system(int n_side, Real density) {
    int n = n_side * n_side * n_side;
    Real box_len = std::cbrt(Real(n) / density);
    Real spacing = box_len / Real(n_side);

    cpu::Frame frame(n);
    int idx = 0;
    for (int ix = 0; ix < n_side; ++ix)
        for (int iy = 0; iy < n_side; ++iy)
            for (int iz = 0; iz < n_side; ++iz) {
                frame.atoms.x[idx] = (Real(ix) + Real(0.5)) * spacing;
                frame.atoms.y[idx] = (Real(iy) + Real(0.5)) * spacing;
                frame.atoms.z[idx] = (Real(iz) + Real(0.5)) * spacing;
                frame.atoms.mass[idx] = Real(1.0);
                frame.atoms.type[idx] = 0;
                idx++;
            }
    frame.box = Mat3::cubic(static_cast<double>(box_len));

    // Small random velocities
    std::srand(42);
    for (int i = 0; i < n; ++i) {
        frame.atoms.vx[i] = Real(0.001) * (Real(std::rand()) / Real(RAND_MAX) - Real(0.5));
        frame.atoms.vy[i] = Real(0.001) * (Real(std::rand()) / Real(RAND_MAX) - Real(0.5));
        frame.atoms.vz[i] = Real(0.001) * (Real(std::rand()) / Real(RAND_MAX) - Real(0.5));
    }

    // Remove center-of-mass velocity
    Real mvx = 0, mvy = 0, mvz = 0;
    for (int i = 0; i < n; ++i) {
        mvx += frame.atoms.vx[i];
        mvy += frame.atoms.vy[i];
        mvz += frame.atoms.vz[i];
    }
    mvx /= Real(n); mvy /= Real(n); mvz /= Real(n);
    for (int i = 0; i < n; ++i) {
        frame.atoms.vx[i] -= mvx;
        frame.atoms.vy[i] -= mvy;
        frame.atoms.vz[i] -= mvz;
    }

    return frame;
}

} // anonymous namespace

TEST(MD, NVEEnergyConservation) {
    auto frame = make_lj_system(4, Real(0.8));  // 64 atoms
    Real cutoff = Real(3.0);

    // First run: get initial total energy
    cpu::Frame frame1 = frame;
    {
        md::MD sim;
        sim.add(md::VelocityVerlet{Real(0.001)});
        sim.on(md::Phase::force)
            .add(NeighborList{cutoff}, -10)
            .add(LJ{Real(1.0), Real(1.0), cutoff});
        sim.build();
        sim.run(frame1, 1);
    }
    Real e_initial = frame1.energy + frame1.kinetic;

    // Run 100 steps
    cpu::Frame frame2 = frame;
    {
        md::MD sim;
        sim.add(md::VelocityVerlet{Real(0.001)});
        sim.on(md::Phase::force)
            .add(NeighborList{cutoff}, -10)
            .add(LJ{Real(1.0), Real(1.0), cutoff});
        sim.build();
        sim.run(frame2, 100);
    }
    Real e_final = frame2.energy + frame2.kinetic;

    // KE should be positive after integration
    EXPECT_GT(frame2.kinetic, Real(0));

    // Total energy should be conserved within tolerance
    // (small dt = 0.001 fs, 100 steps, float precision)
    Real tol = std::abs(e_initial) * Real(0.01) + Real(0.01);
    EXPECT_NEAR(static_cast<double>(e_final),
                static_cast<double>(e_initial),
                static_cast<double>(tol));
}

TEST(MD, ThermoOutput) {
    auto frame = make_lj_system(4, Real(0.8));
    Real cutoff = Real(3.0);

    md::MD sim;
    sim.add(md::VelocityVerlet{Real(0.001)});
    sim.on(md::Phase::force)
        .add(NeighborList{cutoff}, -10)
        .add(LJ{Real(1.0), Real(1.0), cutoff});
    sim.on(md::Phase::end_of_step)
        .add(Thermo{5});
    sim.build();

    // Should print thermo output without crashing
    testing::internal::CaptureStdout();
    sim.run(frame, 10);
    std::string output = testing::internal::GetCapturedStdout();

    // Should contain header and data lines
    EXPECT_NE(output.find("Step"), std::string::npos);
    EXPECT_NE(output.find("Temp"), std::string::npos);
    EXPECT_NE(output.find("PotEng"), std::string::npos);
}

TEST(MD, RecorderXYZ) {
    auto frame = make_lj_system(3, Real(0.8));  // 27 atoms
    Real cutoff = Real(3.0);
    std::string filename = "test_traj.xyz";

    {
        md::MD sim;
        sim.add(md::VelocityVerlet{Real(0.001)});
        sim.on(md::Phase::force)
            .add(NeighborList{cutoff}, -10)
            .add(LJ{Real(1.0), Real(1.0), cutoff});
        sim.on(md::Phase::end_of_step)
            .add(Recorder{filename, 3});
        sim.build();
        sim.run(frame, 9);
    }

    // Check file exists and has content
    std::ifstream file(filename);
    ASSERT_TRUE(file.is_open());

    // Count frames: each frame starts with atom count line
    int frame_count = 0;
    std::string line;
    while (std::getline(file, line)) {
        // Atom count line
        if (line == "27") frame_count++;
    }
    // freq=3, steps 0..8 -> steps 0,3,6 = 3 frames
    EXPECT_EQ(frame_count, 3);

    // Cleanup
    std::remove(filename.c_str());
}

TEST(MD, RDFOutput) {
    auto frame = make_lj_system(4, Real(0.8));
    Real cutoff = Real(3.0);
    std::string filename = "test_rdf.dat";

    {
        md::MD sim;
        sim.add(md::VelocityVerlet{Real(0.001)});
        sim.on(md::Phase::force)
            .add(NeighborList{cutoff}, -10)
            .add(LJ{Real(1.0), Real(1.0), cutoff});
        sim.on(md::Phase::end_of_step)
            .add(RDF{cutoff, 50, filename});
        sim.build();
        sim.run(frame, 5);
    }

    // Check file exists
    std::ifstream file(filename);
    ASSERT_TRUE(file.is_open());

    // Read data: skip header, check we have bins
    std::string line;
    int data_lines = 0;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        data_lines++;
    }
    EXPECT_EQ(data_lines, 50);

    // Cleanup
    std::remove(filename.c_str());
}

TEST(MD, InitVelocity) {
    // Create system with zero velocities
    int n = 64;
    cpu::Frame frame(n);
    Real box_len = std::cbrt(Real(n) / Real(0.8));
    Real spacing = box_len / Real(4);
    int idx = 0;
    for (int ix = 0; ix < 4; ++ix)
        for (int iy = 0; iy < 4; ++iy)
            for (int iz = 0; iz < 4; ++iz) {
                frame.atoms.x[idx] = (Real(ix) + Real(0.5)) * spacing;
                frame.atoms.y[idx] = (Real(iy) + Real(0.5)) * spacing;
                frame.atoms.z[idx] = (Real(iz) + Real(0.5)) * spacing;
                frame.atoms.mass[idx] = Real(1.0);
                idx++;
            }
    frame.box = Mat3::cubic(static_cast<double>(box_len));

    // All velocities should be zero initially
    for (int i = 0; i < n; ++i) {
        ASSERT_EQ(frame.atoms.vx[i], Real(0));
    }

    Real target_temp = Real(300.0);
    Real cutoff = Real(3.0);

    // Run 1 step with InitVelocity — velocities should be set
    {
        md::MD sim;
        sim.add(md::InitVelocity{target_temp, 42});
        sim.add(md::VelocityVerlet{Real(0.001)});
        sim.on(md::Phase::force)
            .add(NeighborList{cutoff}, -10)
            .add(LJ{Real(1.0), Real(1.0), cutoff});
        sim.build();
        sim.run(frame, 1);
    }

    // Velocities should now be non-zero
    bool has_nonzero = false;
    for (int i = 0; i < n; ++i) {
        if (frame.atoms.vx[i] != Real(0) || frame.atoms.vy[i] != Real(0) || frame.atoms.vz[i] != Real(0))
            has_nonzero = true;
    }
    EXPECT_TRUE(has_nonzero);

    // Center-of-mass velocity should be near zero
    Real mvx = 0, mvy = 0, mvz = 0, M = 0;
    for (int i = 0; i < n; ++i) {
        mvx += frame.atoms.mass[i] * frame.atoms.vx[i];
        mvy += frame.atoms.mass[i] * frame.atoms.vy[i];
        mvz += frame.atoms.mass[i] * frame.atoms.vz[i];
        M += frame.atoms.mass[i];
    }
    // After 1 step of VV, COM velocity may drift slightly from forces,
    // but should still be small
    EXPECT_NEAR(static_cast<double>(mvx / M), 0.0, 0.01);
    EXPECT_NEAR(static_cast<double>(mvy / M), 0.0, 0.01);
    EXPECT_NEAR(static_cast<double>(mvz / M), 0.0, 0.01);

    // Kinetic energy should be positive (velocities were generated)
    EXPECT_GT(frame.kinetic, Real(0));

    // Temperature should be in the right ballpark (within factor 3 for 64 atoms)
    Int ndof = 3 * n - 3;
    Real temp = Real(2) * frame.kinetic / (Real(ndof) * kBoltzmann);
    EXPECT_GT(static_cast<double>(temp), static_cast<double>(target_temp) * 0.3);
    EXPECT_LT(static_cast<double>(temp), static_cast<double>(target_temp) * 3.0);
}
