/// @file test_scf.cu
/// @brief Integration tests for the SCF driver and CPU/CUDA parity.

#include <gtest/gtest.h>
#include <cmath>

#include "atomiverse/cuda/scf/scf.cuh"
#include "atomiverse/cuda/recorder/thermo.cuh"
#include "atomiverse/cpu/scf/solve.h"

using namespace ATV;

namespace {

/// Create a minimal cpu::Frame with n atoms (used as matrix dimension).
cpu::Frame make_scf_system(Int n) {
    cpu::Frame frame(n);
    for (Int i = 0; i < n; ++i) {
        frame.x[i] = Real(i);
        frame.y[i] = Real(0);
        frame.z[i] = Real(0);
        frame.mass[i] = Real(1);
        frame.type_id[i] = 0;
    }
    frame.box = Mat3::cubic(100.0);
    return frame;
}

} // anonymous namespace

// ── CUDA Driver tests ─────────────────────────────────────────────────

TEST(SCF, CUDAConvergence) {
    Int n = 8;
    Int n_occ = 3;
    auto frame = make_scf_system(n);

    scf::SCF solver;
    solver.add(scf::Hamiltonian{n_occ, Real(0.1)});
    solver.add(scf::EigenSolver{});
    solver.add(scf::DensityBuilder{});
    solver.add(scf::ConvergenceCheck{Real(1e-6)});
    solver.build();
    solver.run(frame, 200);

    EXPECT_TRUE(std::isfinite(static_cast<double>(frame.energy)));
    EXPECT_GT(frame.energy, Real(0));
}

TEST(SCF, CUDAReproducibility) {
    Int n = 8;
    auto frame1 = make_scf_system(n);
    auto frame2 = make_scf_system(n);

    for (auto* f : {&frame1, &frame2}) {
        scf::SCF solver;
        solver.add(scf::Hamiltonian{2, Real(0.1)});
        solver.add(scf::EigenSolver{});
        solver.add(scf::DensityBuilder{});
        solver.add(scf::ConvergenceCheck{Real(1e-7)});
        solver.build();
        solver.run(*f, 200);
    }

    EXPECT_NEAR(static_cast<double>(frame1.energy),
                static_cast<double>(frame2.energy), 1e-10);
}

TEST(SCF, CUDAWithThermo) {
    Int n = 8;
    auto frame = make_scf_system(n);

    scf::SCF solver;
    solver.add(scf::Hamiltonian{2, Real(0.1)});
    solver.add(scf::EigenSolver{});
    solver.add(scf::DensityBuilder{});
    solver.add(scf::ConvergenceCheck{Real(1e-6)});
    solver.on(scf::Phase::end_of_step).add(Thermo{10});
    solver.build();
    solver.run(frame, 50);

    EXPECT_TRUE(std::isfinite(static_cast<double>(frame.energy)));
}

TEST(SCF, CUDAAlphaEffect) {
    Int n = 8;
    auto frame_a0 = make_scf_system(n);
    auto frame_a1 = make_scf_system(n);

    {
        scf::SCF solver;
        solver.add(scf::Hamiltonian{2, Real(0)});
        solver.add(scf::EigenSolver{});
        solver.add(scf::DensityBuilder{});
        solver.add(scf::ConvergenceCheck{Real(1e-6)});
        solver.build();
        solver.run(frame_a0, 50);
    }
    {
        scf::SCF solver;
        solver.add(scf::Hamiltonian{2, Real(0.5)});
        solver.add(scf::EigenSolver{});
        solver.add(scf::DensityBuilder{});
        solver.add(scf::ConvergenceCheck{Real(1e-6)});
        solver.build();
        solver.run(frame_a1, 200);
    }

    EXPECT_NE(frame_a0.energy, frame_a1.energy);
    EXPECT_TRUE(std::isfinite(static_cast<double>(frame_a0.energy)));
    EXPECT_TRUE(std::isfinite(static_cast<double>(frame_a1.energy)));
}

// ── CPU reference tests ───────────────────────────────────────────────

TEST(SCF, CPUConvergence) {
    Real cpu_tol = (sizeof(Real) == sizeof(float)) ? Real(1e-5) : Real(1e-8);
    auto result = cpu::scf::solve(8, 3, Real(0.1), cpu_tol, 200);
    EXPECT_TRUE(result.converged);
    EXPECT_TRUE(std::isfinite(static_cast<double>(result.energy)));
    EXPECT_GT(result.energy, Real(0));
    EXPECT_GT(result.iterations, 1);
}

TEST(SCF, CPUAlphaZero) {
    // alpha=0: no self-consistency, converges in 2 iterations
    Real cpu_tol = (sizeof(Real) == sizeof(float)) ? Real(1e-5) : Real(1e-10);
    auto result = cpu::scf::solve(8, 3, Real(0), cpu_tol, 50);
    EXPECT_TRUE(result.converged);
    EXPECT_LE(result.iterations, 3);  // 2 for double, possibly 3 for float
}

TEST(SCF, CPUReproducibility) {
    Real cpu_tol = (sizeof(Real) == sizeof(float)) ? Real(1e-5) : Real(1e-8);
    auto r1 = cpu::scf::solve(16, 4, Real(0.1), cpu_tol, 100);
    auto r2 = cpu::scf::solve(16, 4, Real(0.1), cpu_tol, 100);
    EXPECT_DOUBLE_EQ(static_cast<double>(r1.energy),
                     static_cast<double>(r2.energy));
    EXPECT_EQ(r1.iterations, r2.iterations);
}

// ── CPU/CUDA parity ───────────────────────────────────────────────────

TEST(SCF, CPUCUDAParity) {
    Int n = 8;
    Int n_occ = 2;
    Real alpha = Real(0.1);
    Real tol = Real(1e-5);

    // CPU reference
    auto cpu_result = cpu::scf::solve(n, n_occ, alpha, tol, 200);

    // CUDA via Driver
    auto frame = make_scf_system(n);
    {
        scf::SCF solver;
        solver.add(scf::Hamiltonian{n_occ, alpha});
        solver.add(scf::EigenSolver{});
        solver.add(scf::DensityBuilder{});
        solver.add(scf::ConvergenceCheck{tol});
        solver.build();
        solver.run(frame, 200);
    }

    // Energy must agree within float tolerance
    double cpu_e = static_cast<double>(cpu_result.energy);
    double cuda_e = static_cast<double>(frame.energy);
    double rel_tol = (sizeof(Real) == sizeof(float)) ? 1e-4 : 1e-8;

    EXPECT_NEAR(cpu_e, cuda_e, std::abs(cpu_e) * rel_tol + 1e-6)
        << "CPU energy: " << cpu_e << ", CUDA energy: " << cuda_e;
}

TEST(SCF, CPUCUDAParityN16) {
    Int n = 16;
    Int n_occ = 4;
    Real alpha = Real(0.05);
    Real tol = Real(1e-5);

    auto cpu_result = cpu::scf::solve(n, n_occ, alpha, tol, 200);

    auto frame = make_scf_system(n);
    {
        scf::SCF solver;
        solver.add(scf::Hamiltonian{n_occ, alpha});
        solver.add(scf::EigenSolver{});
        solver.add(scf::DensityBuilder{});
        solver.add(scf::ConvergenceCheck{tol});
        solver.build();
        solver.run(frame, 200);
    }

    double cpu_e = static_cast<double>(cpu_result.energy);
    double cuda_e = static_cast<double>(frame.energy);
    double rel_tol = (sizeof(Real) == sizeof(float)) ? 1e-4 : 1e-8;

    EXPECT_NEAR(cpu_e, cuda_e, std::abs(cpu_e) * rel_tol + 1e-6)
        << "CPU energy: " << cpu_e << ", CUDA energy: " << cuda_e;
}
