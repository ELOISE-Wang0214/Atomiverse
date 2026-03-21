/// @file examples/scf.cu
/// @brief SCF self-consistent field example: iterative dense-matrix eigenproblem.
///
/// Demonstrates the Atomiverse SCF workflow through the Driver API:
///   1. Construct a system (n atoms = n basis functions for this demo)
///   2. Run iterative SCF: H(P) = H0 + alpha*P, diagonalize, rebuild density
///   3. Observe convergence via Thermo output
///   4. Compare different coupling strengths
///
/// The Hamiltonian H0 is a deterministic diagonally-dominant symmetric matrix.
/// The SCF cycle finds the self-consistent density matrix P and total energy.
///
/// Usage:
///   ./scf
///
/// Output:
///   Per-iteration thermo: Step, Temp(unused), PotEng(=SCF energy), KinEng, TotEng

#include <cmath>
#include <iostream>
#include <iomanip>

#include "atomiverse/cuda/scf/scf.cuh"
#include "atomiverse/cpu/scf/solve.h"
#include "atomiverse/cuda/recorder/thermo.cuh"

using namespace ATV;

int main() {
    // ── System parameters ─────────────────────────────────────────────
    Int n      = 32;         // matrix dimension (= n_atoms for frame)
    Int n_occ  = 8;          // occupied orbitals
    Real alpha = Real(0.1);  // density coupling strength
    Real tol   = Real(1e-5); // convergence tolerance (float-friendly)
    int max_iter = 50;

    std::cout << "SCF Demo: n=" << n << ", n_occ=" << n_occ
              << ", alpha=" << std::fixed << std::setprecision(2) << static_cast<double>(alpha)
              << ", tol=" << std::scientific << std::setprecision(1) << static_cast<double>(tol)
              << "\n";

    // ── Build frame ───────────────────────────────────────────────────
    // Frame dimension n is used as the matrix size.
    // Positions are just indices -- not physically meaningful here.
    cpu::Frame frame(n);
    for (Int i = 0; i < n; ++i) {
        frame.x[i] = Real(i);
        frame.y[i] = Real(0);
        frame.z[i] = Real(0);
        frame.mass[i] = Real(1);
        frame.type_id[i] = 0;
    }
    frame.box = Mat3::cubic(100.0);

    // ── Run SCF ───────────────────────────────────────────────────────
    std::cout << "\n=== SCF (max " << max_iter << " iterations) ===\n";
    {
        scf::SCF solver;
        solver.add(scf::Hamiltonian{n_occ, alpha});
        solver.add(scf::EigenSolver{});
        solver.add(scf::DensityBuilder{});
        solver.add(scf::ConvergenceCheck{tol});
        solver.on(scf::Phase::end_of_step).add(Thermo{5});

        solver.enable_timing();
        solver.build();
        solver.run(frame, max_iter);

        std::cout << "\n" << solver.summary();
    }

    std::cout << "\n=== Result ===\n"
              << std::fixed << std::setprecision(10)
              << "SCF energy: " << static_cast<double>(frame.energy) << "\n";

    // ── CPU reference (xtensor-blas) ────────────────────────────────────
    std::cout << "\n=== CPU reference (xtensor-blas) ===\n";
    {
        auto cpu_result = cpu::scf::solve(n, n_occ, alpha, tol, max_iter);
        std::cout << std::fixed << std::setprecision(10)
                  << "CPU energy:    " << static_cast<double>(cpu_result.energy)
                  << "  (" << (cpu_result.converged ? "converged" : "not converged")
                  << ", " << cpu_result.iterations << " iterations"
                  << ", residual=" << std::scientific << std::setprecision(2) << static_cast<double>(cpu_result.residual)
                  << ")\n"
                  << std::fixed << std::setprecision(10)
                  << "CUDA energy:   " << static_cast<double>(frame.energy) << "\n"
                  << "Difference:    " << std::scientific << std::setprecision(2)
                  << static_cast<double>(std::abs(frame.energy - cpu_result.energy)) << "\n";
    }

    std::cout << "\n=== Done ===\n";
    return 0;
}
