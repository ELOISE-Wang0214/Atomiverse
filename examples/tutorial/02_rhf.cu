/// @file examples/tutorial/02_rhf.cu
/// @brief CUDA translation of the RHF tutorial example.
///
/// Demonstrates the full Atomiverse RHF workflow:
///   1. Build H2 geometry in cpu::Frame
///   2. Compute molecular integrals on the host
///   3. Run CPU reference RHF
///   4. Run CUDA RHF
///   5. Compare results
///
/// Build:
///   cmake -B build -DATV_BUILD_CUDA=ON -DATV_BUILD_RHF=ON -DATV_BUILD_EXAMPLES=ON
///   cmake --build build --target 02_rhf_cuda
///
/// Run:
///   ./build/examples/02_rhf_cuda
///
/// Expected total energy: -1.1167 Hartree (textbook H2/STO-3G/RHF)

#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>

#include "atomiverse/cpu/molrec.h"
#include "atomiverse/cpu/observer/mulliken.h"
#include "atomiverse/cpu/recorder/xyz.h"
#include "atomiverse/cpu/scf/integrals.h"
#include "atomiverse/cpu/scf/rhf_solve.h"
#include "atomiverse/cuda/scf/scf.cuh"
#include "atomiverse/cuda/scf/rhf_fock.cuh"
#include "atomiverse/cuda/scf/rhf_density.cuh"
#include "atomiverse/cuda/recorder/thermo.cuh"
#include "atomiverse/cuda/observer/mulliken.cuh"

using namespace ATV;

int main() {
    // ── 1. Build H2 geometry ──────────────────────────────────────────
    cpu::Frame frame(2);
    frame.atoms.x[0] = Real(0);     frame.atoms.y[0] = Real(0); frame.atoms.z[0] = Real(0);
    frame.atoms.x[1] = Real(0.74);  frame.atoms.y[1] = Real(0); frame.atoms.z[1] = Real(0);
    frame.atoms.mass[0] = Real(1.008); frame.atoms.mass[1] = Real(1.008);
    frame.atoms.type[0] = 1; frame.atoms.type[1] = 1;  // H: Z=1
    frame.box = Mat3::cubic(100.0);

    std::cout << "=== H2/STO-3G/RHF (CPU and CUDA) ===\n"
              << "Geometry: H-H at 0.74 Angstrom\n\n";

    // ── 2. Compute integrals on CPU (libcint) ─────────────────────────
    auto integrals = cpu::scf::build_integrals_sto3g(frame);
    std::cout << "Basis functions: " << integrals.n_basis << "\n"
              << "Electrons:       " << integrals.n_electrons << "\n"
              << "Occupied MOs:    " << integrals.n_occ() << "\n"
              << std::fixed << std::setprecision(10)
              << "V_nuc:           " << static_cast<double>(integrals.V_nuc) << " Hartree\n";

    // Print overlap and core Hamiltonian
    std::cout << "\nOverlap matrix S:\n" << std::setprecision(8);
    for (Int i = 0; i < integrals.n_basis; ++i) {
        for (Int j = 0; j < integrals.n_basis; ++j)
            std::cout << "  " << std::setw(12) << static_cast<double>(integrals.S(i, j));
        std::cout << "\n";
    }
    std::cout << "\nCore Hamiltonian Hcore:\n";
    for (Int i = 0; i < integrals.n_basis; ++i) {
        for (Int j = 0; j < integrals.n_basis; ++j)
            std::cout << "  " << std::setw(12) << static_cast<double>(integrals.Hcore(i, j));
        std::cout << "\n";
    }

    // ── 3. CPU reference RHF ──────────────────────────────────────────
    std::cout << "\n=== CPU RHF (xtensor-blas) ===\n";
    Real tol = (sizeof(Real) == sizeof(float)) ? Real(1e-5) : Real(1e-10);
    auto cpu_result = cpu::scf::rhf_solve(integrals, tol, 100);
    std::cout << "Converged:   " << (cpu_result.converged ? "yes" : "no")
              << " (" << cpu_result.iterations << " iterations"
              << ", residual=" << std::scientific << std::setprecision(2)
              << static_cast<double>(cpu_result.residual) << ")\n"
              << std::fixed << std::setprecision(10)
              << "E_elec:      " << static_cast<double>(cpu_result.electronic_energy) << " Hartree\n"
              << "E_total:     " << static_cast<double>(cpu_result.energy) << " Hartree\n"
              << "Eigenvalues:";
    for (Int i = 0; i < integrals.n_basis; ++i)
        std::cout << " " << std::setprecision(6) << static_cast<double>(cpu_result.eigenvalues(i));
    std::cout << "\n";

    // Molecular record: Mulliken analysis + extxyz output
    cpu::MolRec rec;
    rec.set_geometry(frame);
    rec.add_string("method", "RHF");
    rec.add_scalar("energy",     cpu_result.energy);
    rec.add_scalar("iterations", Real(cpu_result.iterations));
    rec.add_scalar("n_alpha",    Real(integrals.n_occ()));
    rec.add_scalar("n_beta",     Real(integrals.n_occ()));
    rec.add_scalar("spin_contamination", Real(0));
    rec.add_field("charge",       cpu::mulliken_charge(frame, cpu_result.P, integrals.S));
    rec.add_field("spin_density", std::vector<Real>(frame.atoms.size(), Real(0)));
    rec.commit_frame();
    cpu::print_summary(rec);
    cpu::XYZ("h2_rhf.xyz").write_ext(rec);
    std::cout << "  Wrote h2_rhf.xyz  (open in OVITO / Avogadro 2)\n";

    // ── 4. CUDA RHF through Driver ────────────────────────────────────
    std::cout << "\n=== CUDA RHF (scf::SCF Driver) ===\n";
    MullikenObs mulliken_obs(integrals);
    {
        scf::SCF solver;
        solver.add(scf::RHFFock{integrals});
        solver.add(scf::EigenSolver{});
        solver.add(scf::RHFDensity{integrals.n_occ()});
        solver.add(scf::ConvergenceCheck{tol, false});
        solver.on(scf::Phase::end_of_step).add(Thermo{1});
        solver.on(scf::Phase::end_of_step).add(mulliken_obs);  // async D2H -> CPU Mulliken

        solver.enable_timing();
        solver.build();
        solver.run(frame, 100);

        std::cout << "\n" << solver.summary();
    }

    // ── 5. CUDA Mulliken result (observer -> recorder separation) ──────
    std::cout << "\n=== CUDA Mulliken (from MullikenObs) ===\n";
    {
        cpu::MolRec cuda_rec;
        cuda_rec.set_geometry(frame);
        cuda_rec.add_string("method", "RHF");
        cuda_rec.add_scalar("energy",             mulliken_obs.energy());
        cuda_rec.add_scalar("iterations",         Real(mulliken_obs.iterations()));
        cuda_rec.add_scalar("n_alpha",            Real(mulliken_obs.n_occ()));
        cuda_rec.add_scalar("n_beta",             Real(mulliken_obs.n_occ()));
        cuda_rec.add_scalar("spin_contamination", Real(0));
        cuda_rec.add_field("charge",              mulliken_obs.charges());
        cuda_rec.add_field("spin_density",        std::vector<Real>(frame.atoms.size(), Real(0)));
        cuda_rec.commit_frame();
        cpu::print_summary(cuda_rec);
        cpu::XYZ("h2_rhf_cuda.xyz").write_ext(cuda_rec);
    }
    std::cout << "  Wrote h2_rhf_cuda.xyz\n";

    // ── 6. Compare ────────────────────────────────────────────────────
    std::cout << "\n=== Comparison ===\n"
              << std::fixed << std::setprecision(10)
              << "CPU  energy: " << static_cast<double>(cpu_result.energy) << " Hartree\n"
              << "CUDA energy: " << static_cast<double>(frame.energy) << " Hartree\n"
              << "Difference:  " << std::scientific << std::setprecision(2)
              << static_cast<double>(std::abs(frame.energy - cpu_result.energy)) << "\n"
              << "Textbook:    -1.1167 Hartree (approx)\n";

    std::cout << "\n=== Done ===\n";
    return 0;
}
