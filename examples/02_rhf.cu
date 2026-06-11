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
#include <iomanip>
#include <iostream>
#include <vector>

#include "atomiverse/core/units.h"
#include "atomiverse/cpu/io/molrec.h"
#include "atomiverse/cpu/observer/mulliken.h"
#include "atomiverse/cpu/recorder/xyz.h"
#include "atomiverse/cpu/scf/integrals.h"
#include "atomiverse/cpu/scf/rhf_solve.h"
#include "atomiverse/cuda/scf/scf_evaluator.cuh"
#include "atomiverse/cuda/driver/single_point.cuh"

using namespace ATV;

int main() {
    // ── 1. Build H2 geometry ──────────────────────────────────────────
    cpu::Frame frame(2);
    frame.atoms.x[0] = double(0);
    frame.atoms.y[0] = double(0);
    frame.atoms.z[0] = double(0);
    frame.atoms.x[1] = double(0.74);
    frame.atoms.y[1] = double(0);
    frame.atoms.z[1] = double(0);
    frame.atoms.mass[0] = double(1.008);
    frame.atoms.mass[1] = double(1.008);
    frame.atoms.type[0] = 1;
    frame.atoms.type[1] = 1; // H: Z=1
    frame.box = Mat3::cubic(100.0);

    std::cout << "=== H2/STO-3G/RHF (CPU and CUDA) ===\n"
              << "Geometry: H-H at 0.74 Angstrom\n\n";

    // ── 2. Compute integrals on CPU (libcint) ─────────────────────────
    auto integrals =
        cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    std::cout << "Basis functions: " << integrals.n_basis << "\n"
              << "Electrons:       " << integrals.n_electrons << "\n"
              << "Occupied MOs:    " << integrals.n_occ() << "\n"
              << std::fixed << std::setprecision(10)
              << "V_nuc:           " << static_cast<double>(integrals.V_nuc)
              << " Hartree\n";

    // Print overlap and core Hamiltonian
    std::cout << "\nOverlap matrix S:\n" << std::setprecision(8);
    for (int i = 0; i < integrals.n_basis; ++i) {
        for (int j = 0; j < integrals.n_basis; ++j)
            std::cout << "  " << std::setw(12)
                      << static_cast<double>(integrals.S(i, j));
        std::cout << "\n";
    }
    std::cout << "\nCore Hamiltonian Hcore:\n";
    for (int i = 0; i < integrals.n_basis; ++i) {
        for (int j = 0; j < integrals.n_basis; ++j)
            std::cout << "  " << std::setw(12)
                      << static_cast<double>(integrals.Hcore(i, j));
        std::cout << "\n";
    }

    // ── 3. CPU reference RHF ──────────────────────────────────────────
    std::cout << "\n=== CPU RHF (xtensor-blas) ===\n";
    double tol =
        (sizeof(double) == sizeof(float)) ? double(1e-5) : double(1e-10);
    auto cpu_result = cpu::scf::rhf_solve(integrals, tol, 100);
    std::cout << "Converged:   " << (cpu_result.converged ? "yes" : "no")
              << " (" << cpu_result.iterations << " iterations"
              << ", residual=" << std::scientific << std::setprecision(2)
              << static_cast<double>(cpu_result.residual) << ")\n"
              << std::fixed << std::setprecision(10) << "E_elec:      "
              << static_cast<double>(cpu_result.electronic_energy)
              << " Hartree\n"
              << "E_total:     " << static_cast<double>(cpu_result.energy)
              << " Hartree\n"
              << "Eigenvalues:";
    for (int i = 0; i < integrals.n_basis; ++i)
        std::cout << " " << std::setprecision(6)
                  << static_cast<double>(cpu_result.eigenvalues(i));
    std::cout << "\n";

    // Molecular record: Mulliken analysis + extxyz output
    cpu::MolRec rec;
    rec.set_geometry(frame);
    rec.add_string("method", "RHF");
    rec.add_scalar("energy", cpu_result.energy);
    rec.add_scalar("iterations", double(cpu_result.iterations));
    rec.add_scalar("n_alpha", double(integrals.n_occ()));
    rec.add_scalar("n_beta", double(integrals.n_occ()));
    rec.add_scalar("spin_contamination", double(0));
    rec.add_field("charge", cpu::mulliken_charge(integrals, cpu_result.P));
    rec.add_field("spin_density",
                  std::vector<double>(frame.atoms.size(), double(0)));
    rec.commit_frame();
    cpu::print_summary(rec);
    cpu::XYZ("h2_rhf.xyz").write_ext(rec);
    std::cout << "  Wrote h2_rhf.xyz  (open in OVITO / Avogadro 2)\n";

    // ── 4. CUDA RHF — Evaluate + SCFEvaluator (inner-solver shape) ─────
    // SCF is not a driver category: the CUDA SCF cycle runs inside
    // cuda::scf::SCFSolver, owned by the SCFEvaluator component, driven
    // by the one-pass cuda::Evaluate (here via the single_point helper).
    std::cout << "\n=== CUDA RHF (cuda::Evaluate + scf::SCFEvaluator) ===\n";
    cuda::scf::SCFEvaluator cuda_scf{integrals, tol, /*max_iter=*/100};
    cuda::single_point(frame, cuda_scf); // frame.scalars.energy populated

    // ── 6. Compare ────────────────────────────────────────────────────
    std::cout << "\n=== Comparison ===\n"
              << std::fixed << std::setprecision(10)
              << "CPU  energy: " << static_cast<double>(cpu_result.energy)
              << " Hartree\n"
              << "CUDA energy: "
              << static_cast<double>(frame.scalars.energy * units::kEvInHartree)
              << " Hartree\n"
              << "Difference:  " << std::scientific << std::setprecision(2)
              << static_cast<double>(std::abs(
                     frame.scalars.energy * units::kEvInHartree -
                     cpu_result.energy))
              << "\n"
              << "Textbook:    -1.1167 Hartree (approx)\n";

    std::cout << "\n=== Done ===\n";
    return 0;
}
