/// @file examples/02_rhf.cpp
/// @brief Tutorial example for restricted Hartree-Fock on molecular hydrogen.
///
/// End-to-end RHF workflow on CPU:
///   1. Build H2 geometry in cpu::Frame
///   2. Compute molecular integrals
///   3. Run RHF
///   4. Print results and compare with textbook
///
/// The matching CUDA translation of the same tutorial is
/// examples/02_rhf.cu.
///
/// Build:
///   cmake -B build -DATV_BUILD_RHF=ON -DATV_BUILD_EXAMPLES=ON
///   cmake --build build --target 02_rhf
///
/// Run:
///   ./build/examples/02_rhf
///
/// Expected total energy: about -1.1168 Hartree (Szabo & Ostlund, Table 3.7)

#include <cmath>
#include <iomanip>
#include <iostream>

#include "atomiverse/cpu/scf/integrals.h"
#include "atomiverse/cpu/scf/scf_solver.h"

using namespace ATV;

int main() {
    // ── 1. Build H2 geometry ──────────────────────────────────────────
    //
    // cpu::Frame is the user-facing data container. Positions in Angstrom,
    // type_id = atomic number Z. The library converts to Bohr internally.
    cpu::Frame frame(2);
    frame.span(ATV::kAtom_x)[0] = double(0);
    frame.span(ATV::kAtom_y)[0] = double(0);
    frame.span(ATV::kAtom_z)[0] = double(0);
    frame.span(ATV::kAtom_x)[1] = double(0.74);
    frame.span(ATV::kAtom_y)[1] = double(0);
    frame.span(ATV::kAtom_z)[1] = double(0);
    frame.span(ATV::kAtom_mass)[0] = double(1.008);
    frame.span(ATV::kAtom_mass)[1] = double(1.008);
    frame.span(ATV::kAtom_type)[0] = 1;
    frame.span(ATV::kAtom_type)[1] = 1; // H: Z=1
    frame.simbox = ATV::cpu::Box::cubic(100.0);

    std::cout << "=== H2 / STO-3G / RHF (CPU) ===\n"
              << "Geometry: H-H at 0.74 Angstrom\n\n";

    // ── 2. Compute molecular integrals ────────────────────────────────
    //
    // build_integrals() evaluates, with native C++ integral kernels:
    //   - Overlap S, kinetic T, nuclear attraction V -> Hcore = T + V
    //   - Two-electron repulsion integrals ERI (chemist's notation)
    //   - Nuclear repulsion V_nuc = Z_A * Z_B / R_AB
    //
    // compute_eri=true forces the full O(n^4) two-electron tensor onto
    // integrals.ERI. build_integrals() omits it by default because
    // SCFSolver() builds the Fock matrix directly without the dense tensor;
    // we request it here only so this tutorial can print individual ERI
    // values below. Reading integrals.ERI without this flag is an
    // out-of-bounds access on an empty tensor.
    //
    // All matrices are xtensor arrays. No CUDA involved.
    auto integrals =
        cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G,
                                  /*charge=*/0, /*compute_eri=*/true);

    std::cout << "Basis functions: " << integrals.n_basis << "\n"
              << "Electrons:       " << integrals.n_electrons << "\n"
              << "Occupied MOs:    " << integrals.n_occ()
              << " (RHF closed-shell)\n"
              << std::fixed << std::setprecision(10)
              << "V_nuc:           " << static_cast<double>(integrals.V_nuc)
              << " Hartree\n\n";

    std::cout << "Overlap S:\n" << std::fixed << std::setprecision(6);
    for (int i = 0; i < integrals.n_basis; ++i) {
        for (int j = 0; j < integrals.n_basis; ++j)
            std::cout << "  " << std::setw(10)
                      << static_cast<double>(integrals.S(i, j));
        std::cout << "\n";
    }

    std::cout << "\nCore Hamiltonian Hcore = T + V:\n";
    for (int i = 0; i < integrals.n_basis; ++i) {
        for (int j = 0; j < integrals.n_basis; ++j)
            std::cout << "  " << std::setw(10)
                      << static_cast<double>(integrals.Hcore(i, j));
        std::cout << "\n";
    }

    std::cout << "\nKey ERI (chemist's notation):\n"
              << "  (11|11) = "
              << static_cast<double>(integrals.ERI(0, 0, 0, 0)) << "\n"
              << "  (11|22) = "
              << static_cast<double>(integrals.ERI(0, 0, 1, 1)) << "\n"
              << "  (12|12) = "
              << static_cast<double>(integrals.ERI(0, 1, 0, 1)) << "\n";

    // ── 3. Run RHF ────────────────────────────────────────────────────
    //
    // SCFSolver implements the full RHF algorithm on CPU:
    //   X = S^{-1/2}, then iterate:
    //     F = Hcore + G(P)  ->  F' = X^T F X  ->  diag  ->  C = X C'
    //     P = 2 C_occ C_occ^T  ->  E = 0.5 Tr[(Hcore+F)P] + V_nuc
    //
    // All linear algebra via xtensor-blas (LAPACK).
    double tol =
        (sizeof(double) == sizeof(float)) ? double(1e-5) : double(1e-10);
    auto config = cpu::scf::rhf_config(integrals, tol);
    config.max_iter = 100;
    cpu::scf::SCFSolver solver{config};
    solver.build(integrals);
    const auto& result = solver.run(integrals);

    // ── 4. Results ────────────────────────────────────────────────────
    std::cout << "\n=== RHF Result ===\n"
              << "Converged:        " << (result.converged ? "yes" : "no")
              << " (" << result.iterations << " iterations)\n"
              << std::fixed << std::setprecision(10) << "E_electronic:     "
              << static_cast<double>(result.energy - result.nuclear_repulsion)
              << " Hartree\n"
              << "V_nuclear:        " << static_cast<double>(integrals.V_nuc)
              << " Hartree\n"
              << "E_total:          " << static_cast<double>(result.energy)
              << " Hartree\n\n";

    std::cout << "Orbital energies:\n" << std::setprecision(6);
    for (int i = 0; i < integrals.n_basis; ++i)
        std::cout << "  epsilon[" << i << "] = " << std::setw(12)
                  << static_cast<double>(result.eigenvalues_alpha(i))
                  << " Hartree" << (i < integrals.n_occ() ? "  (occupied)" : "")
                  << "\n";

    // Density check: Tr(P * S) = n_electrons
    double tr_PS = xt::sum(result.P_alpha * integrals.S)();
    std::cout << "\nDensity check: Tr(P*S) = " << static_cast<double>(tr_PS)
              << " (expected: " << integrals.n_electrons << ")\n";

    std::cout << "\n--- Textbook comparison (Szabo & Ostlund) ---\n"
              << "  Reference: -1.1168 Hartree\n"
              << std::setprecision(4)
              << "  This:      " << static_cast<double>(result.energy)
              << " Hartree\n";

    std::cout << "\n=== Done ===\n";
    return 0;
}
