/// @file examples/14_fci.cpp
/// @brief Tutorial example: exact correlation via tiny full CI on H2/STO-3G.
///
/// Full configuration interaction (full CI) gives the EXACT correlation energy
/// within a finite basis. For a tiny two-electron minimal-basis system it is a
/// four-determinant dense-eigenvalue problem, cheap enough to serve as a golden
/// reference for MP2 and future coupled-cluster methods.
///
/// Workflow (all energies in [Hartree]):
///   1. Build H2 geometry in cpu::Frame
///   2. Compute molecular integrals with the two-electron ERI tensor
///   3. Run RHF to get the mean-field reference (E_RHF)
///   4. Run full CI to get the exact energy (E_FCI) and correlation (E_corr)
///
/// Build:
///   cmake -B build -DATV_BUILD_RHF=ON -DATV_BUILD_EXAMPLES=ON
///   cmake --build build --target 14_fci
///
/// Run:
///   ./build/examples/14_fci

#ifdef ATV_HAS_RHF

#include <cmath>
#include <iomanip>
#include <iostream>

#include "atomiverse/core/mat3.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/scf/fci.h"
#include "atomiverse/cpu/scf/integrals.h"
#include "atomiverse/cpu/scf/scf_solver.h"

using namespace ATV;

int main() {
    // 1. H2 geometry (Angstrom); type = atomic number Z.
    cpu::Frame frame(2);
    frame.span(ATV::kAtom_x)[0] = 0.0;
    frame.span(ATV::kAtom_x)[1] = 0.0;
    frame.span(ATV::kAtom_y)[0] = 0.0;
    frame.span(ATV::kAtom_y)[1] = 0.0;
    frame.span(ATV::kAtom_z)[0] = 0.0;
    frame.span(ATV::kAtom_z)[1] = 0.7414;
    frame.span(ATV::kAtom_type)[0] = 1;
    frame.span(ATV::kAtom_type)[1] = 1;
    frame.span(ATV::kAtom_mass)[0] = 1.008;
    frame.span(ATV::kAtom_mass)[1] = 1.008;
    frame.simbox = ATV::cpu::Box::cubic(50.0);

    std::cout << "=== H2 / STO-3G / full CI (CPU) ===\n"
              << "Geometry: H-H at 0.7414 Angstrom\n\n";

    // 2. Integrals -- the two-electron ERI tensor is needed by full CI.
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G,
                                          /*charge=*/0, /*compute_eri=*/true);

    // 3. RHF reference.
    cpu::scf::SCFSolver solver{cpu::scf::rhf_config(ints, 1e-12)};
    solver.build(ints);
    const auto& scf = solver.run(ints);
    if (!scf.converged) {
        std::cerr << "RHF did not converge\n";
        return 1;
    }

    // 4. Full CI: exact correlation energy in this basis.
    auto fci = cpu::scf::fci_energy(ints, scf);

    std::cout << std::fixed << std::setprecision(8) << "E_RHF:   " << scf.energy
              << " Hartree\n"
              << "E_FCI:   " << fci.e_fci << " Hartree\n"
              << "E_corr:  " << fci.e_corr << " Hartree (= E_FCI - E_RHF)\n";

    // Consistency: E_RHF + E_corr == E_FCI to print precision.
    const double recombined = scf.energy + fci.e_corr;
    if (std::abs(recombined - fci.e_fci) > 5e-9) {
        std::cerr << "inconsistent: E_RHF + E_corr != E_FCI\n";
        return 1;
    }

    std::cout << "\n=== Done ===\n";
    return 0;
}

#else // !ATV_HAS_RHF

#include <iostream>
int main() {
    std::cout << "14_fci requires -DATV_BUILD_RHF=ON\n";
    return 0;
}

#endif // ATV_HAS_RHF
