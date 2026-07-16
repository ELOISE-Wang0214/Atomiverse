/// @file examples/13_mp2.cpp
/// @brief Tutorial example: closed-shell MP2 correlation energy on H2O/STO-3G.
///
/// Second-order Moller-Plesset perturbation theory (MP2) is the cheapest
/// post-Hartree-Fock correlation method. It adds the leading dynamic-correlation
/// correction to the RHF reference via the (ia|jb) ovov integral block obtained
/// by transforming the AO two-electron integrals into the MO basis.
///
/// Workflow (all energies in [Hartree]):
///   1. Build H2O geometry in cpu::Frame
///   2. Compute molecular integrals with the two-electron ERI tensor
///   3. Run RHF to get the mean-field reference (E_HF)
///   4. Run MP2 to get the correlation energy (E_corr) and total (E_MP2)
///
/// Build:
///   cmake -B build -DATV_BUILD_RHF=ON -DATV_BUILD_EXAMPLES=ON
///   cmake --build build --target 13_mp2
///
/// Run:
///   ./build/examples/13_mp2

#ifdef ATV_HAS_RHF

#include <cmath>
#include <iomanip>
#include <iostream>

#include "atomiverse/core/mat3.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/scf/integrals.h"
#include "atomiverse/cpu/scf/mp2.h"
#include "atomiverse/cpu/scf/scf_solver.h"

using namespace ATV;

int main() {
    // 1. H2O geometry (Angstrom); type = atomic number Z.
    cpu::Frame frame(3);
    frame.span(ATV::kAtom_x)[0] = 0.0;
    frame.span(ATV::kAtom_x)[1] = 0.0;
    frame.span(ATV::kAtom_x)[2] = 0.0;
    frame.span(ATV::kAtom_y)[0] = 0.0;
    frame.span(ATV::kAtom_y)[1] = 0.7572;
    frame.span(ATV::kAtom_y)[2] = -0.7572;
    frame.span(ATV::kAtom_z)[0] = 0.1173;
    frame.span(ATV::kAtom_z)[1] = -0.4692;
    frame.span(ATV::kAtom_z)[2] = -0.4692;
    frame.span(ATV::kAtom_type)[0] = 8;
    frame.span(ATV::kAtom_type)[1] = 1;
    frame.span(ATV::kAtom_type)[2] = 1;
    frame.span(ATV::kAtom_mass)[0] = 15.999;
    frame.span(ATV::kAtom_mass)[1] = 1.008;
    frame.span(ATV::kAtom_mass)[2] = 1.008;
    frame.simbox = ATV::cpu::Box::cubic(50.0);

    std::cout << "=== H2O / STO-3G / MP2 (CPU) ===\n\n";

    // 2. Integrals -- the two-electron ERI tensor is needed by the AO->MO
    //    transform inside mp2_energy.
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

    // 4. Closed-shell canonical MP2.
    auto mp2 = cpu::scf::mp2_energy(ints, scf);

    std::cout << std::fixed << std::setprecision(8) << "E_HF:    " << scf.energy
              << " Hartree\n"
              << "E_corr:  " << mp2.e_corr
              << " Hartree (MP2 correlation, <= 0)\n"
              << "E_MP2:   " << mp2.e_tot << " Hartree (= E_HF + E_corr)\n";

    // Consistency: E_HF + E_corr == E_MP2 to print precision.
    if (std::abs(scf.energy + mp2.e_corr - mp2.e_tot) > 5e-9) {
        std::cerr << "inconsistent: E_HF + E_corr != E_MP2\n";
        return 1;
    }

    std::cout << "\n=== Done ===\n";
    return 0;
}

#else // !ATV_HAS_RHF

#include <iostream>
int main() {
    std::cout << "13_mp2 requires -DATV_BUILD_RHF=ON\n";
    return 0;
}

#endif // ATV_HAS_RHF
