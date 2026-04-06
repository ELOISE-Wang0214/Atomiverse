/// @file examples/tutorial/04_uks.cpp
/// @brief Tutorial example for unrestricted Kohn-Sham DFT.
///
/// Demonstrates:
///   - UKS(LSDA) on atoms and molecules
///   - Energy decomposition (E_1e + E_J + E_xc + V_nuc = E_total)
///   - UKS vs UHF comparison
///   - Grid-integrated electron counts
///
/// Build:
///   cmake -B build -DATV_BUILD_CUDA=OFF -DATV_BUILD_RHF=ON -DATV_BUILD_EXAMPLES=ON
///   cmake --build build --target 04_uks
///
/// Run:
///   ./build/examples/04_uks

#ifdef ATV_HAS_RHF

#include <iostream>
#include <iomanip>
#include <cmath>

#include "atomiverse/cpu/frame.h"
#include "atomiverse/cpu/scf/integrals.h"
#include "atomiverse/cpu/scf/uhf_solve.h"
#include "atomiverse/cpu/scf/uks_solve.h"

using namespace ATV;

static cpu::Frame make_atom(int Z) {
    cpu::Frame f(1);
    f.atoms.x[0] = 0; f.atoms.y[0] = 0; f.atoms.z[0] = 0;
    f.atoms.type[0] = Z;
    f.box = Mat3::cubic(100.0);
    return f;
}

static cpu::Frame make_h2(double bond_angstrom) {
    cpu::Frame f(2);
    f.atoms.x[0] = 0;                 f.atoms.y[0] = 0; f.atoms.z[0] = 0;
    f.atoms.x[1] = bond_angstrom;     f.atoms.y[1] = 0; f.atoms.z[1] = 0;
    f.atoms.type[0] = 1; f.atoms.type[1] = 1;
    f.box = Mat3::cubic(100.0);
    return f;
}

static void print_result(const char* name, const cpu::scf::UKSResult& r, Real V_nuc) {
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "  " << name << ":\n";
    std::cout << "    converged  = " << (r.converged ? "yes" : "NO") << " (" << r.iterations << " iter)\n";
    std::cout << "    E_total    = " << r.energy << " Ha\n";
    std::cout << "    E_1e       = " << r.one_electron_energy << " Ha\n";
    std::cout << "    E_J        = " << r.coulomb_energy << " Ha\n";
    std::cout << "    E_xc       = " << r.xc_energy << " Ha\n";
    std::cout << "    V_nuc      = " << V_nuc << " Ha\n";
    std::cout << "    sum check  = " << (r.one_electron_energy + r.coulomb_energy + r.xc_energy + V_nuc) << " Ha\n";
    std::cout << "    int(n_a)   = " << r.integrated_n_alpha << "\n";
    std::cout << "    int(n_b)   = " << r.integrated_n_beta << "\n";
    std::cout << "    <S^2> err  = " << r.spin_contamination << "\n";
}

int main() {
    std::cout << "=== Atomiverse UKS Examples ===\n\n";

    // Example 1: He atom (closed-shell singlet)
    {
        std::cout << "--- Example 1: He atom (LSDA, singlet) ---\n";
        auto frame = make_atom(2);
        auto ints = cpu::scf::build_integrals_sto3g(frame);
        auto uks = cpu::scf::uks_solve(ints, 1, 1, cpu::scf::XCType::LSDA, 50, 110);
        auto uhf = cpu::scf::uhf_solve(ints, 1, 1);
        print_result("UKS-LSDA", uks, ints.V_nuc);
        std::cout << "    UHF E      = " << uhf.energy << " Ha\n";
        std::cout << "    diff       = " << (uks.energy - uhf.energy) << " Ha\n\n";
    }

    // Example 2: H atom (doublet)
    {
        std::cout << "--- Example 2: H atom (LSDA, doublet) ---\n";
        auto frame = make_atom(1);
        auto ints = cpu::scf::build_integrals_sto3g(frame);
        auto uks = cpu::scf::uks_solve(ints, 1, 0, cpu::scf::XCType::LSDA, 50, 110);
        print_result("UKS-LSDA", uks, ints.V_nuc);
        std::cout << "\n";
    }

    // Example 3: Li atom (doublet)
    {
        std::cout << "--- Example 3: Li atom (LSDA, doublet) ---\n";
        auto frame = make_atom(3);
        auto ints = cpu::scf::build_integrals_sto3g(frame);
        auto uks = cpu::scf::uks_solve(ints, 2, 1, cpu::scf::XCType::LSDA, 50, 110);
        print_result("UKS-LSDA", uks, ints.V_nuc);
        std::cout << "    eigenvalues alpha: ";
        for (int i = 0; i < ints.n_basis; ++i) std::cout << uks.eigenvalues_alpha[i] << " ";
        std::cout << "\n    eigenvalues beta:  ";
        for (int i = 0; i < ints.n_basis; ++i) std::cout << uks.eigenvalues_beta[i] << " ";
        std::cout << "\n\n";
    }

    // Example 4: H2 dissociation curve (UHF vs UKS-LSDA)
    {
        std::cout << "--- Example 4: H2 dissociation (UHF vs UKS-LSDA) ---\n";
        std::cout << std::setw(8) << "R (A)" << std::setw(14) << "E_UHF"
                  << std::setw(14) << "E_UKS" << std::setw(14) << "diff" << "\n";
        for (double r : {0.5, 0.6, 0.74, 1.0, 1.5, 2.0, 2.5, 3.0}) {
            auto frame = make_h2(r);
            auto ints = cpu::scf::build_integrals_sto3g(frame);
            auto uhf = cpu::scf::uhf_solve(ints, 1, 1, Real(1e-9), 400);
            auto uks = cpu::scf::uks_solve(ints, 1, 1, cpu::scf::XCType::LSDA, 50, 110, Real(1e-7), 100);

            std::cout << std::setw(8) << r
                      << std::setw(14) << uhf.energy
                      << std::setw(14) << uks.energy
                      << std::setw(14) << (uks.energy - uhf.energy);
            if (!uks.converged) std::cout << " *";
            std::cout << "\n";
        }
    }

    return 0;
}

#else
#include <iostream>
int main() { std::cout << "UKS requires -DATV_BUILD_RHF=ON\n"; return 1; }
#endif
