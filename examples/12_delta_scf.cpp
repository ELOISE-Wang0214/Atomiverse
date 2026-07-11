/// @file examples/12_delta_scf.cpp
/// @brief Tutorial example: ΔSCF vertical ionization energies vs Koopmans.
///
/// For each of {He, H₂, H₂O, NH₃, CH₄} at STO-3G and 6-31G this computes the
/// vertical ionization energy two ways at the SAME neutral geometry:
///
///   - Koopmans:  I_K   = −ε_HOMO           [frozen-orbital estimate]
///   - ΔSCF:      I_ΔSCF = E⁺(N−1) − E(N)   [orbital-relaxed estimate]
///
/// The neutral N-electron system is solved with restricted Hartree-Fock
/// (`build_integrals(frame, basis, charge=0)` + `rhf_solve`). The (N−1)-electron
/// doublet cation is solved with unrestricted Hartree-Fock, warm-started from
/// the neutral density via `SCFInitialGuess` (P_alpha = P_beta = P_neutral/2)
/// so the SCF relaxes toward the ground-state ion rather than a hole-shifted
/// excited solution.
///
/// Spin bookkeeping (doublet, M = 2): removing one electron from the highest
/// doubly-occupied level leaves n_alpha = N/2 (the neutral n_occ) and
/// n_beta = N/2 − 1. The cation ⟨S²⟩ (= spin_contamination + 0.75) is printed
/// as a health check — a clean doublet sits at ≈ 0.75.
///
/// All energies are computed in [Hartree] and printed in [eV] using the
/// canonical conversion `ATV::kHartree2eV` (no inline magic numbers). The
/// variational argument guarantees I_ΔSCF ≤ I_K, so the relaxation column
/// (I_K − I_ΔSCF) is non-negative; it vanishes for He and H₂ in a minimal
/// basis, where the single-electron cation has no orbital to relax.
///
/// Build:
///   cmake -B build -DATV_BUILD_CUDA=OFF -DATV_BUILD_RHF=ON -DATV_BUILD_EXAMPLES=ON
///   cmake --build build --target 12_delta_scf
///   ./build/examples/12_delta_scf
///
/// See docs/tutorial/ch03-open-shell.md ("ΔSCF: ionization energies beyond
/// Koopmans") for the physics, and ch02-rhf.md for the Koopmans section.

#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include <xtensor/containers/xtensor.hpp>

#include "atomiverse/core/units.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/scf/integrals.h"
#include "atomiverse/cpu/scf/rhf_solve.h"
#include "atomiverse/cpu/scf/scf_solve.h"

using namespace ATV;

// ── Geometries (Angstrom); atoms.type = atomic number Z ──────────────────────

static cpu::Frame make_he() {
    cpu::Frame f(1);
    f.atoms.x = {0.0};
    f.atoms.y = {0.0};
    f.atoms.z = {0.0};
    f.atoms.type = {2};
    f.box = Mat3::cubic(100.0);
    return f;
}

static cpu::Frame make_h2() {
    cpu::Frame f(2);
    f.atoms.x = {0.0, 0.7414};
    f.atoms.y = {0.0, 0.0};
    f.atoms.z = {0.0, 0.0};
    f.atoms.type = {1, 1};
    f.box = Mat3::cubic(100.0);
    return f;
}

static cpu::Frame make_h2o() {
    cpu::Frame f(3);
    f.atoms.x = {0.0, 0.0, 0.0};
    f.atoms.y = {0.0, 0.7572, -0.7572};
    f.atoms.z = {0.1173, -0.4692, -0.4692};
    f.atoms.type = {8, 1, 1};
    f.box = Mat3::cubic(100.0);
    return f;
}

static cpu::Frame make_nh3() {
    cpu::Frame f(4);
    f.atoms.x = {0.0, 0.0, 0.812328, -0.812328};
    f.atoms.y = {0.0, 0.938035, -0.469017, -0.469017};
    f.atoms.z = {0.113880, -0.265825, -0.265825, -0.265825};
    f.atoms.type = {7, 1, 1, 1};
    f.box = Mat3::cubic(100.0);
    return f;
}

static cpu::Frame make_ch4() {
    cpu::Frame f(5);
    f.atoms.x = {0.0, 0.627640, -0.627640, 0.627640, -0.627640};
    f.atoms.y = {0.0, 0.627640, -0.627640, -0.627640, 0.627640};
    f.atoms.z = {0.0, 0.627640, 0.627640, -0.627640, -0.627640};
    f.atoms.type = {6, 1, 1, 1, 1};
    f.box = Mat3::cubic(100.0);
    return f;
}

struct System {
    const char* name;
    cpu::Frame (*make)();
};

int main() {
    const std::vector<System> systems = {
        {"He", &make_he},   {"H2", &make_h2},   {"H2O", &make_h2o},
        {"NH3", &make_nh3}, {"CH4", &make_ch4},
    };
    const std::vector<std::pair<const char*, cpu::scf::BasisSet>> bases = {
        {"STO-3G", cpu::scf::BasisSet::STO3G},
        {"6-31G", cpu::scf::BasisSet::Pople631G},
    };

    const double tol = 1e-10; // tight: resolve the He/H2 exact-relaxation zero

    std::cout << "ΔSCF vertical ionization vs Koopmans (energies in eV)\n"
              << "  I_K    = -eps_HOMO           (frozen orbitals)\n"
              << "  I_dSCF = E+(N-1) - E(N)      (orbital relaxed, warm-started)\n"
              << "  relax  = I_K - I_dSCF >= 0   (= 0 for He, H2 minimal basis)\n\n";

    std::cout << std::left << std::setw(6) << "sys" << std::setw(8) << "basis"
              << std::right << std::setw(12) << "I_K[eV]" << std::setw(12)
              << "I_dSCF[eV]" << std::setw(11) << "relax[eV]" << std::setw(11)
              << "<S^2>" << std::setw(7) << "conv" << "\n";
    std::cout << std::string(67, '-') << "\n";

    for (const auto& sys : systems) {
        for (const auto& [basis_name, basis] : bases) {
            cpu::Frame frame = sys.make();

            // Neutral N-electron RHF.
            auto ints_n = cpu::scf::build_integrals(frame, basis, /*charge=*/0);
            auto neutral = cpu::scf::rhf_solve(ints_n, tol, 200);
            const int n_occ = ints_n.n_occ();
            const double i_koopmans = -neutral.eigenvalues(n_occ - 1); // [Ha]

            // Doublet cation: (N-1) electrons, n_alpha=N/2, n_beta=N/2-1.
            const int n_alpha = n_occ;
            const int n_beta = n_occ - 1;
            auto ints_c = cpu::scf::build_integrals(frame, basis, /*charge=*/+1);
            auto cfg = cpu::scf::uhf_config(ints_c, n_alpha, n_beta, tol);
            cfg.max_iter = 300;

            // Warm start from the neutral density: P_alpha = P_beta = P/2.
            cpu::scf::SCFInitialGuess guess;
            guess.P_alpha = 0.5 * neutral.P;
            guess.P_beta = 0.5 * neutral.P;
            guess.has_alpha = true;
            guess.has_beta = true;
            auto cation = cpu::scf::scf_solve(ints_c, cfg, &guess);

            const double i_delta = cation.energy - neutral.energy; // [Ha]
            const double relax = i_koopmans - i_delta;             // [Ha] >= 0
            const double s2 = cation.spin_contamination + 0.75;    // <S^2>
            const bool conv = neutral.converged && cation.converged;

            std::cout << std::left << std::setw(6) << sys.name << std::setw(8)
                      << basis_name << std::right << std::fixed
                      << std::setprecision(4) << std::setw(12)
                      << i_koopmans * kHartree2eV << std::setw(12)
                      << i_delta * kHartree2eV << std::setw(11)
                      << relax * kHartree2eV << std::setw(11) << s2
                      << std::setw(7) << (conv ? "yes" : "NO") << "\n";
        }
    }

    std::cout
        << "\nReading the table:\n"
        << "  - relax >= 0 everywhere: orbital relaxation lowers the ionization\n"
        << "    energy below the Koopmans estimate (variational bound).\n"
        << "  - He and H2 show relax = 0 in STO-3G: the single-electron cation\n"
        << "    orbital is symmetry-fixed, so there is nothing to relax. In\n"
        << "    6-31G the lone electron gains relaxation freedom, so relax > 0.\n"
        << "  - <S^2> ~ 0.75 confirms a clean doublet (no spin contamination).\n";
    return 0;
}
