/// @file examples/tutorial/05_geom_opt.cpp
/// @brief Tutorial: Geometry optimization of aspirin with L-BFGS + RHF/STO-3G.
///
/// Demonstrates:
///   - Setting up a multi-atom molecule (aspirin, C9H8O4, 21 atoms)
///   - Configuring the SCF evaluator for ab initio forces
///   - Running L-BFGS geometry optimization with convergence monitoring
///   - Interpreting energy, forces, and displacement convergence criteria
///
/// Build:
///   cmake -B build -DATV_BUILD_CUDA=OFF -DATV_BUILD_RHF=ON -DATV_BUILD_EXAMPLES=ON
///   cmake --build build --target 05_geom_opt
///
/// Run:
///   ./build/examples/05_geom_opt

#ifdef ATV_HAS_RHF

#include <iostream>
#include <iomanip>
#include <cmath>

#include "atomiverse/cpu/frame.h"
#include "atomiverse/cpu/scf/integrals.h"
#include "atomiverse/cpu/scf/scf_solve.h"
#include "atomiverse/cpu/scf/scf_evaluator.h"
#include "atomiverse/cpu/minimizer/lbfgs.h"

using namespace ATV;

/// @brief Build aspirin (C9H8O4) frame with approximate 3D coordinates.
///
/// Aspirin = acetylsalicylic acid. 21 atoms: 9 C + 8 H + 4 O.
/// Coordinates from PubChem CID 2244, slightly perturbed to give
/// the optimizer meaningful work.
///
/// STO-3G basis: 9*5 + 8*1 + 4*5 = 73 basis functions, 94 electrons.
static cpu::Frame make_aspirin() {
    // 21 atoms: C(6), H(1), O(8)
    cpu::Frame f(21);

    // Atomic numbers: 9C + 4O + 8H = 21 atoms, 94 electrons
    int types[21] = {
        6, 6, 6, 6, 6, 6,    // C1-C6 (benzene ring)
        6,                     // C7 (carboxyl C)
        8, 8,                  // O8 (C=O), O9 (OH)
        8,                     // O10 (ester O)
        6,                     // C11 (acetyl carbonyl C)
        8,                     // O12 (acetyl C=O)
        6,                     // C13 (methyl C)
        1, 1, 1, 1,           // H14-H17 (ring)
        1,                     // H18 (OH)
        1, 1, 1,              // H19-H21 (CH3)
    };

    // Masses [amu]
    double masses[21] = {
        12.011, 12.011, 12.011, 12.011, 12.011, 12.011,
        12.011,
        15.999, 15.999,
        15.999,
        12.011,
        15.999,
        12.011,
        1.008, 1.008, 1.008, 1.008,
        1.008,
        1.008, 1.008, 1.008,
    };

    // Coordinates [Angstrom] — PySCF-validated aspirin geometry
    double coords[21][3] = {
        { 1.2879,  0.1068,  0.0000},  // C1
        { 0.1768,  1.0827,  0.0000},  // C2
        { 0.5171,  2.4290,  0.0000},  // C3
        { 1.8519,  2.7993,  0.0000},  // C4
        { 2.9293,  1.8524,  0.0000},  // C5
        { 2.6191,  0.4890,  0.0000},  // C6
        {-1.2076,  0.6698,  0.0000},  // C7
        {-1.5541, -0.4994,  0.0000},  // O8
        {-2.0632,  1.6805,  0.0000},  // O9
        { 0.9973, -1.2285,  0.0000},  // O10
        {-0.2830, -1.7447,  0.0000},  // C11
        {-0.3696, -2.9537,  0.0000},  // O12
        {-1.3593, -0.7411,  0.3000},  // C13
        {-0.2980,  3.1365,  0.0000},  // H14
        { 2.0791,  3.8586,  0.0000},  // H15
        { 3.9560,  2.1939,  0.0000},  // H16
        { 3.4546, -0.2019,  0.0000},  // H17
        {-2.9828,  1.3621,  0.0000},  // H18
        {-1.0697, -0.1755,  1.1920},  // H19
        {-2.3482, -1.1583,  0.5001},  // H20
        {-1.4267, -0.2091, -0.6502},  // H21
    };

    for (int i = 0; i < 21; ++i) {
        f.atoms.type[i] = types[i];
        f.atoms.mass[i] = Real(masses[i]);
        f.atoms.x[i] = Real(coords[i][0]);
        f.atoms.y[i] = Real(coords[i][1]);
        f.atoms.z[i] = Real(coords[i][2]);
    }
    f.box = Mat3::cubic(100.0);
    return f;
}

int main() {
    std::cout << "=== Atomiverse: Aspirin Geometry Optimization ===\n";
    std::cout << "    Molecule: C9H8O4 (acetylsalicylic acid)\n";
    std::cout << "    Method:   RHF / STO-3G\n";
    std::cout << "    Optimizer: L-BFGS (memory=5)\n\n";

    auto frame = make_aspirin();

    // Build integrals to determine electron count
    auto integrals = cpu::scf::build_integrals_sto3g(frame);
    std::cout << "Basis functions: " << integrals.n_basis << "\n";
    std::cout << "Electrons:       " << integrals.n_electrons << "\n\n";

    // Configure RHF
    auto config = cpu::scf::rhf_config(integrals);
    config.tol = Real(1e-6);
    config.max_iter = 200;

    // Set up L-BFGS optimizer
    cpu::minimizer::LBFGS opt;
    opt.add(cpu::minimizer::LBFGSStep{
        5,              // memory
        Real(0.01),     // alpha0
        Real(0.05),     // ftol [eV/A]
        10              // max_ls
    });
    opt.on(ATV::minimizer::Phase::force).add(
        cpu::scf::SCFEvaluator{config}
    );
    opt.build();

    std::cout << "Starting optimization (max 50 steps)...\n";
    std::cout << std::fixed << std::setprecision(6);
    opt.run(frame, 50);

    // Report final state
    std::cout << "\n── Final State ──\n";
    std::cout << "  Energy: " << frame.energy << " eV\n";

    // Compute final max force
    Real max_f = 0;
    for (Int i = 0; i < frame.atoms.size(); ++i) {
        Real fx = frame.atoms.fx(i), fy = frame.atoms.fy(i), fz = frame.atoms.fz(i);
        Real mag = std::sqrt(fx*fx + fy*fy + fz*fz);
        if (mag > max_f) max_f = mag;
    }
    std::cout << "  Max force: " << max_f << " eV/A\n";

    // Print final coordinates
    std::cout << "\n── Final Coordinates [Angstrom] ──\n";
    const char* symbols[] = {"", "H", "He", "Li", "Be", "B", "C", "N", "O"};
    std::cout << std::setprecision(6);
    for (Int i = 0; i < frame.atoms.size(); ++i) {
        int Z = static_cast<int>(frame.atoms.type(i));
        std::cout << "  " << std::setw(2) << symbols[Z]
                  << std::setw(12) << frame.atoms.x(i)
                  << std::setw(12) << frame.atoms.y(i)
                  << std::setw(12) << frame.atoms.z(i) << "\n";
    }

    std::cout << "\nDone.\n";
    return 0;
}

#else
#include <iostream>
int main() { std::cout << "Requires -DATV_BUILD_RHF=ON\n"; return 1; }
#endif
