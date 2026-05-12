/// @file examples/aimd.cpp
/// @brief Ab initio molecular dynamics example — H2 with RHFForce.
///
/// Demonstrates BOMD and second-generation CPMD on a hydrogen molecule.

#include <iostream>
#include <iomanip>

#include "atomiverse/core/units.h"
#include "atomiverse/cpu/frame.h"
#include "atomiverse/cpu/scf/rhf_force.h"

using namespace ATV;

int main() {
    // ── H2 setup ─────────────────────────────────────────────────────
    cpu::Frame h2(2);
    h2.x = {0.0, 0.74};
    h2.y = {0.0, 0.0};
    h2.z = {0.0, 0.0};
    h2.vx = {0.0, 0.001};  // small initial velocity
    h2.vy = {0.0, 0.0};
    h2.vz = {0.0, 0.0};
    h2.fx = {0.0, 0.0};
    h2.fy = {0.0, 0.0};
    h2.fz = {0.0, 0.0};
    h2.mass = {1.008, 1.008};
    h2.type_id = {1, 1};
    h2.energy = 0;
    h2.kinetic = 0;
    h2.box = Mat3::diagonal(100.0, 100.0, 100.0);

    // ── Phase 1: BOMD ────────────────────────────────────────────────
    std::cout << "=== Born-Oppenheimer MD (full SCF each step) ===\n"
              << std::setw(5) << "Step"
              << std::setw(13) << "E_total"
              << std::setw(13) << "E_pot"
              << std::setw(13) << "E_kin"
              << std::setw(13) << "Bond"
              << std::setw(9)  << "SCF_it" << "\n";

    cpu::scf::RHFForce bomd{Real(1e-8), 100};
    bomd.on_step([&](const cpu::scf::AIMDStepRecord& rec) {
        double bond = static_cast<double>(h2.x(1) - h2.x(0));
        std::cout << std::setw(5) << rec.step
                  << std::fixed << std::setprecision(6)
                  << std::setw(13) << static_cast<double>(rec.total_energy)
                  << std::setw(13) << static_cast<double>(rec.potential_energy)
                  << std::setw(13) << static_cast<double>(rec.kinetic_energy)
                  << std::setw(13) << bond
                  << std::setw(9)  << rec.scf_iterations << "\n";
    });

    bomd.run(h2, 10, Real(0.5));  // 10 steps, 0.5 fs

    // ── Phase 2: CPMD2 ──────────────────────────────────────────────
    std::cout << "\n=== Second-generation CPMD (ASPC K=3, 2 corrector steps) ===\n"
              << std::setw(5) << "Step"
              << std::setw(13) << "E_total"
              << std::setw(13) << "E_pot"
              << std::setw(13) << "E_kin"
              << std::setw(13) << "Bond"
              << std::setw(9)  << "SCF_it"
              << std::setw(6)  << "Pred" << "\n";

    // Reset H2
    h2.x = {0.0, 0.74};
    h2.vx = {0.0, 0.001};
    h2.fx = {0.0, 0.0};
    h2.energy = 0;
    h2.kinetic = 0;

    cpu::scf::RHFForce cpmd{Real(1e-8), 100};
    cpmd.predictor(3, 2, Real(1e-4), true);
    cpmd.on_step([&](const cpu::scf::AIMDStepRecord& rec) {
        double bond = static_cast<double>(h2.x(1) - h2.x(0));
        std::cout << std::setw(5) << rec.step
                  << std::fixed << std::setprecision(6)
                  << std::setw(13) << static_cast<double>(rec.total_energy)
                  << std::setw(13) << static_cast<double>(rec.potential_energy)
                  << std::setw(13) << static_cast<double>(rec.kinetic_energy)
                  << std::setw(13) << bond
                  << std::setw(9)  << rec.scf_iterations
                  << std::setw(6)  << (rec.used_predictor ? "yes" : "no") << "\n";
    });

    cpmd.run(h2, 10, Real(0.5));

    std::cout << "\nDone.\n";
    return 0;
}
