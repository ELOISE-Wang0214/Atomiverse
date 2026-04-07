/// @file examples/bench_aimd.cpp
/// @brief AIMD benchmark — BOMD vs CPMD2 wall time and energy drift on H2.
///
/// Standalone benchmark (no Google Benchmark) that compares full-SCF BOMD
/// against ASPC predictor-corrector (CPMD2) on a hydrogen molecule.
/// Reports wall time per step, total energy drift, and SCF iterations.

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <iomanip>

#include "atomiverse/core/units.h"
#include "atomiverse/cpu/frame.h"
#include "atomiverse/cpu/scf/rhf_force.h"

using namespace ATV;

/// @brief Build a fresh H2 frame at 0.74 Angstrom with small initial velocity.
static cpu::Frame make_h2() {
    cpu::Frame h2(2);
    h2.x = {0.0, 0.74};
    h2.y = {0.0, 0.0};
    h2.z = {0.0, 0.0};
    h2.vx = {0.0, 0.001};
    h2.vy = {0.0, 0.0};
    h2.vz = {0.0, 0.0};
    h2.fx = {0.0, 0.0};
    h2.fy = {0.0, 0.0};
    h2.fz = {0.0, 0.0};
    h2.mass = {1.008, 1.008};
    h2.type_id = {1, 1};
    h2.energy = 0;
    h2.kinetic = 0;
    h2.box = Mat3::cubic(100.0);
    return h2;
}

int main() {
    constexpr int n_steps = 10;
    constexpr double dt = double(0.5);  // fs

    // ── BOMD benchmark ──────────────────────────────────────────────
    cpu::Frame h2_bomd = make_h2();
    cpu::scf::RHFForce bomd{double(1e-8), 100};

    auto t0_bomd = std::chrono::high_resolution_clock::now();
    bomd.run(h2_bomd, n_steps, dt);
    auto t1_bomd = std::chrono::high_resolution_clock::now();

    double bomd_wall_ms = std::chrono::duration<double, std::milli>(t1_bomd - t0_bomd).count();

    const auto& bomd_trace = bomd.trace();

    double bomd_e_min = static_cast<double>(bomd_trace[0].total_energy);
    double bomd_e_max = bomd_e_min;
    for (const auto& rec : bomd_trace) {
        double e = static_cast<double>(rec.total_energy);
        bomd_e_min = std::min(bomd_e_min, e);
        bomd_e_max = std::max(bomd_e_max, e);
    }
    double bomd_drift = bomd_e_max - bomd_e_min;

    // ── CPMD2 benchmark ─────────────────────────────────────────────
    cpu::Frame h2_cpmd = make_h2();
    cpu::scf::RHFForce cpmd{double(1e-8), 100};
    cpmd.predictor(3, 2, double(1e-4), true);

    auto t0_cpmd = std::chrono::high_resolution_clock::now();
    cpmd.run(h2_cpmd, n_steps, dt);
    auto t1_cpmd = std::chrono::high_resolution_clock::now();

    double cpmd_wall_ms = std::chrono::duration<double, std::milli>(t1_cpmd - t0_cpmd).count();

    const auto& cpmd_trace = cpmd.trace();

    double cpmd_e_min = static_cast<double>(cpmd_trace[0].total_energy);
    double cpmd_e_max = cpmd_e_min;
    int cpmd_scf_total = 0;
    for (const auto& rec : cpmd_trace) {
        double e = static_cast<double>(rec.total_energy);
        cpmd_e_min = std::min(cpmd_e_min, e);
        cpmd_e_max = std::max(cpmd_e_max, e);
        cpmd_scf_total += rec.scf_iterations;
    }
    double cpmd_drift = cpmd_e_max - cpmd_e_min;

    // ── Summary table ───────────────────────────────────────────────
    std::cout << "\n"
              << "========================================================\n"
              << "  AIMD Benchmark: H2 at 0.74 A, " << n_steps
              << " steps, dt = " << std::fixed << std::setprecision(1)
              << static_cast<double>(dt) << " fs\n"
              << "========================================================\n\n";

    std::cout << std::setw(20) << std::left << ""
              << std::setw(13) << std::right << "BOMD"
              << std::setw(13) << "CPMD2" << "\n"
              << std::setw(20) << std::left << "--------------------"
              << std::setw(13) << std::right << "------------"
              << std::setw(13) << "------------" << "\n";

    std::cout << std::fixed;
    std::cout << std::setw(20) << std::left << "Total wall time"
              << std::setprecision(3) << std::right
              << std::setw(10) << bomd_wall_ms << " ms"
              << std::setw(10) << cpmd_wall_ms << " ms\n";

    std::cout << std::setw(20) << std::left << "Wall time / step"
              << std::setprecision(3) << std::right
              << std::setw(10) << (bomd_wall_ms / n_steps) << " ms"
              << std::setw(10) << (cpmd_wall_ms / n_steps) << " ms\n";

    std::cout << std::scientific << std::setprecision(2);
    std::cout << std::setw(20) << std::left << "Energy drift"
              << std::right
              << std::setw(10) << bomd_drift << " eV"
              << std::setw(10) << cpmd_drift << " eV\n";

    std::cout << std::fixed << std::setprecision(6);
    std::cout << std::setw(20) << std::left << "Final E_total"
              << std::right
              << std::setw(10) << static_cast<double>(bomd_trace.back().total_energy) << " eV"
              << std::setw(10) << static_cast<double>(cpmd_trace.back().total_energy) << " eV\n";

    std::cout << std::setw(20) << std::left << "Avg SCF iter/step"
              << std::right
              << std::setw(10) << "---" << "   "
              << std::fixed << std::setprecision(1)
              << std::setw(10) << (static_cast<double>(cpmd_scf_total) / n_steps) << "\n";

    std::cout << std::setw(20) << std::left << "Speedup vs BOMD"
              << std::right
              << std::setw(10) << "1.00x" << "   "
              << std::fixed << std::setprecision(2)
              << std::setw(10) << (bomd_wall_ms / cpmd_wall_ms) << "x\n";

    // ── Per-step CPMD2 detail ───────────────────────────────────────
    std::cout << "\nCPMD2 per-step SCF iterations:\n"
              << std::setw(5) << "Step"
              << std::setw(9) << "SCF_iter"
              << std::setw(6) << "Pred" << "\n"
              << std::setw(5) << "-----"
              << std::setw(9) << "--------"
              << std::setw(6) << "-----" << "\n";
    for (const auto& rec : cpmd_trace) {
        std::cout << std::setw(5) << rec.step
                  << std::setw(9) << rec.scf_iterations
                  << std::setw(6) << (rec.used_predictor ? "yes" : "no") << "\n";
    }

    std::cout << "\nDone.\n";
    return 0;
}
