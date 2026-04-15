/// @file examples/rhf_multi.cu
/// @brief Multi-system RHF: He, HeH+, LiH, H2 PES scan — CPU + CUDA.
///
/// Demonstrates that Atomiverse RHF works beyond the original H2 demo:
///   1. He atom (1 basis function)
///   2. HeH+ cation (2 basis functions, charge=+1)
///   3. LiH (3 basis functions: Li 1s+2s, H 1s)
///   4. H2 bond-length scan (20 points, 0.4--3.0 A)
///   5. LiH SCF convergence trace
///   6. CPU/CUDA parity on each system
///
/// When built with -DATV_USE_MATPLOT=ON, produces publication-quality figures:
///   h2_pes.png        -- H2 potential energy surface
///   scf_convergence.png -- LiH SCF convergence trace
///
/// Build:
///   cmake -B build -DATV_BUILD_RHF=ON -DATV_BUILD_EXAMPLES=ON [-DATV_USE_MATPLOT=ON]
///   cmake --build build --target rhf_multi

#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

#include "atomiverse/cpu/scf/integrals.h"
#include "atomiverse/cpu/scf/rhf_solve.h"
#include "atomiverse/cuda/scf/scf.cuh"
#include "atomiverse/cuda/scf/rhf_fock.cuh"
#include "atomiverse/cuda/scf/rhf_density.cuh"
#include "atomiverse/cuda/recorder/thermo.cuh"

#ifdef ATV_HAS_MATPLOT
#include <matplot/matplot.h>
#endif

using namespace ATV;

// ── Helpers ──────────────────────────────────────────────────────────

struct SystemResult {
    std::string name;
    int n_basis;
    int n_electrons;
    double cpu_energy;
    double cuda_energy;
    int iterations;
    bool converged;
    std::vector<cpu::scf::IterationRecord> trace;
};

static SystemResult run_system(const char* name, cpu::Frame& frame,
                               const cpu::scf::IntegralPack& integrals) {
    double tol = (sizeof(double) == sizeof(float)) ? double(1e-5) : double(1e-10);

    auto cpu_result = cpu::scf::rhf_solve(integrals, tol, 100);

    scf::SCF solver;
    solver.add(scf::RHFFock{integrals});
    solver.add(scf::EigenSolver{});
    solver.add(scf::RHFDensity{integrals.n_occ()});
    solver.add(scf::ConvergenceCheck{tol, false});
    solver.build();
    solver.run(frame, 100);

    double cpu_e = static_cast<double>(cpu_result.energy);
    double cuda_e = static_cast<double>(frame.energy);

    std::cout << "  " << std::left << std::setw(6) << name
              << std::right
              << "  basis=" << std::setw(2) << integrals.n_basis
              << "  elec=" << std::setw(2) << integrals.n_electrons
              << std::fixed << std::setprecision(8)
              << "  E_cpu=" << std::setw(12) << cpu_e
              << "  E_cuda=" << std::setw(12) << cuda_e
              << std::scientific << std::setprecision(1)
              << "  diff=" << std::abs(cpu_e - cuda_e)
              << "  iter=" << cpu_result.iterations
              << "  DIIS=" << ((!cpu_result.trace.empty() && cpu_result.trace.back().diis_active)
                    ? "on" : "off")
              << "\n";

    return {name, integrals.n_basis, integrals.n_electrons,
            cpu_e, cuda_e, cpu_result.iterations, cpu_result.converged,
            cpu_result.trace};
}

// ── Plotting (ACS style) ────────────────────────────────────────────

#ifdef ATV_HAS_MATPLOT

static void acs_style() {
    // ACS single-column: 3.25 in wide, ~2.5 in tall
    auto fig = matplot::gcf();
    fig->size(780, 600);
    auto ax = matplot::gca();
    ax->font_size(11);
}

/// Fig 1: H2 potential energy surface — the classic diatomic PES
static void plot_h2_pes(const std::vector<double>& r_vals,
                        const std::vector<double>& e_cpu,
                        const std::vector<double>& e_cuda) {
    auto fig = matplot::figure(true);
    fig->size(780, 600);

    // CPU curve: solid black line
    auto p1 = matplot::plot(r_vals, e_cpu);
    p1->line_width(2.0);
    p1->color({0.0f, 0.0f, 0.0f, 1.0f});
    p1->display_name("CPU (DIIS)");

    matplot::hold(matplot::on);

    // CUDA points: red circles
    auto p2 = matplot::plot(r_vals, e_cuda);
    p2->line_width(0);
    p2->marker(matplot::line_spec::marker_style::circle);
    p2->marker_size(6);
    p2->marker_face_color({0.8f, 0.2f, 0.2f, 1.0f});
    p2->color({0.8f, 0.2f, 0.2f, 1.0f});
    p2->display_name("CUDA");

    // Find and annotate minimum
    double e_min = 1e10, r_min = 0;
    for (size_t i = 0; i < r_vals.size(); ++i)
        if (e_cpu[i] < e_min) { e_min = e_cpu[i]; r_min = r_vals[i]; }

    // Axis labels
    auto ax = matplot::gca();
    ax->font_size(11);
    matplot::xlabel("R(H-H) / \\AA");
    matplot::ylabel("E / Hartree");
    matplot::title("H_{2} / STO-3G / RHF Potential Energy Surface");
    matplot::legend()->location(matplot::legend::general_alignment::topright);

    matplot::hold(matplot::off);
    matplot::save("h2_pes.png");
    std::cout << "  Saved h2_pes.png\n";
}

/// Fig 2: SCF convergence trace — shows DIIS acceleration
static void plot_convergence(const std::vector<cpu::scf::IterationRecord>& trace,
                             const char* system_name) {
    if (trace.size() < 2) return;

    std::vector<double> iters, energies, residuals;
    for (auto& rec : trace) {
        iters.push_back(rec.iter);
        energies.push_back(static_cast<double>(rec.energy));
        residuals.push_back(std::max(static_cast<double>(rec.residual), 1e-16));
    }

    auto fig = matplot::figure(true);
    fig->size(780, 600);

    // Residual on log scale
    auto p1 = matplot::semilogy(iters, residuals);
    p1->line_width(2.0);
    p1->color({0.2f, 0.4f, 0.8f, 1.0f});
    p1->marker(matplot::line_spec::marker_style::diamond);
    p1->marker_size(5);
    p1->marker_face_color({0.2f, 0.4f, 0.8f, 1.0f});

    auto ax = matplot::gca();
    ax->font_size(11);
    matplot::xlabel("SCF Iteration");
    matplot::ylabel("||\\DeltaP||_F");
    std::string title = std::string(system_name) + " / STO-3G / RHF Convergence";
    matplot::title(title);

    matplot::save("scf_convergence.png");
    std::cout << "  Saved scf_convergence.png\n";
}

#endif // ATV_HAS_MATPLOT

// ── Main ─────────────────────────────────────────────────────────────

int main() {
    std::cout << "=== Atomiverse Multi-System RHF (STO-3G) ===\n\n";

    // ── Single-point calculations ────────────────────────────────────
    std::cout << "Single-point energies:\n";

    SystemResult he_result, heh_result, lih_result;

    // He
    {
        cpu::Frame f(1);
        f.x[0] = 0; f.y[0] = 0; f.z[0] = 0;
        f.mass[0] = double(4.003); f.type_id[0] = 2;
        f.box = Mat3::cubic(100.0);
        auto integrals = cpu::scf::build_integrals(f, cpu::scf::BasisSet::STO3G);
        he_result = run_system("He", f, integrals);
    }

    // HeH+
    {
        cpu::Frame f(2);
        f.x[0] = double(0);     f.y[0] = 0; f.z[0] = 0;
        f.x[1] = double(0.772); f.y[1] = 0; f.z[1] = 0;
        f.mass[0] = double(4.003); f.mass[1] = double(1.008);
        f.type_id[0] = 2; f.type_id[1] = 1;
        f.box = Mat3::cubic(100.0);
        auto integrals = cpu::scf::build_integrals(f, cpu::scf::BasisSet::STO3G, 1);
        heh_result = run_system("HeH+", f, integrals);
    }

    // LiH
    {
        cpu::Frame f(2);
        f.x[0] = double(0);     f.y[0] = 0; f.z[0] = 0;
        f.x[1] = double(1.595); f.y[1] = 0; f.z[1] = 0;
        f.mass[0] = double(6.941); f.mass[1] = double(1.008);
        f.type_id[0] = 3; f.type_id[1] = 1;
        f.box = Mat3::cubic(100.0);
        auto integrals = cpu::scf::build_integrals(f, cpu::scf::BasisSet::STO3G);
        lih_result = run_system("LiH", f, integrals);
    }

    // ── H2 PES scan ──────────────────────────────────────────────────
    std::cout << "\nH2 bond-length scan:\n"
              << std::setw(8) << "R(A)"
              << std::setw(15) << "E_cpu(Ha)"
              << std::setw(15) << "E_cuda(Ha)"
              << std::setw(5) << "iter" << "\n";

    double tol = (sizeof(double) == sizeof(float)) ? double(1e-5) : double(1e-10);
    std::vector<double> r_vals, e_cpu_vals, e_cuda_vals;

    int n_points = 20;
    for (int i = 0; i < n_points; ++i) {
        double r = 0.4 + i * (3.0 - 0.4) / (n_points - 1);
        cpu::Frame f(2);
        f.x[0] = double(0); f.y[0] = 0; f.z[0] = 0;
        f.x[1] = double(r); f.y[1] = 0; f.z[1] = 0;
        f.mass[0] = double(1.008); f.mass[1] = double(1.008);
        f.type_id[0] = 1; f.type_id[1] = 1;
        f.box = Mat3::cubic(100.0);

        auto integrals = cpu::scf::build_integrals(f, cpu::scf::BasisSet::STO3G);
        auto cpu_result = cpu::scf::rhf_solve(integrals, tol, 100);

        scf::SCF solver;
        solver.add(scf::RHFFock{integrals});
        solver.add(scf::EigenSolver{});
        solver.add(scf::RHFDensity{integrals.n_occ()});
        solver.add(scf::ConvergenceCheck{tol, false});
        solver.build();
        solver.run(f, 100);

        double cpu_e = static_cast<double>(cpu_result.energy);
        double cuda_e = static_cast<double>(f.energy);

        std::cout << std::fixed
                  << std::setw(8) << std::setprecision(3) << r
                  << std::setw(15) << std::setprecision(10) << cpu_e
                  << std::setw(15) << cuda_e
                  << std::setw(5) << cpu_result.iterations << "\n";

        r_vals.push_back(r);
        e_cpu_vals.push_back(cpu_e);
        e_cuda_vals.push_back(cuda_e);
    }

    // ── Plots ────────────────────────────────────────────────────────
#ifdef ATV_HAS_MATPLOT
    std::cout << "\nGenerating figures...\n";
    plot_h2_pes(r_vals, e_cpu_vals, e_cuda_vals);
    plot_convergence(lih_result.trace, "LiH");
#else
    std::cout << "\nPlotting disabled. Rebuild with -DATV_USE_MATPLOT=ON for figures.\n";
#endif

    std::cout << "\n=== Done ===\n";
    return 0;
}
