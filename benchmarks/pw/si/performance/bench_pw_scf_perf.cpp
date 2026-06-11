/// @file benchmarks/pw/si/performance/bench_pw_scf_perf.cpp
/// @brief Plane-wave SCF time-to-solution: Atomiverse CPU/CUDA vs CP2K.
///
/// Times full SCF runs of bulk Si (2-atom fcc primitive cell,
/// ONCVPSP LDA pseudopotential) at increasing plane-wave cutoffs and
/// emits one CSV row per run. The paired CP2K driver
/// (cp2k/run_perf.sh) runs the same physical system with Quickstep
/// and appends rows with backend=cp2k to the same CSV schema, so the
/// plotters can overlay all three backends.
///
/// Protocol: wall time of the complete SCF to the backend's own
/// convergence (tol 1e-7 on the density), single process. Atomiverse
/// pw is single-threaded; record OMP_NUM_THREADS for CP2K rows. The
/// discretizations differ (plane waves vs GPW) — this benchmark
/// compares time-to-converged-total-energy for the same physical
/// system, not FLOP parity. See ../README.md.
///
/// Usage: bench_pw_scf_perf [--backend cpu|cuda] [--ecut-list "20,30,40"]
///                          [--kmesh N] [--a LATTICE] [--out csv]

#include "atomiverse/core/mat3.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/pw/pw_scf.h"
#ifdef ATV_BUILD_CUDA
#include "atomiverse/cuda/pw/pw_scf.cuh"
#endif

#include <chrono>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace ATV;

namespace {

cpu::Frame make_si_fcc(double a) {
    cpu::Frame f(2);
    for (int i = 0; i < 2; ++i) {
        f.atoms.type[i] = 14;
        f.atoms.mass[i] = 28.0855;
    }
    f.atoms.x[0] = 0;
    f.atoms.y[0] = 0;
    f.atoms.z[0] = 0;
    f.atoms.x[1] = a / 4.0;
    f.atoms.y[1] = a / 4.0;
    f.atoms.z[1] = a / 4.0;
    double h = a / 2.0;
    Mat3 cell;
    cell(0, 0) = 0;
    cell(0, 1) = h;
    cell(0, 2) = h;
    cell(1, 0) = h;
    cell(1, 1) = 0;
    cell(1, 2) = h;
    cell(2, 0) = h;
    cell(2, 1) = h;
    cell(2, 2) = 0;
    f.box = cell;
    return f;
}

} // namespace

int main(int argc, char** argv) {
    std::string backend = "cpu";
    std::string ecut_list = "20,30,40";
    std::string out_csv = "pw_scf_perf.csv";
    std::string upf = std::string(ATV_PW_DATA_DIR) + "/Si_lda_nc.upf";
    int kmesh = 2;
    double a_lat = 5.43;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        auto next = [&]() -> std::string {
            return (i + 1 < argc) ? argv[++i] : "";
        };
        if (arg == "--backend")
            backend = next();
        else if (arg == "--ecut-list")
            ecut_list = next();
        else if (arg == "--kmesh")
            kmesh = std::atoi(next().c_str());
        else if (arg == "--a")
            a_lat = std::atof(next().c_str());
        else if (arg == "--out")
            out_csv = next();
        else if (arg == "--upf")
            upf = next();
    }

    std::vector<double> ecuts;
    {
        std::istringstream ss(ecut_list);
        std::string tok;
        while (std::getline(ss, tok, ','))
            ecuts.push_back(std::atof(tok.c_str()));
    }

    // Append so cpu / cuda / cp2k rows accumulate into one file.
    bool fresh = !std::ifstream(out_csv).good();
    std::ofstream out(out_csv, std::ios::app);
    if (fresh)
        out << "backend,n_atoms,a_angstrom,ecut_ry,k_mesh,threads,"
               "iterations,converged,energy_ry,wall_s,s_per_iter\n";

    for (double ecut : ecuts) {
        auto frame = make_si_fcc(a_lat);
        std::cout << "\n### " << backend << " ecut=" << ecut << " Ry"
                  << std::endl;
        auto t0 = std::chrono::steady_clock::now();
        cpu::pw::PWResult r;
        if (backend == "cuda") {
#ifdef ATV_BUILD_CUDA
            r = cuda::pw::pw_scf(frame, upf, ecut, kmesh, 1e-7, 60);
#else
            std::cerr << "FAIL: built without ATV_BUILD_CUDA" << std::endl;
            return 1;
#endif
        } else {
            r = cpu::pw::pw_scf(frame, upf, ecut, kmesh, 1e-7, 60);
        }
        auto t1 = std::chrono::steady_clock::now();
        double wall = std::chrono::duration<double>(t1 - t0).count();

        char row[256];
        std::snprintf(row, sizeof row,
                      "%s,2,%.2f,%.0f,%d,1,%d,%d,%.8f,%.3f,%.3f\n",
                      backend.c_str(), a_lat, ecut, kmesh, r.iterations,
                      r.converged ? 1 : 0, r.energy, wall,
                      r.iterations > 0 ? wall / r.iterations : 0.0);
        out << row;
        out.flush();
        std::printf("=> %s ecut=%.0f Ry: %.3f s (%d iters, %s)\n",
                    backend.c_str(), ecut, wall, r.iterations,
                    r.converged ? "converged" : "NOT CONVERGED");
    }
    std::cout << "\nCSV appended to " << out_csv << std::endl;
    return 0;
}
