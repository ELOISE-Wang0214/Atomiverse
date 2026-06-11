/// @file examples/lj.cu
/// @brief Complete LJ simulation: minimize -> NVE MD -> RDF + trajectory + thermo.
///
/// Demonstrates the full Atomiverse workflow on a Lennard-Jones argon system:
///   1. Build an FCC lattice of LJ particles
///   2. Energy minimization (steepest descent)
///   3. NVE molecular dynamics with Maxwell-Boltzmann initial velocities
///   4. Output: thermo (stdout + thermo.dat), XYZ trajectory (traj.xyz), RDF (rdf.dat)
///   5. Optional: energy + RDF plots (when built with -DATV_USE_MATPLOT=ON)
///
/// Usage:
///   ./lj
///
/// Output files:
///   traj.xyz    -- XYZ trajectory (every 100 steps)
///   rdf.dat     -- Radial distribution function g(r)
///   thermo.dat  -- Thermodynamic log (Step Temp PE KE TotEng)
///   energy.png  -- Energy time series plot (if matplotplusplus enabled)
///   rdf.png     -- RDF g(r) plot (if matplotplusplus enabled)

#include <cmath>
#include <iomanip>
#include <iostream>

#include "atomiverse/cuda/md/init.cuh"
#include "atomiverse/cuda/md/md.cuh"
#include "atomiverse/cuda/md/velocity_verlet.cuh"
#include "atomiverse/cuda/minimizer/sd.cuh"
#include "atomiverse/cuda/neighbor/neighbor_list.cuh"
#include "atomiverse/cuda/observer/rdf.cuh"
#include "atomiverse/cuda/potentials/pair_lj_cut.cuh"
#include "atomiverse/cuda/recorder/recorder.cuh"
#include "atomiverse/cuda/recorder/thermo.cuh"

#ifdef ATV_HAS_MATPLOT
#include <fstream>
#include <matplot/matplot.h>
#include <sstream>
#include <vector>
#endif

using namespace ATV;

/// Build an FCC lattice of n_cells^3 * 4 atoms.
cpu::Frame make_fcc(int n_cells, double lattice_const) {
    int n = n_cells * n_cells * n_cells * 4;
    double box_len = lattice_const * double(n_cells);

    cpu::Frame frame(n);

    // FCC basis: (0,0,0), (0.5,0.5,0), (0.5,0,0.5), (0,0.5,0.5)
    const double basis[4][3] = {
        {double(0), double(0), double(0)},
        {double(0.5), double(0.5), double(0)},
        {double(0.5), double(0), double(0.5)},
        {double(0), double(0.5), double(0.5)},
    };

    int idx = 0;
    for (int ix = 0; ix < n_cells; ++ix)
        for (int iy = 0; iy < n_cells; ++iy)
            for (int iz = 0; iz < n_cells; ++iz)
                for (int b = 0; b < 4; ++b) {
                    frame.atoms.x[idx] =
                        (double(ix) + basis[b][0]) * lattice_const;
                    frame.atoms.y[idx] =
                        (double(iy) + basis[b][1]) * lattice_const;
                    frame.atoms.z[idx] =
                        (double(iz) + basis[b][2]) * lattice_const;
                    frame.atoms.mass[idx] = double(39.948); // Argon [amu]
                    frame.atoms.type[idx] = 0;
                    idx++;
                }

    frame.box = Mat3::cubic(static_cast<double>(box_len));
    return frame;
}

#ifdef ATV_HAS_MATPLOT
/// Read a space-delimited data file, skipping comment lines starting with '#'.
/// Returns columns as vectors of doubles.
struct DataTable {
    std::vector<std::vector<double>> cols;
    int n_cols = 0;
    int n_rows = 0;
};

DataTable read_data(const std::string& filename) {
    DataTable table;
    std::ifstream in(filename);
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#')
            continue;
        std::istringstream iss(line);
        double val;
        int col = 0;
        while (iss >> val) {
            if (col >= static_cast<int>(table.cols.size()))
                table.cols.emplace_back();
            table.cols[col].push_back(val);
            col++;
        }
        if (col > 0) {
            table.n_cols = col;
            table.n_rows++;
        }
    }
    return table;
}

void plot_energy(const std::string& thermo_file) {
    auto data = read_data(thermo_file);
    if (data.n_cols < 5 || data.n_rows < 2) {
        std::cout << "Skipping energy plot: insufficient data in "
                  << thermo_file << "\n";
        return;
    }

    // Columns: Step(0) Temp(1) PotEng(2) KinEng(3) TotEng(4)
    auto& steps = data.cols[0];
    auto& pe = data.cols[2];
    auto& ke = data.cols[3];
    auto& total = data.cols[4];

    auto fig = matplot::figure(true);
    matplot::hold(matplot::on);
    auto p1 = matplot::plot(steps, pe);
    p1->display_name("PE");
    auto p2 = matplot::plot(steps, ke);
    p2->display_name("KE");
    auto p3 = matplot::plot(steps, total);
    p3->display_name("Total");
    matplot::legend();
    matplot::xlabel("Step");
    matplot::ylabel("Energy [eV]");
    matplot::title("Energy Time Series");
    matplot::save("energy.png");
    std::cout << "Saved energy.png\n";
}

void plot_rdf(const std::string& rdf_file) {
    auto data = read_data(rdf_file);
    if (data.n_cols < 2 || data.n_rows < 2) {
        std::cout << "Skipping RDF plot: insufficient data in " << rdf_file
                  << "\n";
        return;
    }

    // Columns: r(0) g(r)(1)
    auto& r = data.cols[0];
    auto& gr = data.cols[1];

    auto fig = matplot::figure(true);
    auto p = matplot::plot(r, gr);
    p->line_width(1.5);
    matplot::xlabel("r [Angstrom]");
    matplot::ylabel("g(r)");
    matplot::title("Radial Distribution Function");
    matplot::save("rdf.png");
    std::cout << "Saved rdf.png\n";
}
#endif // ATV_HAS_MATPLOT

int main() {
    // ── System setup ──────────────────────────────────────────────────
    // LJ argon: epsilon = 0.0104 eV, sigma = 3.4 A
    // Cutoff = 2.5*sigma = 8.5 A
    // FCC lattice constant for Ar ~ 5.26 A
    double epsilon = double(0.0104);
    double sigma = double(3.4);
    double cutoff = double(6.8); // 2.0*sigma (1.3*cutoff < L/2)
    int n_cells = 5;             // 5x5x5 FCC -> 500 atoms
    double a0 = double(5.26);    // lattice constant [A]

    auto frame = make_fcc(n_cells, a0);
    int n = frame.atoms.size();
    std::cout << "System: " << n << " Ar atoms, FCC lattice a=" << std::fixed
              << std::setprecision(2) << static_cast<double>(a0)
              << " A, box=" << static_cast<double>(a0 * double(n_cells))
              << " A\n";

    // ── Phase 1: Energy minimization ──────────────────────────────────
    std::cout
        << "\n=== Energy Minimization (Steepest Descent, 500 steps) ===\n";
    {
        minimizer::SD min;
        min.add(minimizer::SteepestDescentStep{double(0.01)});
        min.on(minimizer::Phase::force)
            .add(NeighborList{cutoff}, -10)
            .add(LJ{epsilon, sigma, cutoff});
        min.on(minimizer::Phase::end_of_step).add(Thermo{100});
        min.build();
        min.run(frame, 500);
    }
    std::cout << std::fixed << std::setprecision(6)
              << "Final PE after minimization: "
              << static_cast<double>(frame.energy) << " eV ("
              << static_cast<double>(frame.energy) / n << " eV/atom)\n";

    // ── Phase 2: NVE Molecular Dynamics ───────────────────────────────
    // NOTE: Using small dt and low T for stability.
    // The cluster-based NeighborList currently produces duplicate periodic
    // pairs, inflating energies. Absolute energy values are not physical
    // until the NeighborList deduplication is fixed. Energy *conservation*
    // (constant TotEng) can still be verified at small dt.
    double dt = double(0.01);         // 0.01 fs timestep
    double temperature = double(5.0); // 5 K initial temperature
    int n_steps = 1000;
    int thermo_freq = 100;
    int dump_freq = 100;

    std::cout << "\n=== NVE MD (" << n_steps
              << " steps, dt=" << std::setprecision(1)
              << static_cast<double>(dt)
              << " fs, T_init=" << static_cast<double>(temperature)
              << " K) ===\n";
    {
        md::MD sim;
        sim.add(md::InitVelocity{temperature, 12345});
        sim.add(md::VelocityVerlet{dt});
        sim.on(md::Phase::force)
            .add(NeighborList{cutoff}, -10)
            .add(LJ{epsilon, sigma, cutoff});
        sim.on(md::Phase::end_of_step)
            .add(Thermo{thermo_freq, "thermo.dat"})
            .add(Recorder{"traj.xyz", dump_freq})
            .add(RDF{cutoff, 200, "rdf.dat"});

        sim.enable_timing();
        sim.build();
        sim.run(frame, n_steps);

        std::cout << "\n" << sim.summary();
    }

    std::cout << "\n=== Done ===\n"
              << std::setprecision(6)
              << "Final PE: " << static_cast<double>(frame.energy)
              << " eV, KE: " << static_cast<double>(frame.kinetic)
              << " eV, Total: "
              << static_cast<double>(frame.energy + frame.kinetic) << " eV\n"
              << "Output: traj.xyz, rdf.dat, thermo.dat\n";

    // ── Plotting (optional, requires -DATV_USE_MATPLOT=ON) ────────────
#ifdef ATV_HAS_MATPLOT
    std::cout << "\n=== Generating plots ===\n";
    plot_energy("thermo.dat");
    plot_rdf("rdf.dat");
#endif

    return 0;
}
