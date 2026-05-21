/// @file benchmarks/md/scaling.cpp
/// @brief Strong / weak scaling benchmark for the parallel CPU MD engine.
///
/// Builds an LJ argon-like FCC system, runs `--steps` MD steps under
/// `cpu::md::MD`, and prints a one-line CSV row containing wall-clock time
/// per step. Designed to be invoked at multiple rank counts via `mpirun`.
///
/// Strong-scaling form (fixed total atoms):
/// @code
/// mpirun -n 1 ./build/benchmarks/md/scaling --atoms 50000 --steps 200
/// mpirun -n 8 ./build/benchmarks/md/scaling --atoms 50000 --steps 200
/// @endcode
///
/// Weak-scaling form (fixed atoms per rank):
/// @code
/// mpirun -n 1 ./build/benchmarks/md/scaling --atoms-per-rank 50000 --steps 200
/// mpirun -n 8 ./build/benchmarks/md/scaling --atoms-per-rank 50000 --steps 200
/// @endcode
///
/// CSV header:
///   `ranks,atoms,steps,time_total_s,time_per_step_ms,atoms_per_rank,mode`
///
/// @see .claude/specs/parallel-md-engine.md  § Phase 4 / Verification §4

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include "atomiverse/core/mat3.h"
#include "atomiverse/cpu/frame.h"
#include "atomiverse/cpu/md/lj.h"
#include "atomiverse/cpu/md/md.h"
#include "atomiverse/cpu/md/mpi/halo_exchange.h"
#include "atomiverse/cpu/md/neighbor_list.h"
#include "atomiverse/cpu/md/velocity_verlet.h"

#ifdef ATV_HAS_MPI
#include <mpi.h>
#else
#include <chrono>
#endif

using namespace ATV;

namespace {

struct Args {
    int atoms = 0;               ///< total atoms (strong scaling); 0 = unset
    int atoms_per_rank = 50'000; ///< atoms per rank (weak scaling default)
    int steps = 200;             ///< MD steps to time
    int warmup = 10;             ///< unrecorded warmup steps before timing
    double cutoff = 6.0;         ///< LJ cutoff [Angstrom]
    double dt = 1.0;             ///< Time step [fs]
    bool csv_header = false;     ///< print CSV header to stdout (rank 0 only)
    bool quiet = false;          ///< suppress one-line per-rank reports
};

void print_help() {
    std::printf(
        "Usage: scaling [options]\n"
        "  --atoms N           Total atoms across all ranks (strong scaling).\n"
        "                      Mutually exclusive with --atoms-per-rank.\n"
        "  --atoms-per-rank N  Atoms per rank (weak scaling). Default 50000.\n"
        "  --steps N           Timed MD steps. Default 200.\n"
        "  --warmup N          Warmup steps before timing. Default 10.\n"
        "  --cutoff R          LJ cutoff in Angstrom. Default 6.0.\n"
        "  --dt T              Time step in fs. Default 1.0.\n"
        "  --csv-header        Print CSV column header before the data row.\n"
        "  --quiet             Suppress per-rank summaries.\n"
        "  --help              Show this message.\n");
}

bool parse_args(int argc, char** argv, Args& a) {
    for (int i = 1; i < argc; ++i) {
        const char* s = argv[i];
        auto take = [&](const char* name, int& dst) {
            if (std::strcmp(s, name) != 0)
                return false;
            if (i + 1 >= argc) {
                std::fprintf(stderr, "%s requires a value\n", name);
                std::exit(2);
            }
            dst = std::atoi(argv[++i]);
            return true;
        };
        auto take_d = [&](const char* name, double& dst) {
            if (std::strcmp(s, name) != 0)
                return false;
            if (i + 1 >= argc) {
                std::fprintf(stderr, "%s requires a value\n", name);
                std::exit(2);
            }
            dst = std::atof(argv[++i]);
            return true;
        };
        if (take("--atoms", a.atoms))
            continue;
        if (take("--atoms-per-rank", a.atoms_per_rank))
            continue;
        if (take("--steps", a.steps))
            continue;
        if (take("--warmup", a.warmup))
            continue;
        if (take_d("--cutoff", a.cutoff))
            continue;
        if (take_d("--dt", a.dt))
            continue;
        if (std::strcmp(s, "--csv-header") == 0) {
            a.csv_header = true;
            continue;
        }
        if (std::strcmp(s, "--quiet") == 0) {
            a.quiet = true;
            continue;
        }
        if (std::strcmp(s, "--help") == 0) {
            print_help();
            std::exit(0);
        }
        std::fprintf(stderr, "unknown argument: %s\n", s);
        print_help();
        return false;
    }
    return true;
}

// Build an FCC argon lattice with `n_target` total atoms (rounded to the
// nearest cube), a = 5.3 Å (matches the LAMMPS-typical Ar density for σ=3.4).
cpu::Frame make_fcc_argon_total(int n_target) {
    const int atoms_per_unit_cell = 4;
    int n_cells_per_axis =
        std::max(1, static_cast<int>(std::cbrt(static_cast<double>(n_target)
                                               / atoms_per_unit_cell)));
    while (atoms_per_unit_cell * n_cells_per_axis * n_cells_per_axis
               * n_cells_per_axis
           < n_target) {
        ++n_cells_per_axis;
    }
    const double a_lat = 5.3;
    const int n = atoms_per_unit_cell * n_cells_per_axis * n_cells_per_axis
                  * n_cells_per_axis;
    cpu::Frame frame(n);
    const double basis[4][3] = {
        {0.0, 0.0, 0.0},
        {0.5, 0.5, 0.0},
        {0.5, 0.0, 0.5},
        {0.0, 0.5, 0.5},
    };
    const double shift = 0.25 * a_lat;
    int idx = 0;
    for (int ix = 0; ix < n_cells_per_axis; ++ix)
        for (int iy = 0; iy < n_cells_per_axis; ++iy)
            for (int iz = 0; iz < n_cells_per_axis; ++iz)
                for (int b = 0; b < 4; ++b) {
                    frame.atoms.x[idx] = (ix + basis[b][0]) * a_lat + shift;
                    frame.atoms.y[idx] = (iy + basis[b][1]) * a_lat + shift;
                    frame.atoms.z[idx] = (iz + basis[b][2]) * a_lat + shift;
                    frame.atoms.mass[idx] = 39.948;
                    frame.atoms.type[idx] = 0;
                    const double s = 1e-4 * ((idx % 2 == 0) ? 1.0 : -1.0);
                    frame.atoms.vx[idx] = s;
                    frame.atoms.vy[idx] = s * 0.5;
                    frame.atoms.vz[idx] = s * 0.25;
                    ++idx;
                }
    frame.box = Mat3::cubic(n_cells_per_axis * a_lat);
    return frame;
}

double now_seconds() {
#ifdef ATV_HAS_MPI
    return MPI_Wtime();
#else
    using clk = std::chrono::steady_clock;
    static const auto t0 = clk::now();
    const auto t = clk::now();
    return std::chrono::duration<double>(t - t0).count();
#endif
}

} // namespace

int main(int argc, char** argv) {
    Args args;
    if (!parse_args(argc, argv, args))
        return 2;

#ifdef ATV_HAS_MPI
    int provided = 0;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    int world_size = 1, my_rank = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#else
    int world_size = 1, my_rank = 0;
#endif

    // Determine atom count and mode.
    const char* mode = "weak";
    int n_total = args.atoms;
    if (n_total <= 0) {
        n_total = args.atoms_per_rank * world_size; // weak scaling
    } else {
        mode = "strong";
    }

    cpu::Frame frame = make_fcc_argon_total(n_total);
    const int n_actual = frame.atoms.size();

    if (my_rank == 0 && args.csv_header) {
        std::printf("ranks,atoms,steps,time_total_s,time_per_step_ms,atoms_per_"
                    "rank,mode\n");
    }

    cpu::md::MD sim;
    sim.add(cpu::md::VelocityVerlet{args.dt});
    sim.on(ATV::md::Phase::pre_force)
        .add(cpu::md::mpi::HaloExchangeForward{args.cutoff}, -20);
    sim.on(ATV::md::Phase::post_force)
        .add(cpu::md::mpi::HaloExchangeReverse{}, -1000);
    sim.on(ATV::md::Phase::force)
        .add(cpu::md::NeighborList<>{args.cutoff}, -10);
    sim.on(ATV::md::Phase::force).add(cpu::md::LJ{0.0104, 3.4, args.cutoff});
    sim.build();

    if (args.warmup > 0)
        sim.run(frame, args.warmup);

#ifdef ATV_HAS_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    const double t0 = now_seconds();
    sim.run(frame, args.steps);
#ifdef ATV_HAS_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    const double t1 = now_seconds();

    const double total_s = t1 - t0;
    const double per_step_ms = 1000.0 * total_s / std::max(1, args.steps);

    if (my_rank == 0) {
        const int per_rank = (n_actual + world_size - 1) / world_size;
        std::printf("%d,%d,%d,%.6f,%.6f,%d,%s\n", world_size, n_actual,
                    args.steps, total_s, per_step_ms, per_rank, mode);
        if (!args.quiet) {
            std::fprintf(stderr,
                         "scaling: %d ranks  %d atoms  %d steps  %.4f s total  "
                         "%.4f ms/step  (%s, ~%d atoms/rank)\n",
                         world_size, n_actual, args.steps, total_s, per_step_ms,
                         mode, per_rank);
        }
    }

#ifdef ATV_HAS_MPI
    MPI_Finalize();
#endif
    return 0;
}
