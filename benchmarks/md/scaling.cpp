/// @file benchmarks/md/scaling.cpp
/// @brief Strong / weak scaling benchmark for the parallel CPU MD engine.
///
/// Builds a simple-cubic LJ system (sigma=2.5, rho*=0.712, cut=2*sigma —
/// identical to the GPU benchmark bench_multi_gpu_scaling so CPU and GPU
/// panels share one system), runs `--steps` MD steps under `cpu::md::MD`, and
/// prints a one-line CSV row with wall-clock time per step. Designed to be
/// invoked at multiple rank counts via `mpirun`.
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
///   `ranks,atoms,steps,time_total_s,time_per_step_ms,atoms_per_rank,mode,peak_rss_total_kb`
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
#include "atomiverse/cpu/comm/atom_migration.h"
#include "atomiverse/cpu/comm/halo_exchange.h"
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
    double cutoff = 5.0;         ///< LJ cutoff (= 2*sigma, matches GPU bench)
    double skin = 0.6;           ///< Verlet skin [A] — matches LAMMPS lj.in
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
        "  --cutoff R          LJ cutoff in Angstrom. Default 5.0.\n"
        "  --skin R            Verlet skin in Angstrom. Default 0.6.\n"
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
        if (take_d("--skin", a.skin))
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

// Build a simple-cubic LJ system with `n_target` total atoms (rounded to the
// nearest cube), IDENTICAL to the GPU benchmark (bench_multi_gpu_scaling):
// sigma = 2.5, spacing = 1.12*sigma -> reduced density rho* = 0.712, so the
// CPU and GPU panels of the 2x2 study run the same LJ system. Mild
// deterministic jitter + thermal velocities (not a frozen lattice).
cpu::Frame make_sc_lj_total(int n_target) {
    const double sigma = 2.5;
    const double spacing = 1.12 * sigma; // rho* = (1/1.12)^3 = 0.712
    int cells = std::max(
        1, static_cast<int>(std::cbrt(static_cast<double>(n_target)) + 0.5));
    const int n = cells * cells * cells;
    cpu::Frame frame(n);
    int idx = 0;
    for (int ix = 0; ix < cells; ++ix)
        for (int iy = 0; iy < cells; ++iy)
            for (int iz = 0; iz < cells; ++iz) {
                frame.atoms.x[idx] =
                    (ix + 0.5) * spacing + 0.1 * std::sin(0.7 * idx);
                frame.atoms.y[idx] =
                    (iy + 0.5) * spacing + 0.1 * std::sin(1.3 * idx);
                frame.atoms.z[idx] =
                    (iz + 0.5) * spacing + 0.1 * std::sin(1.9 * idx);
                frame.atoms.mass[idx] = 1.0;
                frame.atoms.type[idx] = 0;
                frame.atoms.vx[idx] = 0.01 * std::sin(0.5 * idx);
                frame.atoms.vy[idx] = 0.01 * std::sin(0.9 * idx);
                frame.atoms.vz[idx] = 0.01 * std::sin(1.1 * idx);
                ++idx;
            }
    frame.box = Mat3::cubic(cells * spacing);
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

// Peak resident set size of this process so far [kB] (Linux /proc/self/status
// VmHWM). Summed across ranks this is the benchmark's total memory footprint.
long read_vmhwm_kb() {
    FILE* f = std::fopen("/proc/self/status", "r");
    if (!f)
        return 0;
    char line[256];
    long kb = 0;
    while (std::fgets(line, sizeof(line), f))
        if (std::sscanf(line, "VmHWM: %ld kB", &kb) == 1)
            break;
    std::fclose(f);
    return kb;
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

    cpu::Frame frame = make_sc_lj_total(n_total);
    const int n_actual = frame.atoms.size();

    if (my_rank == 0 && args.csv_header) {
        std::printf("ranks,atoms,steps,time_total_s,time_per_step_ms,atoms_per_"
                    "rank,mode,peak_rss_total_kb\n");
    }

    cpu::md::MD sim;
    sim.add(cpu::md::VelocityVerlet{args.dt});
    // Re-home atoms that drift across subdomain boundaries (pre_exchange).
    // Without this the decomposition corrupts over a run and the neighbor
    // list degrades toward O(N^2). Migrate also bumps neighbor_generation,
    // which triggers the halo full-rebuild on migration steps.
    sim.add(cpu::comm::Migrate{20});
    // HaloSend `every=1` forces a full halo rebuild (count exchange) every
    // step. This is a WORKAROUND for a real bug: the halo's position-only
    // refresh path (cached per-face counts between full rebuilds) desyncs the
    // ranks on a bondless run — they drift apart until a Waitall can't be
    // matched and they deadlock (PE-melt masks it; pure LJ exposes it). Forcing
    // per-step full rebuild keeps the counts fresh and the ranks in lockstep.
    // The proper fix belongs in halo_exchange.cpp's position-only branch.
    sim.on(ATV::md::Phase::pre_force)
        .add(cpu::comm::HaloSend{args.cutoff, 1}, -20);
    sim.on(ATV::md::Phase::post_force)
        .add(cpu::comm::HaloRecv{}, -1000);
    sim.on(ATV::md::Phase::force)
        .add(cpu::md::NeighborList<>{args.cutoff, args.skin}, -10);
    // Plain NVE: no barostat / pressure consumer, so skip the pair virial
    // (LAMMPS-style vflag gating) — ~20-25% off the LJ inner loop.
    sim.on(ATV::md::Phase::force)
        .add(cpu::md::LJ{1.0, 2.5, args.cutoff}.set_compute_stress(false));
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

    // Total memory footprint = sum of per-rank peak RSS [kB].
    long my_hwm = read_vmhwm_kb();
    long total_hwm = my_hwm;
#ifdef ATV_HAS_MPI
    MPI_Reduce(&my_hwm, &total_hwm, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
#endif

    if (my_rank == 0) {
        const int per_rank = (n_actual + world_size - 1) / world_size;
        std::printf("%d,%d,%d,%.6f,%.6f,%d,%s,%ld\n", world_size, n_actual,
                    args.steps, total_s, per_step_ms, per_rank, mode, total_hwm);
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
