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
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

#include "atomiverse/cpu/potentials/pair_lj_cut.h"
#include "atomiverse/cpu/md/md.h"
#include "atomiverse/cpu/comm/atom_migration.h"
#include "atomiverse/cpu/comm/halo_exchange.h"
#include "atomiverse/cpu/neighbor/neighbor_list.h"
#include "atomiverse/cpu/md/velocity_verlet.h"

#include "sc_lj_system.h"

#ifdef ATV_HAS_MPI
#include <mpi.h>
#endif

using namespace ATV;
using atv_bench::make_sc_lj_total;
using atv_bench::now_seconds;
using atv_bench::read_vmhwm_kb;
#ifdef ATV_HAS_MPI
using atv_bench::make_sc_lj_local;
#endif

namespace {

struct Args {
    int atoms = 0;               ///< total atoms (strong scaling); 0 = unset
    int atoms_per_rank = 50'000; ///< atoms per rank (weak scaling default)
    int steps = 200;             ///< MD steps to time
    int warmup = 10;             ///< unrecorded warmup steps before timing
    double cutoff = 5.0;         ///< LJ cutoff (= 2*sigma, matches GPU bench)
    double skin = 0.6;           ///< Verlet skin [A] — matches LAMMPS lj.in
    int reneigh = 20;            ///< reneighbor cadence (steps) — matches LAMMPS
    int migrate = 0;             ///< migrate cadence (0 = follow reneigh)
    // Time step. sigma=2.5, eps=1, m=1 → LJ time unit tau = sigma*sqrt(m/eps)
    // = 2.5. LAMMPS lj.in uses dt* = 0.005 (reduced); the matching value here is
    // 0.005*tau = 0.0125. dt=1.0 (dt/tau=0.4) was ~80x too large — the velocity-
    // Verlet integrator was unstable and the system blew up (single-rank too),
    // which is unphysical and also corrupts the throughput comparison (a blown-up
    // system loses neighbors and runs artificially fast).
    double dt = 0.0125;
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
        if (take("--reneigh", a.reneigh))
            continue;
        if (take("--migrate", a.migrate))
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

// System construction + timing/RSS helpers live in sc_lj_system.h (shared
// with minimize_scaling.cpp — one source of truth for the benchmark lattice).
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

    // Global atom count (cube-rounded, identical to make_sc_lj_total).
    const int cube = std::max(
        1, static_cast<int>(std::cbrt(static_cast<double>(n_total)) + 0.5));
    const int n_actual = cube * cube * cube;

    // Distributed construction: each rank builds ONLY its own subdomain's atoms
    // (make_sc_lj_local), so no rank ever materializes the global frame. The
    // old "rank 0 builds global + scatter" peaked at O(global-N) on rank 0
    // (5.9 GB @256 ranks weak, 23 GB @1024 — the dominant per-rank RSS, and an
    // eventual OOM). Single-rank keeps the simple full build.
#ifdef ATV_HAS_MPI
    cpu::Frame frame = (world_size == 1)
                           ? make_sc_lj_total(n_total)
                           : make_sc_lj_local(n_total, my_rank, world_size);
#else
    cpu::Frame frame = make_sc_lj_total(n_total);
#endif

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
    // Re-home drifted atoms + reneighbor every `reneigh` steps, matching the
    // LAMMPS `neigh_modify every 20` cadence the comparison uses. Migrate shares
    // the cadence (a reneighbor step is also an exchange step). Between these,
    // the halo only refreshes ghost positions and the neighbor list is reused —
    // so the per-step cost is the LJ kernel, not a rebuild.
    sim.add(cpu::comm::Migrate{args.migrate > 0 ? args.migrate : args.reneigh});
    // HaloSend cutoff = pair cutoff + Verlet skin (LAMMPS comm rule): ghosts
    // must cover atoms that drift into the cutoff band before the next
    // reneighbor, else a real pair is missed and the force blows up. Explicit
    // halo here means MD::install_defaults skips its auto-inject (which would
    // otherwise add a SECOND HaloRecv — see fix d8ea6a5: the reverse exchange
    // must run exactly once or it double-applies ghost forces and leaks
    // momentum).
    sim.on(ATV::md::Phase::pre_force)
        .add(cpu::comm::HaloSend{args.cutoff + args.skin, args.reneigh}, -20);
    sim.on(ATV::md::Phase::post_force)
        .add(cpu::comm::HaloRecv{}, -1000);
    // Adaptive Verlet rebuild (LAMMPS `neigh_modify check yes`): re-enumerate
    // the list only when an atom has drifted more than skin/2, not on a fixed
    // cadence. The per-step skin scan is far cheaper than a full rebuild, so
    // for systems that move less than skin/2 over `reneigh` steps this does
    // fewer rebuilds than a fixed `every N` would — a legitimate AV speedup
    // over the deck's `every 20 check no`, kept because it is a general,
    // physics-correct optimisation (a stale list is impossible by construction).
    sim.on(ATV::md::Phase::force)
        .add(cpu::md::NeighborList<>{args.cutoff, args.skin}, -10);
    // Plain NVE: no barostat / pressure consumer, so skip the pair virial
    // (LAMMPS-style vflag gating) — ~20-25% off the LJ inner loop.
    sim.on(ATV::md::Phase::force)
        .add(cpu::PairLJCut{1.0, 2.5, args.cutoff}.set_compute_stress(false));
    // Match the LAMMPS deck, which has no `thermo N` → it computes/reduces
    // global energy only at run start + end, NOT every step. AV's default
    // reduces an 18-double MPI_Allreduce EVERY step; that latency-bound
    // collective is invisible at low rank counts but dominates the per-step
    // cost at 1024 ranks (it is the documented weak-scaling-efficiency leak).
    // This run is plain NVE with no per-step consumer of the global scalars
    // (the only reader is the end-of-run correctness oracle, which does its
    // own final reduce below), so align the reduce cadence with the
    // reneighbor stride. General NVE optimisation, not benchmark-specific.
    sim.set_scalar_reduce_every(args.reneigh);
    // Opt-in per-component timing (ATV_TIMER[_VERBOSE]) — splits the per-step
    // wall into force / neighbor / integrate / Comm(halo+migrate) to expose
    // the multi-rank scaling bottleneck. Harmless when the env var is unset.
    const bool do_timing = std::getenv("ATV_TIMER") != nullptr;
    if (do_timing)
        sim.enable_timing(true);
    sim.build();

    if (args.warmup > 0)
        sim.run(frame, args.warmup);
    if (do_timing)
        sim.reset_timing();

#ifdef ATV_HAS_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    const double t0 = now_seconds();
    sim.run(frame, args.steps);
#ifdef ATV_HAS_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    const double t1 = now_seconds();

    if (do_timing)
        sim.emit_timing("bench");

    // Correctness oracle: frame.scalars.energy is the cross-rank-reduced total
    // LJ potential energy (filled by the Reduce component) on rank 0 — the same
    // field tests/mpi/test_lj_energy_bruteforce.cpp checks. Compare across
    // reneigh cadences to confirm the position-only halo conserves energy.
    if (do_timing) {
        // Total linear momentum over OWNED atoms, reduced across ranks. For NVE
        // it must stay constant — a drift exposes a Newton's-3rd-law violation
        // in the multi-rank force path (reverse-halo / dedup). The authoritative
        // owned count is dd->n_local (frame.atoms.nlocal is left at global N).
        int nloc = frame.atoms.size();
#ifdef ATV_HAS_MPI
        auto dd_probe = sim.release_dd();
        if (world_size > 1 && dd_probe)
            nloc = dd_probe->n_local;
#endif
        double p[3] = {0, 0, 0};
        for (int i = 0; i < nloc; ++i) {
            const double m = frame.atoms.mass[i];
            p[0] += m * frame.atoms.vx[i];
            p[1] += m * frame.atoms.vy[i];
            p[2] += m * frame.atoms.vz[i];
        }
        double pg[3] = {p[0], p[1], p[2]};
#ifdef ATV_HAS_MPI
        MPI_Reduce(p, pg, 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#endif
        const double pmag =
            std::sqrt(pg[0] * pg[0] + pg[1] * pg[1] + pg[2] * pg[2]);
        if (my_rank == 0)
            std::fprintf(stderr, "FINAL_ENERGY %.6f  |P| %.6e  nloc0 %d\n",
                         frame.scalars.energy, pmag, nloc);
    }

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
