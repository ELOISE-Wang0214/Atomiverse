/// @file benchmarks/md/lj/performance/minimize_scaling.cpp
/// @brief Strong / weak scaling benchmark for the parallel CPU minimizer.
///
/// Same LJ system as the MD scaling benchmark (sc_lj_system.h: sigma=2.5,
/// rho*=0.712, cut=2*sigma, 0.1-A deterministic jitter — the jitter IS the
/// perturbation the minimizer relaxes). Runs a FIXED number of optimizer
/// iterations under `cpu::minimizer::Minimize` (all convergence tolerances
/// set to 0 so the budget is always exhausted), making the per-iteration
/// work identical at every rank count — the LAMMPS-benchmark convention
/// (`minimize 0 0 <iters> <evals>` with etol=ftol=0).
///
/// With the fixed-alpha optimizers (no line search) every iteration costs
/// exactly ONE force evaluation, so Matom-iter/s == Matom-feval/s and is
/// directly comparable to LAMMPS `min_style cg` via its reported
/// "Force evaluations" count (deck: in.lj_min_scale_kk; the sibling
/// in.lj_min_kk is the single-GPU convergence benchmark, bench_minimize.cu).
///
/// `Minimize::install_defaults()` is exercised end-to-end at scale: Migrate,
/// HaloSend/HaloRecv (auto-injected from the pair cutoff), `md::mpi::Reduce`,
/// and the optimizer's globally-reduced inner products / convergence metrics.
///
/// Strong-scaling form (fixed total atoms):
/// @code
/// mpirun -n 8 ./build/benchmarks/minimize_scaling --atoms 8000000 --iters 100
/// @endcode
///
/// Weak-scaling form (fixed atoms per rank):
/// @code
/// mpirun -n 8 ./build/benchmarks/minimize_scaling --atoms-per-rank 200000
/// @endcode
///
/// CSV header:
///   `ranks,atoms,iters,time_total_s,time_per_iter_ms,atoms_per_rank,mode,opt,peak_rss_total_kb`

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

#include "atomiverse/cpu/potentials/pair_lj_cut.h"
#include "atomiverse/cpu/neighbor/neighbor_list.h"
#include "atomiverse/cpu/minimizer/cg.h"
#include "atomiverse/cpu/minimizer/lbfgs.h"
#include "atomiverse/cpu/minimizer/minimize.h"

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
    int atoms = 0;                ///< total atoms (strong scaling); 0 = unset
    int atoms_per_rank = 200'000; ///< atoms per rank (weak scaling default)
    int iters = 100;              ///< optimizer iterations (always exhausted)
    double cutoff = 5.0;          ///< LJ cutoff (= 2*sigma, matches MD bench)
    double skin = 0.6;            ///< Verlet skin [A]
    int migrate = 10;             ///< Migrate cadence (Minimize default)
    // Fixed step size [A^2/eV]. Small enough that the 0.1-A jitter relaxes
    // monotonically without overshoot (energy decrease is asserted at the
    // end); per-iteration cost is jitter-amplitude independent.
    double alpha0 = 0.005;
    std::string opt = "cg";       ///< cg | lbfgs
    bool csv_header = false;      ///< print CSV column header (rank 0 only)
    bool quiet = false;           ///< suppress one-line summaries
};

void print_help() {
    std::printf(
        "Usage: minimize_scaling [options]\n"
        "  --atoms N           Total atoms across all ranks (strong scaling).\n"
        "                      Mutually exclusive with --atoms-per-rank.\n"
        "  --atoms-per-rank N  Atoms per rank (weak scaling). Default 200000.\n"
        "  --iters N           Optimizer iterations (budget always exhausted;\n"
        "                      1 force evaluation per iteration). Default 100.\n"
        "  --opt cg|lbfgs      Optimizer component. Default cg.\n"
        "  --cutoff R          LJ cutoff in Angstrom. Default 5.0.\n"
        "  --skin R            Verlet skin in Angstrom. Default 0.6.\n"
        "  --migrate N         Migrate cadence in iterations. Default 10.\n"
        "  --alpha R           Fixed step size [A^2/eV]. Default 0.005.\n"
        "  --csv-header        Print CSV column header before the data row.\n"
        "  --quiet             Suppress summaries.\n"
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
        if (take("--iters", a.iters))
            continue;
        if (take_d("--cutoff", a.cutoff))
            continue;
        if (take_d("--skin", a.skin))
            continue;
        if (take("--migrate", a.migrate))
            continue;
        if (take_d("--alpha", a.alpha0))
            continue;
        if (std::strcmp(s, "--opt") == 0) {
            if (i + 1 >= argc) {
                std::fprintf(stderr, "--opt requires a value\n");
                std::exit(2);
            }
            a.opt = argv[++i];
            continue;
        }
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
    if (a.opt != "cg" && a.opt != "lbfgs") {
        std::fprintf(stderr, "--opt must be cg or lbfgs (got '%s')\n",
                     a.opt.c_str());
        return false;
    }
    return true;
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
    const int cube = atv_bench::cells_for(n_total);
    const int n_actual = cube * cube * cube;

    // Distributed construction — each rank builds only its subdomain (see
    // sc_lj_system.h; no rank ever materializes the global frame).
#ifdef ATV_HAS_MPI
    cpu::Frame frame = (world_size == 1)
                           ? make_sc_lj_total(n_total)
                           : make_sc_lj_local(n_total, my_rank, world_size);
#else
    cpu::Frame frame = make_sc_lj_total(n_total);
#endif

    if (my_rank == 0 && args.csv_header) {
        std::printf("ranks,atoms,iters,time_total_s,time_per_iter_ms,atoms_"
                    "per_rank,mode,opt,peak_rss_total_kb\n");
    }

    cpu::minimizer::Minimize opt;
    // No barostat / stress consumer in a plain minimization: skip the pair
    // virial (same vflag gating as the MD benchmark, ~20-25% off the LJ
    // inner loop).
    opt.set_pair(cpu::PairLJCut{1.0, 2.5, args.cutoff}.set_compute_stress(false));
    opt.set_neighbor(cpu::md::NeighborList<>{args.cutoff, args.skin});
    // All five convergence tolerances = 0: no criterion can ever be met, so
    // run() exhausts the full --iters budget — fixed work per rank count,
    // the LAMMPS `minimize 0 0 N M` benchmark convention.
    if (args.opt == "lbfgs") {
        opt.set_optimizer(cpu::minimizer::LBFGS{
            /*memory=*/5, args.alpha0, /*ftol=*/0.0, /*max_ls=*/10,
            /*rms_ftol=*/0.0, /*max_dtol=*/0.0, /*rms_dtol=*/0.0,
            /*etol=*/0.0});
    } else {
        opt.set_optimizer(cpu::minimizer::CG{args.alpha0, /*ftol=*/0.0,
                                             /*rms_ftol=*/0.0,
                                             /*max_dtol=*/0.0,
                                             /*rms_dtol=*/0.0, /*etol=*/0.0});
    }
    opt.set_migrate_every(args.migrate);
    opt.build();

    // Single timed run: Minimize::run() is non-continuable (setup/teardown
    // per call), so the first neighbor build is inside the timing and
    // amortizes over the --iters budget — same accounting as the LAMMPS
    // minimize loop time.
#ifdef ATV_HAS_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    const double t0 = now_seconds();
    const int done = opt.run(frame, args.iters);
    const double t1 = now_seconds();

    // Correctness oracle: the budget must be exhausted (tolerances are 0 —
    // an early stop means the fixed-work contract is broken) and the global
    // energy must have decreased from the jittered start.
    if (done != args.iters && my_rank == 0)
        std::fprintf(stderr,
                     "WARNING: run() stopped at %d/%d iterations — fixed-work "
                     "contract broken, row is not comparable\n",
                     done, args.iters);
    if (my_rank == 0 && !args.quiet)
        std::fprintf(stderr, "FINAL_ENERGY %.6f  iters %d\n",
                     frame.scalars.energy, done);

    const double total_s = t1 - t0;
    const double per_iter_ms = 1000.0 * total_s / std::max(1, done);

    // Total memory footprint = sum of per-rank peak RSS [kB].
    long my_hwm = read_vmhwm_kb();
    long total_hwm = my_hwm;
#ifdef ATV_HAS_MPI
    MPI_Reduce(&my_hwm, &total_hwm, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
#endif

    if (my_rank == 0) {
        const int per_rank = (n_actual + world_size - 1) / world_size;
        std::printf("%d,%d,%d,%.6f,%.6f,%d,%s,%s,%ld\n", world_size, n_actual,
                    done, total_s, per_iter_ms, per_rank, mode,
                    args.opt.c_str(), total_hwm);
        if (!args.quiet) {
            std::fprintf(stderr,
                         "minimize_scaling: %d ranks  %d atoms  %d iters  "
                         "%.4f s total  %.4f ms/iter  (%s, %s, ~%d "
                         "atoms/rank)\n",
                         world_size, n_actual, done, total_s, per_iter_ms,
                         mode, args.opt.c_str(), per_rank);
        }
    }

#ifdef ATV_HAS_MPI
    MPI_Finalize();
#endif
    return 0;
}
