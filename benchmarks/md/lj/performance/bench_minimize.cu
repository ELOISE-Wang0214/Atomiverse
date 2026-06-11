/// @file benchmarks/bench_minimize.cu
/// @brief Single-GPU LJ energy-minimization benchmark (AV vs LAMMPS).
///
/// Builds a jittered fcc reduced-LJ crystal (rho* = 1.0, so the
/// nearest-neighbour distance a/sqrt(2) = 1.1225 sigma sits at the LJ
/// pair minimum; sigma = epsilon = 1, rc = 2.0 shifted, skin = 0.6) and
/// relaxes it back to the crystal with a CUDA Optimizer component (CG
/// or LBFGS) under a force-only stopping rule (the other four
/// convergence tolerances are disabled with huge values). Emits one
/// CSV line. The MD throughput deck's sc rho* = 0.712 lattice is NOT
/// used here: sc is mechanically unstable for LJ, so minimizing it
/// triggers a glassy structural collapse and no engine reaches a tight
/// force tolerance — fcc makes time-to-ftol well-defined.
///
/// Optionally writes the initial configuration as a LAMMPS data file so
/// a LAMMPS `minimize` run (deck in.lj_min_kk, with `pair_modify shift
/// yes` to match the cutoff-shifted cuda::PairLJCut) starts from the IDENTICAL
/// coordinates — initial energies must then agree to output precision.
///
/// Usage:
///   bench_minimize [n_target] [style] [ftol] [maxiter] [jitter] [data_out] [alpha0]
///     n_target  approx atom count, rounded to 4 x cells^3 fcc (default 256000)
///     style     cg | lbfgs (default lbfgs)
///     ftol      max-per-atom-|F| stopping tolerance (default 1e-4)
///     maxiter   iteration budget (default 5000)
///     jitter    per-coordinate displacement amplitude, full width (default 0.1)
///     data_out  LAMMPS data file to write ("" = skip; default skip)
///     alpha0    optimizer step scale (default: cg 0.002, lbfgs 1.0)
///
/// CSV: CSV,AV,<style>,<atoms>,<iters>,<force_evals>,<wall_s>,
///      <E0_total>,<Ef_total>,<maxF_final>,<converged>
/// One driver iteration = one force evaluation (line-search rollbacks
/// consume a driver iteration each), so iters == force_evals.

#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

#include "atomiverse/core/mat3.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cuda/minimizer/cg.cuh"
#include "atomiverse/cuda/minimizer/lbfgs.cuh"
#include "atomiverse/cuda/minimizer/minimize.cuh"
#include "atomiverse/cuda/neighbor/neighbor_list.cuh"
#include "atomiverse/cuda/potentials/pair_lj_cut.cuh"
#include "atomiverse/facade/stage.h"

using namespace ATV;

namespace {

/// splitmix64 — deterministic per-atom hash, identical on every run.
uint64_t splitmix64(uint64_t x) {
    x += 0x9e3779b97f4a7c15ull;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ull;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebull;
    return x ^ (x >> 31);
}

/// Uniform double in [0, 1) from a 64-bit hash.
double u01(uint64_t h) { return double(h >> 11) * (1.0 / 9007199254740992.0); }

/// Jittered fcc lattice, 4-atom basis, lattice constant a = (4/rho)^(1/3).
cpu::Frame make_jittered_fcc(int cells, double a, double jitter) {
    static const double basis[4][3] = {
        {0.0, 0.0, 0.0}, {0.5, 0.5, 0.0}, {0.5, 0.0, 0.5}, {0.0, 0.5, 0.5}};
    const int n = 4 * cells * cells * cells;
    cpu::Frame f(n);
    f.box = Mat3::cubic(cells * a);
    int i = 0;
    for (int ix = 0; ix < cells; ++ix)
        for (int iy = 0; iy < cells; ++iy)
            for (int iz = 0; iz < cells; ++iz)
                for (int b = 0; b < 4; ++b, ++i) {
                    f.atoms.x[i] = (ix + basis[b][0] + 0.25) * a
                                   + jitter * (u01(splitmix64(3 * i + 0)) - 0.5);
                    f.atoms.y[i] = (iy + basis[b][1] + 0.25) * a
                                   + jitter * (u01(splitmix64(3 * i + 1)) - 0.5);
                    f.atoms.z[i] = (iz + basis[b][2] + 0.25) * a
                                   + jitter * (u01(splitmix64(3 * i + 2)) - 0.5);
                    f.atoms.mass[i] = 1.0;
                    f.atoms.type[i] = 0;
                    f.atoms.tag[i] = i;
                }
    return f;
}

/// Write the configuration as a LAMMPS data file (atom_style atomic,
/// reduced units; coordinates at full double precision for an exact
/// energy cross-check).
bool write_lammps_data(const cpu::Frame& f, const char* path) {
    FILE* fp = std::fopen(path, "w");
    if (!fp)
        return false;
    const int n = static_cast<int>(f.atoms.size());
    const double L = f.box(0, 0);
    std::fprintf(fp, "LAMMPS data file - bench_minimize jittered sc lattice\n\n"
                     "%d atoms\n1 atom types\n\n"
                     "0.0 %.17g xlo xhi\n0.0 %.17g ylo yhi\n0.0 %.17g zlo zhi\n\n"
                     "Masses\n\n1 1.0\n\nAtoms # atomic\n\n",
                 n, L, L, L);
    for (int i = 0; i < n; ++i)
        std::fprintf(fp, "%d 1 %.17g %.17g %.17g\n", i + 1,
                     double(f.atoms.x[i]), double(f.atoms.y[i]),
                     double(f.atoms.z[i]));
    std::fclose(fp);
    return true;
}

struct RunResult {
    int iters = 0;
    double wall = 0, e0 = 0, ef = 0, maxf = 0;
    bool converged = false;
};

/// Assemble the Minimize driver around @p opt, run, and harvest metrics.
template <typename Opt>
RunResult run_minimize(Opt& opt, cpu::Frame& frame, double cutoff,
                       int maxiter) {
    cuda::minimizer::Minimize min;
    min.add(opt);
    min.on(ATV::minimizer::Phase::force)
        .add(cuda::NeighborList{cutoff}, -10)
        .add(cuda::PairLJCut{1.0, 1.0, cutoff});
    min.build();

    auto t0 = std::chrono::high_resolution_clock::now();
    int iters = min.run(frame, maxiter);
    auto t1 = std::chrono::high_resolution_clock::now();

    RunResult r;
    r.iters = iters;
    r.wall = std::chrono::duration<double>(t1 - t0).count();
    r.e0 = opt.trace().front().energy;
    r.ef = opt.trace().back().energy;
    r.maxf = opt.max_force();
    r.converged = opt.converged();
    return r;
}

} // namespace

int main(int argc, char** argv) {
    const int n_target = (argc > 1) ? std::atoi(argv[1]) : 256000;
    const std::string style = (argc > 2) ? argv[2] : "lbfgs";
    const double ftol = (argc > 3) ? std::atof(argv[3]) : 1e-4;
    const int maxiter = (argc > 4) ? std::atoi(argv[4]) : 5000;
    const double jitter = (argc > 5) ? std::atof(argv[5]) : 0.1;
    const char* data_out = (argc > 6) ? argv[6] : "";
    const double alpha0 = (argc > 7) ? std::atof(argv[7])
                          : (style == "cg" ? 0.002 : 1.0);

    const double rho = 1.0; // fcc, nn distance = LJ pair minimum
    const double a = std::cbrt(4.0 / rho);
    const double cutoff = 2.0;
    const double BIG = 1e30; // disables the non-force convergence criteria
    int cells = static_cast<int>(std::round(std::cbrt(double(n_target) / 4.0)));
    if (cells < 2)
        cells = 2;
    const int n = 4 * cells * cells * cells;

    cpu::Frame frame = make_jittered_fcc(cells, a, jitter);
    if (data_out[0] != '\0' && !write_lammps_data(frame, data_out)) {
        std::fprintf(stderr, "cannot write %s\n", data_out);
        return 2;
    }

    // Warmup on a tiny system: pays CUDA context + lazy kernel loading.
    {
        cpu::Frame w = make_jittered_fcc(4, a, jitter);
        if (style == "cg") {
            cuda::minimizer::CG opt{alpha0, ftol, BIG, BIG, BIG, BIG};
            run_minimize(opt, w, cutoff, 3);
        } else {
            cuda::minimizer::LBFGS opt{5, alpha0, ftol, 20, BIG, BIG, BIG, BIG};
            opt.line_search();
            run_minimize(opt, w, cutoff, 3);
        }
    }

    RunResult r;
    if (style == "cg") {
        cuda::minimizer::CG opt{alpha0, ftol, BIG, BIG, BIG, BIG};
        r = run_minimize(opt, frame, cutoff, maxiter);
    } else if (style == "lbfgs") {
        cuda::minimizer::LBFGS opt{5, alpha0, ftol, 20, BIG, BIG, BIG, BIG};
        opt.line_search();
        r = run_minimize(opt, frame, cutoff, maxiter);
    } else {
        std::fprintf(stderr, "unknown style '%s' (cg|lbfgs)\n", style.c_str());
        return 2;
    }

    std::printf("CSV,AV,%s,%d,%d,%d,%.4f,%.10g,%.10g,%.3e,%d\n", style.c_str(),
                n, r.iters, r.iters, r.wall, r.e0, r.ef, r.maxf,
                r.converged ? 1 : 0);
    std::printf("# AV %s n=%d: %d iters in %.3f s (%.1f iter/s, %.2f "
                "Matom-iter/s), E %.6f -> %.6f, max|F| %.3e, converged=%d\n",
                style.c_str(), n, r.iters, r.wall, r.iters / r.wall,
                double(n) * r.iters / r.wall / 1e6, r.e0, r.ef, r.maxf,
                r.converged ? 1 : 0);
    return 0;
}
