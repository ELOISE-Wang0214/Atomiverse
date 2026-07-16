/// @file benchmarks/qc/rhf/performance/bench_scf_singlepoint.cpp
/// @brief SCF single-point wall time on real molecules (RHF / 6-31G,
///        production defaults: MINAO guess + DIIS + occupation guard +
///        canonical-orthogonalization linear-dependence screening), plus
///        micro-benchmarks of the per-run() orthogonalizer build
///        (CanonicalOrth eigendecomposition vs the Löwdin S^{-1/2} build).
///
/// Build / run:
///     cmake -B build -DATV_BUILD_BENCHMARKS=ON -DATV_BUILD_RHF=ON ...
///     cmake --build build --target bench_scf_singlepoint
///     ./build/benchmarks/bench_scf_singlepoint
///
/// Geometries come from the shared validation molecule set
/// (benchmarks/qc/rhf/validation/molecules/<stem>.xyz, PySCF-referenced).
#ifdef ATV_HAS_RHF

#include <benchmark/benchmark.h>

#include <cstdio>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>

#include "atomiverse/core/mat3.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/linalg/eigensolve.h"
#include "atomiverse/cpu/scf/basis_set.h"
#include "atomiverse/cpu/scf/integrals.h"
#include "atomiverse/cpu/scf/scf_solver.h"
#include "atomiverse/cpu/scf/scf_types.h"

using namespace ATV;

namespace {

int z_from_symbol(const std::string& s) {
    static const std::pair<const char*, int> TBL[] = {
        {"H", 1}, {"He", 2}, {"Li", 3}, {"Be", 4}, {"B", 5},
        {"C", 6}, {"N", 7},  {"O", 8},  {"F", 9},  {"Ne", 10},
    };
    for (const auto& [sym, z] : TBL)
        if (s == sym)
            return z;
    throw std::runtime_error("z_from_symbol: unknown element " + s);
}

/// Load an XYZ file from the shared validation molecule set.
cpu::Frame load_xyz(const std::string& stem) {
    const std::string path =
        std::string(ATV_QC_DATA_DIR) + "/molecules/" + stem + ".xyz";
    std::ifstream in(path);
    if (!in)
        throw std::runtime_error("cannot open " + path);
    int n = 0;
    std::string line;
    std::getline(in, line);
    n = std::stoi(line);
    std::getline(in, line); // comment
    cpu::Frame f(n);
    for (int i = 0; i < n; ++i) {
        std::getline(in, line);
        std::istringstream ss(line);
        std::string sym;
        double x, y, z;
        ss >> sym >> x >> y >> z;
        f.span(ATV::kAtom_type)[i] = z_from_symbol(sym);
        f.span(ATV::kAtom_x)[i] = x;
        f.span(ATV::kAtom_y)[i] = y;
        f.span(ATV::kAtom_z)[i] = z;
    }
    f.simbox = ATV::cpu::Box::cubic(100.0);
    return f;
}

/// Full production-default RHF single point: integrals prebuilt (measured
/// separately by bench_native); the timed region is SCFSolver::run —
/// orthogonalizer build + MINAO guess + DIIS-accelerated SCF cycle.
void bm_rhf_run(benchmark::State& state, const char* stem) {
    auto frame = load_xyz(stem);
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::Pople631G);
    auto cfg = cpu::scf::rhf_config(ints);

    bool printed = false;
    for (auto _ : state) {
        cpu::scf::SCFSolver solver{cfg};
        solver.build(ints);
        const auto& r = solver.run(ints);
        if (!r.converged)
            state.SkipWithError("SCF did not converge");
        benchmark::DoNotOptimize(r.energy);
        if (!printed) {
            // Parity read-off on stderr (google-benchmark counters mangle
            // negative values); also on the solver's "[scf] done" log line.
            std::fprintf(stderr, "[bench] %s E = %.10f Ha\n", stem, r.energy);
            printed = true;
        }
    }
    state.counters["n_basis"] = ints.n_basis;
}

/// Per-run() orthogonalizer pieces on the real overlap matrix.
void bm_canonical_orth(benchmark::State& state, const char* stem) {
    auto frame = load_xyz(stem);
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::Pople631G);
    for (auto _ : state) {
        cpu::linalg::CanonicalOrth ortho(ints.S, 1e-8);
        benchmark::DoNotOptimize(ortho.X);
    }
    state.counters["n_basis"] = ints.n_basis;
}

void bm_lowdin_build(benchmark::State& state, const char* stem) {
    auto frame = load_xyz(stem);
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::Pople631G);
    for (auto _ : state) {
        auto X = cpu::linalg::build_orthogonalization_matrix(ints.S);
        benchmark::DoNotOptimize(X);
    }
    state.counters["n_basis"] = ints.n_basis;
}

/// What SCFSolver::run actually pays per run() on a healthy basis after the
/// single-eigh optimization: CanonicalOrth screening + the Löwdin matrix
/// from the same eigendecomposition (compare against
/// bm_canonical_orth + bm_lowdin_build, the pre-optimization sum).
void bm_solver_orthogonalizer(benchmark::State& state, const char* stem) {
    auto frame = load_xyz(stem);
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::Pople631G);
    for (auto _ : state) {
        cpu::linalg::CanonicalOrth ortho(ints.S, 1e-8);
        benchmark::DoNotOptimize(ortho.lowdin);
    }
    state.counters["n_basis"] = ints.n_basis;
}

} // namespace

BENCHMARK_CAPTURE(bm_rhf_run, h2o, "h2o")->Unit(benchmark::kMillisecond);
BENCHMARK_CAPTURE(bm_rhf_run, benzene, "benzene")
    ->Unit(benchmark::kMillisecond)
    ->Iterations(3);
BENCHMARK_CAPTURE(bm_rhf_run, naphthalene, "naphthalene")
    ->Unit(benchmark::kMillisecond)
    ->Iterations(3);
BENCHMARK_CAPTURE(bm_rhf_run, aspirin, "aspirin")
    ->Unit(benchmark::kMillisecond)
    ->Iterations(3);
BENCHMARK_CAPTURE(bm_rhf_run, caffeine, "caffeine")
    ->Unit(benchmark::kMillisecond)
    ->Iterations(3);

BENCHMARK_CAPTURE(bm_canonical_orth, aspirin, "aspirin")
    ->Unit(benchmark::kMicrosecond);
BENCHMARK_CAPTURE(bm_lowdin_build, aspirin, "aspirin")
    ->Unit(benchmark::kMicrosecond);
BENCHMARK_CAPTURE(bm_canonical_orth, caffeine, "caffeine")
    ->Unit(benchmark::kMicrosecond);
BENCHMARK_CAPTURE(bm_lowdin_build, caffeine, "caffeine")
    ->Unit(benchmark::kMicrosecond);
BENCHMARK_CAPTURE(bm_solver_orthogonalizer, aspirin, "aspirin")
    ->Unit(benchmark::kMicrosecond);
BENCHMARK_CAPTURE(bm_solver_orthogonalizer, caffeine, "caffeine")
    ->Unit(benchmark::kMicrosecond);

BENCHMARK_MAIN();

#endif // ATV_HAS_RHF
