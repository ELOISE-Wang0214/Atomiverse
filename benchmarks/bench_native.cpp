/// @file bench_native.cpp
/// @brief Speed benchmarks for all native CPU implementations replacing libcint/libxc.
///
/// Coverage map (native function → external library equivalent):
///
/// ┌──────────────────────────────────────────────────────────────────┐
/// │  Native function               │ Replaces (libcint / libxc)     │
/// ├──────────────────────────────────────────────────────────────────┤
/// │  boys_f0, boys_fn              │ libcint: CINTf_0e, cint_boys  │
/// │  evaluate_xc(LSDA)            │ libxc: xc_lda_exc_vxc (S+VWN5)│
/// │  evaluate_xc(PBE)             │ libxc: xc_gga_exc_vxc (PBE)   │
/// │  build_molecular_grid          │ external grid libraries        │
/// │  build_integral_derivatives    │ libcint: int1e_ipovlp, etc.   │
/// │  scf_solve (UHF)              │ (no direct equivalent)         │
/// │  scf_solve (UKS/LSDA)         │ libxc-backed KS solver         │
/// │  scf_solve (UKS/PBE)          │ libxc-backed KS solver         │
/// │  scf_gradient (RHF)           │ libcint derivative integrals   │
/// │  scf_gradient (UHF)           │ libcint derivative integrals   │
/// └──────────────────────────────────────────────────────────────────┘
///
/// Integral primitives (overlap, kinetic, nuclear, ERI) and RHF solve
/// are benchmarked in bench_rhf.cu — not duplicated here.

#ifdef ATV_HAS_RHF

#include <benchmark/benchmark.h>
#include <cmath>
#include <vector>

#include "atomiverse/cpu/scf/native_integrals.h"
#include "atomiverse/cpu/scf/rys_roots.h"
#include "atomiverse/cpu/scf/rys_eri.h"
#include "atomiverse/cpu/scf/integrals.h"
#include "atomiverse/cpu/scf/shell_pair.h"
#include "atomiverse/cpu/scf/operators.h"
#include "atomiverse/cpu/scf/integral_derivatives.h"
#include "atomiverse/cpu/scf/xc_functional.h"
#ifdef ATV_HAS_LIBXC
#include "atomiverse/cpu/scf/xc_libxc.h"
#endif
#include "atomiverse/cpu/scf/grid.h"
#include "atomiverse/cpu/scf/scf_solve.h"
#include "atomiverse/cpu/scf/gradient.h"
#include "atomiverse/core/constants.h"
#include "atomiverse/core/log.h"

using namespace ATV;
using namespace ATV::cpu::scf;

// ── Helper frames ────────────────────────────────────────────────────────

static cpu::Frame make_h2() {
    cpu::Frame frame(2);
    frame.atoms.x = {0.0, 0.74};
    frame.atoms.y = {0.0, 0.0};
    frame.atoms.z = {0.0, 0.0};
    frame.atoms.mass = {1.008, 1.008};
    frame.atoms.type = {1, 1};
    frame.box = Mat3::cubic(100.0);
    return frame;
}

static cpu::Frame make_lih() {
    cpu::Frame frame(2);
    frame.atoms.x = {0.0, 1.595};
    frame.atoms.y = {0.0, 0.0};
    frame.atoms.z = {0.0, 0.0};
    frame.atoms.mass = {6.941, 1.008};
    frame.atoms.type = {3, 1};
    frame.box = Mat3::cubic(100.0);
    return frame;
}

static cpu::Frame make_aspirin() {
    cpu::Frame frame(21);
    int types[21] = {
        6, 6, 6, 6, 6, 6,
        6,
        8, 8,
        8,
        6,
        8,
        6,
        1, 1, 1, 1,
        1,
        1, 1, 1,
    };
    double masses[21] = {
        12.011, 12.011, 12.011, 12.011, 12.011, 12.011,
        12.011,
        15.999, 15.999,
        15.999,
        12.011,
        15.999,
        12.011,
        1.008, 1.008, 1.008, 1.008,
        1.008,
        1.008, 1.008, 1.008,
    };
    double coords[21][3] = {
        { 1.2879,  0.1068,  0.0000},
        { 0.1768,  1.0827,  0.0000},
        { 0.5171,  2.4290,  0.0000},
        { 1.8519,  2.7993,  0.0000},
        { 2.9293,  1.8524,  0.0000},
        { 2.6191,  0.4890,  0.0000},
        {-1.2076,  0.6698,  0.0000},
        {-1.5541, -0.4994,  0.0000},
        {-2.0632,  1.6805,  0.0000},
        { 0.9973, -1.2285,  0.0000},
        {-0.2830, -1.7447,  0.0000},
        {-0.3696, -2.9537,  0.0000},
        {-1.3593, -0.7411,  0.3000},
        {-0.2980,  3.1365,  0.0000},
        { 2.0791,  3.8586,  0.0000},
        { 3.9560,  2.1939,  0.0000},
        { 3.4546, -0.2019,  0.0000},
        {-2.9828,  1.3621,  0.0000},
        {-1.0697, -0.1755,  1.1920},
        {-2.3482, -1.1583,  0.5001},
        {-1.4267, -0.2091, -0.6502},
    };

    for (int i = 0; i < 21; ++i) {
        frame.atoms.type[i] = types[i];
        frame.atoms.mass[i] = masses[i];
        frame.atoms.x[i] = coords[i][0];
        frame.atoms.y[i] = coords[i][1];
        frame.atoms.z[i] = coords[i][2];
    }
    frame.box = Mat3::cubic(100.0);
    return frame;
}

// ── Additional molecule frames for multi-size benchmarks ──────────────────

static cpu::Frame make_water() {
    // H2O: 3 atoms, 7 AOs (STO-3G)
    cpu::Frame frame(3);
    frame.atoms.type = {8, 1, 1};
    frame.atoms.mass = {15.999, 1.008, 1.008};
    frame.atoms.x = {0.0, 0.757, -0.757};
    frame.atoms.y = {0.0, 0.587, 0.587};
    frame.atoms.z = {0.0, 0.0, 0.0};
    frame.box = Mat3::cubic(100.0);
    return frame;
}

static cpu::Frame make_methane() {
    // CH4: 5 atoms, 9 AOs (STO-3G)
    cpu::Frame frame(5);
    frame.atoms.type = {6, 1, 1, 1, 1};
    frame.atoms.mass = {12.011, 1.008, 1.008, 1.008, 1.008};
    frame.atoms.x = {0.0, 0.629, -0.629, 0.629, -0.629};
    frame.atoms.y = {0.0, 0.629, -0.629, -0.629, 0.629};
    frame.atoms.z = {0.0, 0.629, 0.629, -0.629, -0.629};
    frame.box = Mat3::cubic(100.0);
    return frame;
}

static cpu::Frame make_ethanol() {
    // C2H5OH: 9 atoms, 26 AOs (STO-3G)
    cpu::Frame frame(9);
    int types[] = {6, 6, 8, 1, 1, 1, 1, 1, 1};
    double masses[] = {12.011, 12.011, 15.999, 1.008, 1.008, 1.008, 1.008, 1.008, 1.008};
    double coords[][3] = {
        {-0.748, -0.015,  0.024},  // C1
        { 0.558,  0.683, -0.060},  // C2
        { 1.651, -0.217,  0.018},  // O
        {-0.813, -0.543,  0.978},  // H
        {-0.813, -0.726, -0.801},  // H
        {-1.576,  0.695,  0.026},  // H
        { 0.622,  1.395,  0.768},  // H
        { 0.622,  1.214, -1.009},  // H
        { 1.617, -0.694,  0.865},  // H (OH)
    };
    for (int i = 0; i < 9; ++i) {
        frame.atoms.type[i] = types[i];
        frame.atoms.mass[i] = masses[i];
        frame.atoms.x[i] = coords[i][0];
        frame.atoms.y[i] = coords[i][1];
        frame.atoms.z[i] = coords[i][2];
    }
    frame.box = Mat3::cubic(100.0);
    return frame;
}

// ═════════════════════════════════════════════════════════════════════════
// Boys function  (replaces libcint internal Boys evaluator)
// ═════════════════════════════════════════════════════════════════════════

static void BM_BoysF0(benchmark::State& state) {
    double args[] = {0.0, 0.01, 0.1, 0.5, 1.0, 5.0, 10.0, 50.0, 100.0};
    for (auto _ : state) {
        for (double x : args)
            benchmark::DoNotOptimize(native::boys_f0(x));
    }
}
BENCHMARK(BM_BoysF0);

static void BM_BoysFn(benchmark::State& state) {
    int n = static_cast<int>(state.range(0));
    double args[] = {0.01, 0.1, 0.5, 1.0, 5.0, 10.0, 50.0};
    for (auto _ : state) {
        for (double x : args)
            benchmark::DoNotOptimize(native::boys_fn(n, x));
    }
}
BENCHMARK(BM_BoysFn)->Arg(0)->Arg(1)->Arg(2)->Arg(4);

// ═════════════════════════════════════════════════════════════════════════
// XC functional  (replaces libxc: xc_lda_exc_vxc / xc_gga_exc_vxc)
// ═════════════════════════════════════════════════════════════════════════

namespace {

struct LsdaPointCase {
    const char* label;
    double rho_a;
    double rho_b;
};

struct PbePointCase {
    const char* label;
    double rho_a;
    double rho_b;
    double sigma_aa;
    double sigma_ab;
    double sigma_bb;
};

LsdaPointCase make_lsda_point_case(int case_id) {
    static constexpr LsdaPointCase kCases[] = {
        {"low_density", 1e-6, 1e-8},
        {"mid_density", 0.5, 0.5},
        {"high_density", 50.0, 50.0},
    };
    return kCases[case_id];
}

PbePointCase make_pbe_point_case(int case_id) {
    static constexpr PbePointCase kCases[] = {
        {"low_density_small_grad", 1e-3, 1e-3, 1e-6, 1e-7, 1e-6},
        {"mid_density_small_grad", 0.5, 0.5, 0.01, 0.005, 0.01},
        {"high_density_large_grad", 5.0, 3.0, 2.0, 1.0, 1.5},
    };
    return kCases[case_id];
}

} // namespace

static void BM_XC_LSDA_SinglePoint(benchmark::State& state) {
    const auto point = make_lsda_point_case(static_cast<int>(state.range(0)));
    state.SetLabel(point.label);
    for (auto _ : state) {
        auto out = evaluate_xc(XCType::LSDA, point.rho_a, point.rho_b);
        benchmark::DoNotOptimize(out.exc);
    }
}
BENCHMARK(BM_XC_LSDA_SinglePoint)->DenseRange(0, 2);

#ifdef ATV_HAS_LIBXC
static void BM_XC_LSDA_SinglePoint_Libxc(benchmark::State& state) {
    const auto point = make_lsda_point_case(static_cast<int>(state.range(0)));
    state.SetLabel(point.label);
    for (auto _ : state) {
        auto out = evaluate_xc_libxc(XCType::LSDA, point.rho_a, point.rho_b);
        benchmark::DoNotOptimize(out.exc);
    }
}
BENCHMARK(BM_XC_LSDA_SinglePoint_Libxc)->DenseRange(0, 2);
#endif

static void BM_XC_PBE_SinglePoint(benchmark::State& state) {
    const auto point = make_pbe_point_case(static_cast<int>(state.range(0)));
    state.SetLabel(point.label);
    for (auto _ : state) {
        auto out = evaluate_xc(
            XCType::PBE,
            point.rho_a, point.rho_b,
            point.sigma_aa, point.sigma_ab, point.sigma_bb);
        benchmark::DoNotOptimize(out.exc);
    }
}
BENCHMARK(BM_XC_PBE_SinglePoint)->DenseRange(0, 2);

#ifdef ATV_HAS_LIBXC
static void BM_XC_PBE_SinglePoint_Libxc(benchmark::State& state) {
    const auto point = make_pbe_point_case(static_cast<int>(state.range(0)));
    state.SetLabel(point.label);
    for (auto _ : state) {
        auto out = evaluate_xc_libxc(
            XCType::PBE,
            point.rho_a, point.rho_b,
            point.sigma_aa, point.sigma_ab, point.sigma_bb);
        benchmark::DoNotOptimize(out.exc);
    }
}
BENCHMARK(BM_XC_PBE_SinglePoint_Libxc)->DenseRange(0, 2);
#endif

static void BM_XC_LSDA_Batch(benchmark::State& state) {
    int n = static_cast<int>(state.range(0));
    std::vector<double> rho_a(n), rho_b(n);
    for (int i = 0; i < n; ++i) {
        double r = 0.01 + 2.0 * i / n;
        rho_a[i] = 0.5 * std::exp(-r);
        rho_b[i] = rho_a[i];
    }
    for (auto _ : state) {
        double total = 0.0;
        for (int i = 0; i < n; ++i) {
            auto out = evaluate_xc(XCType::LSDA, rho_a[i], rho_b[i]);
            total += out.exc;
        }
        benchmark::DoNotOptimize(total);
    }
    state.SetItemsProcessed(state.iterations() * n);
}
BENCHMARK(BM_XC_LSDA_Batch)->Arg(100)->Arg(1000)->Arg(10000);

#ifdef ATV_HAS_LIBXC
static void BM_XC_LSDA_Batch_Libxc(benchmark::State& state) {
    int n = static_cast<int>(state.range(0));
    std::vector<double> rho_a(n), rho_b(n);
    for (int i = 0; i < n; ++i) {
        double r = 0.01 + 2.0 * i / n;
        rho_a[i] = 0.5 * std::exp(-r);
        rho_b[i] = rho_a[i];
    }
    for (auto _ : state) {
        double total = 0.0;
        for (int i = 0; i < n; ++i) {
            auto out = evaluate_xc_libxc(XCType::LSDA, rho_a[i], rho_b[i]);
            total += out.exc;
        }
        benchmark::DoNotOptimize(total);
    }
    state.SetItemsProcessed(state.iterations() * n);
}
BENCHMARK(BM_XC_LSDA_Batch_Libxc)->Arg(100)->Arg(1000)->Arg(10000);
#endif

static void BM_XC_PBE_Batch(benchmark::State& state) {
    int n = static_cast<int>(state.range(0));
    std::vector<double> rho_a(n), rho_b(n), sigma_aa(n), sigma_ab(n), sigma_bb(n);
    for (int i = 0; i < n; ++i) {
        double r = 0.01 + 2.0 * i / n;
        rho_a[i] = 0.5 * std::exp(-r);
        rho_b[i] = rho_a[i];
        double grad = 0.3 * std::exp(-r);
        sigma_aa[i] = grad * grad;
        sigma_ab[i] = 0.5 * grad * grad;
        sigma_bb[i] = sigma_aa[i];
    }
    for (auto _ : state) {
        double total = 0.0;
        for (int i = 0; i < n; ++i) {
            auto out = evaluate_xc(XCType::PBE, rho_a[i], rho_b[i],
                                   sigma_aa[i], sigma_ab[i], sigma_bb[i]);
            total += out.exc;
        }
        benchmark::DoNotOptimize(total);
    }
    state.SetItemsProcessed(state.iterations() * n);
}
BENCHMARK(BM_XC_PBE_Batch)->Arg(100)->Arg(1000)->Arg(10000);

#ifdef ATV_HAS_LIBXC
static void BM_XC_PBE_Batch_Libxc(benchmark::State& state) {
    int n = static_cast<int>(state.range(0));
    std::vector<double> rho_a(n), rho_b(n), sigma_aa(n), sigma_ab(n), sigma_bb(n);
    for (int i = 0; i < n; ++i) {
        double r = 0.01 + 2.0 * i / n;
        rho_a[i] = 0.5 * std::exp(-r);
        rho_b[i] = rho_a[i];
        double grad = 0.3 * std::exp(-r);
        sigma_aa[i] = grad * grad;
        sigma_ab[i] = 0.5 * grad * grad;
        sigma_bb[i] = sigma_aa[i];
    }
    for (auto _ : state) {
        double total = 0.0;
        for (int i = 0; i < n; ++i) {
            auto out = evaluate_xc_libxc(
                XCType::PBE,
                rho_a[i], rho_b[i],
                sigma_aa[i], sigma_ab[i], sigma_bb[i]);
            total += out.exc;
        }
        benchmark::DoNotOptimize(total);
    }
    state.SetItemsProcessed(state.iterations() * n);
}
BENCHMARK(BM_XC_PBE_Batch_Libxc)->Arg(100)->Arg(1000)->Arg(10000);
#endif

// ═════════════════════════════════════════════════════════════════════════
// Grid construction  (replaces external grid libraries)
// ═════════════════════════════════════════════════════════════════════════

static void BM_Grid_SingleAtom(benchmark::State& state) {
    int n_angular = static_cast<int>(state.range(0));
    native::Atom h = {1, {0.0, 0.0, 0.0}};
    for (auto _ : state) {
        auto grid = build_molecular_grid(&h, 1, 50, n_angular);
        benchmark::DoNotOptimize(grid.points.size());
    }
}
BENCHMARK(BM_Grid_SingleAtom)
    ->Arg(6)->Arg(26)->Arg(110)->Arg(302)->Arg(590)
    ->Unit(benchmark::kMicrosecond);

static void BM_Grid_H2(benchmark::State& state) {
    int n_angular = static_cast<int>(state.range(0));
    double d = 0.74 * kAngstromToBohr;
    native::Atom atoms[] = {{1, {0, 0, 0}}, {1, {d, 0, 0}}};
    for (auto _ : state) {
        auto grid = build_molecular_grid(atoms, 2, 50, n_angular);
        benchmark::DoNotOptimize(grid.points.size());
    }
}
BENCHMARK(BM_Grid_H2)
    ->Arg(6)->Arg(26)->Arg(110)->Arg(302)
    ->Unit(benchmark::kMicrosecond);

// ═════════════════════════════════════════════════════════════════════════
// Integral derivatives  (replaces libcint: int1e_ipovlp, int1e_ipkin, etc.)
// ═════════════════════════════════════════════════════════════════════════

static void BM_IntegralDeriv_H2(benchmark::State& state) {
    auto frame = make_h2();
    for (auto _ : state) {
        auto dpack = build_integral_derivatives_sto3g(frame);
        benchmark::DoNotOptimize(dpack.dV_nuc(0, 0));
    }
}
BENCHMARK(BM_IntegralDeriv_H2)->Unit(benchmark::kMicrosecond);

static void BM_IntegralDeriv_LiH(benchmark::State& state) {
    auto frame = make_lih();
    for (auto _ : state) {
        auto dpack = build_integral_derivatives_sto3g(frame);
        benchmark::DoNotOptimize(dpack.dV_nuc(0, 0));
    }
}
BENCHMARK(BM_IntegralDeriv_LiH)->Unit(benchmark::kMicrosecond);

static void BM_IntegralDeriv_Aspirin(benchmark::State& state) {
    auto frame = make_aspirin();
    for (auto _ : state) {
        auto dpack = build_integral_derivatives_sto3g(frame);
        benchmark::DoNotOptimize(dpack.dV_nuc(0, 0));
    }
}
BENCHMARK(BM_IntegralDeriv_Aspirin)->Unit(benchmark::kMillisecond);

// ═════════════════════════════════════════════════════════════════════════
// UHF solve  (open-shell HF without libcint dependency)
// ═════════════════════════════════════════════════════════════════════════

static void BM_UHFSolve_H2(benchmark::State& state) {
    auto frame = make_h2();
    auto integrals = build_integrals(frame, cpu::scf::BasisSet::STO3G);
    constexpr double tol = 1e-10;
    auto cfg = uhf_config(integrals, 1, 1, tol);
    cfg.max_iter = 100;
    for (auto _ : state) {
        auto result = scf_solve(integrals, cfg);
        benchmark::DoNotOptimize(result.energy);
    }
}
BENCHMARK(BM_UHFSolve_H2)->Unit(benchmark::kMicrosecond);

static void BM_UHFSolve_LiH(benchmark::State& state) {
    auto frame = make_lih();
    auto integrals = build_integrals(frame, cpu::scf::BasisSet::STO3G);
    constexpr double tol = 1e-10;
    auto cfg = uhf_config(integrals, 2, 2, tol);
    cfg.max_iter = 100;
    for (auto _ : state) {
        auto result = scf_solve(integrals, cfg);
        benchmark::DoNotOptimize(result.energy);
    }
}
BENCHMARK(BM_UHFSolve_LiH)->Unit(benchmark::kMicrosecond);

// ═════════════════════════════════════════════════════════════════════════
// UKS solve  (DFT with native XC, replaces libxc-backed KS solver)
// ═════════════════════════════════════════════════════════════════════════

static void BM_UKSSolve_H2_LSDA(benchmark::State& state) {
    auto frame = make_h2();
    auto integrals = build_integrals(frame, cpu::scf::BasisSet::STO3G);
    constexpr double tol = 1e-8;
    auto cfg = uks_config(integrals, 1, 1, XCType::LSDA, tol);
    cfg.max_iter = 100;
    for (auto _ : state) {
        auto result = scf_solve(integrals, cfg);
        benchmark::DoNotOptimize(result.energy);
    }
}
BENCHMARK(BM_UKSSolve_H2_LSDA)->Unit(benchmark::kMillisecond);

static void BM_UKSSolve_H2_PBE(benchmark::State& state) {
    auto frame = make_h2();
    auto integrals = build_integrals(frame, cpu::scf::BasisSet::STO3G);
    constexpr double tol = 1e-8;
    auto cfg = uks_config(integrals, 1, 1, XCType::PBE, tol);
    cfg.max_iter = 100;
    for (auto _ : state) {
        auto result = scf_solve(integrals, cfg);
        benchmark::DoNotOptimize(result.energy);
    }
}
BENCHMARK(BM_UKSSolve_H2_PBE)->Unit(benchmark::kMillisecond);

// ═════════════════════════════════════════════════════════════════════════
// Gradient  (replaces libcint derivative integral API)
// ═════════════════════════════════════════════════════════════════════════

static void BM_RHFGradient_H2(benchmark::State& state) {
    auto frame = make_h2();
    auto integrals = build_integrals(frame, cpu::scf::BasisSet::STO3G);
    constexpr double tol = 1e-10;
    auto cfg = rhf_config(integrals, tol);
    cfg.max_iter = 200;
    auto result = scf_solve(integrals, cfg);
    for (auto _ : state) {
        auto grad = scf_gradient(frame, integrals, result, cfg);
        benchmark::DoNotOptimize(grad.energy);
    }
}
BENCHMARK(BM_RHFGradient_H2)->Unit(benchmark::kMicrosecond);

static void BM_UHFGradient_H2(benchmark::State& state) {
    auto frame = make_h2();
    auto integrals = build_integrals(frame, cpu::scf::BasisSet::STO3G);
    constexpr double tol = 1e-10;
    auto cfg = uhf_config(integrals, 1, 1, tol);
    cfg.max_iter = 200;
    auto result = scf_solve(integrals, cfg);
    for (auto _ : state) {
        auto grad = scf_gradient(frame, integrals, result, cfg);
        benchmark::DoNotOptimize(grad.energy);
    }
}
BENCHMARK(BM_UHFGradient_H2)->Unit(benchmark::kMicrosecond);

static void BM_RHFIntegrals_Aspirin(benchmark::State& state) {
    auto frame = make_aspirin();
    for (auto _ : state) {
        auto integrals = build_integrals(frame, cpu::scf::BasisSet::STO3G);
        benchmark::DoNotOptimize(integrals.Hcore(0, 0));
    }
}
BENCHMARK(BM_RHFIntegrals_Aspirin)->Unit(benchmark::kMillisecond);

static void BM_RHFIntegrals_Aspirin_631Gs(benchmark::State& state) {
    auto frame = make_aspirin();
    for (auto _ : state) {
        auto integrals = build_integrals(frame, cpu::scf::BasisSet::Pople631Gs);
        benchmark::DoNotOptimize(integrals.Hcore(0, 0));
    }
}
BENCHMARK(BM_RHFIntegrals_Aspirin_631Gs)->Unit(benchmark::kMillisecond);

static void BM_RHFIntegrals_Aspirin_def2svp(benchmark::State& state) {
    auto frame = make_aspirin();
    for (auto _ : state) {
        auto integrals = build_integrals(frame, cpu::scf::BasisSet::Def2SVP);
        benchmark::DoNotOptimize(integrals.Hcore(0, 0));
    }
}
BENCHMARK(BM_RHFIntegrals_Aspirin_def2svp)->Unit(benchmark::kMillisecond);

// ── G-matrix builder: AO-level vs shell-block ──────────────────────────────

/// @brief Direct RHF G = J - 0.5*K using AO-level shell pairs (baseline).
static void BM_GBuild_AO_Aspirin(benchmark::State& state) {
    auto frame = make_aspirin();
    auto integrals = build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto sp = build_shell_pair_list(integrals);
    int n = integrals.n_basis;
    RealMatrix P = xt::zeros<double>({(size_t)n, (size_t)n});
    for (int i = 0; i < n; ++i) P(i, i) = 1.0;
    RealMatrix G;
    for (auto _ : state) {
        build_g_rhf_direct(sp, P, G);
        benchmark::DoNotOptimize(G(0, 0));
    }
}
BENCHMARK(BM_GBuild_AO_Aspirin)->Unit(benchmark::kMillisecond);

/// @brief Direct RHF G = J - 0.5*K using shell-block pairs (libcint/PySCF approach).
static void BM_GBuild_Block_Aspirin(benchmark::State& state) {
    auto frame = make_aspirin();
    auto integrals = build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto sp_blk = build_shell_block_pair_list(integrals);
    int n = integrals.n_basis;
    RealMatrix P = xt::zeros<double>({(size_t)n, (size_t)n});
    for (int i = 0; i < n; ++i) P(i, i) = 1.0;
    RealMatrix G;
    for (auto _ : state) {
        build_g_rhf_direct_block(sp_blk, P, G);
        benchmark::DoNotOptimize(G(0, 0));
    }
}
BENCHMARK(BM_GBuild_Block_Aspirin)->Unit(benchmark::kMillisecond);

/// @brief Direct RHF G = J - 0.5*K for aspirin in a d-polarized basis.
///
/// Phase 4 (eri-hpc-opt) benchmark: exercises the d-shell ERI path
/// (`native::eri_block` → `rys_eri_block`, the Rys-quadrature engine
/// algorithmically identical to libcint) inside the shell-block direct
/// Fock builder. 6-31G* puts a d-shell on every heavy atom, so the bulk
/// of the quartet work is (d..|d..) — the cc-pVDZ-class regime the spec
/// targets.
///
/// Baseline reference: the eri-hpc-opt spec records PySCF direct-Fock
/// baselines only for STO-3G aspirin (~952 ms). No PySCF cc-pVDZ /
/// 6-31G* direct-Fock baseline is documented in the spec, and PySCF is
/// not installed in this build environment — so this benchmark records
/// the native median wall time without a fabricated PySCF number.
/// For scale: STO-3G aspirin (73 AOs) is Phase-3 ~500 ms; 6-31G*
/// aspirin has ~280 AOs, an ~4x AO-count increase, and d-shell quartets
/// cost far more per block than s/p, so a multi-second build is
/// expected.
static void BM_GBuild_Block_Aspirin_ccPVDZ(benchmark::State& state) {
    auto frame = make_aspirin();
    // 6-31G* — d-polarized; the d-shell ERI path is the Phase-4 target.
    auto integrals = build_integrals(frame, cpu::scf::BasisSet::Pople631Gs);
    auto sp_blk = build_shell_block_pair_list(integrals);
    int n = integrals.n_basis;
    RealMatrix P = xt::zeros<double>({(size_t)n, (size_t)n});
    for (int i = 0; i < n; ++i) P(i, i) = 1.0;
    RealMatrix G;
    for (auto _ : state) {
        build_g_rhf_direct_block(sp_blk, P, G);
        benchmark::DoNotOptimize(G(0, 0));
    }
    state.counters["n_AOs"] = n;
    state.SetLabel("aspirin 6-31G* (d-polarized)");
}
BENCHMARK(BM_GBuild_Block_Aspirin_ccPVDZ)->Unit(benchmark::kMillisecond);

// ── Multi-size G-build (block): water → methane → ethanol → aspirin ───

static void BM_GBuild_Block_Water(benchmark::State& state) {
    auto frame = make_water();
    auto integrals = build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto sp_blk = build_shell_block_pair_list(integrals);
    int n = integrals.n_basis;
    RealMatrix P = xt::zeros<double>({(size_t)n, (size_t)n});
    for (int i = 0; i < n; ++i) P(i, i) = 1.0;
    RealMatrix G;
    for (auto _ : state) {
        build_g_rhf_direct_block(sp_blk, P, G);
        benchmark::DoNotOptimize(G(0, 0));
    }
    state.SetLabel("7 AOs");
}
BENCHMARK(BM_GBuild_Block_Water)->Unit(benchmark::kMicrosecond);

static void BM_GBuild_Block_Methane(benchmark::State& state) {
    auto frame = make_methane();
    auto integrals = build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto sp_blk = build_shell_block_pair_list(integrals);
    int n = integrals.n_basis;
    RealMatrix P = xt::zeros<double>({(size_t)n, (size_t)n});
    for (int i = 0; i < n; ++i) P(i, i) = 1.0;
    RealMatrix G;
    for (auto _ : state) {
        build_g_rhf_direct_block(sp_blk, P, G);
        benchmark::DoNotOptimize(G(0, 0));
    }
    state.SetLabel("9 AOs");
}
BENCHMARK(BM_GBuild_Block_Methane)->Unit(benchmark::kMicrosecond);

static void BM_GBuild_Block_Ethanol(benchmark::State& state) {
    auto frame = make_ethanol();
    auto integrals = build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto sp_blk = build_shell_block_pair_list(integrals);
    int n = integrals.n_basis;
    RealMatrix P = xt::zeros<double>({(size_t)n, (size_t)n});
    for (int i = 0; i < n; ++i) P(i, i) = 1.0;
    RealMatrix G;
    for (auto _ : state) {
        build_g_rhf_direct_block(sp_blk, P, G);
        benchmark::DoNotOptimize(G(0, 0));
    }
    state.SetLabel("26 AOs");
}
BENCHMARK(BM_GBuild_Block_Ethanol)->Unit(benchmark::kMillisecond);

// ── Multi-size RHF solve: water → ethanol → aspirin ──────────────────

static void BM_RHFSolve_Water(benchmark::State& state) {
    auto frame = make_water();
    auto integrals = build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto cfg = rhf_config(integrals);
    cfg.tol = 1e-6;
    cfg.max_iter = 200;
    for (auto _ : state) {
        auto result = scf_solve(integrals, cfg);
        benchmark::DoNotOptimize(result.energy);
    }
    state.SetLabel("7 AOs");
}
BENCHMARK(BM_RHFSolve_Water)->Unit(benchmark::kMillisecond);

static void BM_RHFSolve_Ethanol(benchmark::State& state) {
    auto frame = make_ethanol();
    auto integrals = build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto cfg = rhf_config(integrals);
    cfg.tol = 1e-6;
    cfg.max_iter = 200;
    for (auto _ : state) {
        auto result = scf_solve(integrals, cfg);
        benchmark::DoNotOptimize(result.energy);
    }
    state.SetLabel("26 AOs");
}
BENCHMARK(BM_RHFSolve_Ethanol)->Unit(benchmark::kMillisecond);

static void BM_RHFSolve_Aspirin(benchmark::State& state) {
    auto frame = make_aspirin();
    // PySCF approach: store ERI once, BLAS contraction per iteration.
    auto integrals = build_integrals(frame, cpu::scf::BasisSet::STO3G, 0, /*compute_eri=*/true);
    auto cfg = rhf_config(integrals);
    cfg.tol = 1e-9;   // match PySCF conv_tol (|dE|<tol + ||err||<sqrt(tol))
    cfg.max_iter = 50;
    for (auto _ : state) {
        auto result = scf_solve(integrals, cfg);
        benchmark::DoNotOptimize(result.energy);
        state.counters["iters"] = result.iterations;
        state.counters["energy_ha"] = result.energy;
        state.counters["converged"] = result.converged ? 1.0 : 0.0;
    }
}
BENCHMARK(BM_RHFSolve_Aspirin)->Unit(benchmark::kMillisecond);

static void BM_RHFGradient_Aspirin(benchmark::State& state) {
    auto frame = make_aspirin();
    auto integrals = build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto cfg = rhf_config(integrals);
    cfg.tol = 1e-6;
    cfg.max_iter = 200;
    auto result = scf_solve(integrals, cfg);
    for (auto _ : state) {
        auto grad = scf_gradient(frame, integrals, result, cfg);
        benchmark::DoNotOptimize(grad.energy);
    }
}
BENCHMARK(BM_RHFGradient_Aspirin)->Unit(benchmark::kMillisecond);

// ═════════════════════════════════════════════════════════════════════════
// Combined pipeline benchmarks
// ═════════════════════════════════════════════════════════════════════════

static void BM_RHF_IntegralsToGradient_H2(benchmark::State& state) {
    // Full pipeline: integrals → solve → gradient
    auto frame = make_h2();
    constexpr double tol = 1e-10;
    for (auto _ : state) {
        auto integrals = build_integrals(frame, cpu::scf::BasisSet::STO3G);
        auto cfg = rhf_config(integrals, tol);
        cfg.max_iter = 200;
        auto result = scf_solve(integrals, cfg);
        auto grad = scf_gradient(frame, integrals, result, cfg);
        benchmark::DoNotOptimize(grad.energy);
    }
}
BENCHMARK(BM_RHF_IntegralsToGradient_H2)->Unit(benchmark::kMicrosecond);

static void BM_RHF_IntegralsToGradient_LiH(benchmark::State& state) {
    auto frame = make_lih();
    constexpr double tol = 1e-10;
    for (auto _ : state) {
        auto integrals = build_integrals(frame, cpu::scf::BasisSet::STO3G);
        auto cfg = rhf_config(integrals, tol);
        cfg.max_iter = 200;
        auto result = scf_solve(integrals, cfg);
        auto grad = scf_gradient(frame, integrals, result, cfg);
        benchmark::DoNotOptimize(grad.energy);
    }
}
BENCHMARK(BM_RHF_IntegralsToGradient_LiH)->Unit(benchmark::kMicrosecond);

static void BM_RHF_IntegralsToGradient_Aspirin(benchmark::State& state) {
    auto frame = make_aspirin();
    for (auto _ : state) {
        auto integrals = build_integrals(frame, cpu::scf::BasisSet::STO3G);
        auto cfg = rhf_config(integrals);
        cfg.tol = 1e-6;
        cfg.max_iter = 200;
        auto result = scf_solve(integrals, cfg);
        auto grad = scf_gradient(frame, integrals, result, cfg);
        benchmark::DoNotOptimize(grad.energy);
    }
}
BENCHMARK(BM_RHF_IntegralsToGradient_Aspirin)->Unit(benchmark::kMillisecond);

// ═════════════════════════════════════════════════════════════════════════
// Rys quadrature roots  (ported from libcint, drives the Rys ERI engine)
// ═════════════════════════════════════════════════════════════════════════

static void BM_RysRoots(benchmark::State& state) {
    const int nroots = static_cast<int>(state.range(0));
    double x_vals[] = {0.01, 0.5, 2.0, 7.5, 15.0, 35.0, 80.0};
    double u[5], w[5];
    for (auto _ : state) {
        for (double x : x_vals)
            native::rys_roots(nroots, x, u, w);
        benchmark::DoNotOptimize(u[0]);
        benchmark::DoNotOptimize(w[0]);
    }
    state.SetItemsProcessed(state.iterations() * 7);
}
BENCHMARK(BM_RysRoots)->DenseRange(1, 5);

// ═════════════════════════════════════════════════════════════════════════
// Rys ERI block  (shell-level: direct comparison with BM_ERIBlock_Libcint
//                 in bench_libcint.cpp)
// ═════════════════════════════════════════════════════════════════════════

namespace {

/// @brief Build a minimal ShellBlock from STO-3G hydrogen (l=0) at origin.
static native::ShellBlock make_s_block() {
    static const double alpha[] = {3.42525091, 0.62391373, 0.16885540};
    static const double coeff[] = {0.15432897, 0.53532814, 0.44463454};
    native::ShellBlock b{};
    b.n_prim = 3; b.alpha = alpha; b.coeff = coeff;
    b.center[0] = b.center[1] = b.center[2] = 0.0;
    b.l = 0; b.n_cart = 1; b.ao_start = 0;
    return b;
}

/// @brief Build a minimal p-ShellBlock at (1.4, 0, 0) Bohr.
static native::ShellBlock make_p_block() {
    // STO-3G carbon 2sp p-coefficients (Z=6, second raw shell)
    static const double alpha[] = {2.9412494, 0.6834831, 0.2222899};
    static const double coeff[] = {-0.09996723, 0.39951283, 0.70011547};
    native::ShellBlock b{};
    b.n_prim = 3; b.alpha = alpha; b.coeff = coeff;
    b.center[0] = 1.4; b.center[1] = 0.0; b.center[2] = 0.0;
    b.l = 1; b.n_cart = 3; b.ao_start = 1;
    return b;
}

} // namespace

static void BM_RysERI_ssss(benchmark::State& state) {
    auto A = make_s_block(); auto B = A; auto C = A; auto D = A;
    double out[1];
    for (auto _ : state) {
        out[0] = 0.0;
        native::rys_eri_block(A, B, C, D, out);
        benchmark::DoNotOptimize(out[0]);
    }
}
BENCHMARK(BM_RysERI_ssss);

static void BM_RysERI_psss(benchmark::State& state) {
    auto A = make_p_block(); auto B = make_s_block();
    auto C = make_s_block(); auto D = make_s_block();
    double out[3];
    for (auto _ : state) {
        std::fill(out, out + 3, 0.0);
        native::rys_eri_block(A, B, C, D, out);
        benchmark::DoNotOptimize(out[0]);
    }
}
BENCHMARK(BM_RysERI_psss);

static void BM_RysERI_ppss(benchmark::State& state) {
    auto A = make_p_block(); auto B = make_p_block();
    auto C = make_s_block(); auto D = make_s_block();
    double out[9];
    for (auto _ : state) {
        std::fill(out, out + 9, 0.0);
        native::rys_eri_block(A, B, C, D, out);
        benchmark::DoNotOptimize(out[0]);
    }
}
BENCHMARK(BM_RysERI_ppss);

static void BM_RysERI_pppp(benchmark::State& state) {
    auto A = make_p_block(); auto B = make_p_block();
    auto C = make_p_block(); auto D = make_p_block();
    double out[81];
    for (auto _ : state) {
        std::fill(out, out + 81, 0.0);
        native::rys_eri_block(A, B, C, D, out);
        benchmark::DoNotOptimize(out[0]);
    }
}
BENCHMARK(BM_RysERI_pppp);

int main(int argc, char** argv) {
    ATV::log_init(spdlog::level::warn);
    benchmark::Initialize(&argc, argv);
    if (benchmark::ReportUnrecognizedArguments(argc, argv))
        return 1;
    benchmark::RunSpecifiedBenchmarks();
    benchmark::Shutdown();
    return 0;
}

#else
int main() { return 0; }
#endif
