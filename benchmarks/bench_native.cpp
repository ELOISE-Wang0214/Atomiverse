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
#include "atomiverse/cpu/scf/integrals.h"
#include "atomiverse/cpu/scf/integral_derivatives.h"
#include "atomiverse/cpu/scf/xc_functional.h"
#include "atomiverse/cpu/scf/grid.h"
#include "atomiverse/cpu/scf/scf_solve.h"
#include "atomiverse/cpu/scf/gradient.h"
#include "atomiverse/core/constants.h"

using namespace ATV;
using namespace ATV::cpu::scf;

// ── Helper frames ────────────────────────────────────────────────────────

static cpu::Frame make_h2() {
    cpu::Frame frame(2);
    frame.atoms.x = {Real(0), Real(0.74)};
    frame.atoms.y = {Real(0), Real(0)};
    frame.atoms.z = {Real(0), Real(0)};
    frame.atoms.mass = {Real(1.008), Real(1.008)};
    frame.atoms.type = {1, 1};
    frame.box = Mat3::cubic(100.0);
    return frame;
}

static cpu::Frame make_lih() {
    cpu::Frame frame(2);
    frame.atoms.x = {Real(0), Real(1.595)};
    frame.atoms.y = {Real(0), Real(0)};
    frame.atoms.z = {Real(0), Real(0)};
    frame.atoms.mass = {Real(6.941), Real(1.008)};
    frame.atoms.type = {3, 1};
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

static void BM_XC_LSDA_SinglePoint(benchmark::State& state) {
    Real rho_a = 0.5, rho_b = 0.5;
    for (auto _ : state) {
        auto out = evaluate_xc(XCType::LSDA, rho_a, rho_b);
        benchmark::DoNotOptimize(out.exc);
    }
}
BENCHMARK(BM_XC_LSDA_SinglePoint);

static void BM_XC_PBE_SinglePoint(benchmark::State& state) {
    Real rho_a = 0.5, rho_b = 0.5;
    Real sigma_aa = 0.1, sigma_ab = 0.05, sigma_bb = 0.1;
    for (auto _ : state) {
        auto out = evaluate_xc(XCType::PBE, rho_a, rho_b, sigma_aa, sigma_ab, sigma_bb);
        benchmark::DoNotOptimize(out.exc);
    }
}
BENCHMARK(BM_XC_PBE_SinglePoint);

static void BM_XC_LSDA_Batch(benchmark::State& state) {
    int n = static_cast<int>(state.range(0));
    std::vector<Real> rho_a(n), rho_b(n);
    for (int i = 0; i < n; ++i) {
        double r = 0.01 + 2.0 * i / n;
        rho_a[i] = Real(0.5 * std::exp(-r));
        rho_b[i] = rho_a[i];
    }
    for (auto _ : state) {
        Real total = 0;
        for (int i = 0; i < n; ++i) {
            auto out = evaluate_xc(XCType::LSDA, rho_a[i], rho_b[i]);
            total += out.exc;
        }
        benchmark::DoNotOptimize(total);
    }
    state.SetItemsProcessed(state.iterations() * n);
}
BENCHMARK(BM_XC_LSDA_Batch)->Arg(100)->Arg(1000)->Arg(10000);

static void BM_XC_PBE_Batch(benchmark::State& state) {
    int n = static_cast<int>(state.range(0));
    std::vector<Real> rho_a(n), rho_b(n), sigma_aa(n), sigma_ab(n), sigma_bb(n);
    for (int i = 0; i < n; ++i) {
        double r = 0.01 + 2.0 * i / n;
        rho_a[i] = Real(0.5 * std::exp(-r));
        rho_b[i] = rho_a[i];
        double grad = 0.3 * std::exp(-r);
        sigma_aa[i] = Real(grad * grad);
        sigma_ab[i] = Real(0.5 * grad * grad);
        sigma_bb[i] = sigma_aa[i];
    }
    for (auto _ : state) {
        Real total = 0;
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

// ═════════════════════════════════════════════════════════════════════════
// UHF solve  (open-shell HF without libcint dependency)
// ═════════════════════════════════════════════════════════════════════════

static void BM_UHFSolve_H2(benchmark::State& state) {
    auto frame = make_h2();
    auto integrals = build_integrals_sto3g(frame);
    Real tol = (sizeof(Real) == sizeof(float)) ? Real(1e-5) : Real(1e-10);
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
    auto integrals = build_integrals_sto3g(frame);
    Real tol = (sizeof(Real) == sizeof(float)) ? Real(1e-5) : Real(1e-10);
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
    auto integrals = build_integrals_sto3g(frame);
    Real tol = (sizeof(Real) == sizeof(float)) ? Real(1e-5) : Real(1e-8);
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
    auto integrals = build_integrals_sto3g(frame);
    Real tol = (sizeof(Real) == sizeof(float)) ? Real(1e-5) : Real(1e-8);
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
    auto integrals = build_integrals_sto3g(frame);
    Real tol = (sizeof(Real) == sizeof(float)) ? Real(1e-5) : Real(1e-10);
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
    auto integrals = build_integrals_sto3g(frame);
    Real tol = (sizeof(Real) == sizeof(float)) ? Real(1e-5) : Real(1e-10);
    auto cfg = uhf_config(integrals, 1, 1, tol);
    cfg.max_iter = 200;
    auto result = scf_solve(integrals, cfg);
    for (auto _ : state) {
        auto grad = scf_gradient(frame, integrals, result, cfg);
        benchmark::DoNotOptimize(grad.energy);
    }
}
BENCHMARK(BM_UHFGradient_H2)->Unit(benchmark::kMicrosecond);

// ═════════════════════════════════════════════════════════════════════════
// Combined pipeline benchmarks
// ═════════════════════════════════════════════════════════════════════════

static void BM_RHF_IntegralsToGradient_H2(benchmark::State& state) {
    // Full pipeline: integrals → solve → gradient
    auto frame = make_h2();
    Real tol = (sizeof(Real) == sizeof(float)) ? Real(1e-5) : Real(1e-10);
    for (auto _ : state) {
        auto integrals = build_integrals_sto3g(frame);
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
    Real tol = (sizeof(Real) == sizeof(float)) ? Real(1e-5) : Real(1e-10);
    for (auto _ : state) {
        auto integrals = build_integrals_sto3g(frame);
        auto cfg = rhf_config(integrals, tol);
        cfg.max_iter = 200;
        auto result = scf_solve(integrals, cfg);
        auto grad = scf_gradient(frame, integrals, result, cfg);
        benchmark::DoNotOptimize(grad.energy);
    }
}
BENCHMARK(BM_RHF_IntegralsToGradient_LiH)->Unit(benchmark::kMicrosecond);

BENCHMARK_MAIN();

#else
int main() { return 0; }
#endif
