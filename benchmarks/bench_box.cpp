/// @file benchmarks/bench_box.cpp
/// @brief Microbench for `cpu::Box::min_image` — scalar + xsimd batch path.
///
/// Reports ns/op for:
///   • scalar `Box::min_image` (ortho fast path)
///   • SIMD `Box::min_image_batch` at the platform's native xsimd batch width
/// Used to gate ac-012: scalar within ±5% of the pre-refactor branchy form,
/// SIMD ≥ 1.5× scalar.

#include <benchmark/benchmark.h>

#include <random>
#include <vector>

#include "atomiverse/cpu/box.h"
#include "atomiverse/cpu/box_simd.h"

#include <xsimd/xsimd.hpp>

namespace {

constexpr int kNumDisp = 1 << 16; // 64k displacements per iteration

std::vector<double> random_disps(std::size_t n, double L, std::uint32_t seed) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> dist(-1.4 * L, 1.4 * L);
    std::vector<double> v;
    v.reserve(n);
    for (std::size_t i = 0; i < n; ++i)
        v.push_back(dist(rng));
    return v;
}

} // namespace

static void BM_BoxMinImageScalar(benchmark::State& state) {
    const auto bx = ATV::cpu::Box::orthorhombic(20.0, 20.0, 20.0);
    auto dx = random_disps(kNumDisp, 20.0, 0xA1);
    auto dy = random_disps(kNumDisp, 20.0, 0xA2);
    auto dz = random_disps(kNumDisp, 20.0, 0xA3);
    for (auto _ : state) {
        for (int i = 0; i < kNumDisp; ++i) {
            double x = dx[i], y = dy[i], z = dz[i];
            bx.min_image(x, y, z);
            benchmark::DoNotOptimize(x);
            benchmark::DoNotOptimize(y);
            benchmark::DoNotOptimize(z);
        }
    }
    state.SetItemsProcessed(state.iterations() * kNumDisp);
}
BENCHMARK(BM_BoxMinImageScalar);

static void BM_BoxMinImageBatch(benchmark::State& state) {
    using batch_t = xsimd::batch<double>;
    constexpr std::size_t W = batch_t::size;
    static_assert(kNumDisp % W == 0,
                  "kNumDisp must be a multiple of the SIMD width");

    const auto bx = ATV::cpu::Box::orthorhombic(20.0, 20.0, 20.0);
    auto dx = random_disps(kNumDisp, 20.0, 0xB1);
    auto dy = random_disps(kNumDisp, 20.0, 0xB2);
    auto dz = random_disps(kNumDisp, 20.0, 0xB3);

    for (auto _ : state) {
        for (int i = 0; i < kNumDisp; i += static_cast<int>(W)) {
            batch_t bx_d = batch_t::load_unaligned(dx.data() + i);
            batch_t by_d = batch_t::load_unaligned(dy.data() + i);
            batch_t bz_d = batch_t::load_unaligned(dz.data() + i);
            bx.min_image_batch(bx_d, by_d, bz_d);
            benchmark::DoNotOptimize(bx_d);
            benchmark::DoNotOptimize(by_d);
            benchmark::DoNotOptimize(bz_d);
        }
    }
    state.SetItemsProcessed(state.iterations() * kNumDisp);
}
BENCHMARK(BM_BoxMinImageBatch);

// Reference: pre-refactor inline branchy form, copied verbatim from
// the original `static inline minimum_image` in lj.cpp.
static void BM_BranchyReference(benchmark::State& state) {
    const double L = 20.0;
    const double half_L = 0.5 * L;
    auto dx = random_disps(kNumDisp, L, 0xC1);
    auto dy = random_disps(kNumDisp, L, 0xC2);
    auto dz = random_disps(kNumDisp, L, 0xC3);
    for (auto _ : state) {
        for (int i = 0; i < kNumDisp; ++i) {
            double x = dx[i], y = dy[i], z = dz[i];
            if (x > half_L)
                x -= L;
            if (x < -half_L)
                x += L;
            if (y > half_L)
                y -= L;
            if (y < -half_L)
                y += L;
            if (z > half_L)
                z -= L;
            if (z < -half_L)
                z += L;
            benchmark::DoNotOptimize(x);
            benchmark::DoNotOptimize(y);
            benchmark::DoNotOptimize(z);
        }
    }
    state.SetItemsProcessed(state.iterations() * kNumDisp);
}
BENCHMARK(BM_BranchyReference);

BENCHMARK_MAIN();
