/// @file tests/cpu/scf/test_scf_evaluator_no_alloc.cpp
/// @brief Phase 4 — assert that `SCFEvaluator::launch()` allocates the
///        same number of times across steady-state launches.
///
/// The MD-side equivalent (`tests/cpu/test_md_no_alloc.cpp`) asserts
/// **zero** allocations after warmup, because LJ + bonded force kernels
/// already use pre-allocated per-thread buffers and direct double arrays.
/// SCF launches still allocate per-step (xtensor temporaries in the Fock
/// build, integral rebuild, gradient evaluation) — driving them to true
/// zero would require rewriting the SCF cycle in raw BLAS/LAPACK.
///
/// What the spec's acceptance gate ac-014 requires is a **stable** alloc
/// count across steady-state launches: bounded spread (≤5%) and no upward
/// trend over a 12-launch window. This catches:
///
///   * Allocations that grow over time (memory leaks, accumulating
///     diagnostic buffers, history that never shrinks).
///   * SCF state that drifts iteration counts from one launch to the
///     next (which would step the per-launch allocation count up).
///
/// Exact launch-to-launch equality is NOT required: the pre-built
/// persistent SCFSolver carries buffer capacities and DIIS history
/// across launches, and OpenMP workers grow per-thread buffers lazily
/// under scheduler control, so the count wobbles ~1% without any leak.
/// (Exact equality only ever held for the historical transient-solver
/// implementation, whose cold-start allocation sequence was identical
/// every launch.)
///
/// It also does **not** require zero. Reaching zero is a separate
/// `rebuild_integrals_inplace` + `scf_gradient_inplace` follow-up
/// (Phase 4.5).
///
/// Linux/glibc only — uses `__libc_malloc/calloc/realloc` interposition
/// (mirrors `tests/cpu/test_md_no_alloc.cpp`). macOS skips with a clear
/// reason.

#ifdef ATV_HAS_RHF

#include <gtest/gtest.h>

#include <atomic>
#include <cstdlib>

#include "atomiverse/core/mat3.h"
#include "atomiverse/cpu/dynamics_driver.h"
#include "atomiverse/cpu/frame.h"
#include "atomiverse/cpu/md/md.h"
#include "atomiverse/cpu/md/velocity_verlet.h"
#include "atomiverse/cpu/scf/scf_evaluator.h"
#include "atomiverse/cpu/scf/scf_solve.h"

#if defined(__GLIBC__) && defined(__linux__)
#define ATV_TEST_HAVE_MALLOC_HOOK 1
#include <malloc.h>
#endif

#ifdef ATV_TEST_HAVE_MALLOC_HOOK

namespace {

std::atomic<uint64_t>& alloc_count() {
    static std::atomic<uint64_t> c{0};
    return c;
}

bool& alloc_count_armed() {
    static bool armed = false;
    return armed;
}

extern "C" void* __libc_malloc(size_t);
extern "C" void* __libc_calloc(size_t, size_t);
extern "C" void* __libc_realloc(void*, size_t);
extern "C" void __libc_free(void*);

extern "C" void* malloc(size_t sz) {
    if (alloc_count_armed())
        alloc_count().fetch_add(1, std::memory_order_relaxed);
    return __libc_malloc(sz);
}
extern "C" void* calloc(size_t n, size_t sz) {
    if (alloc_count_armed())
        alloc_count().fetch_add(1, std::memory_order_relaxed);
    return __libc_calloc(n, sz);
}
extern "C" void* realloc(void* p, size_t sz) {
    if (alloc_count_armed())
        alloc_count().fetch_add(1, std::memory_order_relaxed);
    return __libc_realloc(p, sz);
}
extern "C" void free(void* p) { __libc_free(p); }

} // namespace

#endif // ATV_TEST_HAVE_MALLOC_HOOK

namespace {

ATV::cpu::Frame make_h2() {
    ATV::cpu::Frame f(2);
    f.atoms.x = {0.0, 0.74};
    f.atoms.y = {0.0, 0.0};
    f.atoms.z = {0.0, 0.0};
    f.atoms.mass = {1.008, 1.008};
    f.atoms.type = {1, 1};
    f.box = ATV::Mat3::cubic(50.0);
    return f;
}

uint64_t armed_launch_count(ATV::cpu::md::MD& sim, ATV::cpu::Frame& frame) {
#ifdef ATV_TEST_HAVE_MALLOC_HOOK
    alloc_count().store(0, std::memory_order_relaxed);
    alloc_count_armed() = true;
    sim.run(frame, /*n_steps=*/1);
    alloc_count_armed() = false;
    return alloc_count().load(std::memory_order_relaxed);
#else
    sim.run(frame, /*n_steps=*/1);
    return 0;
#endif
}

} // namespace

TEST(SCFEvaluatorNoAlloc, ConstantAllocCountAcrossSteadyStateLaunches) {
#ifndef ATV_TEST_HAVE_MALLOC_HOOK
    GTEST_SKIP() << "malloc-hook unavailable on this platform "
                 << "(__GLIBC__ + __linux__ required); the constant-count "
                 << "invariant is still verified by code review and by "
                 << "Linux/glibc CI.";
#else
    using namespace ATV;

    auto frame = make_h2();
    auto integrals =
        cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto cfg = cpu::scf::rhf_config(integrals);

    cpu::md::MD sim;
    sim.add(cpu::md::VelocityVerlet{0.5});
    sim.on(ATV::md::Phase::force).add(cpu::scf::SCFEvaluator{cfg});
    sim.build();

    // Warm-up: lets thread pools, xtensor caches, persistent-solver
    // buffer capacities, and any one-shot allocator behaviour settle
    // into steady state.
    sim.run(frame, /*n_steps=*/4);

    // Take a window of armed measurements. The original contract
    // asserted exact equality launch-to-launch; that held only while
    // SCFEvaluator constructed a transient SCFSolver per launch (an
    // identical cold-start allocation sequence every time). The
    // pre-built persistent solver (no-alloc-in-launch rule) carries
    // vector capacities and DIIS history across launches, and OpenMP
    // workers grow their per-thread buffers lazily on whichever thread
    // the scheduler picks, so the per-launch count wobbles by ~1%
    // without any leak. What ac-014 guards against is *unbounded
    // growth* (leaks, accumulating diagnostic buffers, history that
    // never shrinks) — assert that directly: bounded spread plus no
    // upward trend across the window.
    constexpr int kProbes = 12;
    uint64_t counts[kProbes];
    uint64_t lo = ~uint64_t{0}, hi = 0;
    double mean_first = 0.0, mean_last = 0.0;
    for (int i = 0; i < kProbes; ++i) {
        counts[i] = armed_launch_count(sim, frame);
        lo = counts[i] < lo ? counts[i] : lo;
        hi = counts[i] > hi ? counts[i] : hi;
        (i < kProbes / 2 ? mean_first : mean_last) +=
            static_cast<double>(counts[i]);
    }
    mean_first /= kProbes / 2;
    mean_last /= kProbes / 2;

    // A leak of even one buffer per launch grows the count
    // monotonically across 12 probes; scheduler wobble stays within a
    // few allocations of the floor.
    EXPECT_LE(hi - lo, lo / 20)
        << "per-launch alloc spread exceeds 5%: min=" << lo
        << " max=" << hi;
    EXPECT_LE(mean_last, mean_first * 1.02)
        << "per-launch alloc count trends upward: first-half mean="
        << mean_first << " last-half mean=" << mean_last;

    // Surface the steady-state level via gtest's RecordProperty so
    // xUnit reporters can plot drift over time without cluttering
    // stdout. (Driving this toward zero is the Phase 4.5
    // rebuild_integrals_inplace + scf_gradient_inplace follow-up.)
    RecordProperty("steady_state_alloc_count_per_launch_max",
                   std::to_string(hi).c_str());
#endif
}

#endif // ATV_HAS_RHF
