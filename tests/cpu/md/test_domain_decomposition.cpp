/// @file tests/cpu/md/test_domain_decomposition.cpp
/// @brief RED tests for single-rank `DomainDecomposition` state (Phase 1).
///
/// Spec: .claude/specs/parallel-md-engine.md § A.1 + § B (Context fields).
/// The `Backend` must expose a non-null `Context::dd` in single-rank mode with
/// identity topology, zero ghosts, and an identity tag→local map.

#include <gtest/gtest.h>

#include <cstdint>

#include "atomiverse/core/mat3.h"
#include "atomiverse/cpu/driver/backend.h"
#include "atomiverse/cpu/component/context.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/comm/domain.h"

using namespace ATV;

namespace {

cpu::Frame make_small_frame(int n) {
    cpu::Frame frame(n);
    frame.box = Mat3::cubic(10.0);
    for (int i = 0; i < n; ++i) {
        frame.atoms.x[i] = 0.1 * i;
        frame.atoms.y[i] = 0.2 * i;
        frame.atoms.z[i] = 0.3 * i;
        frame.atoms.mass[i] = 1.0;
    }
    return frame;
}

cpu::Backend make_backend_for(cpu::Frame& frame) {
    cpu::Backend backend;
    backend.init();
    backend.set_frame(frame);
    return backend;
}

} // namespace

// ── single-rank defaults ─────────────────────────────────────────────

TEST(DomainDecomposition, SingleRankDefaults) {
    auto frame = make_small_frame(64);
    auto backend = make_backend_for(frame);

    cpu::Context& ctx = backend.context();
    ASSERT_NE(ctx.dd, nullptr) << "Context::dd must be non-null after "
                                  "Backend::init() even in single-rank mode";

    EXPECT_EQ(ctx.dd->dims[0], 1);
    EXPECT_EQ(ctx.dd->dims[1], 1);
    EXPECT_EQ(ctx.dd->dims[2], 1);
    EXPECT_EQ(ctx.dd->coords[0], 0);
    EXPECT_EQ(ctx.dd->coords[1], 0);
    EXPECT_EQ(ctx.dd->coords[2], 0);
    EXPECT_EQ(ctx.dd->periods[0], 1);
    EXPECT_EQ(ctx.dd->periods[1], 1);
    EXPECT_EQ(ctx.dd->periods[2], 1);
    EXPECT_EQ(ctx.dd->n_local, frame.atoms.size());
    EXPECT_EQ(ctx.dd->n_ghost, 0);
    EXPECT_EQ(ctx.dd->capacity, frame.atoms.size());
    EXPECT_EQ(ctx.dd->neighbor_generation, static_cast<std::uint64_t>(0));
}

// ── send lists empty ─────────────────────────────────────────────────

TEST(DomainDecomposition, SendListsEmpty) {
    auto frame = make_small_frame(32);
    auto backend = make_backend_for(frame);

    cpu::Context& ctx = backend.context();
    ASSERT_NE(ctx.dd, nullptr);
    for (int face = 0; face < 6; ++face) {
        EXPECT_TRUE(ctx.dd->send_lists[face].empty())
            << "send_lists[" << face << "] must be empty in single-rank mode";
    }
}

// ── identity tag→local map ───────────────────────────────────────────
//
// The tag→local map exists solely for bonded-tuple resolution, so it is
// maintained ONLY for frames carrying bonded topology. Two cases:
// bonded frame → identity map; bare frame → map stays empty (and the
// per-atom hash-node memory is never paid).

TEST(DomainDecomposition, TagToLocalIdentityMap) {
    const int n = 48;
    auto frame = make_small_frame(n);
    // Give the frame bonded topology so the map is maintained.
    frame.atoms.map_clear();
    frame.atoms.map_set(0, n);
    for (int i = 0; i + 1 < n; i += 2)
        frame.atoms.add_bond(i, i + 1, 0);
    auto backend = make_backend_for(frame);

    cpu::Context& ctx = backend.context();
    ASSERT_NE(ctx.dd, nullptr);

    EXPECT_EQ(static_cast<int>(ctx.dd->tag_to_local.size()), n);
    for (int i = 0; i < n; ++i) {
        auto it = ctx.dd->tag_to_local.find(frame.atoms.tag(i));
        ASSERT_NE(it, ctx.dd->tag_to_local.end())
            << "tag " << frame.atoms.tag(i) << " missing from tag_to_local";
        EXPECT_EQ(it->second, i);
    }
}

TEST(DomainDecomposition, TagToLocalSkippedWithoutBondedTopology) {
    auto frame = make_small_frame(48);
    auto backend = make_backend_for(frame);

    cpu::Context& ctx = backend.context();
    ASSERT_NE(ctx.dd, nullptr);
    EXPECT_TRUE(ctx.dd->tag_to_local.empty())
        << "tag_to_local must not be populated for a frame without "
           "bonded topology";
    EXPECT_TRUE(frame.atoms.map.empty())
        << "atoms.map must be released for a frame without bonded topology";
}

// ── mpi_rank / mpi_size set in Context ───────────────────────────────

TEST(DomainDecomposition, MpiFieldsInContext) {
    auto frame = make_small_frame(8);
    auto backend = make_backend_for(frame);

    cpu::Context& ctx = backend.context();
    EXPECT_EQ(ctx.mpi_rank, 0);
    EXPECT_EQ(ctx.mpi_size, 1);
}
