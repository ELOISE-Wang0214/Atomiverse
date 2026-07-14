/// @file tests/cpu/frame/test_frame_registry.cpp
/// @brief Unit tests for the frame-block registry (frame-block-01): schema
///        completeness, the typed mapping API (span/has/insert), registry-
///        owned extension columns riding grow/copy_row/pack, the typed meta
///        region, and the launch-boundary resolution discipline.

#include <gtest/gtest.h>

#include <set>
#include <string>
#include <vector>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/cpu/component/component.h"
#include "atomiverse/cpu/comm/atom_pack.h"
#include "atomiverse/cpu/driver/evaluate.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/facade/component_role.h"

using namespace ATV;

// ── Schema completeness (ac-001) ────────────────────────────────────────────

TEST(FrameRegistry, BindingVisitsEverySchemaRowExactlyOnce) {
    cpu::Frame frame(4);
    int visited = 0;
    std::set<std::string> keys;
    frame.atoms.for_each_column([&](const FieldSpec& spec, auto& /*arr*/) {
        ++visited;
        EXPECT_TRUE(keys.insert(spec.key).second)
            << "duplicate schema key " << spec.key;
    });
    // One binding line per table row — a drifting macro list breaks this.
    EXPECT_EQ(visited, kAtomFieldCount);
    EXPECT_EQ(static_cast<int>(keys.size()), kAtomFieldCount);
}

TEST(FrameRegistry, ColumnWidthsMatchTheBoundStorage) {
    cpu::Frame frame(4);
    frame.atoms.reserve_topology(2, 1, 1, 6);
    frame.atoms.for_each_column([&](const FieldSpec& spec, auto& arr) {
        const auto cols =
            static_cast<std::size_t>(frame.atoms.cols_of(spec.stride));
        if (arr.dimension() == 1) {
            EXPECT_EQ(cols, 1u) << spec.key;
        } else {
            EXPECT_EQ(arr.shape()[1], cols) << spec.key;
        }
    });
}

// ── Typed spans (ac-002) ────────────────────────────────────────────────────

TEST(FrameRegistry, SpanResolvesTypedColumns) {
    cpu::Frame frame(3);
    auto x = frame.span(kAtom_x);
    x[1] = StorageReal(4.5);
    EXPECT_EQ(frame.atoms.x[1], StorageReal(4.5)); // same storage

    auto img = frame.span(kAtom_image);
    EXPECT_EQ(img.cols, 3);
    img(2, 1) = 7;
    EXPECT_EQ(frame.atoms.image(2, 1), 7);
}

TEST(FrameRegistry, SpanRejectsUnknownKeyAndWrongType) {
    cpu::Frame frame(2);
    EXPECT_THROW(frame.span(ColKey<StorageReal>{"atoms", "nope"}),
                 cpu::column_error);
    EXPECT_THROW(frame.span(ColKey<int>{"atoms", "x"}), cpu::column_error);
    EXPECT_THROW(frame.span(ColKey<StorageReal>{"bonds", "x"}),
                 cpu::column_error);
}

TEST(FrameRegistry, HasTracksTheLazyTopologyGroup) {
    cpu::Frame frame(3);
    EXPECT_TRUE(frame.has("atoms", "x"));
    EXPECT_TRUE(frame.has("atoms", "num_bond")); // counts always present
    EXPECT_FALSE(frame.has("atoms", "bond_atom")); // lazy wide storage
    frame.atoms.map_clear();
    frame.atoms.map_set(0, 3);
    frame.atoms.add_bond(0, 1, 1);
    EXPECT_TRUE(frame.has("atoms", "bond_atom"));
    EXPECT_FALSE(frame.has("atoms", "rho"));
    EXPECT_FALSE(frame.has("bonds", "x"));
}

// ── Extension columns (ac-002) ──────────────────────────────────────────────

TEST(FrameRegistry, ExtensionColumnLivesTheFullLifecycle) {
    cpu::Frame frame(3);
    frame.insert("rho", ColDType::Storage, 1);
    EXPECT_TRUE(frame.has("atoms", "rho"));

    auto rho = frame.span(ColKey<StorageReal>{"atoms", "rho"});
    rho[0] = 1.5;
    rho[2] = 2.5;

    // Survives a chunked reallocation (the probe-4 bug class)…
    frame.atoms.grow(frame.atoms.nmax + 1);
    rho = frame.span(ColKey<StorageReal>{"atoms", "rho"});
    EXPECT_EQ(rho[0], StorageReal(1.5));
    EXPECT_EQ(rho[2], StorageReal(2.5));
    EXPECT_EQ(rho[frame.atoms.nmax - 1], StorageReal(0)); // zeroed tail

    // …and rides copy_row.
    frame.atoms.copy_row(1, 2);
    EXPECT_EQ(rho[1], StorageReal(2.5));
}

TEST(FrameRegistry, InsertRejectsCollisions) {
    cpu::Frame frame(2);
    EXPECT_THROW(frame.insert("x", ColDType::Storage, 1),
                 std::invalid_argument);
    frame.insert("rho", ColDType::Storage, 1);
    EXPECT_THROW(frame.insert("rho", ColDType::Storage, 1),
                 std::invalid_argument);
    EXPECT_THROW(frame.insert("bad", ColDType::Storage, 0),
                 std::invalid_argument);
}

TEST(FrameRegistry, ExchangedExtensionColumnRidesThePayload) {
    cpu::Frame a(2);
    a.insert("rho", ColDType::Storage, 2, /*migrates=*/true,
             /*exchanged=*/true);
    auto rho = a.span(ColKey<StorageReal>{"atoms", "rho"});
    rho(1, 0) = 3.5;
    rho(1, 1) = 4.5;
    a.atoms.x[1] = 9.0;

    cpu::Frame b(2); // same extras layout on the receiving side (SPMD)
    b.insert("rho", ColDType::Storage, 2, true, true);

    std::vector<double> buf(
        static_cast<std::size_t>(cpu::comm::max_exchange_doubles(a.atoms)));
    const int packed = cpu::comm::pack_exchange(a.atoms, 1, buf.data());
    const int unpacked = cpu::comm::unpack_exchange(b.atoms, 0, buf.data());
    EXPECT_EQ(packed, unpacked);
    EXPECT_LE(packed, static_cast<int>(buf.size()));

    EXPECT_EQ(b.atoms.x[0], StorageReal(9.0));
    auto rho_b = b.span(ColKey<StorageReal>{"atoms", "rho"});
    EXPECT_EQ(rho_b(0, 0), StorageReal(3.5));
    EXPECT_EQ(rho_b(0, 1), StorageReal(4.5));
}

// ── Typed meta (ac-002) ─────────────────────────────────────────────────────

TEST(FrameRegistry, MetaViewsAliasTheLiveBookkeeping) {
    cpu::Frame frame(5);
    EXPECT_EQ(std::get<std::int64_t>(frame.meta("nlocal")), 5);
    frame.set_meta("nlocal", cpu::MetaValue{std::int64_t{3}});
    EXPECT_EQ(frame.atoms.nlocal, 3); // the view routes to the member

    frame.scalars.energy = -1.25;
    EXPECT_EQ(std::get<double>(frame.meta("energy")), -1.25);
    frame.set_meta("energy", cpu::MetaValue{-2.5});
    EXPECT_EQ(frame.scalars.energy, -2.5);
}

TEST(FrameRegistry, MetaOwnedEntriesAndErrors) {
    cpu::Frame frame(1);
    EXPECT_FALSE(frame.has_meta("note"));
    frame.set_meta("note", cpu::MetaValue{std::string("equilibrated")});
    EXPECT_TRUE(frame.has_meta("note"));
    EXPECT_EQ(std::get<std::string>(frame.meta("note")), "equilibrated");
    EXPECT_TRUE(frame.has_meta("nlocal"));
    EXPECT_THROW(frame.meta("absent"), std::out_of_range);
}

// ── Resolution discipline (ac-003) ──────────────────────────────────────────

TEST(FrameRegistry, ResolutionWindowCountsAndThrows) {
    cpu::Frame frame(2);
    frame.begin_resolution_window(4);
    (void)frame.span(kAtom_x);
    (void)frame.span(kAtom_vx);
    (void)frame.has("atoms", "x");
    EXPECT_EQ(frame.resolution_count(), 3);
    EXPECT_NO_THROW(frame.end_resolution_window("ok_component"));

    frame.begin_resolution_window(2);
    for (int i = 0; i < 5; ++i)
        (void)frame.span(kAtom_x);
    EXPECT_THROW(frame.end_resolution_window("greedy_component"),
                 std::logic_error);
}

namespace {

/// @brief A component that breaks the 13.7× rule on purpose: resolves a key
///        per "atom" inside its loop.
class KeyLookupInALoop : public cpu::ComponentBase {
  public:
    std::string name() const override { return "key_lookup_in_a_loop"; }
    ComponentRole role() const override { return ComponentRole::Evaluator; }
    Capability capabilities() const override { return Capability::Energy; }
    void launch(const cpu::LaunchContext& ctx) override {
        double acc = 0.0;
        for (int i = 0; i < 100; ++i) {
            auto x = ctx.frame.span(kAtom_x); // WRONG: inside the loop
            acc += x[0];
        }
        ctx.frame.scalars.energy = acc;
    }
};

} // namespace

TEST(FrameRegistry, DriverTripsTheViolatingComponent) {
    cpu::Frame frame(4);
    cpu::Evaluate eval;
    eval.on(evaluate::Phase::evaluate).add(KeyLookupInALoop{});
    eval.build();
    EXPECT_THROW(eval.run(frame), std::logic_error);
}
