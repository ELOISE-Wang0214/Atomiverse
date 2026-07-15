/// @file tests/cpu/test_task_codec.cpp
/// @brief Unit tests for the `cpu/task/` TaskSpec data model and its MPI-free
///        byte codec (`cpu::task::pack` / `cpu::task::unpack`).
///
/// Spec: .claude/specs/mpi-launcher-embedded-python-02-codec.md. Covers the
/// two validating TaskSpec constructors (MD / AM1 SinglePoint), field-exact
/// pack→unpack roundtrips over the full nine-component MD bindable set, and
/// the five wire-tampering rejections (version, task_kind, component_id,
/// truncation, kind-invariant violation). Pure serial, no MPI, no I/O.

#include <gtest/gtest.h>

#include <algorithm>
#include <cctype>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <variant>
#include <vector>

#include "atomiverse/cpu/task/task_codec.h"
#include "atomiverse/cpu/task/task_spec.h"
#include "atomiverse/facade/stage.h"

using namespace ATV;
using namespace ATV::cpu::task;

namespace {

// ── Wire-format offsets ─────────────────────────────────────────────────
//
// Per the spec's wire-format section the stream is a tightly packed,
// native-endianness byte sequence:
//   {magic 'ATVT' (4 B), u32 version, u8 task_kind, u32 n_entries,
//    [entry ...], i64 n_steps}
// so the fixed header offsets are:
constexpr std::size_t kMagicOffset = 0;   ///< 4 bytes: 'A' 'T' 'V' 'T'
constexpr std::size_t kVersionOffset = 4; ///< u32 version, bytes [4, 8)
constexpr std::size_t kTaskKindOffset = 8;  ///< u8 task_kind, byte 8
constexpr std::size_t kNEntriesOffset = 9;  ///< u32 n_entries, bytes [9, 13)
// entry = {u16 component_id, u8 has_binding, (u32 phase_slot, i32 order)?,
//          per-field payload}; the first entry starts right after n_entries:
constexpr std::size_t kFirstComponentIdOffset = 13; ///< u16, bytes [13, 15)

template <typename T>
T read_scalar(const std::vector<std::byte>& buf, std::size_t offset) {
    T value{};
    std::memcpy(&value, buf.data() + offset, sizeof(T));
    return value;
}

template <typename T>
void write_scalar(std::vector<std::byte>& buf, std::size_t offset, T value) {
    std::memcpy(buf.data() + offset, &value, sizeof(T));
}

/// @brief All-nine-MD-component spec: three default-binding entries
///        (the `md.add(comp)` path components of bind_md.cpp — VelocityVerlet,
///        InitVelocity, NoseHoover) and six explicit `{phase, order}` entries
///        covering three distinct md::Phase values and nonzero orders.
TaskSpec make_full_md_spec() {
    std::vector<ComponentEntry> entries;
    entries.push_back(ComponentEntry{
        ComponentParams{NeighborListParams{6.0, 0.3}},
        ExplicitBinding{ATV::md::Phase::force, -10}});
    entries.push_back(ComponentEntry{
        ComponentParams{PairLJCutParams{0.0104, 3.4, 6.0}},
        ExplicitBinding{ATV::md::Phase::force, 0}});
    entries.push_back(ComponentEntry{
        ComponentParams{BondHarmonicParams{{15.2, 2.5}, {1.53, 1.09}}},
        ExplicitBinding{ATV::md::Phase::force, 1}});
    entries.push_back(ComponentEntry{
        ComponentParams{AngleHarmonicParams{{2.7}, {1.911}}},
        ExplicitBinding{ATV::md::Phase::force, 2}});
    entries.push_back(ComponentEntry{
        ComponentParams{
            DihedralMultiHarmonicParams{{{1.73, -4.49, 0.776, 6.99, 0.0}}}},
        ExplicitBinding{ATV::md::Phase::pre_force, 3}});
    entries.push_back(
        ComponentEntry{ComponentParams{VelocityVerletParams{0.5}},
                       std::nullopt});
    entries.push_back(
        ComponentEntry{ComponentParams{InitVelocityParams{300.0, 4242u}},
                       std::nullopt});
    entries.push_back(
        ComponentEntry{ComponentParams{NoseHooverParams{300.0, 100.0, 0.5}},
                       std::nullopt});
    entries.push_back(ComponentEntry{ComponentParams{ThermoParams{100}},
                                     ExplicitBinding{
                                         ATV::md::Phase::end_of_step, 50}});
    return TaskSpec(std::move(entries), /*n_steps=*/7);
}

/// @brief Field-by-field ComponentEntry equality (TaskSpec itself need not
///        provide operator==; the ten param structs must).
void expect_entry_eq(const ComponentEntry& got, const ComponentEntry& want) {
    EXPECT_EQ(got.params.index(), want.params.index());
    EXPECT_TRUE(got.params == want.params);
    ASSERT_EQ(got.binding.has_value(), want.binding.has_value());
    if (want.binding.has_value()) {
        EXPECT_EQ(got.binding->phase, want.binding->phase);
        EXPECT_EQ(got.binding->order, want.binding->order);
    }
}

std::string lowercase(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    return s;
}

/// @brief unpack() must throw std::runtime_error whose what() contains
///        `needle` (the offending value, per the codec error contract).
void expect_unpack_runtime_error_containing(const std::vector<std::byte>& buf,
                                            const std::string& needle) {
    try {
        (void)unpack(buf.data(), buf.size());
        FAIL() << "unpack accepted a tampered buffer (expected "
                  "std::runtime_error containing '"
               << needle << "')";
    } catch (const std::runtime_error& e) {
        EXPECT_NE(std::string(e.what()).find(needle), std::string::npos)
            << "what() does not name the offending value '" << needle
            << "': " << e.what();
    }
}

// ── Variant layout: component_id stability ──────────────────────────────

// The wire component_id is the variant index; the nine MD alternatives'
// indices are frozen and AM1EvaluatorParams is appended LAST (index 9) so
// kWireVersion stays 1.
static_assert(std::variant_size_v<ComponentParams> == 10,
              "ComponentParams must have exactly ten alternatives");
static_assert(
    std::is_same_v<std::variant_alternative_t<9, ComponentParams>,
                   AM1EvaluatorParams>,
    "AM1EvaluatorParams must be the LAST (10th) variant alternative");
static_assert(
    std::is_same_v<std::variant_alternative_t<0, ComponentParams>,
                   PairLJCutParams>,
    "PairLJCutParams must be the first variant alternative (id 0)");

TEST(TaskCodec, Am1ParamsOccupyLastComponentId) {
    const ComponentParams p{AM1EvaluatorParams{}};
    EXPECT_EQ(p.index(), 9u);
}

TEST(TaskCodec, WireVersionIsOne) {
    EXPECT_EQ(TaskSpec::kWireVersion, 1u);
}

// ── (a) MD roundtrip over the full bindable set ─────────────────────────

TEST(TaskCodec, MdRoundtripAllNineComponentTypesFieldExact) {
    const TaskSpec spec = make_full_md_spec();
    ASSERT_EQ(spec.kind(), TaskKind::MD);

    const std::vector<std::byte> buf = pack(spec);
    const TaskSpec back = unpack(buf.data(), buf.size());

    EXPECT_EQ(back.kind(), TaskKind::MD);
    EXPECT_EQ(back.n_steps(), spec.n_steps());
    ASSERT_EQ(back.entries().size(), spec.entries().size());
    for (std::size_t i = 0; i < spec.entries().size(); ++i) {
        SCOPED_TRACE("entry " + std::to_string(i));
        expect_entry_eq(back.entries()[i], spec.entries()[i]);
    }
}

TEST(TaskCodec, MdRoundtripPreservesVectorParamsExactly) {
    const TaskSpec spec = make_full_md_spec();
    const std::vector<std::byte> buf = pack(spec);
    const TaskSpec back = unpack(buf.data(), buf.size());

    // Spot-check the three vector-carrying params through std::get — the
    // u64-length-prefixed wire vectors must come back element-exact.
    const auto& bond = std::get<BondHarmonicParams>(back.entries()[2].params);
    EXPECT_EQ(bond.K, (std::vector<double>{15.2, 2.5}));
    EXPECT_EQ(bond.r0, (std::vector<double>{1.53, 1.09}));
    const auto& angle = std::get<AngleHarmonicParams>(back.entries()[3].params);
    EXPECT_EQ(angle.K, (std::vector<double>{2.7}));
    EXPECT_EQ(angle.theta0, (std::vector<double>{1.911}));
    const auto& dih =
        std::get<DihedralMultiHarmonicParams>(back.entries()[4].params);
    EXPECT_EQ(dih.Cn, (std::vector<std::vector<double>>{
                          {1.73, -4.49, 0.776, 6.99, 0.0}}));
}

TEST(TaskCodec, MdRoundtripPreservesExplicitAndDefaultBindings) {
    const TaskSpec spec = make_full_md_spec();
    const std::vector<std::byte> buf = pack(spec);
    const TaskSpec back = unpack(buf.data(), buf.size());

    // Explicit entries keep their typed phase + order.
    ASSERT_TRUE(back.entries()[0].binding.has_value());
    EXPECT_EQ(back.entries()[0].binding->phase, ATV::md::Phase::force);
    EXPECT_EQ(back.entries()[0].binding->order, -10);
    ASSERT_TRUE(back.entries()[4].binding.has_value());
    EXPECT_EQ(back.entries()[4].binding->phase, ATV::md::Phase::pre_force);
    EXPECT_EQ(back.entries()[4].binding->order, 3);
    ASSERT_TRUE(back.entries()[8].binding.has_value());
    EXPECT_EQ(back.entries()[8].binding->phase, ATV::md::Phase::end_of_step);
    EXPECT_EQ(back.entries()[8].binding->order, 50);
    // Default-binding entries stay nullopt (the md.add(comp) path).
    EXPECT_FALSE(back.entries()[5].binding.has_value());
    EXPECT_FALSE(back.entries()[6].binding.has_value());
    EXPECT_FALSE(back.entries()[7].binding.has_value());
}

// ── Wire header layout (pins the tamper offsets below) ──────────────────

TEST(TaskCodec, WireHeaderLayoutMatchesSpec) {
    const TaskSpec spec = make_full_md_spec();
    const std::vector<std::byte> buf = pack(spec);

    ASSERT_GE(buf.size(), kFirstComponentIdOffset + sizeof(std::uint16_t));
    EXPECT_EQ(static_cast<char>(buf[kMagicOffset + 0]), 'A');
    EXPECT_EQ(static_cast<char>(buf[kMagicOffset + 1]), 'T');
    EXPECT_EQ(static_cast<char>(buf[kMagicOffset + 2]), 'V');
    EXPECT_EQ(static_cast<char>(buf[kMagicOffset + 3]), 'T');
    EXPECT_EQ(read_scalar<std::uint32_t>(buf, kVersionOffset),
              TaskSpec::kWireVersion);
    EXPECT_EQ(static_cast<std::uint8_t>(buf[kTaskKindOffset]),
              static_cast<std::uint8_t>(TaskKind::MD));
    EXPECT_EQ(read_scalar<std::uint32_t>(buf, kNEntriesOffset),
              static_cast<std::uint32_t>(spec.entries().size()));
    // First entry is NeighborListParams = variant index 1.
    EXPECT_EQ(read_scalar<std::uint16_t>(buf, kFirstComponentIdOffset),
              std::uint16_t{1});
}

// ── (b) SinglePoint roundtrip ───────────────────────────────────────────

TEST(TaskCodec, SinglePointRoundtripAllFiveAm1FieldsExact) {
    // All five fields away from the NDDOConfig defaults.
    AM1EvaluatorParams am1{};
    am1.charge = -1;
    am1.multiplicity = 2;
    am1.max_iter = 50;
    am1.density_tol = 1.0e-6;
    am1.energy_tol = 1.0e-8;

    const TaskSpec spec(am1);
    ASSERT_EQ(spec.kind(), TaskKind::SinglePoint);
    ASSERT_EQ(spec.entries().size(), 1u);
    EXPECT_FALSE(spec.entries()[0].binding.has_value());

    const std::vector<std::byte> buf = pack(spec);
    EXPECT_EQ(static_cast<std::uint8_t>(buf[kTaskKindOffset]),
              static_cast<std::uint8_t>(TaskKind::SinglePoint));

    const TaskSpec back = unpack(buf.data(), buf.size());
    EXPECT_EQ(back.kind(), TaskKind::SinglePoint);
    ASSERT_EQ(back.entries().size(), 1u);
    EXPECT_FALSE(back.entries()[0].binding.has_value());
    const auto& got = std::get<AM1EvaluatorParams>(back.entries()[0].params);
    EXPECT_EQ(got.charge, -1);
    EXPECT_EQ(got.multiplicity, 2);
    EXPECT_EQ(got.max_iter, 50);
    EXPECT_EQ(got.density_tol, 1.0e-6);
    EXPECT_EQ(got.energy_tol, 1.0e-8);
    EXPECT_TRUE(got == am1);
}

// ── (c) Version tamper ──────────────────────────────────────────────────

TEST(TaskCodec, VersionTamperThrowsNamingOffendingVersion) {
    std::vector<std::byte> buf = pack(make_full_md_spec());
    const std::uint32_t bad_version = TaskSpec::kWireVersion + 1;
    write_scalar<std::uint32_t>(buf, kVersionOffset, bad_version);
    expect_unpack_runtime_error_containing(buf, std::to_string(bad_version));
}

// ── (d) component_id / task_kind tamper ─────────────────────────────────

TEST(TaskCodec, ComponentIdTamperThrowsNamingOffendingId) {
    std::vector<std::byte> buf = pack(make_full_md_spec());
    // Ten alternatives → valid ids are [0, 9]; 999 is well outside.
    write_scalar<std::uint16_t>(buf, kFirstComponentIdOffset,
                                std::uint16_t{999});
    expect_unpack_runtime_error_containing(buf, "999");
}

TEST(TaskCodec, TaskKindTamperThrowsNamingOffendingValue) {
    std::vector<std::byte> buf = pack(make_full_md_spec());
    // TaskKind is {MD = 0, SinglePoint = 1}; 7 is outside the enum.
    buf[kTaskKindOffset] = std::byte{7};
    expect_unpack_runtime_error_containing(buf, "7");
}

// ── (e) Truncated buffer ────────────────────────────────────────────────

TEST(TaskCodec, TruncatedBufferThrows) {
    const std::vector<std::byte> buf = pack(make_full_md_spec());
    ASSERT_GT(buf.size(), 2u);
    EXPECT_THROW((void)unpack(buf.data(), buf.size() / 2),
                 std::runtime_error);
}

// ── (f) MD ctor validates n_steps ───────────────────────────────────────

TEST(TaskCodec, MdCtorRejectsNStepsBelowOne) {
    std::vector<ComponentEntry> entries;
    entries.push_back(
        ComponentEntry{ComponentParams{VelocityVerletParams{0.5}},
                       std::nullopt});
    EXPECT_THROW(TaskSpec(entries, /*n_steps=*/0), std::invalid_argument);
    EXPECT_THROW(TaskSpec(entries, /*n_steps=*/-5), std::invalid_argument);
}

// ── (g) MD ctor rejects evaluator-only params, naming the component ─────

TEST(TaskCodec, MdCtorRejectsAm1EntryNamingComponent) {
    std::vector<ComponentEntry> entries;
    entries.push_back(
        ComponentEntry{ComponentParams{VelocityVerletParams{0.5}},
                       std::nullopt});
    entries.push_back(ComponentEntry{ComponentParams{AM1EvaluatorParams{}},
                                     std::nullopt});
    try {
        TaskSpec spec(entries, /*n_steps=*/10);
        FAIL() << "MD TaskSpec accepted an AM1EvaluatorParams entry";
    } catch (const std::invalid_argument& e) {
        // The message must name the offending component.
        EXPECT_NE(lowercase(e.what()).find("am1"), std::string::npos)
            << "what() does not name the AM1 component: " << e.what();
    }
}

// ── (h) Kind-invariant violations on the wire ───────────────────────────

TEST(TaskCodec, KindFlipToSinglePointOnMultiEntryMdBufferThrows) {
    // A packed multi-entry MD task whose kind byte is flipped to SinglePoint
    // violates the exactly-one-entry invariant; decode goes through the
    // validating ctors, so unpack must throw.
    std::vector<std::byte> buf = pack(make_full_md_spec());
    buf[kTaskKindOffset] =
        static_cast<std::byte>(static_cast<std::uint8_t>(
            TaskKind::SinglePoint));
    EXPECT_ANY_THROW((void)unpack(buf.data(), buf.size()));
}

TEST(TaskCodec, KindFlipToSinglePointOnNonAm1SingleEntryBufferThrows) {
    // One entry, but not AM1EvaluatorParams: a SinglePoint task cannot carry
    // a non-AM1 component in v1.
    std::vector<ComponentEntry> entries;
    entries.push_back(ComponentEntry{
        ComponentParams{PairLJCutParams{0.0104, 3.4, 6.0}},
        std::nullopt});
    std::vector<std::byte> buf = pack(TaskSpec(std::move(entries), 3));
    buf[kTaskKindOffset] =
        static_cast<std::byte>(static_cast<std::uint8_t>(
            TaskKind::SinglePoint));
    EXPECT_ANY_THROW((void)unpack(buf.data(), buf.size()));
}

// ── (i) Empty-entry MD spec is a codec-level citizen ────────────────────

TEST(TaskCodec, EmptyEntryMdSpecRoundtrips) {
    // Driver-level component validation is the materializer's build() job,
    // not the codec's: an empty MD task packs and unpacks cleanly.
    const TaskSpec spec(std::vector<ComponentEntry>{}, /*n_steps=*/10);
    const std::vector<std::byte> buf = pack(spec);
    const TaskSpec back = unpack(buf.data(), buf.size());
    EXPECT_EQ(back.kind(), TaskKind::MD);
    EXPECT_TRUE(back.entries().empty());
    EXPECT_EQ(back.n_steps(), 10);
}

} // namespace
