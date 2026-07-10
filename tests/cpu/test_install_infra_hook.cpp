/// @file tests/cpu/test_install_infra_hook.cpp
/// @brief Phase 1 (driver-infra-auto-injection) — assert the foundation:
///        (a) `ComponentBase::reads()` defaults to `Capability::None`;
///        (b) `ComponentBase::needs_neighbor_list()` defaults to `false`;
///        (c) a CPU category-driver subclass can override
///            `install_defaults()` to inject a component, and the hook fires
///            *before* `validate_components()` so an injection can satisfy
///            an otherwise-unsatisfied validation rule.
///
/// Acceptance: T1 + T2 of `.claude/specs/driver-infra-auto-injection.md`.
///
/// arch-refactor-02-halo additions (`HaloParticipant` protocol scan):
///        (d) `ComponentBase::halo_participant()` defaults to `nullptr`; the
///            eight legacy halo virtuals are gone from the base vtable;
///        (e) `HaloSend`/`HaloRecv` expose `Endpoint::send`/`Endpoint::recv`
///            through the protocol, and a manually registered pair +
///            `set_pair` builds to exactly ONE endpoint pair (no double
///            injection);
///        (f) a pure field participant (`Endpoint::none`) is not treated as
///            an endpoint and does not suppress the automatic injection.

#include <gtest/gtest.h>

#include <type_traits>
#include <unordered_set>
#include <utility>
#include <vector>

#include "atomiverse/cpu/comm/halo_exchange.h"
#include "atomiverse/cpu/component/component.h"
#include "atomiverse/cpu/component/halo_participant.h"
#include "atomiverse/cpu/driver/dynamics_driver.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/md/md.h"
#include "atomiverse/cpu/neighbor/neighbor_list.h"
#include "atomiverse/cpu/potentials/pair_lj_cut.h"
#include "atomiverse/facade/component_role.h"
#include "atomiverse/facade/plan.h"
#include "atomiverse/facade/stage.h"

#include "cpu/driver_test_helpers.h"

using atv_test::make_trivial_frame;
using atv_test::StubForceEvaluator;
using atv_test::StubIntegrator;

namespace {

// ─── Foundation virtuals: defaults ──────────────────────────────────

TEST(InstallInfraHook, ReadsDefaultsToNone) {
    StubIntegrator c;
    EXPECT_EQ(c.reads(), ATV::Capability::None);
}

TEST(InstallInfraHook, NeedsNeighborListDefaultsToFalse) {
    StubIntegrator c;
    EXPECT_FALSE(c.needs_neighbor_list());
}

// ─── install_defaults() hook fires before validate_components() ────────

/// @brief Counts how often `install_defaults()` is called per `build()`.
class HookCounterDynamics : public ATV::cpu::DynamicsDriver {
  public:
    int install_calls = 0;
    StubIntegrator integ;

  protected:
    void install_defaults() override {
        ++install_calls;
        // Inject the integrator from inside the hook — driver's
        // `validate_components()` must see this component on the
        // pending list, otherwise validation throws "missing Integrator".
        this->on(ATV::md::Phase::initial_integrate).add(integ);
    }
};

TEST(InstallInfraHook, HookFiresBeforeValidateComponents) {
    HookCounterDynamics d;
    // Only register the force evaluator from outside; the integrator
    // arrives via install_defaults(). If validate_components() ran first,
    // the build would throw "missing Integrator".
    StubForceEvaluator force;
    d.on(ATV::md::Phase::force).add(force);

    EXPECT_NO_THROW(d.build());
    EXPECT_EQ(d.install_calls, 1);
}

TEST(InstallInfraHook, HookFiresExactlyOncePerBuild) {
    HookCounterDynamics d;
    StubForceEvaluator force;
    d.on(ATV::md::Phase::force).add(force);

    d.build();
    EXPECT_EQ(d.install_calls, 1);
}

// ─── Default install_defaults is no-op (existing tests still pass) ─────

TEST(InstallInfraHook, DefaultInstallInfraIsNoOp) {
    // The plain `cpu::DynamicsDriver` (no subclass override) must still
    // throw on missing Integrator — proves the base `install_defaults()` does
    // not silently inject anything.
    ATV::cpu::DynamicsDriver d;
    StubForceEvaluator force;
    d.on(ATV::md::Phase::force).add(force);
    EXPECT_THROW(d.build(), std::runtime_error);
}

// ═══ HaloParticipant protocol scan (arch-refactor-02-halo) ══════════════

/// @brief Pure field participant: opts into `HaloParticipant`, keeps the
/// default `Endpoint::none`, and declares one width-1 forward field. Must
/// neither count as a halo endpoint nor suppress endpoint auto-injection.
class PureFieldParticipant : public ATV::cpu::ComponentBase,
                             public ATV::cpu::HaloParticipant {
  public:
    std::string name() const override { return "pure_field_participant"; }
    ATV::cpu::HaloParticipant* halo_participant() override { return this; }
    std::vector<CommField> comm_forward_fields() const override {
        CommField f;
        f.data = nullptr; // storage is only allocated in build(); the
        f.width = 1;      // declaration shape alone drives the driver scan
        return {f};
    }
    void launch(const ATV::cpu::LaunchContext&) override {}
};

/// @brief Count distinct send/recv halo endpoints in a compiled Plan via the
/// protocol probe (`halo_participant()->halo_endpoint()`), deduplicating
/// components that appear under multiple phase bindings. No `dynamic_cast`,
/// no `name()` matching — same discipline as `scan_pending_infra_()`.
std::pair<int, int>
count_halo_endpoints(const ATV::Plan<ATV::cpu::ComponentBase>& plan) {
    std::unordered_set<const ATV::cpu::ComponentBase*> seen;
    int sends = 0;
    int recvs = 0;
    for (const auto& slot : plan.slots) {
        for (const auto& op : slot.ops) {
            ATV::cpu::ComponentBase* c = op.component;
            if (!seen.insert(c).second)
                continue;
            ATV::cpu::HaloParticipant* p = c->halo_participant();
            if (p == nullptr)
                continue;
            switch (p->halo_endpoint()) {
            case ATV::cpu::HaloParticipant::Endpoint::send:
                ++sends;
                break;
            case ATV::cpu::HaloParticipant::Endpoint::recv:
                ++recvs;
                break;
            case ATV::cpu::HaloParticipant::Endpoint::none:
                break;
            }
        }
    }
    return {sends, recvs};
}

// ─── (d) base-vtable contract ───────────────────────────────────────────

TEST(HaloParticipantProtocol, ComponentBaseDefaultsToNullParticipant) {
    StubIntegrator c;
    EXPECT_EQ(c.halo_participant(), nullptr);
}

// Detection idiom: the eight legacy halo virtuals must be GONE from
// `cpu::ComponentBase` (replaced by the single `halo_participant()` query).
#define ATV_TEST_HAS_MEMBER(memb)                                              \
    template <typename T, typename = void>                                     \
    struct HasMember_##memb : std::false_type {};                              \
    template <typename T>                                                      \
    struct HasMember_##                                                        \
        memb<T, std::void_t<decltype(&T::memb)>> : std::true_type {}

ATV_TEST_HAS_MEMBER(provides_halo_send);
ATV_TEST_HAS_MEMBER(provides_halo_recv);
ATV_TEST_HAS_MEMBER(comm_forward);
ATV_TEST_HAS_MEMBER(comm_reverse);
ATV_TEST_HAS_MEMBER(pack_forward);
ATV_TEST_HAS_MEMBER(unpack_forward);
ATV_TEST_HAS_MEMBER(pack_reverse);
ATV_TEST_HAS_MEMBER(unpack_reverse);

#undef ATV_TEST_HAS_MEMBER

TEST(HaloParticipantProtocol, LegacyHaloVirtualsAreGoneFromComponentBase) {
    using CB = ATV::cpu::ComponentBase;
    static_assert(!HasMember_provides_halo_send<CB>::value,
                  "provides_halo_send must be deleted from cpu::ComponentBase");
    static_assert(!HasMember_provides_halo_recv<CB>::value,
                  "provides_halo_recv must be deleted from cpu::ComponentBase");
    static_assert(!HasMember_comm_forward<CB>::value,
                  "comm_forward must be deleted from cpu::ComponentBase");
    static_assert(!HasMember_comm_reverse<CB>::value,
                  "comm_reverse must be deleted from cpu::ComponentBase");
    static_assert(!HasMember_pack_forward<CB>::value,
                  "pack_forward must be deleted from cpu::ComponentBase");
    static_assert(!HasMember_unpack_forward<CB>::value,
                  "unpack_forward must be deleted from cpu::ComponentBase");
    static_assert(!HasMember_pack_reverse<CB>::value,
                  "pack_reverse must be deleted from cpu::ComponentBase");
    static_assert(!HasMember_unpack_reverse<CB>::value,
                  "unpack_reverse must be deleted from cpu::ComponentBase");
    SUCCEED(); // contract is compile-time
}

// ─── (e) endpoint exposure + no double injection ────────────────────────

TEST(HaloParticipantProtocol, HaloSendRecvExposeEndpoints) {
    ATV::cpu::comm::HaloSend send{2.0};
    ATV::cpu::HaloParticipant* ps = send.halo_participant();
    ASSERT_NE(ps, nullptr);
    EXPECT_EQ(ps, static_cast<ATV::cpu::HaloParticipant*>(&send));
    EXPECT_EQ(ps->halo_endpoint(), ATV::cpu::HaloParticipant::Endpoint::send);

    ATV::cpu::comm::HaloRecv recv{};
    ATV::cpu::HaloParticipant* pr = recv.halo_participant();
    ASSERT_NE(pr, nullptr);
    EXPECT_EQ(pr, static_cast<ATV::cpu::HaloParticipant*>(&recv));
    EXPECT_EQ(pr->halo_endpoint(), ATV::cpu::HaloParticipant::Endpoint::recv);
}

TEST(HaloParticipantProtocol, ManualHaloPairIsNotDoubleInjected) {
    // User registers the halo pair explicitly (custom cutoff) AND sets a
    // pair evaluator. install_defaults() must detect the existing endpoints
    // via the protocol probe and inject nothing on top — reverse exchange
    // must run exactly once (see halo_exchange.h).
    ATV::cpu::md::MD sim;
    StubIntegrator integ;
    sim.on(ATV::md::Phase::initial_integrate).add(integ);
    sim.set_pair(ATV::cpu::PairLJCut{0.0104, 3.4, 6.0});
    sim.set_neighbor(ATV::cpu::md::NeighborList<>{6.0});
    sim.add(ATV::cpu::comm::HaloSend{6.6});
    sim.add(ATV::cpu::comm::HaloRecv{});
    sim.build();

    const auto [sends, recvs] = count_halo_endpoints(sim.plan());
    EXPECT_EQ(sends, 1) << "expected exactly one HaloSend endpoint "
                           "(manual registration must suppress auto-inject)";
    EXPECT_EQ(recvs, 1) << "expected exactly one HaloRecv endpoint "
                           "(manual registration must suppress auto-inject)";
}

// ─── (f) Endpoint::none neither endpoint nor injection suppressor ───────

TEST(HaloParticipantProtocol, FieldParticipantIsNotAnEndpoint) {
    PureFieldParticipant fp;
    ATV::cpu::HaloParticipant* p = fp.halo_participant();
    ASSERT_NE(p, nullptr);
    EXPECT_EQ(p->halo_endpoint(), ATV::cpu::HaloParticipant::Endpoint::none);
    EXPECT_TRUE(fp.comm_reverse_fields().empty()); // base default
    ASSERT_EQ(fp.comm_forward_fields().size(), 1u);
    EXPECT_EQ(fp.comm_forward_fields()[0].width, 1);
}

TEST(HaloParticipantProtocol, FieldParticipantDoesNotSuppressAutoInjection) {
    // A pure field participant is registered, but NO halo endpoints are.
    // Under an MPI build, install_defaults() must still auto-inject the
    // HaloSend/HaloRecv pair for the set_pair cutoff: Endpoint::none must
    // not trip the has_halo_send/has_halo_recv scan booleans.
    ATV::cpu::md::MD sim;
    StubIntegrator integ;
    PureFieldParticipant fp;
    sim.on(ATV::md::Phase::initial_integrate).add(integ);
    sim.set_pair(ATV::cpu::PairLJCut{0.0104, 3.4, 6.0});
    sim.set_neighbor(ATV::cpu::md::NeighborList<>{6.0});
    sim.on(ATV::md::Phase::force).add(fp, +100);
    sim.build();

    const auto [sends, recvs] = count_halo_endpoints(sim.plan());
#ifdef ATV_HAS_MPI
    EXPECT_EQ(sends, 1) << "Endpoint::none participant suppressed the "
                           "HaloSend auto-injection";
    EXPECT_EQ(recvs, 1) << "Endpoint::none participant suppressed the "
                           "HaloRecv auto-injection";
#else
    EXPECT_EQ(sends, 0) << "no MPI build: nothing may inject endpoints";
    EXPECT_EQ(recvs, 0) << "no MPI build: nothing may inject endpoints";
#endif
}

} // namespace
