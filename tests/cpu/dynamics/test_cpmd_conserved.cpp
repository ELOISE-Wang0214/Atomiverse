/// @file tests/cpu/dynamics/test_cpmd_conserved.cpp
/// @brief esmm-cpmd-03-observability — `cpu::md::CpmdConserved` assembles the
///        Car–Parrinello conserved energy onto `frame.meta[kConserved]` [eV]:
///
///            E_cons = kEnergy + kKinetic
///                   + kEnergyElectronKinetic + kEnergyElectronThermostat
///
/// Unit scope: the assembler only. Nothing here runs GPW, an SCF, or an MD
/// step — the four inputs are pre-seeded on a one-atom Frame and `launch()` is
/// driven through a hand-built `LaunchContext` (the `test_multipole.cpp`
/// pattern), so the test pins arithmetic and contract, not physics. The
/// producer side (every CPMD branch writing both electronic keys each launch,
/// T_e matching `last_electronic_ke_ev()`) is pinned in
/// tests/cpu/gpw/test_gpw_cpmd.cpp; the drift envelope is
/// regressions/dynamics/cpmd-water.
///
/// Why the throw matters (AC-003): with no CPMD-enabled evaluator in the plan
/// the two electronic keys are simply absent, and an assembler that treated
/// "absent" as 0 would publish a *plausible* E_cons that silently omits the
/// fictitious electron energy — i.e. a conserved quantity that looks conserved
/// because the drifting term was dropped. Fail fast, no fallback.
///
/// Determinism: no RNG, no wall clock, no filesystem, no network. Every seeded
/// value is exactly representable in binary64 and the expected sum
/// (−10 + 2.5 + 0.25 + 0.125 = −7.125) is exact, so `EXPECT_DOUBLE_EQ` is the
/// right comparison for this energy — it pins the assembly, not a tolerance.
///
/// Spec: `.claude/specs/esmm-cpmd-03-observability/spec.md` (Design).

#include <gtest/gtest.h>

#include <algorithm>
#include <cctype>
#include <exception>
#include <stdexcept>
#include <string>

#include "atomiverse/core/binding.h"
#include "atomiverse/core/frame_schema.h"
#include "atomiverse/core/stage.h"
#include "atomiverse/cpu/component/component.h"
#include "atomiverse/cpu/dynamics/cpmd_conserved.h"
#include "atomiverse/cpu/frame/frame.h"

using namespace ATV;

namespace {

// The four seeded terms [eV]. Deliberately distinct magnitudes and signs so a
// dropped, doubled or swapped term cannot land on the same total.
constexpr double kEks = -10.0;     ///< kEnergy — Kohn–Sham potential energy
constexpr double kTNuc = 2.5;      ///< kKinetic — nuclear kinetic energy
constexpr double kTe = 0.25;       ///< kEnergyElectronKinetic — fictitious T_e
constexpr double kEThermo = 0.125; ///< kEnergyElectronThermostat

bool message_mentions(const std::exception& e, const std::string& needle) {
    return std::string(e.what()).find(needle) != std::string::npos;
}

/// Case-insensitive substring probe: the diagnostic may spell the method
/// "CPMD", "Cpmd" or "cpmd_conserved" — what is pinned is that it names the
/// method at all, not its capitalization.
bool message_mentions_ci(const std::exception& e, std::string needle) {
    std::string haystack(e.what());
    auto lower = [](std::string& s) {
        std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) {
            return static_cast<char>(std::tolower(c));
        });
    };
    lower(haystack);
    lower(needle);
    return haystack.find(needle) != std::string::npos;
}

/// @brief One-atom frame with the two *pre-inserted* energy terms seeded.
///
/// `Frame(int)` pre-inserts kEnergy / kKinetic / kConserved (frame.cpp:484-492)
/// but NOT the two electronic keys — those are open-schema, insert-guarded by
/// their producer (the `kDipole` pattern). Whether they exist is exactly what
/// separates the two cases below, so this helper never inserts them.
cpu::Frame seeded_frame() {
    cpu::Frame f(1);
    f.meta.set(kEnergy, kEks);
    f.meta.set(kKinetic, kTNuc);
    f.meta.set(kConserved, 0.0); // must be overwritten, not accumulated into
    return f;
}

/// A LaunchContext at the component's own stage, so an `at(end_of_step)` guard
/// inside `launch()` (if the implementation adds one) is satisfied rather than
/// silently skipping the work under test.
struct EndOfStepLaunch {
    cpu::Context ctx{};
    cpu::LaunchContext lc;

    explicit EndOfStepLaunch(cpu::Frame& f)
        : lc{ctx, f, ATV::slot_of(ATV::md::Phase::end_of_step), /*step=*/0} {}
};

} // namespace

// ─── The assembly ───────────────────────────────────────────────────────

/// All four terms, once each, onto kConserved. The pre-existing kConserved
/// value (0 here) is replaced: the component publishes E_cons, it does not
/// accumulate into whatever an earlier owner left behind.
TEST(CpmdConserved, SumsAllFourTermsOntoConserved) {
    auto f = seeded_frame();
    f.meta.insert(kEnergyElectronKinetic, kTe);
    f.meta.insert(kEnergyElectronThermostat, kEThermo);

    cpu::md::CpmdConserved c;
    EndOfStepLaunch launch{f};
    c.launch(launch.lc);

    EXPECT_DOUBLE_EQ(f.meta.get(kConserved), kEks + kTNuc + kTe + kEThermo);
    // The inputs are read-only to this component.
    EXPECT_DOUBLE_EQ(f.meta.get(kEnergy), kEks);
    EXPECT_DOUBLE_EQ(f.meta.get(kKinetic), kTNuc);
    EXPECT_DOUBLE_EQ(f.meta.get(kEnergyElectronKinetic), kTe);
    EXPECT_DOUBLE_EQ(f.meta.get(kEnergyElectronThermostat), kEThermo);
}

// ─── No fallback when the CPMD channel never ran (AC-003) ───────────────

/// Neither electronic key present — no CPMD-enabled evaluator in the plan.
/// The message must name the missing key (so the user knows *what* is absent)
/// AND the method (so they know *who* wanted it): a bare
/// "Frame meta has no key 'energy_electron_kinetic'" escaping from `meta.get`
/// leaves the reader with no route back to the mis-wired plan.
TEST(CpmdConserved, ThrowsWhenElectronicKeysAbsent) {
    auto f = seeded_frame();
    ASSERT_FALSE(f.meta.contains(kEnergyElectronKinetic.key));
    ASSERT_FALSE(f.meta.contains(kEnergyElectronThermostat.key));

    cpu::md::CpmdConserved c;
    EndOfStepLaunch launch{f};

    EXPECT_THROW(c.launch(launch.lc), std::runtime_error);
    try {
        c.launch(launch.lc);
        FAIL() << "expected CpmdConserved::launch to refuse a frame with no "
                  "electronic energy keys";
    } catch (const std::runtime_error& e) {
        EXPECT_TRUE(message_mentions(e, "energy_electron_kinetic")) << e.what();
        EXPECT_TRUE(message_mentions_ci(e, "cpmd")) << e.what();
    }
}

/// The second key is the easy one to forget. Seeding only T_e must still
/// throw — a thermostat term treated as an implicit 0 is the NVE-electrons
/// answer silently substituted for a thermostatted run.
TEST(CpmdConserved, ThrowsWhenThermostatKeyAloneAbsent) {
    auto f = seeded_frame();
    f.meta.insert(kEnergyElectronKinetic, kTe);
    ASSERT_FALSE(f.meta.contains(kEnergyElectronThermostat.key));

    cpu::md::CpmdConserved c;
    EndOfStepLaunch launch{f};

    EXPECT_THROW(c.launch(launch.lc), std::runtime_error);
    try {
        c.launch(launch.lc);
        FAIL() << "expected CpmdConserved::launch to refuse a frame missing "
                  "the electronic thermostat energy";
    } catch (const std::runtime_error& e) {
        EXPECT_TRUE(message_mentions(e, "energy_electron_thermostat"))
            << e.what();
        EXPECT_TRUE(message_mentions_ci(e, "cpmd")) << e.what();
    }
    // A refused launch publishes nothing.
    EXPECT_DOUBLE_EQ(f.meta.get(kConserved), 0.0);
}

// ─── Wiring contract ────────────────────────────────────────────────────

/// end_of_step @ -90: after `KineticEnergy` (@ -100), whose half-step-fresh
/// kKinetic this component consumes, and before `mpi::Reduce` (@ -50) /
/// Thermo — the only slot where all four terms are current. Exactly one
/// binding: this is an MD-phase assembler, it has no business in another
/// phase.
TEST(CpmdConserved, DefaultBindingIsEndOfStepMinus90) {
    cpu::md::CpmdConserved c;
    const auto b = c.default_bindings();
    ASSERT_EQ(b.size(), 1u);
    EXPECT_EQ(b[0].kind, ATV::PhaseKind::MD);
    EXPECT_EQ(b[0].phase, ATV::slot_of(ATV::md::Phase::end_of_step));
    EXPECT_EQ(b[0].order, -90);
}

/// House convention: the component name is the snake_case class name
/// (`kinetic_energy`, `nose_hoover`, `fix_csvr`). It reaches users through
/// Thermo/timing labels and plan-validation messages, so it is a pinned
/// string, not an incidental one.
TEST(CpmdConserved, NameIsStable) {
    cpu::md::CpmdConserved c;
    EXPECT_EQ(c.name(), "cpmd_conserved");
}
