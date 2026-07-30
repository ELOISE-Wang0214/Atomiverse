/// @file tests/cpu/semiempirical/test_am1_thermochemistry.cpp
/// @brief AM1 heat-of-formation constants vs MOPAC (one constant per test).
///
/// `ΔHf = (E_tot − Σ_A EISOL_A)·23.061 + Σ_A ΔHf_atom(A)` [kcal/mol]
/// (openmopac.net/Manual/calc_of_hof.html). The end-to-end value is pinned by
/// `regressions/semiempirical/am1-water-mopac-energy/`; these unit tests pin
/// the three *constant* inputs to that formula independently, so a ΔHf
/// regression failure separates "the SCF energy moved" from "a thermochemistry
/// constant was mis-transcribed".
///
/// MOPAC values: `src/models/parameters_for_AM1_C.F90` (EISOL) and the
/// `atheat` table; AM1 itself is Dewar, Zoebisch, Healy & Stewart,
/// J. Am. Chem. Soc. 107, 3902 (1985).

#include <gtest/gtest.h>

#include <stdexcept>

#include "atomiverse/cpu/semiempirical/am1_parameters.h"

namespace {

using ATV::cpu::semiempirical::am1_atomic_heat_kcal_mol;
using ATV::cpu::semiempirical::am1_parameter;
using ATV::cpu::semiempirical::kAm1EvToKcalMol;

TEST(TestAm1Thermochemistry, EvToKcalMolMatchesMopacFactor) {
    // MOPAC's default eV→kcal/mol factor is 2018 CODATA e·N_A·10⁻³/4.184
    // (`conref_C.F90` `fpcref(1,9)`); the legacy 23.061 needs `OLDCON`.
    EXPECT_DOUBLE_EQ(kAm1EvToKcalMol, 23.060547830619029);
}

TEST(TestAm1Thermochemistry, HydrogenAtomicHeatMatchesMopac) {
    EXPECT_DOUBLE_EQ(am1_atomic_heat_kcal_mol(1), 52.102);
}

TEST(TestAm1Thermochemistry, OxygenAtomicHeatMatchesMopac) {
    EXPECT_DOUBLE_EQ(am1_atomic_heat_kcal_mol(8), 59.559);
}

TEST(TestAm1Thermochemistry, AtomicHeatThrowsForUntabulatedElement) {
    // No silent fallback: an element with no tabulated ΔHf must throw rather
    // than contribute 0 kcal/mol to a heat of formation.
    EXPECT_THROW((void)am1_atomic_heat_kcal_mol(43), std::runtime_error);
}

TEST(TestAm1Thermochemistry, HydrogenEisolMatchesMopac) {
    const auto* h = am1_parameter(1);
    ASSERT_NE(h, nullptr);
    EXPECT_DOUBLE_EQ(h->eisol, -11.396427);
}

TEST(TestAm1Thermochemistry, OxygenEisolMatchesMopac) {
    const auto* o = am1_parameter(8);
    ASSERT_NE(o, nullptr);
    EXPECT_DOUBLE_EQ(o->eisol, -316.09952);
}

TEST(TestAm1Thermochemistry, WaterHeatOfFormationFormulaFromMopacTotalEnergy) {
    // Feed MOPAC's own printed TOTAL ENERGY through ATV's constants: the
    // formula plus the constants must land on MOPAC's printed ΔHf. This
    // isolates the thermochemistry conversion from the SCF energy — the
    // regression covers the energy itself.
    const auto* h = am1_parameter(1);
    const auto* o = am1_parameter(8);
    ASSERT_NE(h, nullptr);
    ASSERT_NE(o, nullptr);

    constexpr double kMopacTotalEnergyEv = -348.56254;
    const double sum_eisol = o->eisol + 2.0 * h->eisol;
    const double sum_atheat =
        am1_atomic_heat_kcal_mol(8) + 2.0 * am1_atomic_heat_kcal_mol(1);
    const double hof =
        (kMopacTotalEnergyEv - sum_eisol) * kAm1EvToKcalMol + sum_atheat;

    // Tolerance is the transcript's own precision: TOTAL ENERGY is printed to
    // 5 decimals (±5×10⁻⁶ eV ⇒ ±1.2×10⁻⁴ kcal/mol), and rounding it back
    // through the formula is what the residual here is.
    EXPECT_NEAR(hof, -59.24079, 5.0e-3);
}

} // namespace
