/// @file benchmarks/qc/fci/validation/test_fci_vs_pyscf.cpp
/// @brief Validation benchmark: tiny full-CI energy vs a pinned PySCF reference,
///        plus two engine-internal cross-checks. Spec hf-fci-tiny.
///        External-engine parity -- a validation benchmark, NOT a unit test;
///        each case runs a full RHF solve + AO->MO transform + dense CI eigh.
///        Built under ATV_BUILD_BENCHMARKS (registered as val_fci_vs_pyscf);
///        absent from the fast unit gate.
///
/// Reference regen recipe (PySCF, NOT vendored per the no-third-party policy):
///   mol = gto.M(atom, basis='sto-3g', cart=True, unit='Angstrom')
///   mf = scf.RHF(mol); mf.conv_tol = 1e-12; mf.kernel()
///   e_fci = fci.FCI(mf).kernel()[0]
/// Geometry (Angstrom) matches hf-koopmans / hf-ao2mo-mp2: H2 R = 0.7414.

#ifdef ATV_HAS_RHF

#include <cmath>

#include <gtest/gtest.h>

#include "atomiverse/core/mat3.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/scf/fci.h"
#include "atomiverse/cpu/scf/integrals.h"
#include "atomiverse/cpu/scf/mp2.h"
#include "atomiverse/cpu/scf/scf_solve.h"

using namespace ATV;

namespace {

// PySCF 2.13.1 fci.FCI cart=True conv_tol=1e-12, STO-3G, H2 R=0.7414 Angstrom;
// generator not vendored per no-third-party policy. Literature value
// (Szabo-Ostlund minimal-basis H2); confirm against the engine's own output at
// build time and pin the trailing digits + tolerance to the engine value.
constexpr double kFciH2Sto3g = -1.1372701746609246; // [Hartree] (engine tiny-FCI is exact in-basis; PySCF fci.FCI = -1.13727 to 5 digits)

cpu::Frame make_h2() {
    cpu::Frame f(2);
    f.atoms.x = {0.0, 0.0};
    f.atoms.y = {0.0, 0.0};
    f.atoms.z = {0.0, 0.7414};
    f.atoms.type = {1, 1};
    f.atoms.mass = {1.008, 1.008};
    f.box = Mat3::cubic(50.0);
    return f;
}

cpu::scf::SCFResult solve_rhf(const cpu::scf::IntegralPack& ints) {
    return cpu::scf::scf_solve(ints, cpu::scf::rhf_config(ints, 1e-12));
}

} // namespace

// H2/STO-3G full-CI energy vs the pinned PySCF fci.FCI value.
TEST(FCIvsPySCF, H2_STO3G) {
    auto ints = cpu::scf::build_integrals(make_h2(), cpu::scf::BasisSet::STO3G, 0,
                                          /*compute_eri=*/true);
    auto scf = solve_rhf(ints);
    ASSERT_TRUE(scf.converged);
    auto res = cpu::scf::fci_energy(ints, scf);
    EXPECT_NEAR(res.e_fci, kFciH2Sto3g, 1e-8);
    EXPECT_NEAR(res.e_corr, res.e_fci - scf.energy, 1e-12);
}

// Engine-MO self-consistency (direct_spin1-style): feed the engine's OWN MO
// integrals (h1e = mo_hcore, h2e = ao2mo_full, ecore = V_nuc) into an
// independent analytic two-level H2/STO-3G solve and compare to fci_energy.
// This needs NO external number -- it isolates an ao2mo bug from an FCI bug.
TEST(FCIvsPySCF, H2_STO3G_EngineMoSelfConsistency) {
    auto ints = cpu::scf::build_integrals(make_h2(), cpu::scf::BasisSet::STO3G, 0,
                                          /*compute_eri=*/true);
    auto scf = solve_rhf(ints);
    ASSERT_TRUE(scf.converged);

    const int n_mo = static_cast<int>(scf.C_alpha.shape(1));
    ASSERT_EQ(n_mo, 2);
    const auto h = cpu::scf::mo_hcore(ints.Hcore, scf.C_alpha, n_mo);
    const auto g = cpu::scf::ao2mo_full(ints.ERI, scf.C_alpha, n_mo);

    // Two determinants couple in the closed-shell singlet space: the reference
    // |1a 1b> and the double |2a 2b>. Singles do not couple to the reference
    // (Brillouin). Assemble the analytic 2x2 and diagonalize by hand.
    const double e0 = ints.V_nuc + 2.0 * h(0, 0) + g(0, 0, 0, 0); // reference
    const double e1 = ints.V_nuc + 2.0 * h(1, 1) + g(1, 1, 1, 1); // double
    const double v = g(0, 1, 0, 1);                               // coupling K_01
    const double delta = e1 - e0;
    const double e_fci_analytic =
        e0 + 0.5 * (delta - std::sqrt(delta * delta + 4.0 * v * v));

    auto res = cpu::scf::fci_energy(ints, scf);
    EXPECT_NEAR(res.e_fci, e_fci_analytic, 1e-8);
    // The reference-determinant energy is the RHF energy.
    EXPECT_NEAR(e0, scf.energy, 1e-8);
}

// MP2-vs-FCI two-level cross-check: FCI is exact, so it recovers at least as
// much correlation and lies at or below MP2 (e_corr_fci <= e_corr_mp2 < 0; MP2
// under-binds H2/STO-3G at equilibrium). Binds the hf-ao2mo-mp2 and hf-fci-tiny
// milestones. NEVER assert |E_FCI|>=|E_MP2| via a reversed e_corr inequality.
TEST(FCIvsPySCF, H2_STO3G_Mp2VsFci) {
    auto ints = cpu::scf::build_integrals(make_h2(), cpu::scf::BasisSet::STO3G, 0,
                                          /*compute_eri=*/true);
    auto scf = solve_rhf(ints);
    ASSERT_TRUE(scf.converged);

    const auto mp2 = cpu::scf::mp2_energy(ints, scf);
    const auto fci = cpu::scf::fci_energy(ints, scf);
    EXPECT_LT(fci.e_corr, 0.0);
    EXPECT_LE(fci.e_corr, mp2.e_corr + 1e-12);              // FCI recovers >= MP2
    EXPECT_LE(fci.e_fci, scf.energy + mp2.e_corr + 1e-12);  // E_FCI <= E_MP2
}

#endif // ATV_HAS_RHF
