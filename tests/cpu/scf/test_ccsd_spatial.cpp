/// @file tests/cpu/scf/test_ccsd_spatial.cpp
/// @brief Unit checks for the closed-shell spatial-orbital CCSD kernel
///        (internal/ccsd_spatial.h, ATV::cpu::scf::ccsd_energy_spatial).
///
/// scf-review-08-ccsd-spatial-diis, deliverable B. `ccsd_energy_spatial` solves
/// closed-shell CCSD in the SPATIAL MO basis (n_o²·n_v² storage, no spin
/// expansion), reusing `ao2mo_full` + `CcAmplitudeDiis`. It is the internal
/// computing function the closed-shell energy entry routes to.
///
/// This test builds the same private inputs the public path builds
/// (`internal::IntegralWorkspace{frame, …, compute_eri=true}` + the converged RHF
/// determinant + the Hartree SCF energy) and calls the kernel in isolation.
/// Written RED-first: internal/ccsd_spatial.h does not exist yet.
///
/// Reference: Scuseria, Janssen, Schaefer, J. Chem. Phys. 89, 7382 (1988).

#ifdef ATV_HAS_HF

#include <gtest/gtest.h>

#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/scf/internal/ccsd_spatial.h"
#include "atomiverse/cpu/scf/internal/integrals.h"
#include "atomiverse/cpu/scf/methods.h"

using namespace ATV;

namespace {

cpu::Frame make_he() {
    cpu::Frame f(1);
    f.span(kAtom_x)[0] = 0.0;
    f.span(kAtom_y)[0] = 0.0;
    f.span(kAtom_z)[0] = 0.0;
    f.span(kAtom_type)[0] = 2;
    f.span(kAtom_mass)[0] = 4.0026;
    f.box = cpu::Box::cubic(20.0);
    return f;
}

cpu::Frame make_h2() {
    cpu::Frame f(2);
    f.span(kAtom_x)[0] = 0.0;
    f.span(kAtom_x)[1] = 0.74;
    f.span(kAtom_y)[0] = 0.0;
    f.span(kAtom_y)[1] = 0.0;
    f.span(kAtom_z)[0] = 0.0;
    f.span(kAtom_z)[1] = 0.0;
    f.span(kAtom_type)[0] = 1;
    f.span(kAtom_type)[1] = 1;
    f.span(kAtom_mass)[0] = 1.008;
    f.span(kAtom_mass)[1] = 1.008;
    f.box = cpu::Box::cubic(20.0);
    return f;
}

} // namespace

// (a) domain: H2/STO-3G closed-shell CCSD correlation. The spatial kernel must
//     reproduce the same PySCF 2.14.0 pin as regressions/scf/ccsd-h2-sto3g
//     (conv_tol 1e-12). Hard-coded golden, ±1e-9 Ha; converged == true.
TEST(TestCcsdSpatial, H2Sto3gClosedShellCorrelation) {
    auto frame = make_h2();
    cpu::scf::RHF rhf;
    rhf.basis = cpu::scf::BasisSet::STO3G;
    rhf.tol = 1e-12;
    rhf.max_iter = 100;
    rhf.run(frame);
    ASSERT_TRUE(rhf.converged());

    const cpu::scf::internal::IntegralWorkspace ints{frame, cpu::scf::BasisSet::STO3G,
                                      /*charge=*/0, /*compute_eri=*/true};
    const cpu::scf::CCSDResult r = cpu::scf::ccsd_energy_spatial(
        ints, rhf.determinant(), rhf.energy(), /*tol=*/1e-10, /*max_iter=*/100);

    EXPECT_TRUE(r.converged);
    EXPECT_NEAR(r.e_corr, -0.020524527092106, 1e-9); // hard-coded golden
    EXPECT_NEAR(r.e_tot, rhf.energy() + r.e_corr, 1e-12);
}

// (b) edge: He/STO-3G has a single basis function → n_vir = 0. With no virtual
//     orbitals there are no excitations, so E_corr must be exactly 0 and the
//     kernel must report converged == true.
TEST(TestCcsdSpatial, HeSto3gZeroVirtualsZeroCorrelation) {
    auto frame = make_he();
    cpu::scf::RHF rhf;
    rhf.basis = cpu::scf::BasisSet::STO3G;
    rhf.tol = 1e-12;
    rhf.run(frame);
    ASSERT_TRUE(rhf.converged());

    const cpu::scf::internal::IntegralWorkspace ints{frame, cpu::scf::BasisSet::STO3G,
                                      /*charge=*/0, /*compute_eri=*/true};
    const cpu::scf::CCSDResult r = cpu::scf::ccsd_energy_spatial(
        ints, rhf.determinant(), rhf.energy(), /*tol=*/1e-10, /*max_iter=*/100);

    EXPECT_TRUE(r.converged);
    EXPECT_NEAR(r.e_corr, 0.0, 1e-12);
    EXPECT_NEAR(r.e_tot, rhf.energy(), 1e-12);
}

#endif // ATV_HAS_HF
