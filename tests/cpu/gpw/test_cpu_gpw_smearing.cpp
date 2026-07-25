/// @file test_cpu_gpw_smearing.cpp
/// @brief Unit tests for GPW metallic occupation smearing API + fill math.
///
/// Full multi-k SCF is the regression `reg_gpw-metallic-smearing`. Here: public
/// members (no Config bag), build() fail-fast, and reuse of pw::internal
/// occupation fill (FD recovers Ne).

#include "atomiverse/cpu/gpw/gpw_solver.h"
#include "atomiverse/cpu/planewave/internal/smearing.h"

#include <cmath>
#include <gtest/gtest.h>
#include <stdexcept>
#include <vector>

using namespace ATV;
using namespace ATV::cpu::gpw;
using ATV::cpu::pw::internal::SmearKind;
using ATV::cpu::pw::internal::fill_occupations;
using ATV::cpu::pw::internal::find_fermi_level;
using ATV::cpu::pw::internal::smear_f_band;

TEST(GpwSmearing, PublicMembersNoConfigBag) {
    GpwDft m;
    m.smearing = Smearing::fermi_dirac;
    m.degauss = 0.02;
    m.mp_order = 1;
    m.extra_bands = 4;
    m.n_occ = 9;
    EXPECT_EQ(m.smearing, Smearing::fermi_dirac);
    EXPECT_DOUBLE_EQ(m.degauss, 0.02);
    EXPECT_EQ(m.extra_bands, 4);
    GpwResult r;
    EXPECT_DOUBLE_EQ(r.e_fermi, 0.0);
    EXPECT_DOUBLE_EQ(r.e_entropy, 0.0);
    EXPECT_DOUBLE_EQ(r.n_elec_occ, 0.0);
}

TEST(GpwSmearing, BuildRejectsUksAndOtAndBadDegauss) {
    std::vector<GpwAtom> atoms; // empty → empty-atom reject first
    GpwDft empty;
    empty.n_occ = 1;
    empty.smearing = Smearing::fermi_dirac;
    empty.degauss = 0.02;
    EXPECT_THROW(empty.build(atoms), std::invalid_argument);

    // Minimal non-empty atom for later gates (geometry irrelevant for build).
    // Use explicit GpwAtom fields without full basis when only build() is
    // exercised — still need non-empty vector; GpwAtom default is enough for
    // emptiness check only. For UKS/OT we need at least one atom.
    GpwAtom a;
    a.z_val = 1.0;
    a.z = 1;
    atoms.push_back(a);

    GpwDft uks;
    uks.n_alpha = 1;
    uks.n_beta = 0;
    uks.cell[0] = uks.cell[1] = uks.cell[2] = 10.0;
    uks.smearing = Smearing::fermi_dirac;
    uks.degauss = 0.02;
    EXPECT_THROW(uks.build(atoms), std::invalid_argument);

    GpwDft ot;
    ot.n_occ = 1;
    ot.cell[0] = ot.cell[1] = ot.cell[2] = 10.0;
    ot.use_ot = true;
    ot.smearing = Smearing::fermi_dirac;
    ot.degauss = 0.02;
    EXPECT_THROW(ot.build(atoms), std::invalid_argument);

    GpwDft bad_deg;
    bad_deg.n_occ = 1;
    bad_deg.cell[0] = bad_deg.cell[1] = bad_deg.cell[2] = 10.0;
    bad_deg.smearing = Smearing::fermi_dirac;
    bad_deg.degauss = 0.0;
    EXPECT_THROW(bad_deg.build(atoms), std::invalid_argument);
}

TEST(GpwSmearing, FermiDiracFillRecoversNelec) {
    // Shared with PW unit tests: three levels, want Ne=3, band form.
    std::vector<AccReal> evals = {AccReal(-1), AccReal(0), AccReal(1),
                                  AccReal(2)};
    std::vector<AccReal> w(evals.size(), AccReal(1));
    const AccReal ne = AccReal(3);
    const AccReal deg = AccReal(0.1);
    const AccReal mu =
        find_fermi_level(evals, w, ne, deg, SmearKind::fermi_dirac, 1);
    auto fill = fill_occupations(evals, w, mu, deg, SmearKind::fermi_dirac, 1);
    EXPECT_NEAR(static_cast<double>(fill.n_occ), 3.0, 1e-8);
    EXPECT_NEAR(static_cast<double>(
                    smear_f_band(AccReal(0), AccReal(0), deg,
                                 SmearKind::fermi_dirac)),
                1.0, 1e-12);
}

TEST(GpwSmearing, MethfesselPaxtonMemberAccepted) {
    GpwDft m;
    m.smearing = Smearing::methfessel_paxton;
    m.mp_order = 1;
    m.degauss = 0.02;
    m.n_occ = 1;
    m.cell[0] = m.cell[1] = m.cell[2] = 8.0;
    GpwAtom a;
    a.z_val = 1.0;
    a.z = 1;
    std::vector<GpwAtom> atoms{a};
    // build validates configuration; does not require a complete basis set.
    EXPECT_NO_THROW(m.build(atoms));
}
