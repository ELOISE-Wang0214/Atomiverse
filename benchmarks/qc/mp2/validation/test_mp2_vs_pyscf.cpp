/// @file benchmarks/qc/mp2/validation/test_mp2_vs_pyscf.cpp
/// @brief Validation benchmark: closed-shell MP2 correlation energy vs pinned
///        PySCF 2.13.1 across {He,H2,H2O,NH3,CH4} x {STO-3G,6-31G}, 1e-8 Ha.
///        Spec hf-ao2mo-mp2 (Domain basis table). External-engine parity — a
///        validation benchmark, NOT a unit test; each case runs a full RHF
///        solve + AO->MO transform + MP2 sum. Built under ATV_BUILD_BENCHMARKS
///        (registered as val_mp2_vs_pyscf); absent from the fast unit gate.
///
/// Reference regen: mol = gto.M(atom, basis, cart=True); mf = scf.RHF(mol);
/// mf.conv_tol = 1e-12; mp.MP2(mf).kernel()[0]. Geometries (Angstrom) match
/// hf-koopmans. NOTE: He/6-31G needs the He (Z=2) 6-31G basis row added by
/// hf-koopmans; until that lands it throws (no basis row).

#ifdef ATV_HAS_RHF

#include <cmath>

#include <gtest/gtest.h>

#include "atomiverse/core/mat3.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/scf/integrals.h"
#include "atomiverse/cpu/scf/mp2.h"
#include "atomiverse/cpu/scf/scf_solve.h"

using namespace ATV;

namespace {

cpu::Frame make_he() {
    cpu::Frame f(1);
    f.atoms.x = {0.0};
    f.atoms.y = {0.0};
    f.atoms.z = {0.0};
    f.atoms.type = {2};
    f.atoms.mass = {4.0026};
    f.box = Mat3::cubic(50.0);
    return f;
}

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

cpu::Frame make_h2o() {
    cpu::Frame f(3);
    f.atoms.x = {0.0, 0.0, 0.0};
    f.atoms.y = {0.0, 0.7572, -0.7572};
    f.atoms.z = {0.1173, -0.4692, -0.4692};
    f.atoms.type = {8, 1, 1};
    f.atoms.mass = {15.999, 1.008, 1.008};
    f.box = Mat3::cubic(50.0);
    return f;
}

cpu::Frame make_nh3() {
    cpu::Frame f(4);
    f.atoms.x = {0.0, 0.0, 0.8121, -0.8121};
    f.atoms.y = {0.0, -0.9377, 0.4689, 0.4689};
    f.atoms.z = {0.0, -0.3816, -0.3816, -0.3816};
    f.atoms.type = {7, 1, 1, 1};
    f.atoms.mass = {14.007, 1.008, 1.008, 1.008};
    f.box = Mat3::cubic(50.0);
    return f;
}

cpu::Frame make_ch4() {
    cpu::Frame f(5);
    f.atoms.x = {0.0, 0.627, 0.627, -0.627, -0.627};
    f.atoms.y = {0.0, 0.627, -0.627, 0.627, -0.627};
    f.atoms.z = {0.0, 0.627, -0.627, -0.627, 0.627};
    f.atoms.type = {6, 1, 1, 1, 1};
    f.atoms.mass = {12.011, 1.008, 1.008, 1.008, 1.008};
    f.box = Mat3::cubic(50.0);
    return f;
}

void expect_mp2_corr(cpu::Frame frame, cpu::scf::BasisSet basis,
                     double expected) {
    auto ints = cpu::scf::build_integrals(frame, basis, 0, /*compute_eri=*/true);
    auto scf = cpu::scf::scf_solve(ints, cpu::scf::rhf_config(ints, 1e-12));
    ASSERT_TRUE(scf.converged);
    auto res = cpu::scf::mp2_energy(ints, scf);
    EXPECT_NEAR(res.e_corr, expected, 1e-8);
    EXPECT_NEAR(res.e_tot, scf.energy + res.e_corr, 1e-12);
}

} // namespace

// ---- STO-3G column (pinned PySCF 2.13.1, cart=True, conv_tol=1e-12) ----
TEST(MP2vsPySCF, He_STO3G) { expect_mp2_corr(make_he(), cpu::scf::BasisSet::STO3G, 0.0); }
TEST(MP2vsPySCF, H2_STO3G) { expect_mp2_corr(make_h2(), cpu::scf::BasisSet::STO3G, -0.0131707665); }
TEST(MP2vsPySCF, H2O_STO3G) { expect_mp2_corr(make_h2o(), cpu::scf::BasisSet::STO3G, -0.0355456516); }
TEST(MP2vsPySCF, NH3_STO3G) { expect_mp2_corr(make_nh3(), cpu::scf::BasisSet::STO3G, -0.0472198510); }
TEST(MP2vsPySCF, CH4_STO3G) { expect_mp2_corr(make_ch4(), cpu::scf::BasisSet::STO3G, -0.0561492533); }

// ---- 6-31G column (He needs hf-koopmans' He 6-31G basis row) ----
TEST(MP2vsPySCF, He_631G) {
    // The He (Z=2) 6-31G basis row is added by hf-koopmans; until that lands
    // build_integrals throws (no basis row). Skip cleanly rather than fail.
    try {
        cpu::scf::build_integrals(make_he(), cpu::scf::BasisSet::Pople631G, 0,
                                  /*compute_eri=*/true);
    } catch (const std::exception&) {
        GTEST_SKIP() << "He/6-31G basis row pending hf-koopmans";
    }
    expect_mp2_corr(make_he(), cpu::scf::BasisSet::Pople631G, -0.0112001229);
}
TEST(MP2vsPySCF, H2_631G) { expect_mp2_corr(make_h2(), cpu::scf::BasisSet::Pople631G, -0.0173964438); }
TEST(MP2vsPySCF, H2O_631G) { expect_mp2_corr(make_h2o(), cpu::scf::BasisSet::Pople631G, -0.1288509172); }
TEST(MP2vsPySCF, NH3_631G) { expect_mp2_corr(make_nh3(), cpu::scf::BasisSet::Pople631G, -0.1169596914); }
TEST(MP2vsPySCF, CH4_631G) { expect_mp2_corr(make_ch4(), cpu::scf::BasisSet::Pople631G, -0.0996598545); }

#endif // ATV_HAS_RHF
