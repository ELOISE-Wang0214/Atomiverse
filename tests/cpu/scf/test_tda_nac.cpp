#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/scf/methods.h"
#include <gtest/gtest.h>
using namespace ATV;
TEST(TdaNac, H2StrengthFinite) {
    // Need ≥2 roots: H2 STO-3G has only 1 virtual → n_ov = n_occ*n_virt = 1*1 = 1
    // Use larger system: He? He has no virt. Use LiH or H2 with better basis.
    // Pople 6-31G H2: more virtuals.
    cpu::Frame f(2);
    f.span(kAtom_x)[0] = 0;
    f.span(kAtom_y)[0] = 0;
    f.span(kAtom_z)[0] = 0;
    f.span(kAtom_x)[1] = 0.85;
    f.span(kAtom_y)[1] = 0;
    f.span(kAtom_z)[1] = 0;
    f.span(kAtom_type)[0] = f.span(kAtom_type)[1] = 1;
    f.box = cpu::Box::cubic(15);
    cpu::scf::RHF rhf;
    rhf.basis = cpu::scf::BasisSet::Pople631G;
    rhf.tol = 1e-8;
    rhf.run(f);
    ASSERT_TRUE(rhf.converged());
    auto es = rhf.tda(2);
    ASSERT_GE(es.energy_ev.size(), 2u);
    // NAC strength is expensive (many SCFs) — just call with large h
    double s = rhf.tda_nac_strength(f, 0, 1, 5e-3);
    EXPECT_TRUE(std::isfinite(s));
    EXPECT_GT(s, 0.0);
}
