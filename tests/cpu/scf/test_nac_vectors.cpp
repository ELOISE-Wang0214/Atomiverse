/// @file tests/cpu/scf/test_nac_vectors.cpp
/// @brief True TDA NAC vectors: shape, antisymmetry, norm vs strength.
#include <cmath>
#include <gtest/gtest.h>

#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/scf/methods.h"

using namespace ATV;

TEST(NacVectors, ShapeAndAntisymmetry) {
    // H2 / 6-31G has ≥2 TDA roots.
    cpu::Frame f(2);
    f.span(kAtom_x)[0] = 0;
    f.span(kAtom_y)[0] = 0;
    f.span(kAtom_z)[0] = 0;
    f.span(kAtom_x)[1] = 0.85;
    f.span(kAtom_y)[1] = 0;
    f.span(kAtom_z)[1] = 0;
    f.span(kAtom_type)[0] = f.span(kAtom_type)[1] = 1;
    f.simbox = cpu::Box::cubic(15);

    cpu::scf::RHF rhf;
    rhf.basis = cpu::scf::BasisSet::Pople631G;
    rhf.tol = 1e-8;
    rhf.run(f);
    ASSERT_TRUE(rhf.converged());
    auto es = rhf.tda(2);
    ASSERT_GE(es.energy_ev.size(), 2u);

    auto dij = rhf.tda_nac_vectors(f, 0, 1);
    auto dji = rhf.tda_nac_vectors(f, 1, 0);
    ASSERT_EQ(dij.size(), 6u);
    ASSERT_EQ(dji.size(), 6u);
    for (std::size_t k = 0; k < dij.size(); ++k) {
        EXPECT_TRUE(std::isfinite(dij[k]));
        // Antisymmetry: d_ji = −d_ij. Independent FD force runs inject noise.
        EXPECT_NEAR(dji[k], -dij[k], 0.05);
    }
    // Phase: first large component of d_01 ≥ 0
    for (double v : dij) {
        if (std::fabs(v) > 1e-8) {
            EXPECT_GE(v, 0.0);
            break;
        }
    }
    double n2 = 0;
    for (double v : dij)
        n2 += v * v;
    EXPECT_GT(std::sqrt(n2), 0.0);
    // strength is ‖d‖ — same formula; allow FD noise across independent runs.
    const double strength = rhf.tda_nac_strength(f, 0, 1);
    EXPECT_NEAR(std::sqrt(n2), strength, 0.1 * std::max(1.0, std::sqrt(n2)));
}

TEST(NacVectors, ThrowsOnSameState) {
    cpu::Frame f(2);
    f.span(kAtom_x)[0] = 0;
    f.span(kAtom_y)[0] = 0;
    f.span(kAtom_z)[0] = 0;
    f.span(kAtom_x)[1] = 0.85;
    f.span(kAtom_y)[1] = 0;
    f.span(kAtom_z)[1] = 0;
    f.span(kAtom_type)[0] = f.span(kAtom_type)[1] = 1;
    f.simbox = cpu::Box::cubic(15);
    cpu::scf::RHF rhf;
    rhf.basis = cpu::scf::BasisSet::Pople631G;
    rhf.tol = 1e-8;
    rhf.run(f);
    EXPECT_THROW(rhf.tda_nac_vectors(f, 0, 0), std::invalid_argument);
}
