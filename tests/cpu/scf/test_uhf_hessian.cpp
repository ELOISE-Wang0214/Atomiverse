/// @file tests/cpu/scf/test_uhf_hessian.cpp
/// @brief UHF::hessian micro first-gate: OH radical / STO-3G analytic
///        Hessian vs central FD of the analytic forces — both in-house.
///        Structural invariants (symmetry, translational sum rules).

#ifdef ATV_HAS_HF

#include <cmath>
#include <stdexcept>

#include <gtest/gtest.h>
#include <xtensor/containers/xtensor.hpp>

#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/scf/uhf.h"
#include "atomiverse/cpu/vibration/hessian_result.h"

using namespace ATV;
using cpu::scf::BasisSet;

namespace {

/// OH radical, collinear along a skew axis so all Cartesian blocks engage.
cpu::Frame make_oh() {
    cpu::Frame f(2);
    f.span(kAtom_type)[0] = 8; // O
    f.span(kAtom_x)[0] = 0.0;
    f.span(kAtom_y)[0] = 0.0;
    f.span(kAtom_z)[0] = 0.0;
    f.span(kAtom_type)[1] = 1; // H
    f.span(kAtom_x)[1] = 0.20;
    f.span(kAtom_y)[1] = 0.30;
    f.span(kAtom_z)[1] = 0.90;
    f.simbox = cpu::Box::cubic(20.0);
    return f;
}

cpu::Frame displaced(const cpu::Frame& frame, int coord, double dx) {
    cpu::Frame f = frame;
    const int A = coord / 3, d = coord % 3;
    if (d == 0)
        f.span(kAtom_x)[A] += static_cast<StorageReal>(dx);
    else if (d == 1)
        f.span(kAtom_y)[A] += static_cast<StorageReal>(dx);
    else
        f.span(kAtom_z)[A] += static_cast<StorageReal>(dx);
    return f;
}

xt::xtensor<double, 1> forces_at(const cpu::Frame& frame) {
    cpu::scf::UHF uhf(5, 4); // OH: 9 electrons, doublet
    uhf.basis = BasisSet::STO3G;
    uhf.tol = 1e-13;
    cpu::Frame f = frame;
    const auto g = uhf.forces(f);
    const int natm = static_cast<int>(g.forces.shape()[0]);
    xt::xtensor<double, 1> out =
        xt::zeros<double>({static_cast<std::size_t>(3 * natm)});
    for (int A = 0; A < natm; ++A)
        for (int d = 0; d < 3; ++d)
            out(3 * A + d) = g.forces(A, d);
    return out;
}

} // namespace

TEST(UhfAnalyticHessian, OhMatchesFdOfAnalyticForces) {
    const auto frame = make_oh();

    cpu::scf::UHF uhf(5, 4);
    uhf.basis = BasisSet::STO3G;
    uhf.tol = 1e-13;
    cpu::Frame f = frame;
    const cpu::vib::HessianResult res = uhf.hessian(f);
    ASSERT_EQ(res.n_atoms, 2);
    ASSERT_EQ(res.hessian.shape()[0], 6u);
    ASSERT_EQ(res.units, "eV/Angstrom^2");

    // H(x, y) = −∂F_y/∂x — 5-point central FD per coordinate.
    const double h = 2e-3; // Å
    double max_dev = 0.0;
    for (int x = 0; x < 6; ++x) {
        const auto fp1 = forces_at(displaced(frame, x, h));
        const auto fm1 = forces_at(displaced(frame, x, -h));
        const auto fp2 = forces_at(displaced(frame, x, 2.0 * h));
        const auto fm2 = forces_at(displaced(frame, x, -2.0 * h));
        for (int y = 0; y < 6; ++y) {
            const double dF =
                (8.0 * (fp1(y) - fm1(y)) - (fp2(y) - fm2(y))) / (12.0 * h);
            const double dev = std::abs(res.hessian(x, y) + dF);
            max_dev = std::max(max_dev, dev);
            EXPECT_NEAR(res.hessian(x, y), -dF, 5e-5)
                << "element (" << x << ", " << y << ")";
        }
    }
    EXPECT_LT(max_dev, 5e-5);
}

TEST(UhfAnalyticHessian, OhSymmetricAndTranslationClean) {
    auto frame = make_oh();
    cpu::scf::UHF uhf(5, 4);
    uhf.basis = BasisSet::STO3G;
    uhf.tol = 1e-13;
    const cpu::vib::HessianResult res = uhf.hessian(frame);

    for (int x = 0; x < 6; ++x)
        for (int y = 0; y < 6; ++y)
            EXPECT_NEAR(res.hessian(x, y), res.hessian(y, x), 1e-9);

    for (int x = 0; x < 6; ++x)
        for (int d = 0; d < 3; ++d) {
            double sum = 0.0;
            for (int B = 0; B < 2; ++B)
                sum += res.hessian(x, 3 * B + d);
            EXPECT_NEAR(sum, 0.0, 1e-5) << "row " << x << " axis " << d;
        }
}

#endif // ATV_HAS_HF
