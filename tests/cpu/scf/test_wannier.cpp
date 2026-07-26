/// @file test_wannier.cpp
/// @brief Unit tests for Foster–Boys pair angle math and WannierCenters
///        on synthetic / minimal determinants. No full multi-orbital SCF
///        e2e (that lives in regressions/scf/wannier-centers).

#include <array>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>
#include <xtensor/containers/xtensor.hpp>

#include "atomiverse/core/units.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/scf/slater_determinant.h"
#include "atomiverse/cpu/scf/wannier.h"

using namespace ATV;
using cpu::scf::boys_pair_angle;
using cpu::scf::MeanFieldOutputs;
using cpu::scf::SlaterDeterminant;
using cpu::scf::WannierCenters;
using cpu::scf::native::Shell;

namespace {

constexpr double kPi = 3.14159265358979323846;

/// One normalized s primitive at `center_bohr`.
struct OneS {
    std::vector<double> storage;
    std::shared_ptr<std::vector<Shell>> shells;
    std::shared_ptr<SlaterDeterminant> det;
    std::shared_ptr<bool> ready;

    explicit OneS(const double center_bohr[3]) {
        const double alpha = 1.0;
        const double norm = std::pow(2.0 * alpha / kPi, 0.75);
        storage = {alpha, norm};
        shells = std::make_shared<std::vector<Shell>>();
        shells->push_back(
            Shell{1,
                  &storage[0],
                  &storage[1],
                  {center_bohr[0], center_bohr[1], center_bohr[2]},
                  0,
                  0});
        cpu::linalg::RealMatrix C =
            xt::ones<double>({std::size_t(1), std::size_t(1)});
        cpu::linalg::RealVector eps = xt::zeros<double>({std::size_t(1)});
        cpu::linalg::RealVector occ = xt::ones<double>({std::size_t(1)});
        det = std::make_shared<SlaterDeterminant>(C, eps, occ, 1);
        ready = std::make_shared<bool>(true);
    }

    MeanFieldOutputs handle() const {
        return {det, shells, ready,
                std::make_shared<const xt::xtensor<double, 2>>()};
    }
};

/// Two non-orthogonal s AOs at ±d on z; two occupied MOs = identity in AO
/// basis (not orthonormal — only for exercise of multi-orbital path when
/// multipole is well-defined). Uses a restricted det with n_e=2, f=(1,1).
struct TwoS {
    std::vector<double> storage;
    std::shared_ptr<std::vector<Shell>> shells;
    std::shared_ptr<SlaterDeterminant> det;
    std::shared_ptr<bool> ready;

    TwoS() {
        const double alpha = 1.0;
        const double norm = std::pow(2.0 * alpha / kPi, 0.75);
        storage = {alpha, norm};
        shells = std::make_shared<std::vector<Shell>>();
        const double d = 1.0; // Bohr
        shells->push_back(
            Shell{1, &storage[0], &storage[1], {0.0, 0.0, -d}, 0, 0});
        shells->push_back(
            Shell{1, &storage[0], &storage[1], {0.0, 0.0, +d}, 0, 0});
        // Two orthonormal-looking columns (AO basis identity): not
        // S-orthogonal, but WannierCenters does not re-orthogonalize —
        // moments remain well-defined for the path.
        cpu::linalg::RealMatrix C =
            xt::eye<double>({std::size_t(2), std::size_t(2)});
        cpu::linalg::RealVector eps = xt::zeros<double>({std::size_t(2)});
        cpu::linalg::RealVector occ = xt::ones<double>({std::size_t(2)});
        det = std::make_shared<SlaterDeterminant>(C, eps, occ, 2);
        ready = std::make_shared<bool>(true);
    }

    MeanFieldOutputs handle() const {
        return {det, shells, ready,
                std::make_shared<const xt::xtensor<double, 2>>()};
    }
};

cpu::Frame dummy_frame(int n) {
    cpu::Frame f(n);
    for (int i = 0; i < n; ++i) {
        f.span(kAtom_type)[i] = 1;
        f.span(kAtom_x)[i] = 0.0;
        f.span(kAtom_y)[i] = 0.0;
        f.span(kAtom_z)[i] = 0.0;
    }
    f.box = cpu::Box::cubic(20.0);
    return f;
}

} // namespace

TEST(BoysPairAngle, Quadrants) {
    EXPECT_NEAR(boys_pair_angle(1.0, 0.0), 0.0, 1e-14);
    EXPECT_NEAR(boys_pair_angle(0.0, 1.0), kPi / 8.0, 1e-14);
    EXPECT_NEAR(boys_pair_angle(0.0, -1.0), -kPi / 8.0, 1e-14);
    // A < 0, B = 0 → atan2(0, -1) = π → γ = π/4
    EXPECT_NEAR(boys_pair_angle(-1.0, 0.0), kPi / 4.0, 1e-14);
}

TEST(WannierCenters, OneOrbitalIdentity) {
    const double c[3] = {0.5, -0.25, 1.0}; // Bohr
    OneS sys{c};
    WannierCenters wb;
    auto frame = dummy_frame(1);
    auto r = wb.compute(sys.handle(), frame);
    ASSERT_EQ(r.centers.size(), 1u);
    ASSERT_EQ(r.spreads.size(), 1u);
    EXPECT_TRUE(r.converged);
    EXPECT_EQ(r.sweeps, 0);
    EXPECT_NEAR(r.total_spread, r.total_spread_initial, 1e-12);
    // Centroid of a single s AO is its center.
    EXPECT_NEAR(r.centers[0][0], c[0] * units::kBohrInAngstrom, 1e-8);
    EXPECT_NEAR(r.centers[0][1], c[1] * units::kBohrInAngstrom, 1e-8);
    EXPECT_NEAR(r.centers[0][2], c[2] * units::kBohrInAngstrom, 1e-8);
    EXPECT_GT(r.spreads[0], 0.0);
    EXPECT_TRUE(std::isfinite(r.spreads[0]));
}

TEST(WannierCenters, TwoOrbitalSpreadDoesNotIncrease) {
    TwoS sys;
    WannierCenters wb;
    wb.tol = 1e-12;
    wb.max_sweeps = 50;
    auto frame = dummy_frame(2);
    auto r = wb.compute(sys.handle(), frame);
    ASSERT_EQ(r.centers.size(), 2u);
    EXPECT_LE(r.total_spread, r.total_spread_initial + 1e-10);
    for (double s : r.spreads) {
        EXPECT_GT(s, 0.0);
        EXPECT_TRUE(std::isfinite(s));
    }
}

TEST(WannierCenters, NotReadyThrows) {
    const double c[3] = {0.0, 0.0, 0.0};
    OneS sys{c};
    *sys.ready = false;
    WannierCenters wb;
    auto frame = dummy_frame(1);
    EXPECT_THROW(wb.compute(sys.handle(), frame), std::invalid_argument);
}

TEST(WannierCenters, BadParamsThrow) {
    const double c[3] = {0.0, 0.0, 0.0};
    OneS sys{c};
    WannierCenters wb;
    wb.max_sweeps = 0;
    auto frame = dummy_frame(1);
    EXPECT_THROW(wb.compute(sys.handle(), frame), std::invalid_argument);
    wb.max_sweeps = 10;
    wb.tol = 0.0;
    EXPECT_THROW(wb.compute(sys.handle(), frame), std::invalid_argument);
}

TEST(WannierCenters, OpenShellThrows) {
    // Unrestricted det with α only occupied in a trivial way.
    cpu::linalg::RealMatrix Ca =
        xt::ones<double>({std::size_t(1), std::size_t(1)});
    cpu::linalg::RealMatrix Cb =
        xt::ones<double>({std::size_t(1), std::size_t(1)});
    cpu::linalg::RealVector eps = xt::zeros<double>({std::size_t(1)});
    cpu::linalg::RealVector fa = xt::ones<double>({std::size_t(1)});
    cpu::linalg::RealVector fb = xt::zeros<double>({std::size_t(1)});
    auto det =
        std::make_shared<SlaterDeterminant>(Ca, eps, fa, 1, Cb, eps, fb, 0);
    const double alpha = 1.0;
    const double norm = std::pow(2.0 * alpha / kPi, 0.75);
    auto storage =
        std::make_shared<std::vector<double>>(std::vector<double>{alpha, norm});
    auto shells = std::make_shared<std::vector<Shell>>();
    shells->push_back(
        Shell{1, &(*storage)[0], &(*storage)[1], {0.0, 0.0, 0.0}, 0, 0});
    // Keep storage alive via shells? Shell holds raw pointers — hold storage
    // in a static... Use a longer-lived vector on the stack carefully:
    // actually storage is shared_ptr so it lives. But Shell points into it.
    MeanFieldOutputs outs{det, shells, std::make_shared<bool>(true),
                          std::make_shared<const xt::xtensor<double, 2>>()};
    // Keep storage alive for Shell pointers.
    auto keep = storage;
    WannierCenters wb;
    auto frame = dummy_frame(1);
    EXPECT_THROW(wb.compute(outs, frame), std::invalid_argument);
}
