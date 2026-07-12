/// @file tests/cpu/potentials/test_eam.cpp
/// @brief Unit: analytic EAM (Sutton-Chen) — lookup, density, and the FD force
///        oracle.
///
/// The FD oracle is the hard gate. Energy-only tests cannot see a force that is
/// the gradient of a DIFFERENT functional than the one the energy evaluates —
/// that exact bug hid in the GPW path for 17 commits. Every force this engine
/// grows gets differenced against its own energy.

#include <cmath>
#include <vector>

#include <gtest/gtest.h>

#include "atomiverse/core/mat3.h"
#include "atomiverse/cpu/driver/evaluate.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/neighbor/neighbor_list.h"
#include "atomiverse/cpu/potentials/eam.h"
#include "atomiverse/cpu/potentials/eam_potential.h"

using namespace ATV;

namespace {

/// Distorted fcc copper — several call sites below.
cpu::Frame copper(double a0, int n_cell, double jitter) {
    const int n = 4 * n_cell * n_cell * n_cell;
    cpu::Frame f(n);
    const double b[4][3] = {{0, 0, 0}, {0.5, 0.5, 0}, {0.5, 0, 0.5}, {0, 0.5, 0.5}};
    int idx = 0;
    for (int ix = 0; ix < n_cell; ++ix)
        for (int iy = 0; iy < n_cell; ++iy)
            for (int iz = 0; iz < n_cell; ++iz)
                for (int k = 0; k < 4; ++k) {
                    const double j = jitter * std::sin(double(13 * idx + 3));
                    f.atoms.x[idx] = (double(ix) + b[k][0]) * a0 + j;
                    f.atoms.y[idx] = (double(iy) + b[k][1]) * a0 - 0.7 * j;
                    f.atoms.z[idx] = (double(iz) + b[k][2]) * a0 + 0.3 * j;
                    f.atoms.mass[idx] = 63.546;
                    f.atoms.type[idx] = 29;
                    ++idx;
                }
    f.box = Mat3::cubic(a0 * n_cell);
    return f;
}

/// Single-point energy through the public Evaluate driver — called once per FD
/// displacement, i.e. many times.
double energy_of(cpu::Frame& f, double cutoff) {
    cpu::Evaluate eval;
    eval.set_neighbor(cpu::md::NeighborList<>{cutoff, 0.3, 1, 0, true});
    eval.on(evaluate::Phase::evaluate)
        .add(cpu::EAM{cpu::EAMPotential{Element::Cu, "sutton-chen"}, cutoff});
    eval.build();
    eval.run(f);
    return f.scalars.energy;
}

} // namespace

/// The lookup constructor throws rather than inventing a fallback.
TEST(EAM, LookupThrowsOnUnsupported) {
    EXPECT_NO_THROW(cpu::EAMPotential(Element::Cu, "sutton-chen"));
    EXPECT_THROW(cpu::EAMPotential(Element::Cu, "no-such-parameterization"),
                 std::out_of_range);
    EXPECT_THROW(cpu::EAMPotential(Element::C, "sutton-chen"), std::out_of_range);
}

/// The embedding function is nonlinear — the property that makes EAM many-body.
/// If F were linear in rho, EAM would collapse to a pair potential.
TEST(EAM, EmbeddingIsNonlinear) {
    const cpu::EAMPotential p(Element::Cu, "sutton-chen");
    // F(2 rho) != 2 F(rho) for a sqrt embedding.
    const double rho = 3.0;
    EXPECT_GT(std::fabs(p.F(2.0 * rho) - 2.0 * p.F(rho)), 1e-6);
    // F is negative (binding) and dF is negative and finite away from zero.
    EXPECT_LT(p.F(rho), 0.0);
    EXPECT_LT(p.dF(rho), 0.0);
    // rho = 0 must not produce a NaN even though sqrt has infinite slope there.
    EXPECT_TRUE(std::isfinite(p.dF(0.0)));
}

/// Analytic derivatives of phi and f match finite differences.
TEST(EAM, PairAndDensityDerivativesMatchFD) {
    const cpu::EAMPotential p(Element::Cu, "sutton-chen");
    const double h = 1e-6;
    for (double r : {2.2, 2.556, 3.0, 4.5, 6.0}) {
        const double dphi_fd = (p.phi(r + h) - p.phi(r - h)) / (2 * h);
        const double df_fd = (p.f(r + h) - p.f(r - h)) / (2 * h);
        EXPECT_NEAR(p.dphi(r), dphi_fd, 1e-6 * std::fabs(dphi_fd) + 1e-9);
        EXPECT_NEAR(p.df(r), df_fd, 1e-6 * std::fabs(df_fd) + 1e-9);
    }
}

/// THE GATE: the analytic many-body force is the gradient of the energy the
/// component actually reports.
TEST(EAM, AnalyticForceMatchesFiniteDifference) {
    const double cutoff = 7.0;
    cpu::Frame f = copper(3.61, 2, /*jitter=*/0.12); // 32 atoms, distorted
    const int n = f.atoms.size();

    // Analytic forces.
    cpu::Evaluate eval;
    eval.set_neighbor(cpu::md::NeighborList<>{cutoff, 0.3, 1, 0, true});
    eval.on(evaluate::Phase::evaluate)
        .add(cpu::EAM{cpu::EAMPotential{Element::Cu, "sutton-chen"}, cutoff});
    eval.build();
    eval.run(f);

    std::vector<double> fx(n), fy(n), fz(n);
    for (int i = 0; i < n; ++i) {
        fx[i] = f.atoms.fx[i];
        fy[i] = f.atoms.fy[i];
        fz[i] = f.atoms.fz[i];
    }

    // Central-difference the energy for a handful of atoms (all 3N would be slow
    // and buys nothing — a sign or factor error shows up on the first atom).
    const double h = 1e-5;
    double worst = 0.0;
    for (int i : {0, 5, 17, 31}) {
        for (int d = 0; d < 3; ++d) {
            auto& coord = (d == 0) ? f.atoms.x : (d == 1) ? f.atoms.y : f.atoms.z;
            const double x0 = coord[i];

            coord[i] = x0 + h;
            cpu::Frame plus = f;
            const double e_plus = energy_of(plus, cutoff);

            coord[i] = x0 - h;
            cpu::Frame minus = f;
            const double e_minus = energy_of(minus, cutoff);

            coord[i] = x0;

            const double f_fd = -(e_plus - e_minus) / (2 * h);
            const double f_an = (d == 0) ? fx[i] : (d == 1) ? fy[i] : fz[i];
            worst = std::max(worst, std::fabs(f_an - f_fd));
            EXPECT_NEAR(f_an, f_fd, 1e-6)
                << "atom " << i << " dim " << d << ": analytic " << f_an
                << " vs FD " << f_fd;
        }
    }
    EXPECT_LT(worst, 1e-6);
}

/// Newton's third law: the many-body forces still sum to zero.
TEST(EAM, ForcesSumToZero) {
    const double cutoff = 7.0;
    cpu::Frame f = copper(3.61, 2, 0.12);
    cpu::Evaluate eval;
    eval.set_neighbor(cpu::md::NeighborList<>{cutoff, 0.3, 1, 0, true});
    eval.on(evaluate::Phase::evaluate)
        .add(cpu::EAM{cpu::EAMPotential{Element::Cu, "sutton-chen"}, cutoff});
    eval.build();
    eval.run(f);

    double sx = 0, sy = 0, sz = 0;
    for (int i = 0; i < f.atoms.size(); ++i) {
        sx += f.atoms.fx[i];
        sy += f.atoms.fy[i];
        sz += f.atoms.fz[i];
    }
    EXPECT_NEAR(sx, 0.0, 1e-9);
    EXPECT_NEAR(sy, 0.0, 1e-9);
    EXPECT_NEAR(sz, 0.0, 1e-9);
}
