/// @file tests/cpu/kspace/test_pme_coul.cpp
/// @brief Validation of CPU Particle-Mesh Ewald (CoulLong + KSpacePME).
///
/// Same checks as the direct-Ewald test, plus a direct PME-vs-Ewald total
/// energy comparison. Built only with ATV_BUILD_KSPACE (FFTW3).

#include <gtest/gtest.h>

#include <cmath>
#include <vector>

#include "atomiverse/core/mat3.h"
#include "atomiverse/core/units.h"
#include "atomiverse/cpu/component/component.h"
#include "atomiverse/cpu/component/context.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/kspace/ewald.h"
#include "atomiverse/cpu/kspace/pme.h"
#include "atomiverse/cpu/neighbor/pair_list.h"
#include "atomiverse/cpu/potentials/coul_long.h"

using namespace ATV;

namespace {

constexpr double kMadelungNaCl = 1.747564594633182;

cpu::Frame make_rocksalt(double d, int n) {
    const int N = n * n * n;
    cpu::Frame frame(N);
    int idx = 0;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            for (int k = 0; k < n; ++k) {
                frame.atoms.x[idx] = i * d;
                frame.atoms.y[idx] = j * d;
                frame.atoms.z[idx] = k * d;
                frame.atoms.charge[idx] = ((i + j + k) % 2 == 0) ? 1.0 : -1.0;
                frame.atoms.mass[idx] = 1.0;
                frame.atoms.type[idx] = 0;
                frame.atoms.tag[idx] = idx + 1;
                ++idx;
            }
    frame.box = Mat3::cubic(n * d);
    return frame;
}

void build_full_halflist(const cpu::Frame& f, cpu::md::PairList& pairs) {
    const int N = f.atoms.size();
    pairs.inum = N;
    pairs.ilist.resize(N);
    pairs.numneigh.resize(N);
    pairs.firstneigh.resize(N);
    pairs.neighbors.clear();
    int off = 0;
    for (int i = 0; i < N; ++i) {
        pairs.ilist[i] = i;
        pairs.firstneigh[i] = off;
        for (int j = i + 1; j < N; ++j)
            pairs.neighbors.push_back(j);
        pairs.numneigh[i] = N - 1 - i;
        off += N - 1 - i;
    }
    pairs.total_pairs = off;
    pairs.has_shifts = false;
    pairs.version = 1;
}

// Combined real-space + a reciprocal component (templated on the kspace type).
template <typename KSpace>
double electrostatic_energy(cpu::Frame& f, const cpu::md::PairList& pairs,
                            double alpha, double cutoff, KSpace&& recip) {
    for (int i = 0; i < f.atoms.size(); ++i)
        f.atoms.fx[i] = f.atoms.fy[i] = f.atoms.fz[i] = 0;
    f.scalars.energy = 0;
    f.scalars.energy_pair = 0;

    cpu::Context ctx;
    ctx.pair_list = &pairs;
    cpu::BuildContext bctx{ctx, f};
    cpu::LaunchContext lctx{ctx, f, 0, 0};

    cpu::CoulLong real(alpha, cutoff);
    real.build(bctx);
    recip.build(bctx);
    real.launch(lctx);
    recip.launch(lctx);
    real.destroy({ctx, f});
    recip.destroy({ctx, f});
    return f.scalars.energy;
}

} // namespace

TEST(PmeCoul, MadelungNaCl) {
    const double d = 2.82;
    const int n = 4;
    auto frame = make_rocksalt(d, n);
    cpu::md::PairList pairs;
    build_full_halflist(frame, pairs);

    const double cutoff = 1.9 * d;
    const double alpha = 3.5 / cutoff; // moderate α so the grid resolves k-space
    const double e = electrostatic_energy(frame, pairs, alpha, cutoff,
                                          cpu::KSpacePME(alpha, 0.5, 6));

    const int N = n * n * n;
    const double e_expected = -0.5 * N * kMadelungNaCl * ATV::kCoulomb / d;
    EXPECT_NEAR(e / e_expected, 1.0, 5e-4)
        << "E=" << e << " expected=" << e_expected;
}

TEST(PmeCoul, AgreesWithEwald) {
    const double d = 2.82;
    const int n = 4;
    auto frame = make_rocksalt(d, n);
    cpu::md::PairList pairs;
    build_full_halflist(frame, pairs);

    const double cutoff = 1.9 * d;
    const double alpha = 3.5 / cutoff;
    const double e_ewald = electrostatic_energy(
        frame, pairs, alpha, cutoff, cpu::KSpaceEwald(alpha, 1e-8));
    const double e_pme = electrostatic_energy(frame, pairs, alpha, cutoff,
                                              cpu::KSpacePME(alpha, 0.5, 6));
    EXPECT_NEAR(e_pme / e_ewald, 1.0, 5e-4)
        << "ewald=" << e_ewald << " pme=" << e_pme;
}

TEST(PmeCoul, FiniteDifferenceForces) {
    cpu::Frame frame(4);
    const double coords[4][3] = {
        {1.0, 1.0, 1.0}, {3.2, 1.4, 1.1}, {1.3, 3.1, 2.4}, {2.7, 2.2, 3.6}};
    const double q[4] = {0.6, -0.4, 0.5, -0.7};
    for (int i = 0; i < 4; ++i) {
        frame.atoms.x[i] = coords[i][0];
        frame.atoms.y[i] = coords[i][1];
        frame.atoms.z[i] = coords[i][2];
        frame.atoms.charge[i] = q[i];
        frame.atoms.mass[i] = 1.0;
        frame.atoms.type[i] = 0;
        frame.atoms.tag[i] = i + 1;
    }
    frame.box = Mat3::cubic(8.0);
    cpu::md::PairList pairs;
    build_full_halflist(frame, pairs);

    const double cutoff = 3.9;
    const double alpha = 0.45;
    auto run = [&]() {
        return electrostatic_energy(frame, pairs, alpha, cutoff,
                                    cpu::KSpacePME(alpha, 0.7, 6));
    };
    auto force = [&](int a, int axis) {
        run();
        return axis == 0 ? frame.atoms.fx[a]
               : axis == 1 ? frame.atoms.fy[a]
                           : frame.atoms.fz[a];
    };
    auto energy_at = [&](int a, int axis, double delta) {
        auto& pos = axis == 0 ? frame.atoms.x
                    : axis == 1 ? frame.atoms.y
                                : frame.atoms.z;
        const double orig = pos[a];
        pos[a] = orig + delta;
        const double e = run();
        pos[a] = orig;
        return e;
    };

    const double h = 1e-5;
    for (int a = 0; a < 4; ++a)
        for (int axis = 0; axis < 3; ++axis) {
            const double f_analytic = force(a, axis);
            const double grad =
                (energy_at(a, axis, h) - energy_at(a, axis, -h)) / (2 * h);
            EXPECT_NEAR(f_analytic, -grad, 2e-3)
                << "atom " << a << " axis " << axis;
        }
}
