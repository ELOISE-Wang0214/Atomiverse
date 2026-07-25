/// @file tests/cpu/gpw/test_cpu_gpw_kpoints.cpp
/// @brief Unit tests for the GPW k-point machinery (internal/kpoints.h):
///        Monkhorst-Pack mesh generation, phased Bloch lattice sums (Γ limit,
///        Hermiticity, k ↔ −k conjugation), the zheevd Hermitian eigensolver,
///        the per-k pre-summed basis-on-grid Γ limit, and the grid
///        commensuration the supercell-folding regression relies on.
///
/// Unit-scale: no SCF is run — each test exercises one function on a toy
/// basis (single H-like s Gaussians in a small cubic cell).

#include <gtest/gtest.h>

#include <cmath>
#include <complex>
#include <vector>

#include "atomiverse/cpu/gpw/internal/collocate.h"
#include "atomiverse/cpu/gpw/internal/kpoints.h"
#include "atomiverse/cpu/gpw/internal/one_electron.h"
#include "atomiverse/cpu/gpw/internal/periodic.h"

using namespace ATV::cpu::gpw;
using internal::CplxMatrix;
using internal::KPoint;
using Cplx = std::complex<double>;

namespace {

// Two-function toy basis: normalized-ish s Gaussians at two centers.
std::vector<GridBasisFunction> toy_basis(double L) {
    GridPrimitive p1;
    p1.alpha = 0.35;
    p1.center[0] = 0.2 * L;
    p1.center[1] = 0.3 * L;
    p1.center[2] = 0.1 * L;
    p1.coeff = 1.0;
    GridPrimitive p2 = p1;
    p2.alpha = 0.22;
    p2.center[0] = 0.7 * L;
    p2.center[1] = 0.6 * L;
    p2.center[2] = 0.5 * L;
    p2.coeff = 0.8;
    return {GridBasisFunction{p1}, GridBasisFunction{p2}};
}

} // namespace

TEST(GpwKpoints, MonkhorstPackMeshCp2kConvention) {
    const int mesh[3] = {2, 2, 2};
    const auto kpts = internal::monkhorst_pack(mesh, /*time_reversal=*/false);
    ASSERT_EQ(kpts.size(), 8u);
    double wsum = 0.0;
    for (const auto& kp : kpts) {
        wsum += kp.weight;
        for (int d = 0; d < 3; ++d) // CP2K MP N=2: kappa = ±1/4
            EXPECT_NEAR(std::abs(kp.frac[d]), 0.25, 1e-14);
    }
    EXPECT_NEAR(wsum, 1.0, 1e-14);
}

TEST(GpwKpoints, MonkhorstPackTimeReversalFold) {
    const int mesh[3] = {2, 2, 2};
    const auto kpts = internal::monkhorst_pack(mesh, /*time_reversal=*/true);
    ASSERT_EQ(kpts.size(), 4u); // 8 points fold into 4 ±k pairs
    double wsum = 0.0;
    for (const auto& kp : kpts) {
        EXPECT_NEAR(kp.weight, 0.25, 1e-14); // each pair carries 2/8
        wsum += kp.weight;
    }
    EXPECT_NEAR(wsum, 1.0, 1e-14);
}

TEST(GpwKpoints, MonkhorstPackGammaMesh) {
    const int mesh[3] = {1, 1, 1};
    const auto kpts = internal::monkhorst_pack(mesh, true);
    ASSERT_EQ(kpts.size(), 1u);
    for (int d = 0; d < 3; ++d)
        EXPECT_NEAR(kpts[0].frac[d], 0.0, 1e-15);
    EXPECT_NEAR(kpts[0].weight, 1.0, 1e-15);
}

TEST(GpwKpoints, OverlapKineticGammaLimitMatchesRealLatticeSum) {
    const double L = 6.0;
    const double cell[3] = {L, L, L};
    const auto basis = toy_basis(L);
    std::vector<double> radii;
    for (const auto& bf : basis)
        radii.push_back(basis_function_radius(bf, 1e-8));
    const double reach = radii[0] + radii[1];
    int m[3];
    pair_image_counts(cell, true, 2.0 * std::max(radii[0], radii[1]), m);

    std::vector<KPoint> kpts(1); // Γ
    std::vector<CplxMatrix> S_k, T_k;
    internal::build_overlap_kinetic_k(basis, radii, cell, m, 1e-12, kpts, S_k,
                                      T_k);
    // Real Γ lattice sums from one_electron.h — must agree to round-off.
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j) {
            const double s_ref =
                overlap(basis[i], basis[j], cell, m,
                        radii[static_cast<std::size_t>(i)]
                            + radii[static_cast<std::size_t>(j)]);
            const double t_ref =
                kinetic(basis[i], basis[j], cell, m,
                        radii[static_cast<std::size_t>(i)]
                            + radii[static_cast<std::size_t>(j)]);
            EXPECT_NEAR(S_k[0](i, j).real(), s_ref, 1e-12);
            EXPECT_NEAR(S_k[0](i, j).imag(), 0.0, 1e-14);
            EXPECT_NEAR(T_k[0](i, j).real(), t_ref, 1e-12);
            EXPECT_NEAR(T_k[0](i, j).imag(), 0.0, 1e-14);
        }
    (void)reach;
}

TEST(GpwKpoints, OverlapHermitianAndConjugateAtMinusK) {
    const double L = 6.0;
    const double cell[3] = {L, L, L};
    const auto basis = toy_basis(L);
    std::vector<double> radii;
    for (const auto& bf : basis)
        radii.push_back(basis_function_radius(bf, 1e-8));
    int m[3];
    pair_image_counts(cell, true, 2.0 * std::max(radii[0], radii[1]), m);

    std::vector<KPoint> kpts(2);
    kpts[0].frac[0] = 0.25;
    kpts[0].frac[1] = -0.10;
    kpts[0].frac[2] = 0.40;
    kpts[1].frac[0] = -0.25; // −k
    kpts[1].frac[1] = 0.10;
    kpts[1].frac[2] = -0.40;
    std::vector<CplxMatrix> S_k, T_k;
    internal::build_overlap_kinetic_k(basis, radii, cell, m, 1e-12, kpts, S_k,
                                      T_k);
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j) {
            // Hermiticity at each k.
            EXPECT_NEAR(S_k[0](i, j).real(), S_k[0](j, i).real(), 1e-12);
            EXPECT_NEAR(S_k[0](i, j).imag(), -S_k[0](j, i).imag(), 1e-12);
            // S_{−k} = conj(S_k).
            EXPECT_NEAR(S_k[1](i, j).real(), S_k[0](i, j).real(), 1e-12);
            EXPECT_NEAR(S_k[1](i, j).imag(), -S_k[0](i, j).imag(), 1e-12);
        }
    // A genuine off-Γ k must produce a nonzero imaginary part somewhere
    // (two inequivalent centers → complex Bloch phase survives).
    double imax = 0.0;
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
            imax = std::max(imax, std::abs(S_k[0](i, j).imag()));
    EXPECT_GT(imax, 1e-10);
}

TEST(GpwKpoints, EighHermitianAnalytic2x2) {
    // A = [[2, i],[−i, 2]] has eigenvalues 1 and 3.
    CplxMatrix A = xt::zeros<Cplx>({2, 2});
    A(0, 0) = Cplx(2, 0);
    A(0, 1) = Cplx(0, 1);
    A(1, 0) = Cplx(0, -1);
    A(1, 1) = Cplx(2, 0);
    std::vector<double> w;
    CplxMatrix V;
    internal::eigh_hermitian(A, w, V);
    ASSERT_EQ(w.size(), 2u);
    EXPECT_NEAR(w[0], 1.0, 1e-12);
    EXPECT_NEAR(w[1], 3.0, 1e-12);
    // Columns are eigenvectors: A v = w v.
    for (int col = 0; col < 2; ++col)
        for (int row = 0; row < 2; ++row) {
            Cplx av(0, 0);
            for (int t = 0; t < 2; ++t)
                av += A(row, t) * V(t, col);
            const Cplx wv = w[static_cast<std::size_t>(col)] * V(row, col);
            EXPECT_NEAR(std::abs(av - wv), 0.0, 1e-12);
        }
}

TEST(GpwKpoints, KpointBasisGridGammaMatchesCompactBasisGrid) {
    const double L = 6.0;
    const double cell[3] = {L, L, L};
    const auto basis = toy_basis(L);
    const UniformGrid grid = make_uniform_grid(cell, 60.0);

    std::vector<KPoint> kpts(1); // Γ
    const auto kbg =
        internal::build_kpoint_basis_grid(grid, basis, kpts, 1e-10);

    // Same symmetric P through both paths; the Γ CompactBasisGrid double
    // loop over duplicated image triples must equal the pre-summed complex
    // path at k = 0.
    const int nb = 2;
    std::vector<double> P = {1.2, 0.3, 0.3, 0.7};
    const CompactBasisGrid cbg =
        build_compact_basis_grid(grid, basis, /*periodic=*/true, 1e-10);
    const std::vector<double> rho_ref = collocate_density(cbg, nb, P);

    CplxMatrix Pk = xt::zeros<Cplx>({nb, nb});
    Pk(0, 0) = 1.2;
    Pk(0, 1) = 0.3;
    Pk(1, 0) = 0.3;
    Pk(1, 1) = 0.7;
    std::vector<double> rho(static_cast<std::size_t>(grid.npts()), 0.0);
    internal::collocate_density_k(kbg, 0, Pk, 1.0, rho);

    double q_ref = 0.0, q = 0.0, max_diff = 0.0;
    for (long g = 0; g < grid.npts(); ++g) {
        q_ref += rho_ref[static_cast<std::size_t>(g)];
        q += rho[static_cast<std::size_t>(g)];
        max_diff = std::max(max_diff,
                            std::abs(rho[static_cast<std::size_t>(g)]
                                     - rho_ref[static_cast<std::size_t>(g)]));
    }
    // Same eps drop threshold semantics differ per-triple vs pre-summed, so
    // compare at physical accuracy, not bit identity.
    EXPECT_NEAR(max_diff, 0.0, 1e-8);
    EXPECT_NEAR(q * grid.weight(), q_ref * grid.weight(), 1e-8);
}

TEST(GpwKpoints, SupercellGridsCommensurateAtRegressionCutoffs) {
    // The supercell-folding regression (acceptance A2) doubles the bct Si
    // cell along z and demands IDENTICAL grid spacing: n_z(2c) == 2 n_z(c)
    // and unchanged n_x/n_y. Pin the cutoff (95 Ry) that guarantees it for
    // a = 5.43 Å (10.26120 Bohr), a_t = a/√2.
    const double a = 10.26120;          // Bohr
    const double at = a / std::sqrt(2.0); // 7.25577
    const double cell1[3] = {at, at, a};
    const double cell2[3] = {at, at, 2.0 * a};
    const UniformGrid g1 = make_uniform_grid(cell1, 95.0);
    const UniformGrid g2 = make_uniform_grid(cell2, 95.0);
    EXPECT_EQ(g2.n[0], g1.n[0]);
    EXPECT_EQ(g2.n[1], g1.n[1]);
    EXPECT_EQ(g2.n[2], 2 * g1.n[2]);
}
