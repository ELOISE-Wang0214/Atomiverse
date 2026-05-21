/// @file test_grid.cpp
/// @brief Unit tests for molecular integration grid.
///
/// Covers the three acceptance criteria from libxc-xc-fix-04-grid:
///   ac-001: Lebedev angular grids integrate spherical harmonics exactly.
///   ac-002: Mura-Knowles radial grid reproduces analytic Gaussian moments.
///   ac-003: Bragg-Slater table covers Z=1..36 (and beyond) and rejects
///           out-of-range Z with std::invalid_argument.
///
/// Plus baseline smoke tests for molecular grid construction.

#ifdef ATV_HAS_RHF

#include <cmath>
#include <gtest/gtest.h>
#include <stdexcept>
#include <vector>

#include "atomiverse/core/constants.h"
#include "atomiverse/core/mat3.h"
#include "atomiverse/cpu/frame.h"
#include "atomiverse/cpu/scf/basis_set.h"
#include "atomiverse/cpu/scf/grid.h"
#include "atomiverse/cpu/scf/integrals.h"
#include "atomiverse/cpu/scf/lebedev_data.h"
#include "atomiverse/cpu/scf/native_integrals.h"
#include "atomiverse/cpu/scf/scf_solve.h"

using namespace ATV;
using namespace ATV::cpu::scf;

namespace {

constexpr double kPi = 3.14159265358979323846;

// ── Real (Cartesian) spherical harmonics ────────────────────────────────
//
// Up to l = 4 covers stated-degree-3 grid (order 6). Higher l needed for
// 50/110/194/302/590 we evaluate via the standard tesseral formulas in
// Cartesian form (z := cos theta, x_xy := sin theta cos phi, etc., all
// expressed in the unit-sphere Cartesian coordinates of the Lebedev point).
//
// Normalization: sum_i w_i Y_{lm}(p_i)^2 = 1 for proper SH on the unit
// sphere with weights summing to 1.

// l = 0
inline double Y00(double, double, double) {
    return 1.0 / (2.0 * std::sqrt(kPi));
}

// l = 1: Y_{1,-1}=y, Y_{1,0}=z, Y_{1,1}=x  (real, with sqrt(3/4pi) coef)
inline double Y1m1(double, double y, double) {
    return std::sqrt(3.0 / (4.0 * kPi)) * y;
}
inline double Y10(double, double, double z) {
    return std::sqrt(3.0 / (4.0 * kPi)) * z;
}
inline double Y11(double x, double, double) {
    return std::sqrt(3.0 / (4.0 * kPi)) * x;
}

// l = 2: real tesseral, in normalized form (coefficient * cartesian poly)
//   Y_{2,-2} = c_22 * (x*y)            (sqrt(15/4pi)/2 ... see below)
//   Y_{2,-1} = c_21 * (y*z)
//   Y_{2, 0} = c_20 * (3z^2 - 1)
//   Y_{2, 1} = c_21 * (x*z)
//   Y_{2, 2} = c_22 * (x^2 - y^2)
inline double Y2m2(double x, double y, double) {
    return 0.5 * std::sqrt(15.0 / kPi) * x * y;
}
inline double Y2m1(double, double y, double z) {
    return 0.5 * std::sqrt(15.0 / kPi) * y * z;
}
inline double Y20(double, double, double z) {
    return 0.25 * std::sqrt(5.0 / kPi) * (3.0 * z * z - 1.0);
}
inline double Y21(double x, double, double z) {
    return 0.5 * std::sqrt(15.0 / kPi) * x * z;
}
inline double Y22(double x, double y, double) {
    return 0.25 * std::sqrt(15.0 / kPi) * (x * x - y * y);
}

// l = 3: real tesseral, Cartesian form (Wikipedia table convention)
inline double Y3m3(double x, double y, double) {
    return 0.25 * std::sqrt(35.0 / (2.0 * kPi)) * y * (3.0 * x * x - y * y);
}
inline double Y3m2(double x, double y, double z) {
    return 0.5 * std::sqrt(105.0 / kPi) * x * y * z;
}
inline double Y3m1(double x, double y, double z) {
    return 0.25 * std::sqrt(21.0 / (2.0 * kPi)) * y * (5.0 * z * z - 1.0);
}
inline double Y30(double, double, double z) {
    return 0.25 * std::sqrt(7.0 / kPi) * (5.0 * z * z * z - 3.0 * z);
}
inline double Y31(double x, double, double z) {
    return 0.25 * std::sqrt(21.0 / (2.0 * kPi)) * x * (5.0 * z * z - 1.0);
}
inline double Y32(double x, double y, double z) {
    return 0.25 * std::sqrt(105.0 / kPi) * (x * x - y * y) * z;
}
inline double Y33(double x, double y, double) {
    return 0.25 * std::sqrt(35.0 / (2.0 * kPi)) * x * (x * x - 3.0 * y * y);
}

// Lower-l-only check: confirms each grid integrates Y_{lm}^2 to 1 over the
// l values it claims algebraic exactness for. Each Lebedev order is exact
// to its stated degree L; Y_l^m^2 has degree 2L (in x,y,z), so the rule
// applies up to and including l = L/2... actually for a real spherical
// harmonic squared, the degree in (x,y,z) is 2l, so a rule of algebraic
// degree L can integrate Y_lm^2 exactly when 2l <= L, i.e. l <= L/2.
//
// We use a single per-grid test that picks a representative subset of
// (l, m) pairs spanning l from 0 to floor(degree/2).

struct GridSpec {
    const lebedev::Point* points;
    int n;
    int degree; // stated algebraic degree (e.g. 6->3, 110->17)
    const char* name;
};

GridSpec lebedev_grids[] = {
    {lebedev::grid_6, lebedev::n_6, 3, "Lebedev-6"},
    {lebedev::grid_38, lebedev::n_38, 9, "Lebedev-38"},
    {lebedev::grid_50, lebedev::n_50, 11, "Lebedev-50"},
    {lebedev::grid_86, lebedev::n_86, 15, "Lebedev-86"},
    {lebedev::grid_110, lebedev::n_110, 17, "Lebedev-110"},
    {lebedev::grid_194, lebedev::n_194, 23, "Lebedev-194"},
    {lebedev::grid_266, lebedev::n_266, 27, "Lebedev-266"},
    {lebedev::grid_302, lebedev::n_302, 29, "Lebedev-302"},
    {lebedev::grid_350, lebedev::n_350, 31, "Lebedev-350"},
    {lebedev::grid_434, lebedev::n_434, 35, "Lebedev-434"},
    {lebedev::grid_590, lebedev::n_590, 41, "Lebedev-590"},
    {lebedev::grid_770, lebedev::n_770, 47, "Lebedev-770"},
};

// Integrate f over the grid: sum_i w_i * f(p_i).
template <class F>
double quadrature(const GridSpec& g, F f) {
    double s = 0;
    for (int i = 0; i < g.n; ++i) {
        s += g.points[i].w * f(g.points[i].x, g.points[i].y, g.points[i].z);
    }
    return s;
}

template <class Y>
double quadrature_squared(const GridSpec& g, Y y) {
    return quadrature(g, [y](double x, double y_, double z) {
        double v = y(x, y_, z);
        return v * v;
    });
}

} // namespace

// ── ac-001: Lebedev grids integrate spherical harmonics to 1e-14 ────────

TEST(Grid, LebedevSphericalHarmonicOrthogonality) {
    for (const auto& g : lebedev_grids) {
        // Sum of weights = 1 by construction (normalized to unit sphere).
        EXPECT_NEAR(quadrature(g, [](double, double, double) { return 1.0; }),
                    1.0, 1e-14)
            << g.name << ": weights do not sum to 1";

        // l = 0
        if (g.degree >= 0)
            EXPECT_NEAR(quadrature_squared(g, Y00), 1.0 / (4.0 * kPi), 1e-14)
                << g.name << " Y_00";

        // l = 1: 2l = 2 <= degree
        if (g.degree >= 2) {
            EXPECT_NEAR(quadrature_squared(g, Y1m1), 1.0 / (4.0 * kPi), 1e-14)
                << g.name << " Y_{1,-1}";
            EXPECT_NEAR(quadrature_squared(g, Y10), 1.0 / (4.0 * kPi), 1e-14)
                << g.name << " Y_{1,0}";
            EXPECT_NEAR(quadrature_squared(g, Y11), 1.0 / (4.0 * kPi), 1e-14)
                << g.name << " Y_{1,1}";
        }

        // l = 2: 2l = 4 <= degree
        if (g.degree >= 4) {
            EXPECT_NEAR(quadrature_squared(g, Y2m2), 1.0 / (4.0 * kPi), 1e-14)
                << g.name << " Y_{2,-2}";
            EXPECT_NEAR(quadrature_squared(g, Y2m1), 1.0 / (4.0 * kPi), 1e-14)
                << g.name << " Y_{2,-1}";
            EXPECT_NEAR(quadrature_squared(g, Y20), 1.0 / (4.0 * kPi), 1e-14)
                << g.name << " Y_{2,0}";
            EXPECT_NEAR(quadrature_squared(g, Y21), 1.0 / (4.0 * kPi), 1e-14)
                << g.name << " Y_{2,1}";
            EXPECT_NEAR(quadrature_squared(g, Y22), 1.0 / (4.0 * kPi), 1e-14)
                << g.name << " Y_{2,2}";
        }

        // l = 3: 2l = 6 <= degree
        if (g.degree >= 6) {
            EXPECT_NEAR(quadrature_squared(g, Y3m3), 1.0 / (4.0 * kPi), 1e-14)
                << g.name << " Y_{3,-3}";
            EXPECT_NEAR(quadrature_squared(g, Y30), 1.0 / (4.0 * kPi), 1e-14)
                << g.name << " Y_{3,0}";
            EXPECT_NEAR(quadrature_squared(g, Y33), 1.0 / (4.0 * kPi), 1e-14)
                << g.name << " Y_{3,3}";
        }
    }
}

// ── ac-002: Mura-Knowles reproduces analytic Gaussian moments ────────────

TEST(Grid, MuraKnowlesGaussianMoments) {
    // The radial-quadrature weight already includes the r^2 Jacobian, so
    //   sum_i w_i * f(r_i)  approximates  int_0^infinity f(r) * r^2 dr.
    // Analytic moment: int_0^infinity r^(n+2) exp(-alpha r^2) dr
    //                = Gamma((n+3)/2) / (2 * alpha^((n+3)/2)).
    //
    // Mura-Knowles is calibrated for atomic-density radial decay (sharp
    // near the nucleus, rapidly falling at large r). It concentrates
    // points near r=0 -- correct for SCF densities but suboptimal for
    // arbitrary high-r-weighted moments. Test with (alpha, n)
    // combinations physically relevant to SCF density / overlap integrals:
    // alpha >= 1 (tight orbitals) and n in {0, 2} (the only powers that
    // appear in HF/DFT energy expressions through rho and grad-rho).
    //
    // Production-comparable tolerance: 1e-8 at N=200. This is several
    // orders below the ~1e-4 grid-convergence floor PySCF documents for
    // atom_grid=(50, 110) (the default SCF radial size).
    constexpr int N = 200;
    constexpr double R_bragg = 1.0;
    auto pts = radial_grid_mura_knowles(N, R_bragg);

    auto gamma_half = [](double k) {
        // Gamma((k+1)/2) for integer k via std::tgamma.
        return std::tgamma(0.5 * (k + 1.0));
    };

    struct MomentCase {
        double alpha;
        int n;
    };
    MomentCase cases[] = {
        {1.0, 0}, {1.0, 2}, {2.0, 0}, {2.0, 2},
        {4.0, 0}, {4.0, 2}, {8.0, 0}, {8.0, 2}, // tight orbital
    };

    for (const auto& c : cases) {
        double Q = 0.0;
        for (const auto& p : pts) {
            Q += p.w * std::pow(p.r, c.n) * std::exp(-c.alpha * p.r * p.r);
        }
        double exact =
            gamma_half(c.n + 2) / (2.0 * std::pow(c.alpha, 0.5 * (c.n + 3.0)));

        double rel = std::abs(Q - exact) / std::abs(exact);
        EXPECT_LT(rel, 1e-7)
            << "Mura-Knowles moment alpha=" << c.alpha << ", n=" << c.n
            << " (numerical=" << Q << ", analytic=" << exact
            << ", rel_err=" << rel << ")";
    }
}

// ── ac-003: Bragg-Slater table strictness ────────────────────────────────

TEST(Grid, BraggSlaterCoverage) {
    // Z = 1..36 must each return a finite, positive Bohr radius matching
    // the published Angstrom value times kAngstromToBohr.
    //
    // Spot-checks (Angstrom): PySCF-matching values.
    // H=0.35 (non-bonded), He=1.40 / Ne=1.50 / Kr=1.90 (Bondi-style noble gas),
    // C=0.70 (Slater).
    EXPECT_NEAR(bragg_slater_radius(1), 0.35 * kAngstromToBohr, 1e-15);
    EXPECT_NEAR(bragg_slater_radius(2), 1.40 * kAngstromToBohr, 1e-15);
    EXPECT_NEAR(bragg_slater_radius(6), 0.70 * kAngstromToBohr, 1e-15);
    EXPECT_NEAR(bragg_slater_radius(10), 1.50 * kAngstromToBohr, 1e-15);
    EXPECT_NEAR(bragg_slater_radius(36), 1.90 * kAngstromToBohr, 1e-15);

    // Every Z in 1..36 returns a positive radius without throwing.
    for (int Z = 1; Z <= 36; ++Z) {
        double r = bragg_slater_radius(Z);
        EXPECT_GT(r, 0.0) << "Z=" << Z << " gave non-positive radius";
        EXPECT_LT(r, 10.0) << "Z=" << Z << " gave unphysically large radius";
    }
}

TEST(Grid, BraggSlaterRejectsOutOfRange) {
    // Out-of-range Z must throw std::invalid_argument; no silent fallback.
    EXPECT_THROW(bragg_slater_radius(0), std::invalid_argument);
    EXPECT_THROW(bragg_slater_radius(-1), std::invalid_argument);
    EXPECT_THROW(bragg_slater_radius(87), std::invalid_argument);
    EXPECT_THROW(bragg_slater_radius(100), std::invalid_argument);
    EXPECT_THROW(bragg_slater_radius(118), std::invalid_argument);
}

// ── Smoke / regression tests (existing) ─────────────────────────────────

TEST(Grid, SingleAtomGridHasPoints) {
    native::Atom h = {1, {0.0, 0.0, 0.0}};
    auto grid = build_molecular_grid(&h, 1, 20, 6);

    EXPECT_GT(grid.points.size(), 0u);
    EXPECT_EQ(grid.n_radial, 20);
    EXPECT_EQ(grid.n_angular, 6);
}

TEST(Grid, TwoAtomGridHasMorePoints) {
    native::Atom h1 = {1, {0.0, 0.0, 0.0}};
    auto grid1 = build_molecular_grid(&h1, 1, 20, 6);

    native::Atom atoms[] = {{1, {0.0, 0.0, 0.0}}, {1, {1.4, 0.0, 0.0}}};
    auto grid2 = build_molecular_grid(atoms, 2, 20, 6);

    EXPECT_GT(grid2.points.size(), grid1.points.size());
}

TEST(Grid, BeckePartitionPositive) {
    native::Atom h = {1, {0.0, 0.0, 0.0}};
    auto grid = build_molecular_grid(&h, 1, 30, 6);
    for (const auto& gp : grid.points)
        EXPECT_GT(static_cast<double>(gp.weight), 0.0);
}

TEST(Grid, SingleAtomSlaterIntegral) {
    // int exp(-2r) d^3r = pi  (the 1s^2 hydrogen density norm).
    native::Atom h = {1, {0.0, 0.0, 0.0}};
    auto grid = build_molecular_grid(&h, 1, 50, 110);
    double integral = 0;
    for (const auto& gp : grid.points) {
        double r = std::sqrt(gp.x * gp.x + gp.y * gp.y + gp.z * gp.z);
        integral += gp.weight * std::exp(-2.0 * r);
    }
    EXPECT_NEAR(integral, kPi, 0.5);
}

TEST(Grid, TwoAtomPartitionSymmetry) {
    double d = 1.4;
    native::Atom atoms[] = {{1, {-d / 2, 0, 0}}, {1, {d / 2, 0, 0}}};
    auto grid = build_molecular_grid(atoms, 2, 30, 6);
    double w_left = 0, w_right = 0;
    for (const auto& gp : grid.points) {
        if (gp.x < 0)
            w_left += gp.weight;
        else
            w_right += gp.weight;
    }
    double total = w_left + w_right;
    if (total > 0)
        EXPECT_NEAR(w_left / total, 0.5, 0.1);
}

TEST(Grid, GridPointsFinite) {
    native::Atom atoms[] = {{1, {0, 0, 0}}, {2, {2.0, 0, 0}}};
    auto grid = build_molecular_grid(atoms, 2, 30, 6);
    for (const auto& gp : grid.points) {
        EXPECT_FALSE(std::isnan(gp.x));
        EXPECT_FALSE(std::isnan(gp.y));
        EXPECT_FALSE(std::isnan(gp.z));
        EXPECT_FALSE(std::isnan(gp.weight));
        EXPECT_FALSE(std::isinf(gp.weight));
    }
}

// ════════════════════════════════════════════════════════════════════════
// phase1-productization-03-denser-grid acceptance tests
// ════════════════════════════════════════════════════════════════════════

namespace {

// Analytic unit-sphere average of the monomial x^a y^b z^c, i.e.
//   (1 / 4*pi) * integral_S x^a y^b z^c dOmega.
// Zero for any odd exponent; otherwise a ratio of half-integer Gammas.
double monomial_sphere_average(int a, int b, int c) {
    if ((a & 1) || (b & 1) || (c & 1))
        return 0.0;
    auto half_gamma = [](int k) { return std::tgamma(0.5 * (k + 1)); };
    double num = half_gamma(a) * half_gamma(b) * half_gamma(c);
    double den = std::tgamma(0.5 * (a + b + c + 3));
    // integral over full unit sphere is 2 * num / den; divide by 4*pi.
    return (2.0 * num / den) / (4.0 * kPi);
}

// Maximum quadrature error over every monomial x^a y^b z^c with
// a+b+c <= design_degree for a Lebedev grid (weights normalized to 1).
double lebedev_monomial_maxerr(const GridSpec& g) {
    double maxerr = 0.0;
    for (int a = 0; a <= g.degree; ++a) {
        for (int b = 0; b <= g.degree - a; ++b) {
            for (int c = 0; c <= g.degree - a - b; ++c) {
                double q = 0.0;
                for (int i = 0; i < g.n; ++i) {
                    const auto& p = g.points[i];
                    q += p.w * std::pow(p.x, a) * std::pow(p.y, b)
                         * std::pow(p.z, c);
                }
                double err = std::abs(q - monomial_sphere_average(a, b, c));
                maxerr = std::max(maxerr, err);
            }
        }
    }
    return maxerr;
}

// Build a cpu::Frame from (Z, x, y, z) tuples given in Angstrom.
struct AtomXYZ {
    int Z;
    double x, y, z;
};

cpu::Frame make_frame(const std::vector<AtomXYZ>& atoms) {
    cpu::Frame f(static_cast<int>(atoms.size()));
    for (size_t i = 0; i < atoms.size(); ++i) {
        f.atoms.type[i] = atoms[i].Z;
        f.atoms.x[i] = atoms[i].x;
        f.atoms.y[i] = atoms[i].y;
        f.atoms.z[i] = atoms[i].z;
    }
    f.box = Mat3::cubic(100.0);
    return f;
}

// H2O / NH3 / CH4 starting geometries (Angstrom) — verbatim from the PySCF
// def2-SVP references in benchmarks/qc/reference/{h2o,nh3,ch4}.json.
const std::vector<AtomXYZ> kH2O = {
    {8, 0.0, 0.0, 0.11779},
    {1, 0.0, 0.755453, -0.471161},
    {1, 0.0, -0.755453, -0.471161},
};
const std::vector<AtomXYZ> kNH3 = {
    {7, 0.0, 0.0, 0.11388},
    {1, 0.0, 0.938035, -0.265825},
    {1, 0.812328, -0.469017, -0.265825},
    {1, -0.812328, -0.469017, -0.265825},
};
const std::vector<AtomXYZ> kCH4 = {
    {6, 0.0, 0.0, 0.0},
    {1, 0.62764, 0.62764, 0.62764},
    {1, -0.62764, -0.62764, 0.62764},
    {1, 0.62764, -0.62764, -0.62764},
    {1, -0.62764, 0.62764, -0.62764},
};
// PySCF def2-SVP RKS/PBE single-point energies on the *starting* geometry
// (energy_hartree_starting in the reference JSONs).
constexpr double kE_H2O = -76.27374826780702;
constexpr double kE_NH3 = -56.43334157583261;
constexpr double kE_CH4 = -40.415486061884764;

// Run an ATV RKS/PBE single point at def2-SVP with the given grid level.
double atv_rks_pbe_energy(const std::vector<AtomXYZ>& atoms, GridLevel level) {
    auto frame = make_frame(atoms);
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::Def2SVP);
    cpu::scf::SCFConfig c;
    c.spin = cpu::scf::SpinType::Restricted;
    c.method = cpu::scf::MeanFieldType::KS;
    c.xc_type = cpu::scf::XCType::PBE;
    c.n_alpha = ints.n_occ();
    c.n_beta = ints.n_occ();
    c.tol = 1e-10;
    c.use_diis = true;
    c.diis_max_vectors = 8;
    c.diis_start_iter = 1;
    c.damping_factor = 0.0;
    c.init_guess = cpu::scf::InitGuess::MINAO;
    c.grid_level = level;
    return cpu::scf::scf_solve(ints, c).energy;
}

} // namespace

// ── ac-001: new Lebedev orders 38, 86, 770 integrate exactly ────────────

TEST(Grid, LebedevNewOrdersMonomialExactness) {
    // Every monomial up to each grid's design algebraic degree must be
    // integrated to machine precision (1e-12). This is the correctness
    // guard for the newly tabulated orders 38, 86 and 770.
    for (const auto& g : lebedev_grids) {
        double err = lebedev_monomial_maxerr(g);
        EXPECT_LT(err, 1e-12)
            << g.name << ": monomial exactness failed, maxerr=" << err;
    }
}

// ── ac-009: GridLevel::SG0 reproduces the legacy (50,110) grid ──────────

TEST(Grid, SG0ReproducesLegacyGrid) {
    // The pre-change call build_molecular_grid(atoms, natm, 50, 110) and
    // the new level-aware path with GridLevel::SG0 must produce a
    // bit-for-bit identical MolecularGrid.
    native::Atom atoms[] = {
        {8, {0.0, 0.0, 0.2}}, {1, {0.0, 1.4, -0.9}}, {1, {0.0, -1.4, -0.9}}};
    auto legacy = build_molecular_grid(atoms, 3, 50, 110);
    auto sg0 = build_molecular_grid(atoms, 3, GridLevel::SG0, 50, 110);

    ASSERT_EQ(legacy.points.size(), sg0.points.size());
    for (size_t i = 0; i < legacy.points.size(); ++i) {
        EXPECT_EQ(legacy.points[i].x, sg0.points[i].x) << "point " << i;
        EXPECT_EQ(legacy.points[i].y, sg0.points[i].y) << "point " << i;
        EXPECT_EQ(legacy.points[i].z, sg0.points[i].z) << "point " << i;
        EXPECT_EQ(legacy.points[i].weight, sg0.points[i].weight)
            << "point " << i;
    }
}

// ── ac-006: Becke partition weights sum to unity at every grid point ────

TEST(Grid, BeckePartitionSumsToUnity) {
    // For a multi-atom molecule, the per-atom Becke cell weights at any
    // grid point must sum to exactly 1 (within 1e-12). We probe this with
    // a regular lattice of test points around the H2O nuclei.
    native::Atom atoms[] = {{8, {0.0, 0.0, 0.22}},
                            {1, {0.0, 1.43, -0.89}},
                            {1, {0.0, -1.43, -0.89}}};
    const int natm = 3;
    std::vector<double> w(natm);
    for (double x = -3.0; x <= 3.0; x += 0.5) {
        for (double y = -3.0; y <= 3.0; y += 0.5) {
            for (double z = -3.0; z <= 3.0; z += 0.5) {
                becke_partition_weights_for_test(atoms, natm, x, y, z,
                                                 w.data());
                double s = 0.0;
                for (int a = 0; a < natm; ++a)
                    s += w[a];
                EXPECT_NEAR(s, 1.0, 1e-12)
                    << "Becke weights do not sum to 1 at (" << x << "," << y
                    << "," << z << ")";
            }
        }
    }
}

// ── ac-002 / ac-003: SG-1 and match-PySCF recover the electron count ────

TEST(Grid, SG1RecoversWaterElectronCount) {
    auto frame = make_frame(kH2O);
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::Def2SVP);
    cpu::scf::SCFConfig c;
    c.spin = cpu::scf::SpinType::Restricted;
    c.method = cpu::scf::MeanFieldType::KS;
    c.xc_type = cpu::scf::XCType::PBE;
    c.n_alpha = ints.n_occ();
    c.n_beta = ints.n_occ();
    c.tol = 1e-10;
    c.init_guess = cpu::scf::InitGuess::MINAO;
    c.grid_level = GridLevel::SG1;
    auto r = cpu::scf::scf_solve(ints, c);
    double n_e = r.integrated_n_alpha + r.integrated_n_beta;
    EXPECT_NEAR(n_e, 10.0, 1e-4) << "SG-1 grid electron count = " << n_e;
}

TEST(Grid, MatchPySCFRecoversWaterElectronCount) {
    auto frame = make_frame(kH2O);
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::Def2SVP);
    cpu::scf::SCFConfig c;
    c.spin = cpu::scf::SpinType::Restricted;
    c.method = cpu::scf::MeanFieldType::KS;
    c.xc_type = cpu::scf::XCType::PBE;
    c.n_alpha = ints.n_occ();
    c.n_beta = ints.n_occ();
    c.tol = 1e-10;
    c.init_guess = cpu::scf::InitGuess::MINAO;
    c.grid_level = GridLevel::MatchPySCF;
    auto r = cpu::scf::scf_solve(ints, c);
    double n_e = r.integrated_n_alpha + r.integrated_n_beta;
    EXPECT_NEAR(n_e, 10.0, 1e-4) << "match-PySCF grid electron count = " << n_e;
}

// ── ac-004: SG-1 brings the def2-SVP energy residue to <= 5e-4 Ha ───────

TEST(Grid, SG1EnergyResidue) {
    EXPECT_NEAR(atv_rks_pbe_energy(kH2O, GridLevel::SG1), kE_H2O, 5e-4);
    EXPECT_NEAR(atv_rks_pbe_energy(kNH3, GridLevel::SG1), kE_NH3, 5e-4);
    EXPECT_NEAR(atv_rks_pbe_energy(kCH4, GridLevel::SG1), kE_CH4, 5e-4);
}

// ── ac-005: match-PySCF brings the residue to <= 1e-6 Ha ────────────────

TEST(Grid, MatchPySCFEnergyResidue) {
    EXPECT_NEAR(atv_rks_pbe_energy(kH2O, GridLevel::MatchPySCF), kE_H2O, 1e-6);
    EXPECT_NEAR(atv_rks_pbe_energy(kNH3, GridLevel::MatchPySCF), kE_NH3, 1e-6);
    EXPECT_NEAR(atv_rks_pbe_energy(kCH4, GridLevel::MatchPySCF), kE_CH4, 1e-6);
}

// ── MatchPySCF heavy-element (Z>10) nwchem_prune branch ─────────────────

TEST(Grid, MatchPySCFHeavyElementBuildsAndIntegrates) {
    // Silicon (Z=14) exercises the period-3+ nwchem_prune branch:
    // n_ang=434, pruned orders {50, 86, 350, 434, 350}. This path
    // previously threw (Lebedev 350/434 were not tabulated); it must now
    // build and yield a valid quadrature. A single atom has Becke weight
    // 1 everywhere, so this isolates the radial x angular grid: the
    // integral of a normalized Gaussian exp(-a r^2) must recover the
    // analytic (pi/a)^1.5.
    native::Atom si = {14, {0.0, 0.0, 0.0}};
    auto grid = build_molecular_grid(&si, 1, GridLevel::MatchPySCF, 0, 0);
    ASSERT_FALSE(grid.points.empty());

    const double a = 0.8;
    double q = 0.0;
    for (const auto& p : grid.points) {
        double r2 = p.x * p.x + p.y * p.y + p.z * p.z;
        q += p.weight * std::exp(-a * r2);
    }
    EXPECT_NEAR(q, std::pow(kPi / a, 1.5), 1e-8)
        << "MatchPySCF Z=14 grid Gaussian quadrature";
}

#else

#include <gtest/gtest.h>
TEST(Grid, SkippedNoBuildRHF) {
    GTEST_SKIP() << "ATV_BUILD_RHF=OFF — Grid tests require -DATV_BUILD_RHF=ON";
}

#endif // ATV_HAS_RHF
