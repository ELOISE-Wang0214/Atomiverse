/// @file tests/cpu/scf/test_rys_eri_fshell.cpp
/// @brief f-shell (l=3) Rys-quadrature ERI correctness tests.
///
/// Exercises the engine after the l<=2 -> l<=3 bump: f-shell quartets
/// push the total angular momentum to Ltot = 12 ((ff|ff)) and Ltot = 13
/// for the first derivative of an (ff|ff) quartet, which need the
/// iterative 6/7-root Rys finder.
///
/// Three checks:
///   1. Finiteness / non-zero (ac-005): rys_eri_block produces finite,
///      non-trivial blocks for (f,s|s,s), (f,p|s,s), (f,d|s,s),
///      (f,f|f,f).
///   2. FD self-consistency (ac-006): rys_eri_block_deriv agrees with a
///      central finite difference of rys_eri_block (h=1e-4) — this
///      exercises the +1 derivative angular shift and the 7-root budget.
///   3. libcint cross-check (ac-007): one f-quartet ERI matches
///      build_eri_libcint_from_pack's l=3 path (gated on ATV_HAS_LIBCINT).

#ifdef ATV_HAS_RHF

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <vector>

#include <xtensor/containers/xtensor.hpp>

#include "atomiverse/cpu/scf/integrals.h"
#include "atomiverse/cpu/scf/native_integrals.h"
#include "atomiverse/cpu/scf/rys_eri.h"
#include "atomiverse/cpu/scf/shell_pair.h"

#ifdef ATV_HAS_LIBCINT
#include "atomiverse/core/units.h"
#include "libcint_reference.cuh"
#endif

namespace {

using namespace ATV::cpu::scf::native;

/// @brief Minimal contracted-shell descriptor for the Rys engine.
struct PrimShell {
    int l;
    double alpha[3];
    double coeff[3];
    int n_prim;
    double center[3];
};

/// @brief Convert a PrimShell into a ShellBlock (pointers alias the
/// caller-owned PrimShell, which must outlive the ShellBlock).
ShellBlock make_shell_block(const PrimShell& s) {
    ShellBlock b;
    b.n_prim = s.n_prim;
    b.alpha = s.alpha;
    b.coeff = s.coeff;
    b.center[0] = s.center[0];
    b.center[1] = s.center[1];
    b.center[2] = s.center[2];
    b.l = s.l;
    b.n_cart = (s.l + 1) * (s.l + 2) / 2;
    b.ao_start = 0;
    return b;
}

/// @brief Build a test shell of angular momentum `l` with 3 primitives.
/// Exponents/coefficients resemble real contracted Gaussians but do not
/// depend on any basis table.
PrimShell make_shell(int l, const double center[3]) {
    PrimShell s;
    s.l = l;
    s.n_prim = 3;
    s.alpha[0] = 3.0;
    s.coeff[0] = 0.4;
    s.alpha[1] = 1.0;
    s.coeff[1] = 0.5;
    s.alpha[2] = 0.3;
    s.coeff[2] = 0.3;
    s.center[0] = center[0];
    s.center[1] = center[1];
    s.center[2] = center[2];
    return s;
}

/// @brief Compute the full ERI block (A B | C D) via rys_eri_block.
std::vector<double> eri_block(const PrimShell& A, const PrimShell& B,
                              const PrimShell& C, const PrimShell& D) {
    const auto Ab = make_shell_block(A);
    const auto Bb = make_shell_block(B);
    const auto Cb = make_shell_block(C);
    const auto Db = make_shell_block(D);
    std::vector<double> out(static_cast<size_t>(Ab.n_cart) * Bb.n_cart
                                * Cb.n_cart * Db.n_cart,
                            0.0);
    rys_eri_block(Ab, Bb, Cb, Db, out.data());
    return out;
}

/// @brief Compute the 3-axis derivative block via rys_eri_block_deriv.
std::vector<double> deriv_block(const PrimShell& A, const PrimShell& B,
                                const PrimShell& C, const PrimShell& D,
                                int deriv_center) {
    const auto Ab = make_shell_block(A);
    const auto Bb = make_shell_block(B);
    const auto Cb = make_shell_block(C);
    const auto Db = make_shell_block(D);
    const size_t n =
        static_cast<size_t>(Ab.n_cart) * Bb.n_cart * Cb.n_cart * Db.n_cart;
    std::vector<double> out(3 * n, 0.0);
    rys_eri_block_deriv(Ab, Bb, Cb, Db, deriv_center, out.data());
    return out;
}

/// @brief Central finite difference of the ERI block w.r.t. one
/// Cartesian component of one shell centre.
std::vector<double> fd_deriv_block(PrimShell A, PrimShell B, PrimShell C,
                                   PrimShell D, int deriv_center, int axis,
                                   double h) {
    PrimShell& target = (deriv_center == 0)   ? A
                        : (deriv_center == 1) ? B
                        : (deriv_center == 2) ? C
                                              : D;
    const double saved = target.center[axis];
    target.center[axis] = saved + h;
    auto plus = eri_block(A, B, C, D);
    target.center[axis] = saved - h;
    auto minus = eri_block(A, B, C, D);
    target.center[axis] = saved;
    std::vector<double> fd(plus.size(), 0.0);
    for (size_t i = 0; i < plus.size(); ++i)
        fd[i] = (plus[i] - minus[i]) / (2.0 * h);
    return fd;
}

/// @brief Assert every entry of `block` is finite, and at least one is
/// non-trivially non-zero.
void expect_finite_nonzero(const std::vector<double>& block,
                           const char* label) {
    bool any_nonzero = false;
    for (double v : block) {
        ASSERT_TRUE(std::isfinite(v)) << label << ": non-finite ERI entry";
        if (std::abs(v) > 1e-12)
            any_nonzero = true;
    }
    EXPECT_TRUE(any_nonzero) << label << ": all ERI entries are ~zero";
}

/// @brief Run the analytic-deriv-vs-finite-difference comparison for one
/// quartet and one deriv-centre, over all 3 axes and all block entries.
void check_deriv_vs_fd(const PrimShell& A, const PrimShell& B,
                       const PrimShell& C, const PrimShell& D, int deriv_center,
                       const char* label, double tol = 1e-6) {
    const auto analytic = deriv_block(A, B, C, D, deriv_center);
    const double h = 1e-4;
    const std::array<std::vector<double>, 3> fd = {
        fd_deriv_block(A, B, C, D, deriv_center, 0, h),
        fd_deriv_block(A, B, C, D, deriv_center, 1, h),
        fd_deriv_block(A, B, C, D, deriv_center, 2, h),
    };
    const size_t n_block = fd[0].size();
    double max_diff = 0.0;
    int worst_axis = -1;
    long worst_idx = -1;
    for (int axis = 0; axis < 3; ++axis) {
        for (size_t i = 0; i < n_block; ++i) {
            const double a = analytic[axis * n_block + i];
            const double f = fd[axis][i];
            // Mixed relative/absolute: |a-f| <= tol*(1+|a|) — the FD
            // noise floor scales with the integral magnitude.
            const double d = std::abs(a - f);
            const double scaled = d / (1.0 + std::abs(a));
            if (scaled > max_diff) {
                max_diff = scaled;
                worst_axis = axis;
                worst_idx = static_cast<long>(i);
            }
        }
    }
    EXPECT_LT(max_diff, tol)
        << label << ": max scaled |analytic - FD| = " << max_diff
        << " at axis=" << worst_axis << " idx=" << worst_idx;
}

// ── Shared centres (non-collinear, generic) ──────────────────────────
const double RA[3] = {0.0, 0.0, 0.0};
const double RB[3] = {0.9, 0.0, 0.0};
const double RC[3] = {0.3, 0.7, 0.0};
const double RD[3] = {-0.2, 0.1, 0.5};

/// @brief Build a single methane (CH4) frame — carbon carries the
/// def2-TZVP f-shell, so its integral pack exercises f-quartets.
ATV::cpu::Frame make_ch4() {
    using ATV::Mat3;
    using ATV::cpu::Frame;
    Frame f(5);
    // Tetrahedral geometry, C-H ~ 1.09 Angstrom.
    const double d = 0.6291;
    const double cx[5] = {0.0, d, d, -d, -d};
    const double cy[5] = {0.0, d, -d, d, -d};
    const double cz[5] = {0.0, d, -d, -d, d};
    const int z[5] = {6, 1, 1, 1, 1};
    for (int i = 0; i < 5; ++i) {
        f.span(ATV::kAtom_x)[i] = cx[i];
        f.span(ATV::kAtom_y)[i] = cy[i];
        f.span(ATV::kAtom_z)[i] = cz[i];
        f.span(ATV::kAtom_type)[i] = z[i];
        f.span(ATV::kAtom_mass)[i] = (z[i] == 6) ? 12.011 : 1.008;
    }
    f.simbox = ATV::cpu::Box::cubic(100.0);
    return f;
}

} // namespace

// ── ac-005: f-quartet ERIs are finite and non-zero ───────────────────

TEST(RysEriFShell, FsSs_FiniteNonZero) {
    const auto A = make_shell(3, RA); // f
    const auto B = make_shell(0, RB); // s
    const auto C = make_shell(0, RC); // s
    const auto D = make_shell(0, RD); // s
    expect_finite_nonzero(eri_block(A, B, C, D), "(f,s|s,s)");
}

TEST(RysEriFShell, FpSs_FiniteNonZero) {
    const auto A = make_shell(3, RA); // f
    const auto B = make_shell(1, RB); // p
    const auto C = make_shell(0, RC); // s
    const auto D = make_shell(0, RD); // s
    expect_finite_nonzero(eri_block(A, B, C, D), "(f,p|s,s)");
}

TEST(RysEriFShell, FdSs_FiniteNonZero) {
    const auto A = make_shell(3, RA); // f
    const auto B = make_shell(2, RB); // d
    const auto C = make_shell(0, RC); // s
    const auto D = make_shell(0, RD); // s
    expect_finite_nonzero(eri_block(A, B, C, D), "(f,d|s,s)");
}

TEST(RysEriFShell, FfFf_FiniteNonZero) {
    const auto A = make_shell(3, RA); // f
    const auto B = make_shell(3, RB); // f
    const auto C = make_shell(3, RC); // f
    const auto D = make_shell(3, RD); // f
    auto block = eri_block(A, B, C, D);
    ASSERT_EQ(block.size(), 10u * 10u * 10u * 10u)
        << "(f,f|f,f) block must be 10^4 Cartesian entries";
    expect_finite_nonzero(block, "(f,f|f,f)");
}

// ── ac-006: f-quartet ERI derivatives pass FD self-consistency ───────

TEST(RysEriFShell, FsSs_DerivVsFD) {
    const auto A = make_shell(3, RA);
    const auto B = make_shell(0, RB);
    const auto C = make_shell(0, RC);
    const auto D = make_shell(0, RD);
    for (int dc = 0; dc < 4; ++dc)
        check_deriv_vs_fd(A, B, C, D, dc, "(f,s|s,s)");
}

TEST(RysEriFShell, FpSs_DerivVsFD) {
    const auto A = make_shell(3, RA);
    const auto B = make_shell(1, RB);
    const auto C = make_shell(0, RC);
    const auto D = make_shell(0, RD);
    for (int dc = 0; dc < 4; ++dc)
        check_deriv_vs_fd(A, B, C, D, dc, "(f,p|s,s)");
}

TEST(RysEriFShell, FdSs_DerivVsFD) {
    const auto A = make_shell(3, RA);
    const auto B = make_shell(2, RB);
    const auto C = make_shell(0, RC);
    const auto D = make_shell(0, RD);
    for (int dc = 0; dc < 4; ++dc)
        check_deriv_vs_fd(A, B, C, D, dc, "(f,d|s,s)");
}

TEST(RysEriFShell, FfFf_DerivVsFD) {
    const auto A = make_shell(3, RA);
    const auto B = make_shell(3, RB);
    const auto C = make_shell(3, RC);
    const auto D = make_shell(3, RD);
    // The (f,f|f,f) derivative reaches Ltot = 13 -> 7-root quadrature.
    for (int dc = 0; dc < 4; ++dc)
        check_deriv_vs_fd(A, B, C, D, dc, "(f,f|f,f)");
}

// ── Edge case: non-canonical bra ordering with an f-shell ────────────
//
// rys_eri_block requires li >= lj. The deriv wrapper swaps shells to
// enforce canonical order internally; verify a (p,f|s,s) request (lj >
// li on the bra) still produces a finite derivative agreeing with FD.
TEST(RysEriFShell, PfSs_NonCanonicalDerivVsFD) {
    const auto A = make_shell(1, RA); // p  (lj > li on the bra)
    const auto B = make_shell(3, RB); // f
    const auto C = make_shell(0, RC); // s
    const auto D = make_shell(0, RD); // s
    for (int dc = 0; dc < 4; ++dc)
        check_deriv_vs_fd(A, B, C, D, dc, "(p,f|s,s)");
}

// ── ac-008: def2-TZVP build_integrals succeeds end-to-end ────────────
//
// A real f-shell basis (def2-TZVP places an f-shell on carbon) must run
// through build_integrals to completion and yield a fully finite,
// non-zero ERI tensor — the end-user-visible outcome of f-shell support.
TEST(RysEriFShell, Def2TZVPBuildIntegralsCH4) {
    auto frame = make_ch4();
    auto pack = ATV::cpu::scf::build_integrals(
        frame, ATV::cpu::scf::BasisSet::Def2TZVP, 0, /*need_eri=*/true);

    ASSERT_GT(pack.n_basis, 0) << "def2-TZVP CH4 must resolve a basis";

    // The shell list must actually contain an l=3 (f) shell, otherwise
    // the test would silently not exercise f-shells.
    auto sp = ATV::cpu::scf::build_shell_block_pair_list(pack, 1e-14);
    bool has_f = false;
    for (const auto& b : sp.blocks)
        if (b.l == 3)
            has_f = true;
    ASSERT_TRUE(has_f) << "def2-TZVP CH4 must carry an f-shell on carbon";

    // ERI tensor: every entry finite, the diagonal strictly positive
    // (a self-repulsion integral is always > 0), tensor non-trivial.
    const int n = pack.n_basis;
    ASSERT_EQ(static_cast<int>(pack.ERI.shape()[0]), n);
    bool any_nonzero = false;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            for (int k = 0; k < n; ++k)
                for (int l = 0; l < n; ++l) {
                    const double v = pack.ERI(i, j, k, l);
                    ASSERT_TRUE(std::isfinite(v))
                        << "non-finite ERI at (" << i << "," << j << "," << k
                        << "," << l << ")";
                    if (std::abs(v) > 1e-12)
                        any_nonzero = true;
                }
    EXPECT_TRUE(any_nonzero) << "def2-TZVP ERI tensor is all-zero";

    for (int i = 0; i < n; ++i) {
        EXPECT_GT(pack.ERI(i, i, i, i), 0.0)
            << "diagonal self-repulsion (ii|ii) must be positive at i=" << i;
    }
}

// ── ac-007: f-quartet ERI matches the libcint reference ──────────────

#ifdef ATV_HAS_LIBCINT

namespace {

using ATV::cpu::scf::ShellBlockPairList;

/// @brief Build the full (ij|kl) ERI tensor by looping shell-block
/// quartets through rys_eri_block, exploiting 8-fold permutational
/// symmetry to fill non-canonical orderings the li>=lj / lk>=ll guards
/// skip. Mirrors build_eri_via_rys in test_rys_eri.cpp.
xt::xtensor<double, 4> build_eri_via_rys_f(const ShellBlockPairList& sp) {
    const int n = sp.n_basis;
    xt::xtensor<double, 4> ERI = xt::zeros<double>(
        {static_cast<std::size_t>(n), static_cast<std::size_t>(n),
         static_cast<std::size_t>(n), static_cast<std::size_t>(n)});
    const int nb = static_cast<int>(sp.blocks.size());
    for (int bi = 0; bi < nb; ++bi) {
        const auto& Bi = sp.blocks[bi];
        for (int bj = 0; bj < nb; ++bj) {
            const auto& Bj = sp.blocks[bj];
            if (Bi.l < Bj.l)
                continue; // rys_eri_block needs li>=lj
            for (int bk = 0; bk < nb; ++bk) {
                const auto& Bk = sp.blocks[bk];
                for (int bl = 0; bl < nb; ++bl) {
                    const auto& Bl = sp.blocks[bl];
                    if (Bk.l < Bl.l)
                        continue;
                    const int na = Bi.n_cart, nbc = Bj.n_cart;
                    const int nc = Bk.n_cart, nd = Bl.n_cart;
                    std::vector<double> buf(
                        static_cast<size_t>(na) * nbc * nc * nd, 0.0);
                    ATV::cpu::scf::native::rys_eri_block(Bi, Bj, Bk, Bl,
                                                         buf.data());
                    for (int ia = 0; ia < na; ++ia)
                        for (int ib = 0; ib < nbc; ++ib)
                            for (int ic = 0; ic < nc; ++ic)
                                for (int id = 0; id < nd; ++id) {
                                    const double v =
                                        buf[((ia * nbc + ib) * nc + ic) * nd
                                            + id];
                                    const int ai = Bi.ao_start + ia;
                                    const int aj = Bj.ao_start + ib;
                                    const int ak = Bk.ao_start + ic;
                                    const int al = Bl.ao_start + id;
                                    ERI(ai, aj, ak, al) = v;
                                    ERI(aj, ai, ak, al) = v;
                                    ERI(ai, aj, al, ak) = v;
                                    ERI(aj, ai, al, ak) = v;
                                }
                }
            }
        }
    }
    return ERI;
}

} // namespace

TEST(RysEriFShell, FQuartetMatchesLibcint_Def2TZVP_CH4) {
    auto frame = make_ch4();
    // def2-TZVP places an f-shell on carbon -> f-quartets in the tensor.
    auto pack = ATV::cpu::scf::build_integrals(
        frame, ATV::cpu::scf::BasisSet::Def2TZVP, 0, /*need_eri=*/true);
    auto sp = ATV::cpu::scf::build_shell_block_pair_list(pack, 1e-14);

    // Confirm an f-shell is actually present (def2-TZVP C has l=3).
    bool has_f = false;
    for (const auto& b : sp.blocks)
        if (b.l == 3)
            has_f = true;
    ASSERT_TRUE(has_f) << "def2-TZVP CH4 must carry an f-shell on carbon";

    auto eri_rys = build_eri_via_rys_f(sp);
    auto eri_ref = ATV::test::build_eri_libcint_from_pack(frame, pack);

    const int n = static_cast<int>(eri_rys.shape()[0]);
    double worst = 0.0;
    int wi = -1, wj = -1, wk = -1, wl = -1;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            for (int k = 0; k < n; ++k)
                for (int l = 0; l < n; ++l) {
                    const double diff =
                        std::abs(eri_rys(i, j, k, l) - eri_ref(i, j, k, l));
                    if (diff > worst) {
                        worst = diff;
                        wi = i;
                        wj = j;
                        wk = k;
                        wl = l;
                    }
                }
    EXPECT_LT(worst, 1e-9)
        << "def2-TZVP CH4 f-quartet ERI vs libcint: worst diff " << worst
        << " at (" << wi << "," << wj << "," << wk << "," << wl << ")";
}

#endif // ATV_HAS_LIBCINT

#endif // ATV_HAS_RHF
