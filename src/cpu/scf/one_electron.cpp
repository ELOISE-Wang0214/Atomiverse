/// @file cpu/scf/one_electron.cpp
/// @brief Native-engine implementations of the one-electron integral kernels
///        (IntegralPack private static members) + nuclear repulsion energy.

#include "atomiverse/cpu/scf/internal/one_electron.h"

#include <cmath>

#include "atomiverse/cpu/scf/cart.h"
#include "atomiverse/cpu/scf/internal/integrals.h"
#include "atomiverse/cpu/scf/internal/shell_pair.h"
#include "atomiverse/cpu/scf/native_integrals.h"

namespace ATV {
namespace cpu {
namespace scf {

namespace {

/// Max buffer needed for a Cartesian (AB) block of 1e integrals: the
/// largest case is d×d = 6×6 = 36 doubles.
constexpr int kBlockBuf1e = kMaxNCart * kMaxNCart;

/// @brief Scatter a block's result into a symmetric matrix.
///
/// `buf` is the dense A.n_cart × B.n_cart output from a block kernel.
/// For i != j block pairs we also populate the transposed entries; for
/// diagonal block pairs the kernel already produces the whole block.
inline void scatter_block(xt::xtensor<double, 2>& M,
                          const native::ShellBlock& A,
                          const native::ShellBlock& B, const double* buf,
                          bool diag) {
    const int na = A.n_cart, nb = B.n_cart;
    if (diag) {
        for (int ia = 0; ia < na; ++ia) {
            for (int ib = 0; ib < nb; ++ib) {
                M(A.ao_start + ia, B.ao_start + ib) = buf[ia * nb + ib];
            }
        }
    } else {
        for (int ia = 0; ia < na; ++ia) {
            const int i = A.ao_start + ia;
            for (int ib = 0; ib < nb; ++ib) {
                const int j = B.ao_start + ib;
                const double v = buf[ia * nb + ib];
                M(i, j) = v;
                M(j, i) = v;
            }
        }
    }
}

/// @brief Canonical block-pair loop: apply `kernel(A, B, buf)` to every
///        lower-triangular (i ≥ j) block pair and scatter into `M`.
template <typename Kernel>
void for_each_canonical_block_pair(
    xt::xtensor<double, 2>& M, const std::vector<native::ShellBlock>& blocks,
    Kernel&& kernel) {
    const int nb = static_cast<int>(blocks.size());
    double buf[kBlockBuf1e];
    for (int i = 0; i < nb; ++i) {
        for (int j = 0; j <= i; ++j) {
            const auto& A = blocks[static_cast<std::size_t>(i)];
            const auto& B = blocks[static_cast<std::size_t>(j)];
            const int nab = A.n_cart * B.n_cart;
            for (int k = 0; k < nab; ++k)
                buf[k] = 0.0;
            kernel(A, B, buf);
            scatter_block(M, A, B, buf, /*diag=*/(i == j));
        }
    }
}

} // anonymous namespace

xt::xtensor<double, 2> IntegralPack::overlap(const ShellList& shells) {
    const int n = shells.n_basis;
    xt::xtensor<double, 2> S = xt::zeros<double>({n, n});
    const auto blocks = blocks_from_shells(shells.shells);
    for_each_canonical_block_pair(
        S, blocks,
        [](const native::ShellBlock& A, const native::ShellBlock& B,
           double* buf) { native::overlap_block(A, B, buf); });
    return S;
}

xt::xtensor<double, 2> IntegralPack::kinetic(const ShellList& shells) {
    const int n = shells.n_basis;
    xt::xtensor<double, 2> T = xt::zeros<double>({n, n});
    const auto blocks = blocks_from_shells(shells.shells);
    for_each_canonical_block_pair(
        T, blocks,
        [](const native::ShellBlock& A, const native::ShellBlock& B,
           double* buf) { native::kinetic_block(A, B, buf); });
    return T;
}

xt::xtensor<double, 2>
IntegralPack::nuclear_attraction(const ShellList& shells,
                                 const AtomList& atoms) {
    const int n = shells.n_basis;
    const int natm = static_cast<int>(atoms.atoms.size());
    xt::xtensor<double, 2> V = xt::zeros<double>({n, n});
    const auto blocks = blocks_from_shells(shells.shells);
    for_each_canonical_block_pair(
        V, blocks,
        [&](const native::ShellBlock& A, const native::ShellBlock& B,
            double* buf) {
            native::nuclear_block(A, B, atoms.atoms.data(), natm, buf);
        });
    return V;
}

double nuclear_repulsion_energy(const AtomList& atoms) {
    const int natm = static_cast<int>(atoms.atoms.size());
    double e = 0.0;
    for (int a = 0; a < natm; ++a) {
        for (int b = a + 1; b < natm; ++b) {
            const double dx = atoms.atoms[a].pos[0] - atoms.atoms[b].pos[0];
            const double dy = atoms.atoms[a].pos[1] - atoms.atoms[b].pos[1];
            const double dz = atoms.atoms[a].pos[2] - atoms.atoms[b].pos[2];
            const double r = std::sqrt(dx * dx + dy * dy + dz * dz);
            e += static_cast<double>(atoms.atoms[a].Z * atoms.atoms[b].Z) / r;
        }
    }
    return e;
}

} // namespace scf
} // namespace cpu
} // namespace ATV

namespace ATV {
namespace cpu {
namespace scf {

namespace {

/// 1D integral I(la, lb, e) = ∫ (x−A)^la (x−B)^lb (x−C)^e exp(−p(x−P)²) dx
/// via binomial expansion about the Gaussian product centre P and the
/// even-moment closed form ∫ x^n e^{−px²} dx = (n−1)!! / (2p)^{n/2} √(π/p).
double moment_1d(int la, int lb, int e, double pa, double pb, double pc,
                 double p) {
    auto binom = [](int n, int k) {
        double r = 1.0;
        for (int i = 1; i <= k; ++i)
            r = r * (n - k + i) / i;
        return r;
    };
    const double spi = std::sqrt(M_PI / p);
    double acc = 0.0;
    for (int i = 0; i <= la; ++i) {
        for (int j = 0; j <= lb; ++j) {
            for (int k = 0; k <= e; ++k) {
                const int n = i + j + k;
                if (n % 2)
                    continue;
                double dfac = 1.0; // (n−1)!!
                for (int m = n - 1; m > 1; m -= 2)
                    dfac *= m;
                acc += binom(la, i) * binom(lb, j) * binom(e, k)
                       * std::pow(pa, la - i) * std::pow(pb, lb - j)
                       * std::pow(pc, e - k) * dfac
                       / std::pow(2.0 * p, n / 2.0) * spi;
            }
        }
    }
    return acc;
}

} // namespace

MultipoleIntegrals::MultipoleIntegrals(const native::Shell* shells, int n_ao,
                                       const double* origin_bohr) {
    static const int kPow[9][3] = {// m1: x,y,z; m2: xx,xy,xz,yy,yz,zz
                                   {1, 0, 0}, {0, 1, 0}, {0, 0, 1},
                                   {2, 0, 0}, {1, 1, 0}, {1, 0, 1},
                                   {0, 2, 0}, {0, 1, 1}, {0, 0, 2}};
    xt::xtensor<double, 2>* out[9];
    for (int a = 0; a < 3; ++a) {
        m1_[a] = xt::zeros<double>({std::size_t(n_ao), std::size_t(n_ao)});
        out[a] = &m1_[a];
    }
    for (int k = 0; k < 6; ++k) {
        m2_[k] = xt::zeros<double>({std::size_t(n_ao), std::size_t(n_ao)});
        out[3 + k] = &m2_[k];
    }

    int ca[3], cb[3];
    for (int i = 0; i < n_ao; ++i) {
        const auto& A = shells[i];
        int compsA[kMaxNCart][3];
        fill_cart_comps(A.l, compsA);
        for (int c = 0; c < 3; ++c)
            ca[c] = compsA[A.cart][c];
        for (int j = 0; j <= i; ++j) {
            const auto& B = shells[j];
            int compsB[kMaxNCart][3];
            fill_cart_comps(B.l, compsB);
            for (int c = 0; c < 3; ++c)
                cb[c] = compsB[B.cart][c];

            double acc[9] = {};
            for (int pa = 0; pa < A.n_prim; ++pa) {
                for (int pb = 0; pb < B.n_prim; ++pb) {
                    const double a = A.alpha[pa];
                    const double b = B.alpha[pb];
                    const double p = a + b;
                    const double cc = A.coeff[pa] * B.coeff[pb];
                    double pre = cc;
                    double P[3], PA[3], PB[3], PC[3];
                    double r2 = 0.0;
                    for (int c = 0; c < 3; ++c) {
                        const double d = A.center[c] - B.center[c];
                        r2 += d * d;
                        P[c] = (a * A.center[c] + b * B.center[c]) / p;
                        PA[c] = P[c] - A.center[c];
                        PB[c] = P[c] - B.center[c];
                        PC[c] = P[c] - origin_bohr[c];
                    }
                    pre *= std::exp(-a * b / p * r2);
                    // Per-axis 1D moments up to e = 2 (reused across the 9
                    // operators for this primitive pair).
                    double I[3][3]; // [axis][e]
                    for (int c = 0; c < 3; ++c)
                        for (int e = 0; e < 3; ++e)
                            I[c][e] = moment_1d(ca[c], cb[c], e, PA[c],
                                                PB[c], PC[c], p);
                    for (int op = 0; op < 9; ++op) {
                        double v = pre;
                        for (int c = 0; c < 3; ++c)
                            v *= I[c][kPow[op][c]];
                        acc[op] += v;
                    }
                }
            }
            for (int op = 0; op < 9; ++op) {
                (*out[op])(i, j) = acc[op];
                (*out[op])(j, i) = acc[op];
            }
        }
    }
}

} // namespace scf
} // namespace cpu
} // namespace ATV
