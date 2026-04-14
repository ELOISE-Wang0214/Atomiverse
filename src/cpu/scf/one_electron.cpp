/// @file cpu/scf/one_electron.cpp
/// @brief Native-engine implementations of the one-electron integral builders.

#include "atomiverse/cpu/scf/one_electron.h"

#include <cmath>
#include <vector>

#include "atomiverse/cpu/scf/native_integrals.h"

namespace ATV {
namespace cpu {
namespace scf {

namespace {

constexpr int kMaxCart1eSq = 36;  // (dd) block = 6 * 6

/// @brief Group shells in a ShellList into libcint-style l-blocks.
///
/// Consecutive native::Shells sharing the same l (and same primitive data,
/// by construction of build_shell_list) collapse into one ShellBlock with
/// n_cart = 2l+1-for-spherical ... no, we use Cartesian: 1/3/6 for s/p/d.
/// Each block carries its global `ao_start` so callers can scatter.
std::vector<native::ShellBlock> blocks_from_shells(const ShellList& list) {
    std::vector<native::ShellBlock> blocks;
    blocks.reserve(static_cast<std::size_t>(list.n_basis));

    int shell_idx = 0, ao_idx = 0;
    while (shell_idx < list.n_basis) {
        const auto& sh = list.shells[static_cast<std::size_t>(shell_idx)];
        const int n_cart = (sh.l == 0) ? 1 : (sh.l == 1) ? 3 : 6;

        native::ShellBlock blk;
        blk.n_prim    = sh.n_prim;
        blk.alpha     = sh.alpha;
        blk.coeff     = sh.coeff;
        blk.center[0] = sh.center[0];
        blk.center[1] = sh.center[1];
        blk.center[2] = sh.center[2];
        blk.l         = sh.l;
        blk.n_cart    = n_cart;
        blk.ao_start  = ao_idx;
        blocks.push_back(blk);

        shell_idx += n_cart;
        ao_idx    += n_cart;
    }
    return blocks;
}

/// @brief Scatter a block's upper-triangular result into a symmetric matrix.
///
/// For A != B (block-off-diagonal): fill both M[A,B] and M[B,A] from `buf`.
/// For A == B (block-diagonal): fill all n_cart×n_cart entries; the native
/// kernel produces a symmetric block directly, so copying both halves is safe.
inline void scatter_block(
    xt::xtensor<double, 2>& M,
    const native::ShellBlock& A, const native::ShellBlock& B,
    const double* buf, bool diag)
{
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

} // anonymous namespace

xt::xtensor<double, 2> build_overlap(const ShellList& shells) {
    const int n = shells.n_basis;
    xt::xtensor<double, 2> S = xt::zeros<double>({n, n});

    const auto blocks = blocks_from_shells(shells);
    const int nb = static_cast<int>(blocks.size());

    double buf[kMaxCart1eSq];  // max 6*6 = 36 for (dd) block

    for (int i = 0; i < nb; ++i) {
        for (int j = 0; j <= i; ++j) {
            const auto& A = blocks[static_cast<std::size_t>(i)];
            const auto& B = blocks[static_cast<std::size_t>(j)];
            const int nab = A.n_cart * B.n_cart;
            for (int k = 0; k < nab; ++k) buf[k] = 0.0;
            native::overlap_block(A, B, buf);
            scatter_block(S, A, B, buf, /*diag=*/(i == j));
        }
    }
    return S;
}

xt::xtensor<double, 2> build_kinetic(const ShellList& shells) {
    const int n = shells.n_basis;
    xt::xtensor<double, 2> T = xt::zeros<double>({n, n});

    const auto blocks = blocks_from_shells(shells);
    const int nb = static_cast<int>(blocks.size());

    double buf[kMaxCart1eSq];

    for (int i = 0; i < nb; ++i) {
        for (int j = 0; j <= i; ++j) {
            const auto& A = blocks[static_cast<std::size_t>(i)];
            const auto& B = blocks[static_cast<std::size_t>(j)];
            const int nab = A.n_cart * B.n_cart;
            for (int k = 0; k < nab; ++k) buf[k] = 0.0;
            native::kinetic_block(A, B, buf);
            scatter_block(T, A, B, buf, /*diag=*/(i == j));
        }
    }
    return T;
}

xt::xtensor<double, 2> build_nuclear_attraction(
    const ShellList& shells, const AtomList& atoms)
{
    const int n = shells.n_basis;
    const int natm = static_cast<int>(atoms.atoms.size());
    xt::xtensor<double, 2> V = xt::zeros<double>({n, n});

    const auto blocks = blocks_from_shells(shells);
    const int nb = static_cast<int>(blocks.size());

    double buf[kMaxCart1eSq];

    for (int i = 0; i < nb; ++i) {
        for (int j = 0; j <= i; ++j) {
            const auto& A = blocks[static_cast<std::size_t>(i)];
            const auto& B = blocks[static_cast<std::size_t>(j)];
            const int nab = A.n_cart * B.n_cart;
            for (int k = 0; k < nab; ++k) buf[k] = 0.0;
            native::nuclear_block(A, B, atoms.atoms.data(), natm, buf);
            scatter_block(V, A, B, buf, /*diag=*/(i == j));
        }
    }
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
