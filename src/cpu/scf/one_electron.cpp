/// @file cpu/scf/one_electron.cpp
/// @brief Native-engine implementations of the one-electron integral builders.

#include "atomiverse/cpu/scf/one_electron.h"

#include <cmath>

#include "atomiverse/cpu/scf/cart.h"
#include "atomiverse/cpu/scf/native_integrals.h"
#include "atomiverse/cpu/scf/shell_pair.h"

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

/// @brief Canonical block-pair loop: apply `kernel(A, B, buf)` to every
///        lower-triangular (i ≥ j) block pair and scatter into `M`.
template <typename Kernel>
void for_each_canonical_block_pair(
    xt::xtensor<double, 2>& M,
    const std::vector<native::ShellBlock>& blocks,
    Kernel&& kernel)
{
    const int nb = static_cast<int>(blocks.size());
    double buf[kBlockBuf1e];
    for (int i = 0; i < nb; ++i) {
        for (int j = 0; j <= i; ++j) {
            const auto& A = blocks[static_cast<std::size_t>(i)];
            const auto& B = blocks[static_cast<std::size_t>(j)];
            const int nab = A.n_cart * B.n_cart;
            for (int k = 0; k < nab; ++k) buf[k] = 0.0;
            kernel(A, B, buf);
            scatter_block(M, A, B, buf, /*diag=*/(i == j));
        }
    }
}

} // anonymous namespace

xt::xtensor<double, 2> build_overlap(const ShellList& shells) {
    const int n = shells.n_basis;
    xt::xtensor<double, 2> S = xt::zeros<double>({n, n});
    const auto blocks = blocks_from_shells(shells.shells);
    for_each_canonical_block_pair(S, blocks,
        [](const native::ShellBlock& A, const native::ShellBlock& B, double* buf) {
            native::overlap_block(A, B, buf);
        });
    return S;
}

xt::xtensor<double, 2> build_kinetic(const ShellList& shells) {
    const int n = shells.n_basis;
    xt::xtensor<double, 2> T = xt::zeros<double>({n, n});
    const auto blocks = blocks_from_shells(shells.shells);
    for_each_canonical_block_pair(T, blocks,
        [](const native::ShellBlock& A, const native::ShellBlock& B, double* buf) {
            native::kinetic_block(A, B, buf);
        });
    return T;
}

xt::xtensor<double, 2> build_nuclear_attraction(
    const ShellList& shells, const AtomList& atoms)
{
    const int n = shells.n_basis;
    const int natm = static_cast<int>(atoms.atoms.size());
    xt::xtensor<double, 2> V = xt::zeros<double>({n, n});
    const auto blocks = blocks_from_shells(shells.shells);
    for_each_canonical_block_pair(V, blocks,
        [&](const native::ShellBlock& A, const native::ShellBlock& B, double* buf) {
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
