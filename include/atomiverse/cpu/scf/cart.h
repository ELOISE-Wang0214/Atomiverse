/// @file cpu/scf/cart.h
/// @brief Cartesian-Gaussian enumeration helpers.
///
/// Every shell-level integral kernel needs the same two facts:
///   * how many Cartesian AOs a shell of angular momentum `l` carries
///     (1 for s, 3 for p, 6 for d, 10 for f, ...),
///   * the per-component `(l_x, l_y, l_z)` tuple for each AO.
///
/// Before this header, four separate copies of that logic coexisted
/// (ang_table / cart_exponents / fill_cart_comps / fill_cart_block_1e),
/// each with subtly different signatures and several written for one
/// specific l value. Consolidating them here keeps callers honest:
/// extending to higher l = editing one function.
#pragma once

namespace ATV {
namespace cpu {
namespace scf {

/// @brief Maximum angular momentum currently supported by the CPU
///        integral engines (Rys ERI, MD 1e/2e). l = 3 → f.
inline constexpr int kMaxL = 3;

/// @brief Number of Cartesian AOs carried by a shell of angular
///        momentum `l`: `(l+1)(l+2)/2`.
///
/// `n_cart_of(0) = 1`, `n_cart_of(1) = 3`, `n_cart_of(2) = 6`,
/// `n_cart_of(3) = 10`.
[[nodiscard]] inline constexpr int n_cart_of(int l) {
    return (l + 1) * (l + 2) / 2;
}

/// @brief Maximum Cartesian AO count at the largest supported l.
inline constexpr int kMaxNCart = n_cart_of(kMaxL);

/// @brief Fill the `(l_x, l_y, l_z)` tuple for each Cartesian AO of
///        a shell with angular momentum `l`.
///
/// Enumeration order matches libcint / PySCF:
///   `nx` decreases from `l` to 0; for each `nx`, `ny` decreases from
///   `l - nx` to 0; `nz = l - nx - ny`.
///
/// For `l = 2`: xx (2,0,0), xy (1,1,0), xz (1,0,1), yy (0,2,0),
///              yz (0,1,1), zz (0,0,2).
/// For `l = 3` the same rule gives the 10 f components in order
///              xxx, xxy, xxz, xyy, xyz, xzz, yyy, yyz, yzz, zzz.
///
/// The caller must supply `out` with at least `n_cart_of(l)` rows.
inline void fill_cart_comps(int l, int out[][3]) {
    int c = 0;
    for (int nx = l; nx >= 0; --nx) {
        for (int ny = l - nx; ny >= 0; --ny) {
            out[c][0] = nx;
            out[c][1] = ny;
            out[c][2] = l - nx - ny;
            ++c;
        }
    }
}

} // namespace scf
} // namespace cpu
} // namespace ATV
