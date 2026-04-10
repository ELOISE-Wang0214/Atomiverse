/// @file tests/libcint_reference.cuh
/// @brief libcint ground-truth reference for integral validation tests.
///
/// Isolates ALL libcint C-API interaction into a single test-only header.
/// Gated by ATV_HAS_LIBCINT — only available when -DATV_BUILD_LIBCINT=ON.
///
/// Supports s+p type shells for elements Z=1..8 (STO-3G basis).
#pragma once

#ifdef ATV_HAS_LIBCINT

extern "C" {
#include <cint.h>
#include <cint_funcs.h>
}
#undef atm
#undef bas

#include "atomiverse/cpu/scf/integrals.h"
#include "atomiverse/cpu/scf/sto3g_data.h"
#include "atomiverse/core/constants.h"

#include <cmath>
#include <vector>
#include <numeric>

namespace ATV {
namespace test {

/// @brief Wrap CINTgto_norm for test comparison.
inline double libcint_gto_norm(int l, double alpha) {
    return CINTgto_norm(l, alpha);
}

/// @brief Normalize contracted AO coefficients (PySCF convention).
/// After primitive normalization (coeff *= gto_norm), divide by the
/// contracted orbital norm so that <φ|φ> = 1. This matches the
/// normalization applied by libcint_bridge.cpp and PySCF's make_env.
inline void normalize_contracted_coeffs(int l, const double* exps,
                                        double* coeffs, int nprim) {
    double norm2 = 0.0;
    for (int p = 0; p < nprim; ++p)
        for (int q = 0; q < nprim; ++q) {
            double n1 = 0.5 * (l * 2 + 2 + 1);
            norm2 += coeffs[p] * coeffs[q]
                   * 0.5 * std::tgamma(n1) / std::pow(exps[p] + exps[q], n1);
        }
    double scale = 1.0 / std::sqrt(norm2);
    for (int p = 0; p < nprim; ++p)
        coeffs[p] *= scale;
}

/// @brief Compute STO-3G integrals via libcint (ground truth).
///
/// Mirrors the logic in build_integrals_sto3g but uses libcint as the
/// integral engine. Supports elements Z=1..8 with s and p type shells.
/// Used exclusively for validation — never in production code.
inline cpu::scf::IntegralPack build_integrals_libcint(const cpu::Frame& frame) {
    int natm = static_cast<int>(frame.atoms.size());
    int n_electrons = 0;
    for (int i = 0; i < natm; ++i)
        n_electrons += frame.atoms.type[i];

    // Count basis functions and shells using shared STO-3G data
    int n_basis = 0;
    int n_shells = 0;
    for (int i = 0; i < natm; ++i) {
        int Z = static_cast<int>(frame.atoms.type[i]);
        n_basis += cpu::scf::sto3g_n_basis(Z);
        const auto& basis = cpu::scf::sto3g_data[Z];
        for (int s = 0; s < basis.n_raw_shells; ++s) {
            if (basis.raw[s].type == 0)
                n_shells += 1; // pure s -> 1 shell
            else
                n_shells += 2; // sp -> 1 s-shell + 1 p-shell (libcint treats p as one shell with 3 components)
        }
    }

    cpu::scf::IntegralPack pack;
    pack.n_basis = n_basis;
    pack.n_electrons = n_electrons;
    pack.S = xt::zeros<double>({static_cast<int>(n_basis), static_cast<int>(n_basis)});
    pack.T = xt::zeros<double>({static_cast<int>(n_basis), static_cast<int>(n_basis)});
    pack.V = xt::zeros<double>({static_cast<int>(n_basis), static_cast<int>(n_basis)});
    pack.ERI = xt::zeros<double>({static_cast<int>(n_basis), static_cast<int>(n_basis),
                                static_cast<int>(n_basis), static_cast<int>(n_basis)});

    // Build libcint atm/bas/env arrays
    int off = PTR_ENV_START;
    std::vector<FINT> atm_arr(natm * ATM_SLOTS, 0);
    // Generous env allocation
    std::vector<double> env(off + natm * 3 + n_shells * 10, 0.0);

    for (int i = 0; i < natm; ++i) {
        atm_arr[CHARGE_OF + ATM_SLOTS * i] = static_cast<FINT>(frame.atoms.type[i]);
        atm_arr[PTR_COORD + ATM_SLOTS * i] = off;
        env[off + 0] = static_cast<double>(frame.atoms.x[i]) * kAngstromToBohr;
        env[off + 1] = static_cast<double>(frame.atoms.y[i]) * kAngstromToBohr;
        env[off + 2] = static_cast<double>(frame.atoms.z[i]) * kAngstromToBohr;
        off += 3;
    }

    // Build basis arrays: expand sp shells into separate s and p libcint shells
    std::vector<FINT> bas_arr(n_shells * BAS_SLOTS, 0);
    // Track which basis functions each libcint shell maps to
    std::vector<int> shell_bf_offset(n_shells); // first basis function index
    std::vector<int> shell_bf_count(n_shells);  // number of basis functions

    int shell_idx = 0;
    int bf_idx = 0;
    for (int a = 0; a < natm; ++a) {
        int Z = static_cast<int>(frame.atoms.type[a]);
        const auto& basis = cpu::scf::sto3g_data[Z];
        for (int s = 0; s < basis.n_raw_shells; ++s) {
            const auto& raw = basis.raw[s];

            if (raw.type == 0) {
                // Pure s shell
                int ptr_exp = off;
                for (int p = 0; p < raw.n_prim; ++p) env[off++] = raw.alpha[p];
                int ptr_coeff = off;
                for (int p = 0; p < raw.n_prim; ++p)
                    env[off++] = raw.coeff_s[p] * CINTgto_norm(0, raw.alpha[p]);
                normalize_contracted_coeffs(0, &env[ptr_exp], &env[ptr_coeff], raw.n_prim);

                bas_arr[ATOM_OF   + BAS_SLOTS * shell_idx] = a;
                bas_arr[ANG_OF    + BAS_SLOTS * shell_idx] = 0;
                bas_arr[NPRIM_OF  + BAS_SLOTS * shell_idx] = raw.n_prim;
                bas_arr[NCTR_OF   + BAS_SLOTS * shell_idx] = 1;
                bas_arr[PTR_EXP   + BAS_SLOTS * shell_idx] = ptr_exp;
                bas_arr[PTR_COEFF + BAS_SLOTS * shell_idx] = ptr_coeff;
                shell_bf_offset[shell_idx] = bf_idx;
                shell_bf_count[shell_idx] = 1;
                ++shell_idx;
                bf_idx += 1;
            } else {
                // SP shell: split into s-part and p-part for libcint
                // s-part
                int ptr_exp_s = off;
                for (int p = 0; p < raw.n_prim; ++p) env[off++] = raw.alpha[p];
                int ptr_coeff_s = off;
                for (int p = 0; p < raw.n_prim; ++p)
                    env[off++] = raw.coeff_s[p] * CINTgto_norm(0, raw.alpha[p]);
                normalize_contracted_coeffs(0, &env[ptr_exp_s], &env[ptr_coeff_s], raw.n_prim);

                bas_arr[ATOM_OF   + BAS_SLOTS * shell_idx] = a;
                bas_arr[ANG_OF    + BAS_SLOTS * shell_idx] = 0;
                bas_arr[NPRIM_OF  + BAS_SLOTS * shell_idx] = raw.n_prim;
                bas_arr[NCTR_OF   + BAS_SLOTS * shell_idx] = 1;
                bas_arr[PTR_EXP   + BAS_SLOTS * shell_idx] = ptr_exp_s;
                bas_arr[PTR_COEFF + BAS_SLOTS * shell_idx] = ptr_coeff_s;
                shell_bf_offset[shell_idx] = bf_idx;
                shell_bf_count[shell_idx] = 1;
                ++shell_idx;
                bf_idx += 1;

                // p-part (3 components: px, py, pz)
                int ptr_exp_p = off;
                for (int p = 0; p < raw.n_prim; ++p) env[off++] = raw.alpha[p];
                int ptr_coeff_p = off;
                for (int p = 0; p < raw.n_prim; ++p)
                    env[off++] = raw.coeff_p[p] * CINTgto_norm(1, raw.alpha[p]);
                normalize_contracted_coeffs(1, &env[ptr_exp_p], &env[ptr_coeff_p], raw.n_prim);

                bas_arr[ATOM_OF   + BAS_SLOTS * shell_idx] = a;
                bas_arr[ANG_OF    + BAS_SLOTS * shell_idx] = 1;
                bas_arr[NPRIM_OF  + BAS_SLOTS * shell_idx] = raw.n_prim;
                bas_arr[NCTR_OF   + BAS_SLOTS * shell_idx] = 1;
                bas_arr[PTR_EXP   + BAS_SLOTS * shell_idx] = ptr_exp_p;
                bas_arr[PTR_COEFF + BAS_SLOTS * shell_idx] = ptr_coeff_p;
                shell_bf_offset[shell_idx] = bf_idx;
                shell_bf_count[shell_idx] = 3; // px, py, pz
                ++shell_idx;
                bf_idx += 3;
            }
        }
    }
    env.resize(off);

    FINT* atm_p = atm_arr.data();
    FINT* bas_p = bas_arr.data();
    double* env_p = env.data();

    // Use Cartesian integrals — for l<=1, Cartesian components are (x,y,z)
    // matching our convention (cart=0,1,2) directly. No reordering needed.

    // 1e integrals
    for (int i = 0; i < n_shells; ++i) {
        int ni = shell_bf_count[i];
        int oi = shell_bf_offset[i];
        for (int j = 0; j <= i; ++j) {
            int nj = shell_bf_count[j];
            int oj = shell_bf_offset[j];
            int buf_size = ni * nj;
            std::vector<double> buf(buf_size, 0.0);
            FINT shls[2] = {static_cast<FINT>(i), static_cast<FINT>(j)};

            // libcint Cartesian output: column-major buf[i + di * j]
            int1e_ovlp_cart(buf.data(), nullptr, shls, atm_p, natm, bas_p, n_shells, env_p, nullptr, nullptr);
            for (int ii = 0; ii < ni; ++ii)
                for (int jj = 0; jj < nj; ++jj) {
                    pack.S(oi+ii, oj+jj) = static_cast<double>(buf[ii + ni * jj]);
                    pack.S(oj+jj, oi+ii) = pack.S(oi+ii, oj+jj);
                }

            std::fill(buf.begin(), buf.end(), 0.0);
            int1e_kin_cart(buf.data(), nullptr, shls, atm_p, natm, bas_p, n_shells, env_p, nullptr, nullptr);
            for (int ii = 0; ii < ni; ++ii)
                for (int jj = 0; jj < nj; ++jj) {
                    pack.T(oi+ii, oj+jj) = static_cast<double>(buf[ii + ni * jj]);
                    pack.T(oj+jj, oi+ii) = pack.T(oi+ii, oj+jj);
                }

            std::fill(buf.begin(), buf.end(), 0.0);
            int1e_nuc_cart(buf.data(), nullptr, shls, atm_p, natm, bas_p, n_shells, env_p, nullptr, nullptr);
            for (int ii = 0; ii < ni; ++ii)
                for (int jj = 0; jj < nj; ++jj) {
                    pack.V(oi+ii, oj+jj) = static_cast<double>(buf[ii + ni * jj]);
                    pack.V(oj+jj, oi+ii) = pack.V(oi+ii, oj+jj);
                }
        }
    }

    // 2e integrals
    for (int i = 0; i < n_shells; ++i) {
        int ni = shell_bf_count[i]; int oi = shell_bf_offset[i];
        for (int j = 0; j < n_shells; ++j) {
            int nj = shell_bf_count[j]; int oj = shell_bf_offset[j];
            for (int k = 0; k < n_shells; ++k) {
                int nk = shell_bf_count[k]; int ok = shell_bf_offset[k];
                for (int l = 0; l < n_shells; ++l) {
                    int nl = shell_bf_count[l]; int ol = shell_bf_offset[l];
                    int buf_size = ni * nj * nk * nl;
                    std::vector<double> buf(buf_size, 0.0);
                    FINT shls4[4] = {static_cast<FINT>(i), static_cast<FINT>(j),
                                     static_cast<FINT>(k), static_cast<FINT>(l)};
                    cint2e_cart(buf.data(), shls4, atm_p, natm, bas_p, n_shells, env_p, nullptr);
                    // Column-major: i fastest
                    for (int ii = 0; ii < ni; ++ii)
                      for (int jj = 0; jj < nj; ++jj)
                        for (int kk = 0; kk < nk; ++kk)
                          for (int ll = 0; ll < nl; ++ll)
                            pack.ERI(oi+ii, oj+jj, ok+kk, ol+ll) = static_cast<double>(
                                buf[ii + ni * (jj + nj * (kk + nk * ll))]);
                }
            }
        }
    }

    pack.Hcore = pack.T + pack.V;

    // Nuclear repulsion
    pack.V_nuc = 0;
    for (int a = 0; a < natm; ++a)
        for (int b = a + 1; b < natm; ++b) {
            double dx = (static_cast<double>(frame.atoms.x[a]) - static_cast<double>(frame.atoms.x[b])) * kAngstromToBohr;
            double dy = (static_cast<double>(frame.atoms.y[a]) - static_cast<double>(frame.atoms.y[b])) * kAngstromToBohr;
            double dz = (static_cast<double>(frame.atoms.z[a]) - static_cast<double>(frame.atoms.z[b])) * kAngstromToBohr;
            double r = std::sqrt(dx*dx + dy*dy + dz*dz);
            pack.V_nuc += static_cast<double>(frame.atoms.type[a] * frame.atoms.type[b] / r);
        }

    return pack;
}

} // namespace test
} // namespace ATV

#endif // ATV_HAS_LIBCINT
