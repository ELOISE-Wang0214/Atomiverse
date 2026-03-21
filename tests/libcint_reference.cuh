/// @file tests/libcint_reference.cuh
/// @brief libcint ground-truth reference for integral validation tests.
///
/// Isolates ALL libcint C-API interaction into a single test-only header.
/// Gated by ATV_HAS_LIBCINT — only available when -DATV_BUILD_LIBCINT=ON.
#pragma once

#ifdef ATV_HAS_LIBCINT

extern "C" {
#include <cint.h>
#include <cint_funcs.h>
}
#undef atm
#undef bas

#include "atomiverse/cpu/scf/integrals.h"
#include "atomiverse/core/constants.h"

#include <cmath>
#include <vector>

namespace ATV {
namespace test {

/// @brief Wrap CINTgto_norm for test comparison.
inline double libcint_gto_norm(int l, double alpha) {
    return CINTgto_norm(l, alpha);
}

/// @brief Compute STO-3G integrals via libcint (ground truth).
///
/// Mirrors the logic that was previously in integrals.cpp's libcint path.
/// Used exclusively for validation — never in production code.
inline cpu::scf::IntegralPack build_integrals_libcint(const cpu::Frame& frame) {
    int natm = static_cast<int>(frame.n_atoms());
    Int n_electrons = 0;
    for (Int i = 0; i < natm; ++i)
        n_electrons += frame.type_id[i];

    // STO-3G basis data (same tables as integrals.cpp)
    struct STO3GShell { int l, n_prim; double alpha[3], coeff[3]; };
    struct STO3GBasis { int n_shells; STO3GShell shells[2]; };
    static const STO3GBasis sto3g_basis[] = {
        {0, {}},
        {1, {{0, 3, {3.42525091, 0.62391373, 0.16885540}, {0.15432897, 0.53532814, 0.44463454}}}},
        {1, {{0, 3, {6.36242139, 1.15892300, 0.31364979}, {0.15432897, 0.53532814, 0.44463454}}}},
        {2, {
            {0, 3, {16.1195750, 2.93620070, 0.79465370}, { 0.15432897,  0.53532814,  0.44463454}},
            {0, 3, { 0.63628970, 0.14786010, 0.04808870}, {-0.09996723,  0.39951283,  0.70011547}},
        }},
    };

    Int n_basis = 0;
    for (Int i = 0; i < natm; ++i)
        n_basis += sto3g_basis[static_cast<int>(frame.type_id[i])].n_shells;

    cpu::scf::IntegralPack pack;
    pack.n_basis = n_basis;
    pack.n_electrons = n_electrons;
    pack.S = xt::zeros<Real>({n_basis, n_basis});
    pack.T = xt::zeros<Real>({n_basis, n_basis});
    pack.V = xt::zeros<Real>({n_basis, n_basis});
    pack.ERI = xt::zeros<Real>({n_basis, n_basis, n_basis, n_basis});

    int n_shells = static_cast<int>(n_basis);
    int off = PTR_ENV_START;

    std::vector<FINT> atm_arr(natm * ATM_SLOTS, 0);
    std::vector<FINT> bas_arr(n_shells * BAS_SLOTS, 0);
    std::vector<double> env(off + natm * 3 + n_shells * 6, 0.0);

    for (int i = 0; i < natm; ++i) {
        atm_arr[CHARGE_OF + ATM_SLOTS * i] = static_cast<FINT>(frame.type_id[i]);
        atm_arr[PTR_COORD + ATM_SLOTS * i] = off;
        env[off + 0] = static_cast<double>(frame.x[i]) * kAngstromToBohr;
        env[off + 1] = static_cast<double>(frame.y[i]) * kAngstromToBohr;
        env[off + 2] = static_cast<double>(frame.z[i]) * kAngstromToBohr;
        off += 3;
    }

    int shell_idx = 0;
    for (int a = 0; a < natm; ++a) {
        int Z = static_cast<int>(frame.type_id[a]);
        const auto& basis = sto3g_basis[Z];
        for (int s = 0; s < basis.n_shells; ++s) {
            const auto& sh = basis.shells[s];
            int ptr_exp = off;
            for (int p = 0; p < sh.n_prim; ++p) env[off++] = sh.alpha[p];
            int ptr_coeff = off;
            for (int p = 0; p < sh.n_prim; ++p)
                env[off++] = sh.coeff[p] * CINTgto_norm(sh.l, sh.alpha[p]);
            bas_arr[ATOM_OF   + BAS_SLOTS * shell_idx] = a;
            bas_arr[ANG_OF    + BAS_SLOTS * shell_idx] = sh.l;
            bas_arr[NPRIM_OF  + BAS_SLOTS * shell_idx] = sh.n_prim;
            bas_arr[NCTR_OF   + BAS_SLOTS * shell_idx] = 1;
            bas_arr[PTR_EXP   + BAS_SLOTS * shell_idx] = ptr_exp;
            bas_arr[PTR_COEFF + BAS_SLOTS * shell_idx] = ptr_coeff;
            ++shell_idx;
        }
    }
    env.resize(off);

    FINT* atm_p = atm_arr.data();
    FINT* bas_p = bas_arr.data();
    double* env_p = env.data();

    FINT shls[2];
    for (int i = 0; i < n_shells; ++i)
        for (int j = 0; j <= i; ++j) {
            shls[0] = i; shls[1] = j;
            double buf[1];
            int1e_ovlp_sph(buf, nullptr, shls, atm_p, natm, bas_p, n_shells, env_p, nullptr, nullptr);
            pack.S(i, j) = static_cast<Real>(buf[0]); pack.S(j, i) = pack.S(i, j);
            int1e_kin_sph(buf, nullptr, shls, atm_p, natm, bas_p, n_shells, env_p, nullptr, nullptr);
            pack.T(i, j) = static_cast<Real>(buf[0]); pack.T(j, i) = pack.T(i, j);
            int1e_nuc_sph(buf, nullptr, shls, atm_p, natm, bas_p, n_shells, env_p, nullptr, nullptr);
            pack.V(i, j) = static_cast<Real>(buf[0]); pack.V(j, i) = pack.V(i, j);
        }

    FINT shls4[4];
    for (int i = 0; i < n_shells; ++i)
      for (int j = 0; j < n_shells; ++j)
        for (int k = 0; k < n_shells; ++k)
          for (int l = 0; l < n_shells; ++l) {
            shls4[0] = i; shls4[1] = j; shls4[2] = k; shls4[3] = l;
            double buf[1];
            cint2e_sph(buf, shls4, atm_p, natm, bas_p, n_shells, env_p, nullptr);
            pack.ERI(i, j, k, l) = static_cast<Real>(buf[0]);
          }

    pack.Hcore = pack.T + pack.V;

    pack.V_nuc = 0;
    for (Int a = 0; a < natm; ++a)
        for (Int b = a + 1; b < natm; ++b) {
            double dx = (static_cast<double>(frame.x[a]) - static_cast<double>(frame.x[b])) * kAngstromToBohr;
            double dy = (static_cast<double>(frame.y[a]) - static_cast<double>(frame.y[b])) * kAngstromToBohr;
            double dz = (static_cast<double>(frame.z[a]) - static_cast<double>(frame.z[b])) * kAngstromToBohr;
            double r = std::sqrt(dx*dx + dy*dy + dz*dz);
            pack.V_nuc += static_cast<Real>(frame.type_id[a] * frame.type_id[b] / r);
        }

    return pack;
}

} // namespace test
} // namespace ATV

#endif // ATV_HAS_LIBCINT
