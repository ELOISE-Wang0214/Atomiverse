/// @file cpu/scf/libcint_bridge.cpp
/// @brief STO-3G -> libcint atm/bas/env bridge matching PySCF normalization.

#include "atomiverse/cpu/scf/libcint_bridge.h"

#ifdef ATV_HAS_LIBCINT

extern "C" {
#include <cint.h>
}

#include "atomiverse/cpu/scf/sto3g_data.h"
#include "atomiverse/core/constants.h"
#include "atomiverse/core/gto.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <tuple>
#include <vector>

namespace ATV {
namespace cpu {
namespace scf {

namespace {

inline int cart_dim(int l) {
    return (l + 1) * (l + 2) / 2;
}

double gaussian_int(int n, double alpha) {
    const double p = 0.5 * (static_cast<double>(n) + 1.0);
    return 0.5 * std::tgamma(p) / std::pow(alpha, p);
}

std::array<double, 3> normalize_contracted_ao(
    int l,
    const std::array<double, 3>& exps,
    const std::array<double, 3>& coeffs_with_prim_norm)
{
    double norm2 = 0.0;
    for (int p = 0; p < 3; ++p)
        for (int q = 0; q < 3; ++q)
            norm2 += coeffs_with_prim_norm[p] * coeffs_with_prim_norm[q]
                  * gaussian_int(l * 2 + 2, exps[p] + exps[q]);

    const double scale = 1.0 / std::sqrt(norm2);
    return {
        coeffs_with_prim_norm[0] * scale,
        coeffs_with_prim_norm[1] * scale,
        coeffs_with_prim_norm[2] * scale,
    };
}

void append_shell(
    LibcintSto3gContext& ctx,
    int atom_index,
    int l,
    const std::array<double, 3>& alpha_in,
    const std::array<double, 3>& coeff_in)
{
    std::array<std::tuple<double, double>, 3> pairs = {{
        {alpha_in[0], coeff_in[0]},
        {alpha_in[1], coeff_in[1]},
        {alpha_in[2], coeff_in[2]},
    }};
    std::sort(pairs.begin(), pairs.end(),
              [](const auto& a, const auto& b) {
                  return std::get<0>(a) > std::get<0>(b);
              });

    std::array<double, 3> exps{};
    std::array<double, 3> coeffs{};
    for (int p = 0; p < 3; ++p) {
        exps[p] = std::get<0>(pairs[p]);
        coeffs[p] = std::get<1>(pairs[p]) * ATV::gto_norm(l, exps[p]);
    }
    coeffs = normalize_contracted_ao(l, exps, coeffs);

    const int ptr_exp = static_cast<int>(ctx.env.size());
    ctx.env.insert(ctx.env.end(), exps.begin(), exps.end());
    const int ptr_coeff = static_cast<int>(ctx.env.size());
    ctx.env.insert(ctx.env.end(), coeffs.begin(), coeffs.end());

    ctx.shell_offsets.push_back(ctx.nao);
    ctx.bas.push_back(atom_index);
    ctx.bas.push_back(l);
    ctx.bas.push_back(3);
    ctx.bas.push_back(1);
    ctx.bas.push_back(0);
    ctx.bas.push_back(ptr_exp);
    ctx.bas.push_back(ptr_coeff);
    ctx.bas.push_back(0);

    ctx.nao += cart_dim(l);
    ctx.nshells += 1;
}

} // namespace

LibcintSto3gContext build_libcint_sto3g_context(const Frame& frame) {
    LibcintSto3gContext ctx;
    ctx.natm = static_cast<int>(frame.atoms.size());
    ctx.atm.resize(ctx.natm * ATM_SLOTS, 0);
    ctx.env.resize(PTR_ENV_START, 0.0);

    for (int a = 0; a < ctx.natm; ++a) {
        const int Z = static_cast<int>(frame.atoms.type[a]);
        const int coord_ptr = static_cast<int>(ctx.env.size());
        ctx.env.push_back(static_cast<double>(frame.atoms.x[a]) * kAngstromToBohr);
        ctx.env.push_back(static_cast<double>(frame.atoms.y[a]) * kAngstromToBohr);
        ctx.env.push_back(static_cast<double>(frame.atoms.z[a]) * kAngstromToBohr);

        ctx.atm[a * ATM_SLOTS + CHARGE_OF] = Z;
        ctx.atm[a * ATM_SLOTS + PTR_COORD] = coord_ptr;
        ctx.atm[a * ATM_SLOTS + NUC_MOD_OF] = POINT_NUC;

        const auto& basis = sto3g_data[Z];
        for (int s = 0; s < basis.n_raw_shells; ++s) {
            const auto& raw = basis.raw[s];
            const std::array<double, 3> exps = {
                raw.alpha[0], raw.alpha[1], raw.alpha[2]
            };

            if (raw.type == 0) {
                append_shell(ctx, a, 0, exps, {
                    raw.coeff_s[0], raw.coeff_s[1], raw.coeff_s[2]
                });
            } else {
                append_shell(ctx, a, 0, exps, {
                    raw.coeff_s[0], raw.coeff_s[1], raw.coeff_s[2]
                });
                append_shell(ctx, a, 1, exps, {
                    raw.coeff_p[0], raw.coeff_p[1], raw.coeff_p[2]
                });
            }
        }
    }

    return ctx;
}

} // namespace scf
} // namespace cpu
} // namespace ATV

#endif
