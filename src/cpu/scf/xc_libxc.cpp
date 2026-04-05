/// @file cpu/scf/xc_libxc.cpp
/// @brief libxc-based XC evaluator implementation.
///
/// Wraps libxc's C API to produce XCOutput matching the native evaluate_xc() interface.
/// Used as ground truth: native XC must match this to 1e-12.

#ifdef ATV_HAS_LIBXC

#include "atomiverse/cpu/scf/xc_libxc.h"

#include <xc.h>
#include <stdexcept>
#include <string>
#include <cstring>

namespace ATV {
namespace cpu {
namespace scf {

/// @brief Evaluate a single LDA functional (exchange or correlation) via libxc.
///
/// @param func_id  libxc functional ID (e.g. XC_LDA_X, XC_LDA_C_VWN).
/// @param rho_a    Alpha density.
/// @param rho_b    Beta density.
/// @param[out] exc   Energy per electron (added to, not overwritten).
/// @param[out] vrho  Potential derivatives [vrho_alpha, vrho_beta] (added to).
static void eval_lda(int func_id, double rho_a, double rho_b,
                     double& exc, double* vrho) {
    xc_func_type func;
    if (xc_func_init(&func, func_id, XC_POLARIZED) != 0)
        throw std::runtime_error("libxc: failed to init functional " + std::to_string(func_id));

    double rho[2] = {rho_a, rho_b};
    double zk = 0;
    double vr[2] = {0, 0};

    xc_lda_exc_vxc(&func, 1, rho, &zk, vr);
    xc_func_end(&func);

    exc += zk;
    vrho[0] += vr[0];
    vrho[1] += vr[1];
}

/// @brief Evaluate a single GGA functional (exchange or correlation) via libxc.
static void eval_gga(int func_id, double rho_a, double rho_b,
                     double sigma_aa, double sigma_ab, double sigma_bb,
                     double& exc, double* vrho, double* vsigma) {
    xc_func_type func;
    if (xc_func_init(&func, func_id, XC_POLARIZED) != 0)
        throw std::runtime_error("libxc: failed to init functional " + std::to_string(func_id));

    double rho[2] = {rho_a, rho_b};
    double sigma[3] = {sigma_aa, sigma_ab, sigma_bb};
    double zk = 0;
    double vr[2] = {0, 0};
    double vs[3] = {0, 0, 0};

    xc_gga_exc_vxc(&func, 1, rho, sigma, &zk, vr, vs);
    xc_func_end(&func);

    exc += zk;
    vrho[0] += vr[0];
    vrho[1] += vr[1];
    vsigma[0] += vs[0];
    vsigma[1] += vs[1];
    vsigma[2] += vs[2];
}

XCOutput evaluate_xc_libxc(XCType type, Real rho_a, Real rho_b,
                           Real sigma_aa, Real sigma_ab, Real sigma_bb) {
    XCOutput out{};

    double ra = static_cast<double>(rho_a);
    double rb = static_cast<double>(rho_b);
    double rho_total = ra + rb;

    if (rho_total < 1e-12) return out;

    double exc = 0;
    double vrho[2] = {0, 0};
    double vsigma[3] = {0, 0, 0};

    switch (type) {
    case XCType::LSDA:
        // Slater exchange + VWN5 correlation
        // XC_LDA_C_VWN = 7 (VWN5, parameterization V — NOT VWN3=id 5)
        eval_lda(XC_LDA_X, ra, rb, exc, vrho);
        eval_lda(XC_LDA_C_VWN, ra, rb, exc, vrho);
        break;

    case XCType::PBE:
        // PBE exchange is universal
        eval_gga(XC_GGA_X_PBE, ra, rb,
                 static_cast<double>(sigma_aa),
                 static_cast<double>(sigma_ab),
                 static_cast<double>(sigma_bb),
                 exc, vrho, vsigma);
        // PBE correlation: use VWN5-based variant (id 216) to match native code
        // which uses VWN5 as LDA base. Standard XC_GGA_C_PBE uses PW92.
        eval_gga(XC_GGA_C_PBE_VWN, ra, rb,
                 static_cast<double>(sigma_aa),
                 static_cast<double>(sigma_ab),
                 static_cast<double>(sigma_bb),
                 exc, vrho, vsigma);
        break;
    }

    out.exc = static_cast<Real>(exc);
    out.vrho_alpha = static_cast<Real>(vrho[0]);
    out.vrho_beta  = static_cast<Real>(vrho[1]);
    out.vsigma_aa  = static_cast<Real>(vsigma[0]);
    out.vsigma_ab  = static_cast<Real>(vsigma[1]);
    out.vsigma_bb  = static_cast<Real>(vsigma[2]);

    return out;
}

} // namespace scf
} // namespace cpu
} // namespace ATV

#endif // ATV_HAS_LIBXC
