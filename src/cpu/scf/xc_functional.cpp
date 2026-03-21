/// @file cpu/scf/xc_functional.cpp
/// @brief LSDA (Slater + VWN5) and PBE GGA exchange-correlation functionals.
///
/// All formulas use atomic units (Hartree, Bohr).
///
/// LSDA exchange: Dirac-Slater, exact for uniform electron gas.
/// LSDA correlation: VWN5 (Vosko, Wilk, Nusair, 1980, parameterization V).
///   Must use VWN5 (not VWN3). VWN3 uses RPA data with known pathologies.
///
/// PBE exchange: Perdew-Burke-Ernzerhof enhancement over LSDA.
/// PBE correlation: gradient correction H(r_s, zeta, t) to LSDA correlation.

#include "atomiverse/cpu/scf/xc_functional.h"

#include <cmath>
#include <algorithm>

namespace ATV {
namespace cpu {
namespace scf {

// ── Constants ────────────────────────────────────────────────────────────

static constexpr double kPi = 3.14159265358979323846;
static constexpr double kRhoMin = 1e-12;   // density cutoff
static constexpr double kSigmaMin = 1e-20; // gradient cutoff for GGA

// Slater exchange constant: C_x = (3/4)(3/π)^{1/3}
static constexpr double kCx = 0.7385587663820224; // (3.0/4.0) * std::cbrt(3.0/kPi)

// ── VWN5 parameters (Table 4.4, Vosko et al. 1980) ──────────────────────

struct VWN5Params {
    double A, x0, b, c;
};

// Paramagnetic (ζ=0)
static constexpr VWN5Params kVWN5_P = {0.0621814, -0.10498, 3.72744, 12.9352};
// Ferromagnetic (ζ=1)
static constexpr VWN5Params kVWN5_F = {0.0310907, -0.32500, 7.06042, 18.0578};
// Spin stiffness α_c
static constexpr VWN5Params kVWN5_A = {
    -1.0 / (6.0 * kPi * kPi), -0.0047584, 1.13107, 13.0045
};

// ── VWN5 Pade approximant F(x) ──────────────────────────────────────────

// F(x) = A { ln(x²/X(x)) + 2b/Q arctan(Q/(2x+b))
//          - bx₀/X(x₀) [ln((x-x₀)²/X(x)) + 2(b+2x₀)/Q arctan(Q/(2x+b))] }
// where X(x) = x² + bx + c,  Q = √(4c - b²)
static double vwn5_F(double x, const VWN5Params& p) {
    double X = x * x + p.b * x + p.c;
    double X0 = p.x0 * p.x0 + p.b * p.x0 + p.c;
    double Q = std::sqrt(4.0 * p.c - p.b * p.b);
    double atanQ = std::atan2(Q, 2.0 * x + p.b);

    double term1 = std::log(x * x / X);
    double term2 = 2.0 * p.b / Q * atanQ;
    double term3 = -p.b * p.x0 / X0 * (
        std::log((x - p.x0) * (x - p.x0) / X) +
        2.0 * (p.b + 2.0 * p.x0) / Q * atanQ
    );

    return p.A * (term1 + term2 + term3);
}

// Derivative dF/dx for VWN5
static double vwn5_dFdx(double x, const VWN5Params& p) {
    double X = x * x + p.b * x + p.c;
    double X0 = p.x0 * p.x0 + p.b * p.x0 + p.c;
    double Q = std::sqrt(4.0 * p.c - p.b * p.b);

    // d/dx [ln(x²/X)] = 2/x - (2x+b)/X
    double dterm1 = 2.0 / x - (2.0 * x + p.b) / X;

    // d/dx [2b/Q arctan(Q/(2x+b))] = -2b/(X) * Q/(Q² + (2x+b)²) * 2
    // Simplify: = -4b / (Q * X) ... wait, let's be precise.
    // arctan'(u) = 1/(1+u²), u = Q/(2x+b), u' = -2Q/(2x+b)²
    // d/dx = 2b/Q * 1/(1 + Q²/(2x+b)²) * (-2Q/(2x+b)²)
    //       = 2b/Q * (2x+b)² / ((2x+b)² + Q²) * (-2Q/(2x+b)²)
    //       = 2b * (-2) / ((2x+b)² + Q²)
    // Note: (2x+b)² + Q² = 4x² + 4bx + b² + 4c - b² = 4(x²+bx+c) = 4X
    double dterm2 = -4.0 * p.b / (4.0 * X);  // = -b/X

    // d/dx [-bx₀/X₀ * {ln((x-x₀)²/X) + 2(b+2x₀)/Q arctan(...)}]
    double fac = -p.b * p.x0 / X0;
    double dlog = 2.0 / (x - p.x0) - (2.0 * x + p.b) / X;
    double datan = -4.0 * (p.b + 2.0 * p.x0) / (4.0 * X); // same simplification
    double dterm3 = fac * (dlog + datan);

    return p.A * (dterm1 + dterm2 + dterm3);
}

// ── Spin interpolation function f(ζ) ────────────────────────────────────

// f(ζ) = [(1+ζ)^{4/3} + (1-ζ)^{4/3} - 2] / [2(2^{1/3} - 1)]
static const double kFppZero = 4.0 / (9.0 * (std::pow(2.0, 1.0/3.0) - 1.0));
// f''(0) = 4 / (9 * (2^{1/3} - 1))

static double spin_f(double zeta) {
    double denom = 2.0 * (std::pow(2.0, 1.0/3.0) - 1.0);
    return (std::pow(1.0 + zeta, 4.0/3.0) + std::pow(1.0 - zeta, 4.0/3.0) - 2.0) / denom;
}

// df/dζ
static double spin_df(double zeta) {
    double denom = 2.0 * (std::pow(2.0, 1.0/3.0) - 1.0);
    return (4.0/3.0) * (std::pow(1.0 + zeta, 1.0/3.0) - std::pow(1.0 - zeta, 1.0/3.0)) / denom;
}

// ── Slater exchange ─────────────────────────────────────────────────────

// Spin-polarized Slater exchange:
//   ε_x = -2^{1/3} C_x [ρ_α^{4/3} + ρ_β^{4/3}] / ρ
//   v_x^σ = -(3/π)^{1/3} 2^{1/3} ρ_σ^{1/3}  (= -(6/π)^{1/3} ρ_σ^{1/3} * (4/3) ... )
//
// Actually, the standard form:
//   E_x = -C_x_spin * ∫ [ρ_α^{4/3} + ρ_β^{4/3}] dτ
//   where C_x_spin = 2^{1/3} * C_x = 2^{1/3} * (3/4)(3/π)^{1/3}
//
//   Energy density (per volume): f_x = -C_x_spin * [ρ_α^{4/3} + ρ_β^{4/3}]
//   ε_x = f_x / ρ  (energy per electron, but we track f_x = exc * rho)
//
//   v_x^α = ∂f_x/∂ρ_α = -(4/3) C_x_spin ρ_α^{1/3}
//   v_x^β = ∂f_x/∂ρ_β = -(4/3) C_x_spin ρ_β^{1/3}

static void slater_exchange(double rho_a, double rho_b,
                            double& exc, double& vx_a, double& vx_b) {
    double Cx_spin = std::pow(2.0, 1.0/3.0) * kCx;

    double ra43 = std::pow(rho_a, 4.0/3.0);
    double rb43 = std::pow(rho_b, 4.0/3.0);

    // f_x = -Cx_spin * (ρ_α^{4/3} + ρ_β^{4/3})  (energy per volume)
    double fx = -Cx_spin * (ra43 + rb43);
    double rho = rho_a + rho_b;
    exc = fx / rho;  // energy per electron

    vx_a = -(4.0/3.0) * Cx_spin * std::pow(rho_a, 1.0/3.0);
    vx_b = -(4.0/3.0) * Cx_spin * std::pow(rho_b, 1.0/3.0);
}

// ── VWN5 correlation ────────────────────────────────────────────────────

// ε_c(r_s, ζ) = ε_c^P(r_s) + α_c(r_s) f(ζ)/f''(0) (1 - ζ⁴)
//             + [ε_c^F(r_s) - ε_c^P(r_s)] f(ζ) ζ⁴
//
// where r_s = (3/(4πρ))^{1/3}, ζ = (ρ_α - ρ_β)/ρ

static void vwn5_correlation(double rho_a, double rho_b,
                             double& ec, double& vc_a, double& vc_b) {
    double rho = rho_a + rho_b;
    double rs = std::pow(3.0 / (4.0 * kPi * rho), 1.0/3.0);
    double x = std::sqrt(rs);
    double zeta = (rho_a - rho_b) / rho;

    // Clamp zeta to avoid numerical issues at full polarization
    zeta = std::max(-1.0 + 1e-14, std::min(1.0 - 1e-14, zeta));

    // Evaluate the three VWN5 functions
    double ecP = vwn5_F(x, kVWN5_P);
    double ecF = vwn5_F(x, kVWN5_F);
    double ac  = vwn5_F(x, kVWN5_A);

    double fz = spin_f(zeta);
    double z4 = zeta * zeta * zeta * zeta;

    // ε_c = ε_c^P + α_c f(ζ)/f''(0) (1 - ζ⁴) + [ε_c^F - ε_c^P] f(ζ) ζ⁴
    ec = ecP + ac * fz / kFppZero * (1.0 - z4) + (ecF - ecP) * fz * z4;

    // For the potential v_c^σ = ε_c - (r_s/3) dε_c/dr_s + (∓ζ - 1)(dε_c/dζ)
    // where ∓ = - for α, + for β (i.e., -(ζ-1) for α, -(ζ+1) for β... let's be precise)
    //
    // v_c^α = ε_c - (r_s/3) ∂ε_c/∂r_s - (ζ - 1) ∂ε_c/∂ζ
    // v_c^β = ε_c - (r_s/3) ∂ε_c/∂r_s - (ζ + 1) ∂ε_c/∂ζ
    //
    // Since x = √r_s, ∂/∂r_s = (1/(2x)) ∂/∂x
    // So (r_s/3) ∂ε_c/∂r_s = (x²/3)(1/(2x)) ∂ε_c/∂x = x/(6) ∂ε_c/∂x

    double decP_dx = vwn5_dFdx(x, kVWN5_P);
    double decF_dx = vwn5_dFdx(x, kVWN5_F);
    double dac_dx  = vwn5_dFdx(x, kVWN5_A);

    double dec_dx = decP_dx
                  + dac_dx * fz / kFppZero * (1.0 - z4)
                  + (decF_dx - decP_dx) * fz * z4;

    double rs3_dec = x / 6.0 * dec_dx;

    // ∂ε_c/∂ζ
    double dfz = spin_df(zeta);
    double dz4 = 4.0 * zeta * zeta * zeta;

    double dec_dz = ac / kFppZero * (dfz * (1.0 - z4) - fz * dz4)
                  + (ecF - ecP) * (dfz * z4 + fz * dz4);

    vc_a = ec - rs3_dec - (zeta - 1.0) * dec_dz;
    vc_b = ec - rs3_dec - (zeta + 1.0) * dec_dz;
}

// ── PBE exchange ────────────────────────────────────────────────────────

// PBE exchange enhancement factor:
//   F_x(s) = 1 + κ - κ / (1 + μ s² / κ)
//   s = |∇ρ_σ| / (2 (6π²)^{1/3} ρ_σ^{4/3})
//   κ = 0.804,  μ = 0.21951
//
// The total exchange per spin channel:
//   E_x^PBE = E_x^LSDA * (1 + enhancement from F_x)
//
// Implementation: compute per-spin exchange with enhancement.

static constexpr double kPBE_kappa = 0.804;
static constexpr double kPBE_mu = 0.21951;

// Compute PBE exchange for a single spin channel.
// Returns: (energy density contribution, v_rho, v_sigma) for this spin.
static void pbe_exchange_spin(double rho_s, double sigma_ss,
                              double& ex_s, double& vrho_s, double& vsigma_ss) {
    if (rho_s < kRhoMin * 0.5) {
        ex_s = 0; vrho_s = 0; vsigma_ss = 0;
        return;
    }

    double Cx_spin = std::pow(2.0, 1.0/3.0) * kCx;

    // LSDA exchange for this spin (energy per volume = -Cx_spin * ρ_σ^{4/3})
    double rs43 = std::pow(rho_s, 4.0/3.0);
    double ex_lda = -Cx_spin * rs43;  // per volume

    // Reduced gradient s
    double kf = std::pow(6.0 * kPi * kPi, 1.0/3.0) * std::pow(rho_s, 1.0/3.0);
    double grad = std::sqrt(std::max(sigma_ss, 0.0));
    double s = grad / (2.0 * kf * rho_s);
    double s2 = s * s;

    // Enhancement factor F_x(s)
    double denom = 1.0 + kPBE_mu * s2 / kPBE_kappa;
    double Fx = 1.0 + kPBE_kappa - kPBE_kappa / denom;

    // Energy density (per volume) for this spin
    ex_s = ex_lda * Fx;

    // dF_x/ds² = μ / (1 + μs²/κ)²
    double dFx_ds2 = kPBE_mu / (denom * denom);

    // ds²/dρ_σ: s² = σ / (4 k_f² ρ²), k_f = (6π²)^{1/3} ρ^{1/3}
    // s² = σ / (4 (6π²)^{2/3} ρ^{2/3} ρ²) = σ / (4 (6π²)^{2/3} ρ^{8/3})
    // ds²/dρ = -(8/3) s² / ρ
    double ds2_drho = -(8.0/3.0) * s2 / rho_s;

    // ds²/dσ = 1 / (4 (6π²)^{2/3} ρ^{8/3})
    double kf2_fac = 4.0 * std::pow(6.0 * kPi * kPi, 2.0/3.0) * std::pow(rho_s, 8.0/3.0);
    double ds2_dsigma = 1.0 / kf2_fac;

    // v_rho = d(ex_lda * Fx)/dρ = (d ex_lda/dρ) Fx + ex_lda * dFx/dρ
    double dex_lda_drho = -(4.0/3.0) * Cx_spin * std::pow(rho_s, 1.0/3.0);
    vrho_s = dex_lda_drho * Fx + ex_lda * dFx_ds2 * ds2_drho;

    // v_sigma = ex_lda * dFx/dσ
    vsigma_ss = ex_lda * dFx_ds2 * ds2_dsigma;
}

// ── PBE correlation ─────────────────────────────────────────────────────

// PBE correlation gradient correction:
//   H(r_s, ζ, t) = γ φ³ ln{1 + (β/γ) t² [(1 + At²)/(1 + At² + A²t⁴)]}
//   A = (β/γ) / [exp(-ε_c^unif/(γφ³)) - 1]
//   t = |∇ρ| / (2 φ k_s ρ)
//   φ(ζ) = [(1+ζ)^{2/3} + (1-ζ)^{2/3}] / 2
//   k_s = √(4k_F/π),  k_F = (3π²ρ)^{1/3}
//   β = 0.066725,  γ = (1 - ln2)/π² ≈ 0.031091

static constexpr double kPBE_beta = 0.066725;
static constexpr double kPBE_gamma = 0.031091;  // (1 - ln2)/π²

static double phi_pbe(double zeta) {
    return 0.5 * (std::pow(1.0 + zeta, 2.0/3.0) + std::pow(1.0 - zeta, 2.0/3.0));
}

static void pbe_correlation(double rho_a, double rho_b,
                            double ec_lda, double vc_lda_a, double vc_lda_b,
                            double sigma_aa, double sigma_ab, double sigma_bb,
                            double& ec_pbe, double& vc_a, double& vc_b,
                            double& vsigma_aa, double& vsigma_ab, double& vsigma_bb) {
    double rho = rho_a + rho_b;
    double zeta = (rho_a - rho_b) / rho;
    zeta = std::max(-1.0 + 1e-14, std::min(1.0 - 1e-14, zeta));

    double phi = phi_pbe(zeta);
    double phi3 = phi * phi * phi;

    // k_F and k_s
    double kF = std::pow(3.0 * kPi * kPi * rho, 1.0/3.0);
    double ks = std::sqrt(4.0 * kF / kPi);

    // t = |∇ρ| / (2 φ k_s ρ)
    double sigma_total = sigma_aa + 2.0 * sigma_ab + sigma_bb;
    double grad_rho = std::sqrt(std::max(sigma_total, 0.0));
    double t = grad_rho / (2.0 * phi * ks * rho);
    double t2 = t * t;

    // A coefficient
    double exp_arg = -ec_lda / (kPBE_gamma * phi3);
    double exp_val = std::exp(exp_arg);
    double A = (kPBE_beta / kPBE_gamma) / (exp_val - 1.0);

    // H(r_s, ζ, t)
    double At2 = A * t2;
    double numer = 1.0 + At2;
    double denom = 1.0 + At2 + A * A * t2 * t2;
    double ratio = (kPBE_beta / kPBE_gamma) * t2 * numer / denom;
    double H = kPBE_gamma * phi3 * std::log(1.0 + ratio);

    ec_pbe = ec_lda + H;

    // For the derivatives, we use the chain rule.
    // dH/dt² = γφ³ * 1/(1+ratio) * d(ratio)/dt²
    // This is complex; for now we use numerical-like analytic forms.

    // d(ratio)/dt²:
    // ratio = (β/γ) t² (1+At²) / (1+At²+A²t⁴)
    // Let u = At², then ratio = (β/γ) t² (1+u) / (1+u+u²)
    // d(ratio)/dt² = (β/γ) * [(1+u)/(1+u+u²) + t² * d/dt² {(1+u)/(1+u+u²)}]
    // where du/dt² = A
    double u = At2;
    double one_u = 1.0 + u;
    double one_u_u2 = 1.0 + u + u * u;
    double d_frac_du = (1.0 * one_u_u2 - one_u * (1.0 + 2.0 * u)) / (one_u_u2 * one_u_u2);
    // = (1 + u + u² - (1+u)(1+2u)) / (1+u+u²)²
    // = (1 + u + u² - 1 - 2u - u - 2u²) / ...
    // = (-2u - u²) / (1+u+u²)²
    // = -u(2+u) / (1+u+u²)²

    double dratio_dt2 = (kPBE_beta / kPBE_gamma) * (
        one_u / one_u_u2 + t2 * A * d_frac_du
    );

    double dH_dt2 = kPBE_gamma * phi3 / (1.0 + ratio) * dratio_dt2;

    // dt²/dσ_total: t² = σ_total / (4 φ² k_s² ρ²)
    double dt2_dsigma = 1.0 / (4.0 * phi * phi * ks * ks * rho * rho);

    // dt²/dρ = -(7/3) t²/ρ  (t ~ ρ^{-7/6}, t² ~ ρ^{-7/3})
    // t² = σ/(4φ²ks²ρ²), ks² = 4kF/π, kF ~ ρ^{1/3}
    // t² ~ σ / (φ² ρ^{1/3} ρ²) = σ / (φ² ρ^{7/3})
    // dt²/dρ = -(7/3) t²/ρ  (at fixed σ and ζ)

    // For v_sigma: H depends on σ_total only through t²
    double dH_dsigma_total = dH_dt2 * dt2_dsigma;

    // σ_total = σ_αα + 2σ_αβ + σ_ββ
    vsigma_aa = dH_dsigma_total;
    vsigma_ab = 2.0 * dH_dsigma_total;
    vsigma_bb = dH_dsigma_total;

    // For v_rho: more complex — H depends on ρ through t, k_s, A, ec_lda, phi
    // Use the identity for uniform scaling:
    // v_c^σ = (ec_lda + H) - (r_s/3) d(ec_lda+H)/dr_s ∓ (ζ∓1) d(ec_lda+H)/dζ
    // But this gets very involved. A simpler approach for correctness:
    //
    // v^σ_c = vc_lda_σ + H + dH/dρ_σ (at fixed σ's)
    //
    // dH/dρ_α = (∂H/∂ρ)|_ζ * (∂ρ/∂ρ_α) + (∂H/∂ζ)|_ρ * (∂ζ/∂ρ_α)
    //         + (∂H/∂ec_lda) * (∂ec_lda/∂ρ_α)
    //
    // This is still complex. For a robust first implementation:
    // v_c^α = vc_lda_a + H - (σ_total / rho) * dH_dsigma_total
    //         + correction terms...
    //
    // Actually, the simplest correct form for the KS potential contribution:
    // The PBE correlation to the Fock matrix is handled via the grid integration.
    // The key outputs for the grid assembly are:
    //   vrho_α, vrho_β: partial derivatives of f_c w.r.t. ρ_α, ρ_β
    //   vsigma_αα, vsigma_αβ, vsigma_ββ: partial derivatives w.r.t. σ's
    //
    // where f_c = ρ * ε_c(ρ_α, ρ_β, σ's).
    //
    // For the rho derivatives, we need d(ρ·(ec_lda + H))/dρ_α.
    // = (ec_lda + H) + ρ * d(ec_lda + H)/dρ_α
    // The LDA part gives vc_lda_a = d(ρ·ec_lda)/dρ_α.
    // So we need: vc_lda_a + H + ρ * dH/dρ_α.
    //
    // dH/dρ_α = dH/dt² * dt²/dρ_α + dH/dA * dA/dρ_α + dH/dφ * dφ/dρ_α
    //
    // For numerical stability and correctness, compute all partial derivatives:

    // ∂t²/∂ρ_α = ∂t²/∂ρ * ∂ρ/∂ρ_α + ∂t²/∂ζ * ∂ζ/∂ρ_α
    // ∂ρ/∂ρ_α = 1, ∂ζ/∂ρ_α = (1 - ζ)/ρ
    // ∂t²/∂ρ (at fixed ζ) = -(7/3) t²/ρ
    // ∂t²/∂ζ (at fixed ρ): t² ~ 1/φ², so ∂t²/∂ζ = -2 t²/φ * dφ/dζ

    double dphi_dzeta = (1.0/3.0) * (std::pow(1.0 + zeta, -1.0/3.0) - std::pow(1.0 - zeta, -1.0/3.0));
    double dt2_drho = -(7.0/3.0) * t2 / rho;
    double dt2_dzeta = -2.0 * t2 / phi * dphi_dzeta;

    double dzeta_drho_a = (1.0 - zeta) / rho;
    double dzeta_drho_b = -(1.0 + zeta) / rho;

    double dt2_drho_a = dt2_drho + dt2_dzeta * dzeta_drho_a;
    double dt2_drho_b = dt2_drho + dt2_dzeta * dzeta_drho_b;

    // ∂A/∂ρ_α: A depends on ec_lda and φ
    // A = (β/γ) / (exp(-ec_lda/(γφ³)) - 1)
    // dA/dec_lda = A² * exp_val / (β/γ) * 1/(γφ³)
    // dA/dφ = A² * exp_val / (β/γ) * (-3 ec_lda/(γφ⁴))
    double dA_dfactor = A * A * exp_val / (kPBE_beta / kPBE_gamma);
    double dA_dec_lda = dA_dfactor / (kPBE_gamma * phi3);
    double dA_dphi = dA_dfactor * (-3.0 * ec_lda / (kPBE_gamma * phi3 * phi));

    // dec_lda/dρ_α = vc_lda_a - ec_lda (since d(ρ ec_lda)/dρ_α = vc_lda_a, so dec_lda/dρ_α = (vc_lda_a - ec_lda)/ρ ... no)
    // Actually: vc_lda_α = ε_c + ρ ∂ε_c/∂ρ_α, so ∂ε_c/∂ρ_α = (vc_lda_α - ε_c)/ρ ... no, that's not right either.
    // vc_lda_α = ∂(ρ ε_c)/∂ρ_α = ε_c + ρ ∂ε_c/∂ρ_α
    // So ∂ε_c/∂ρ_α = (vc_lda_α - ε_c) ... no, ρ ∂ε_c/∂ρ_α = vc_lda_α - ε_c
    // ∂ε_c/∂ρ_α = (vc_lda_α - ε_c) / ρ ... hmm this seems off.
    // Actually from the chain rule: d(ρ·ec)/dρ_α = ec + ρ · dec/dρ_α = vc_lda_a
    // So dec/dρ_α = (vc_lda_a - ec_lda) / ρ ... wait, no.
    // d(ρ·ec)/dρ_α = ec · dρ/dρ_α + ρ · dec/dρ_α = ec + ρ · dec/dρ_α
    // = vc_lda_a
    // So dec/dρ_α = (vc_lda_a - ec_lda) / ρ
    double dec_lda_drho_a = (vc_lda_a - ec_lda) / rho;
    double dec_lda_drho_b = (vc_lda_b - ec_lda) / rho;

    double dA_drho_a = dA_dec_lda * dec_lda_drho_a + dA_dphi * dphi_dzeta * dzeta_drho_a;
    double dA_drho_b = dA_dec_lda * dec_lda_drho_b + dA_dphi * dphi_dzeta * dzeta_drho_b;

    // ∂H/∂A: from the ratio
    // ratio = (β/γ) t² (1+At²)/(1+At²+A²t⁴)
    // d(ratio)/dA = (β/γ) t² * d/dA[(1+At²)/(1+At²+A²t⁴)]
    // = (β/γ) t² * [t²(1+At²+A²t⁴) - (1+At²)(t²+2At⁴)] / (1+At²+A²t⁴)²
    // = (β/γ) t² * t² [1+At²+A²t⁴ - (1+At²)(1+2At²)] / (...)²  ... wait let me be more careful
    // Let N = 1+u, D = 1+u+u², du/dA = t²
    // d(N/D)/dA = (dN/dA * D - N * dD/dA) / D²
    //           = (t² * D - N * (t² + 2u t²)) / D²
    //           = t² * (D - N(1+2u)) / D²
    //           = t² * (1+u+u² - (1+u)(1+2u)) / D²
    //           = t² * (1+u+u² - 1 - 2u - u - 2u²) / D²
    //           = t² * (-2u - u²) / D²
    //           = -t² * u * (2+u) / D²
    double d_frac_dA = -t2 * u * (2.0 + u) / (one_u_u2 * one_u_u2);
    double dratio_dA = (kPBE_beta / kPBE_gamma) * t2 * d_frac_dA;
    double dH_dA = kPBE_gamma * phi3 / (1.0 + ratio) * dratio_dA;

    // ∂H/∂φ: H = γ φ³ ln(1+ratio)
    // dH/dφ = 3γφ² ln(1+ratio) + γφ³ * d(ln(1+ratio))/dφ
    // The ratio depends on φ only through t (at fixed A):
    // t² ~ 1/φ², so this is captured by dt²/dφ = -2t²/φ
    // But A also depends on φ, which we handle separately.
    // For the explicit φ³ prefactor:
    double dH_dphi_explicit = 3.0 * kPBE_gamma * phi * phi * std::log(1.0 + ratio);
    // Plus the implicit φ dependence through t²:
    double dH_dphi_implicit = dH_dt2 * (-2.0 * t2 / phi);
    // And through A:
    double dH_dphi_A = dH_dA * dA_dphi;
    double dH_dphi_total = dH_dphi_explicit + dH_dphi_implicit + dH_dphi_A;

    // Now assemble dH/dρ_α:
    // dH/dρ_α = dH/dt² * dt²/dρ_α + dH/dA * dA/dρ_α + dH/dφ_explicit * dφ/dρ_α
    //   where dH/dφ_explicit includes only the prefactor + implicit t² dependence
    //   and dA/dρ_α already includes the φ contribution to A
    // Wait, let me reorganize. Total dH/dρ_α:
    // = (∂H/∂t²)(dt²/dρ_α) + (∂H/∂A)(dA/dρ_α) + (∂H/∂φ)_direct * (dφ/dρ_α)
    // where (∂H/∂φ)_direct = 3γφ² ln(1+ratio) only (the explicit φ³ prefactor)
    // and dt²/dρ_α already includes the φ dependence of t²
    // and dA/dρ_α already includes the φ dependence of A
    //
    // Hmm, actually let's think about it differently.
    // H depends on ρ_α through: ρ, ζ (which determines φ), ec_lda, and σ (fixed).
    // The chain is: ρ_α → (ρ, ζ) → (φ, ec_lda, t², A) → H
    //
    // Let me just compute dH/dρ_α directly:
    double dH_drho_a = dH_dt2 * dt2_drho_a
                     + dH_dA * (dA_dec_lda * dec_lda_drho_a)
                     + (dH_dphi_explicit + dH_dphi_implicit + dH_dA * dA_dphi) * dphi_dzeta * dzeta_drho_a;

    double dH_drho_b = dH_dt2 * dt2_drho_b
                     + dH_dA * (dA_dec_lda * dec_lda_drho_b)
                     + (dH_dphi_explicit + dH_dphi_implicit + dH_dA * dA_dphi) * dphi_dzeta * dzeta_drho_b;

    // Final: v_c^σ = d(ρ·(ec_lda + H))/dρ_σ = vc_lda_σ + H + ρ * dH/dρ_σ
    vc_a = vc_lda_a + H + rho * dH_drho_a;
    vc_b = vc_lda_b + H + rho * dH_drho_b;
}

// ── Public API ──────────────────────────────────────────────────────────

XCOutput evaluate_xc(XCType type, Real rho_a, Real rho_b,
                     Real sigma_aa, Real sigma_ab, Real sigma_bb) {
    XCOutput out{};
    double ra = static_cast<double>(rho_a);
    double rb = static_cast<double>(rho_b);
    double rho = ra + rb;

    if (rho < kRhoMin) return out;

    if (type == XCType::LSDA) {
        double ex, vx_a, vx_b;
        slater_exchange(ra, rb, ex, vx_a, vx_b);

        double ec, vc_a, vc_b;
        vwn5_correlation(ra, rb, ec, vc_a, vc_b);

        out.exc = static_cast<Real>(ex + ec);
        out.vrho_alpha = static_cast<Real>(vx_a + vc_a);
        out.vrho_beta  = static_cast<Real>(vx_b + vc_b);
    } else if (type == XCType::PBE) {
        // PBE = LSDA + gradient corrections

        // LSDA part (Slater + VWN5)
        double ex_lda, vx_lda_a, vx_lda_b;
        slater_exchange(ra, rb, ex_lda, vx_lda_a, vx_lda_b);

        double ec_lda, vc_lda_a, vc_lda_b;
        vwn5_correlation(ra, rb, ec_lda, vc_lda_a, vc_lda_b);

        // PBE exchange (per-spin)
        double sa = static_cast<double>(sigma_aa);
        double sab = static_cast<double>(sigma_ab);
        double sb = static_cast<double>(sigma_bb);

        double ex_pbe_a, vrho_x_a, vsig_x_aa;
        double ex_pbe_b, vrho_x_b, vsig_x_bb;
        pbe_exchange_spin(ra, sa, ex_pbe_a, vrho_x_a, vsig_x_aa);
        pbe_exchange_spin(rb, sb, ex_pbe_b, vrho_x_b, vsig_x_bb);

        // Total PBE exchange energy density (per electron)
        double ex_pbe = (ex_pbe_a + ex_pbe_b) / rho;

        // PBE correlation
        double ec_pbe, vc_pbe_a, vc_pbe_b;
        double vsig_c_aa, vsig_c_ab, vsig_c_bb;
        pbe_correlation(ra, rb, ec_lda, vc_lda_a, vc_lda_b,
                       sa, sab, sb,
                       ec_pbe, vc_pbe_a, vc_pbe_b,
                       vsig_c_aa, vsig_c_ab, vsig_c_bb);

        out.exc = static_cast<Real>(ex_pbe + ec_pbe);
        out.vrho_alpha = static_cast<Real>(vrho_x_a + vc_pbe_a);
        out.vrho_beta  = static_cast<Real>(vrho_x_b + vc_pbe_b);
        out.vsigma_aa = static_cast<Real>(vsig_x_aa + vsig_c_aa);
        out.vsigma_ab = static_cast<Real>(vsig_c_ab);
        out.vsigma_bb = static_cast<Real>(vsig_x_bb + vsig_c_bb);
    }

    return out;
}

} // namespace scf
} // namespace cpu
} // namespace ATV
