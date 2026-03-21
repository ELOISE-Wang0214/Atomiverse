/// @file cpu/scf/native_integrals.cpp
/// @brief Native C++ molecular integral implementations for s-type Gaussians.
///
/// Drop-in replacements for libcint functions. All formulas from:
///   Szabo & Ostlund, "Modern Quantum Chemistry" (1996), Ch. 3 (Appendix A).
///   Obara-Saika recurrence for general case; direct formulas for s-type.
///
/// Validated element-by-element against libcint in test_native_integrals.cu.

#include "atomiverse/cpu/scf/native_integrals.h"
#include "atomiverse/core/gto.h"

#include <cmath>

namespace ATV {
namespace cpu {
namespace scf {
namespace native {

static constexpr double PI = 3.14159265358979323846;

// ── Helpers ───────────────────────────────────────────────────────────

/// Squared distance between two 3D points.
static inline double dist2(const double a[3], const double b[3]) {
    double dx = a[0]-b[0], dy = a[1]-b[1], dz = a[2]-b[2];
    return dx*dx + dy*dy + dz*dz;
}

// ── Phase 1: GTO normalization ────────────────────────────────────────

double gto_norm(int l, double alpha) { return ATV::gto_norm(l, alpha); }
double phys_norm(double alpha) { return ATV::phys_norm_s(alpha); }

// ── Phase 4 prerequisite: Boys function ───────────────────────────────

double boys_f0(double x) { return ATV::boys_f0(x); }
double boys_fn(int n, double x) { return ATV::boys_fn(n, x); }

// ── Phase 2: Overlap integral (s|s) ───────────────────────────────────

/// Primitive overlap for s-type Gaussians.
/// <a|b> = (pi/gamma)^{3/2} * exp(-mu * R^2)
/// where gamma = alpha_a + alpha_b, mu = alpha_a * alpha_b / gamma
static inline double overlap_primitive(double alpha_a, const double A[3],
                                       double alpha_b, const double B[3]) {
    double gamma = alpha_a + alpha_b;
    double mu = alpha_a * alpha_b / gamma;
    double R2 = dist2(A, B);
    return std::pow(PI / gamma, 1.5) * std::exp(-mu * R2);
}

double overlap_ss(const Shell& a, const Shell& b) {
    // S = sum_pq c_p * c_q * <prim_p | prim_q>
    // PERF: unchecked() candidate — 9 iterations for STO-3G
    double S = 0;
    for (int p = 0; p < a.n_prim; ++p)
        for (int q = 0; q < b.n_prim; ++q)
            S += a.coeff[p] * b.coeff[q]
               * overlap_primitive(a.alpha[p], a.center, b.alpha[q], b.center);
    return S;
}

// ── Phase 3: Kinetic energy integral (s|T|s) ──────────────────────────

double kinetic_ss(const Shell& a, const Shell& b) {
    // T_pq = mu * (3 - 2*mu*R^2) * S_pq  for each primitive pair
    // PERF: unchecked() candidate — 9 iterations for STO-3G
    double T = 0;
    double R2 = dist2(a.center, b.center);
    for (int p = 0; p < a.n_prim; ++p)
        for (int q = 0; q < b.n_prim; ++q) {
            double gamma = a.alpha[p] + b.alpha[q];
            double mu = a.alpha[p] * b.alpha[q] / gamma;
            double S_pq = std::pow(PI / gamma, 1.5) * std::exp(-mu * R2);
            T += a.coeff[p] * b.coeff[q] * mu * (3.0 - 2.0 * mu * R2) * S_pq;
        }
    return T;
}

// ── Phase 4: Nuclear attraction integral (s|V|s) ──────────────────────

double nuclear_ss(const Shell& a, const Shell& b,
                  const Atom* atoms, int natm) {
    // V = sum_C (-Z_C) * sum_pq c_p c_q * (2*pi/gamma) * exp(-mu*R^2) * F_0(gamma*|P-C|^2)
    // P = Gaussian product center = (alpha_a * A + alpha_b * B) / gamma
    // PERF: unchecked() candidate — 9 * natm iterations
    double V = 0;
    double R2_AB = dist2(a.center, b.center);

    for (int p = 0; p < a.n_prim; ++p)
        for (int q = 0; q < b.n_prim; ++q) {
            double alpha_a = a.alpha[p];
            double alpha_b = b.alpha[q];
            double gamma = alpha_a + alpha_b;
            double mu = alpha_a * alpha_b / gamma;

            // Gaussian product center P
            double P[3];
            for (int d = 0; d < 3; ++d)
                P[d] = (alpha_a * a.center[d] + alpha_b * b.center[d]) / gamma;

            double prefactor = a.coeff[p] * b.coeff[q]
                             * (2.0 * PI / gamma) * std::exp(-mu * R2_AB);

            // Sum over all nuclei
            for (int c = 0; c < natm; ++c) {
                double PC2 = dist2(P, atoms[c].pos);
                V += prefactor * (-atoms[c].Z) * boys_f0(gamma * PC2);
            }
        }
    return V;
}

// ── Phase 5: Two-electron ERI (ss|ss) ─────────────────────────────────

double eri_ssss(const Shell& a, const Shell& b,
                const Shell& c, const Shell& d) {
    // (ab|cd) = sum_pqrs c_p c_q c_r c_s
    //   * 2*pi^{5/2} / (gamma1 * gamma2 * sqrt(gamma1 + gamma2))
    //   * exp(-mu1*R_AB^2 - mu2*R_CD^2)
    //   * F_0(rho * |PQ|^2)
    //
    // gamma1 = alpha_p + alpha_q, gamma2 = alpha_r + alpha_s
    // mu1 = alpha_p*alpha_q/gamma1, mu2 = alpha_r*alpha_s/gamma2
    // rho = gamma1 * gamma2 / (gamma1 + gamma2)
    // P = (alpha_p*A + alpha_q*B)/gamma1, Q = (alpha_r*C + alpha_s*D)/gamma2
    //
    // PERF: unchecked() candidate — 81 iterations for STO-3G
    double result = 0;
    double R2_AB = dist2(a.center, b.center);
    double R2_CD = dist2(c.center, d.center);

    for (int p = 0; p < a.n_prim; ++p)
      for (int q = 0; q < b.n_prim; ++q) {
        double gamma1 = a.alpha[p] + b.alpha[q];
        double mu1 = a.alpha[p] * b.alpha[q] / gamma1;
        double K_AB = std::exp(-mu1 * R2_AB);

        double P[3];
        for (int dd = 0; dd < 3; ++dd)
            P[dd] = (a.alpha[p] * a.center[dd] + b.alpha[q] * b.center[dd]) / gamma1;

        for (int r = 0; r < c.n_prim; ++r)
          for (int s = 0; s < d.n_prim; ++s) {
            double gamma2 = c.alpha[r] + d.alpha[s];
            double mu2 = c.alpha[r] * d.alpha[s] / gamma2;
            double K_CD = std::exp(-mu2 * R2_CD);

            double Q[3];
            for (int dd = 0; dd < 3; ++dd)
                Q[dd] = (c.alpha[r] * c.center[dd] + d.alpha[s] * d.center[dd]) / gamma2;

            double rho = gamma1 * gamma2 / (gamma1 + gamma2);
            double PQ2 = dist2(P, Q);

            result += a.coeff[p] * b.coeff[q] * c.coeff[r] * d.coeff[s]
                    * 2.0 * std::pow(PI, 2.5)
                    / (gamma1 * gamma2 * std::sqrt(gamma1 + gamma2))
                    * K_AB * K_CD
                    * boys_f0(rho * PQ2);
          }
      }
    return result;
}

} // namespace native
} // namespace scf
} // namespace cpu
} // namespace ATV
