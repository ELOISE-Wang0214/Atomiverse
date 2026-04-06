/// @file cpu/scf/native_integrals.h
/// @brief Native C++ molecular integral functions for s+p type Gaussians.
///
/// Primary integral engine. Validated element-by-element against libcint
/// as ground truth (enable -DATV_BUILD_LIBCINT=ON for validation tests).
///
/// All functions operate in atomic units (Bohr, Hartree).
/// Uses physical normalization: (2*alpha/pi)^{3/4} for s-type.
/// Shell.l selects angular momentum; Shell.cart selects Cartesian component.
///
/// @see cpu::scf::IntegralPack for the output data structure
/// @see cpu::scf::build_integrals_sto3g for the entry point
#pragma once

namespace ATV {
namespace cpu {
namespace scf {
namespace native {

/// @brief Contracted Gaussian shell (s and p type).
struct Shell {
    int n_prim;            ///< Number of primitive Gaussians.
    const double* alpha;   ///< Exponents [n_prim].
    const double* coeff;   ///< Normalized contraction coefficients [n_prim].
    double center[3];      ///< Shell center in Bohr [x, y, z].
    int l;                 ///< Angular momentum (0=s, 1=p).
    int cart;              ///< Cartesian component: s->0, p->0=x,1=y,2=z.
};

/// @brief Atom descriptor for nuclear attraction integrals.
struct Atom {
    int Z;                 ///< Atomic number (nuclear charge).
    double pos[3];         ///< Position in Bohr [x, y, z].
};

// ── Utility ───────────────────────────────────────────────────────────

/// @brief GTO normalization (libcint convention). Delegates to ATV::gto_norm.
/// @param l     Angular momentum quantum number.
/// @param alpha Gaussian exponent.
/// @return Normalization factor.
double gto_norm(int l, double alpha);

/// @brief Physical GTO normalization for s-type. Delegates to ATV::phys_norm_s.
/// @param alpha Gaussian exponent.
/// @return Normalization factor (2*alpha/pi)^{3/4}.
double phys_norm(double alpha);

/// @brief Physical normalization for p-type (l=1) Gaussians.
/// @param alpha Gaussian exponent.
/// @return Normalization factor.
double phys_norm_p(double alpha);

/// @brief Boys function F_0(x) = integral_0^1 exp(-x*t^2) dt.
/// @details F_0(0) = 1. For x > 0: F_0(x) = 0.5 * sqrt(pi/x) * erf(sqrt(x)).
/// @param x Non-negative argument.
/// @return F_0(x).
double boys_f0(double x);

/// @brief General Boys function F_n(x) = integral_0^1 t^{2n} exp(-x*t^2) dt.
/// @details Uses upward recurrence from F_0. F_n(0) = 1/(2n+1).
/// @param n Order (n >= 0).
/// @param x Non-negative argument.
/// @return F_n(x).
double boys_fn(int n, double x);

// ── One-electron integrals ────────────────────────────────────────────

/// @brief Overlap integral between two contracted s-type shells.
/// @details S = sum_pq c_p c_q (pi/gamma)^{3/2} exp(-mu * R^2)
///   where gamma = alpha_p + alpha_q, mu = alpha_p*alpha_q/gamma.
/// @param a First shell.
/// @param b Second shell.
/// @return Overlap integral <a|b>.
double overlap_ss(const Shell& a, const Shell& b);

/// @brief Kinetic energy integral between two contracted s-type shells.
/// @details T_pq = mu * (3 - 2*mu*R^2) * S_pq for each primitive pair.
/// @param a First shell.
/// @param b Second shell.
/// @return Kinetic integral <a|-nabla^2/2|b>.
double kinetic_ss(const Shell& a, const Shell& b);

/// @brief Nuclear attraction integral between two contracted s-type shells.
/// @details V = sum_C (-Z_C) sum_pq c_p c_q (2pi/gamma) exp(-mu*R^2) F_0(gamma*|P-C|^2)
///   summed over all nuclei C.
/// @param a First shell.
/// @param b Second shell.
/// @param atoms Array of atoms (nuclei).
/// @param natm Number of atoms.
/// @return Nuclear attraction integral <a|sum(-Z/r)|b>.
double nuclear_ss(const Shell& a, const Shell& b,
                  const Atom* atoms, int natm);

// ── Two-electron integrals ────────────────────────────────────────────

/// @brief Two-electron repulsion integral (ab|cd) for s-type shells.
/// @details (ab|cd) = sum c_p c_q c_r c_s * 2*pi^{5/2} / (g1*g2*sqrt(g1+g2))
///            * exp(-mu1*R_AB^2 - mu2*R_CD^2) * F_0(rho*|PQ|^2)
/// @param a Shell for electron 1, bra.
/// @param b Shell for electron 1, ket.
/// @param c Shell for electron 2, bra.
/// @param d Shell for electron 2, ket.
/// @return ERI in chemist's notation (ab|cd).
double eri_ssss(const Shell& a, const Shell& b,
                const Shell& c, const Shell& d);

// ── General integrals (s+p type) ─────────────────────────────────────

/// @brief General overlap integral for shells with l in {0,1}.
/// @details Dispatches to Obara-Saika recurrence based on a.l, b.l, a.cart, b.cart.
/// @param a First shell.
/// @param b Second shell.
/// @return Overlap integral <a|b>.
double overlap(const Shell& a, const Shell& b);

/// @brief General kinetic energy integral for shells with l in {0,1}.
/// @param a First shell.
/// @param b Second shell.
/// @return Kinetic integral <a|-nabla^2/2|b>.
double kinetic(const Shell& a, const Shell& b);

/// @brief General nuclear attraction integral for shells with l in {0,1}.
/// @param a First shell.
/// @param b Second shell.
/// @param atoms Array of atoms (nuclei).
/// @param natm Number of atoms.
/// @return Nuclear attraction integral <a|sum(-Z/r)|b>.
double nuclear(const Shell& a, const Shell& b,
               const Atom* atoms, int natm);

/// @brief General two-electron repulsion integral for shells with l in {0,1}.
/// @param a Shell for electron 1, bra.
/// @param b Shell for electron 1, ket.
/// @param c Shell for electron 2, bra.
/// @param d Shell for electron 2, ket.
/// @return ERI in chemist's notation (ab|cd).
double eri(const Shell& a, const Shell& b,
           const Shell& c, const Shell& d);

// ── Primitive-level functions (exposed for derivative code) ──────────

/// @brief Primitive overlap integral for arbitrary Cartesian angular momenta.
/// @details Computes the 3D overlap integral for a single primitive pair using
///   the Obara-Saika recurrence. Supports angular momenta up to 2 per component.
/// @param a_ang Cartesian angular momenta on center A [3] (each 0..2).
/// @param b_ang Cartesian angular momenta on center B [3] (each 0..2).
/// @param alpha_a Exponent on center A.
/// @param A Center A position [3].
/// @param alpha_b Exponent on center B.
/// @param B Center B position [3].
/// @return Primitive overlap integral value.
double overlap_prim_3d(const int a_ang[3], const int b_ang[3],
                       double alpha_a, const double A[3],
                       double alpha_b, const double B[3]);

/// @brief Primitive kinetic energy integral for arbitrary Cartesian angular momenta.
/// @details Computes the 3D kinetic integral for a single primitive pair.
///   Supports angular momenta up to 2 per component.
/// @param a_ang Cartesian angular momenta on center A [3] (each 0..2).
/// @param b_ang Cartesian angular momenta on center B [3] (each 0..2).
/// @param alpha_a Exponent on center A.
/// @param A Center A position [3].
/// @param alpha_b Exponent on center B.
/// @param B Center B position [3].
/// @return Primitive kinetic integral value.
double kinetic_prim_3d(const int a_ang[3], const int b_ang[3],
                       double alpha_a, const double A[3],
                       double alpha_b, const double B[3]);

/// @brief Primitive nuclear attraction integral for one nucleus.
/// @details McMurchie-Davidson scheme with Hermite Coulomb integrals.
///   Supports angular momenta up to 2 per component (Ltot <= 6).
///   Returns the integral without the -Z_C prefactor.
/// @param a_ang Angular momenta on center A [3] (each 0..2).
/// @param b_ang Angular momenta on center B [3] (each 0..2).
/// @param gamma Sum of exponents (alpha_a + alpha_b).
/// @param PA P - A vector [3].
/// @param PB P - B vector [3].
/// @param PC P - C vector [3].
/// @param K_AB exp(-mu * |A-B|^2).
/// @return Primitive nuclear attraction integral value (without -Z).
double nuclear_prim(const int a_ang[3], const int b_ang[3],
                    double gamma, const double PA[3],
                    const double PB[3], const double PC[3],
                    double K_AB);

/// @brief Primitive two-electron repulsion integral (McMurchie-Davidson).
/// @details Supports angular momenta up to 2 per component (Ltot <= 8).
/// @param a_ang Angular momenta on center A [3] (each 0..2).
/// @param b_ang Angular momenta on center B [3] (each 0..2).
/// @param c_ang Angular momenta on center C [3] (each 0..2).
/// @param d_ang Angular momenta on center D [3] (each 0..2).
/// @param alpha_a Exponent on A.
/// @param A Center A position [3].
/// @param alpha_b Exponent on B.
/// @param B Center B position [3].
/// @param alpha_c Exponent on C.
/// @param C Center C position [3].
/// @param alpha_d Exponent on D.
/// @param D Center D position [3].
/// @return Primitive ERI value.
double eri_prim(const int a_ang[3], const int b_ang[3],
                const int c_ang[3], const int d_ang[3],
                double alpha_a, const double A[3],
                double alpha_b, const double B[3],
                double alpha_c, const double C[3],
                double alpha_d, const double D[3]);

} // namespace native
} // namespace scf
} // namespace cpu
} // namespace ATV
