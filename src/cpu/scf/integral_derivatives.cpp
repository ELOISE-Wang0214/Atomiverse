/// @file cpu/scf/integral_derivatives.cpp
/// @brief Analytic derivatives of STO-3G molecular integrals w.r.t. nuclear positions.
///
/// Uses the angular momentum raising formula for derivatives:
///   d/dA_d <a|O|b> = 2*alpha_a * <a+1_d|O|b> - a_d * <a-1_d|O|b>
///
/// where <a+1_d| bumps angular momentum on center A in direction d by 1,
/// and <a-1_d| lowers it by 1 (vanishes if a_d == 0).
///
/// For two-center integrals, translational invariance gives d/dB = -d/dA.
/// Nuclear attraction has an additional term from the nuclear position.
///
/// Supports s-type (l=0) and p-type (l=1) basis functions.
///
/// References:
///   Obara & Saika, JCP 84, 3963 (1986) -- recurrence relations
///   Szabo & Ostlund, Appendix C -- RHF gradient formulas
///   Helgaker, Jorgensen, Olsen, Ch. 9 -- molecular gradients

#include "atomiverse/cpu/scf/integral_derivatives.h"
#include "atomiverse/cpu/scf/native_integrals.h"
#include "atomiverse/cpu/scf/sto3g_data.h"
#include "atomiverse/core/constants.h"
#include "atomiverse/core/gto.h"

#include <cmath>
#include <vector>
#include <array>

namespace ATV {
namespace cpu {
namespace scf {

/// @brief Get Cartesian angular momentum exponents for a shell.
/// @details s-type (l=0): (0,0,0). p-type (l=1, cart=c): unit vector along c.
static inline void cart_exponents(const native::Shell& sh, int a_out[3]) {
    a_out[0] = a_out[1] = a_out[2] = 0;
    if (sh.l == 1) a_out[sh.cart] = 1;
}

// ── Primitive-level derivative helpers ───────────────────────────────

/// @brief Derivative of primitive overlap w.r.t. center A in all 3 directions.
///
/// d/dA_d S(a,b) = 2*alpha_a * S(a+1_d, b) - a_d * S(a-1_d, b)
///
/// @param a_ang Cartesian angular momenta on A [3].
/// @param b_ang Cartesian angular momenta on B [3].
/// @param alpha_a Exponent on center A.
/// @param A Center A position [3].
/// @param alpha_b Exponent on center B.
/// @param B Center B position [3].
/// @param dS_dA Output: dS/dA in each direction [3].
static void overlap_deriv_prim_general(
    const int a_ang[3], const int b_ang[3],
    double alpha_a, const double A[3],
    double alpha_b, const double B[3],
    double dS_dA[3])
{
    for (int d = 0; d < 3; ++d) {
        // Raising term: 2*alpha_a * S(a+1_d, b)
        int a_plus[3] = {a_ang[0], a_ang[1], a_ang[2]};
        a_plus[d] += 1;
        double S_plus = native::overlap_prim_3d(a_plus, b_ang, alpha_a, A, alpha_b, B);

        double result = 2.0 * alpha_a * S_plus;

        // Lowering term: -a_d * S(a-1_d, b) -- only if a_d > 0
        if (a_ang[d] > 0) {
            int a_minus[3] = {a_ang[0], a_ang[1], a_ang[2]};
            a_minus[d] -= 1;
            double S_minus = native::overlap_prim_3d(a_minus, b_ang, alpha_a, A, alpha_b, B);
            result -= a_ang[d] * S_minus;
        }

        dS_dA[d] = result;
    }
}

/// @brief Derivative of primitive kinetic integral w.r.t. center A in all 3 directions.
///
/// d/dA_d T(a,b) = 2*alpha_a * T(a+1_d, b) - a_d * T(a-1_d, b)
///
/// @param a_ang Cartesian angular momenta on A [3].
/// @param b_ang Cartesian angular momenta on B [3].
/// @param alpha_a Exponent on center A.
/// @param A Center A position [3].
/// @param alpha_b Exponent on center B.
/// @param B Center B position [3].
/// @param dT_dA Output: dT/dA in each direction [3].
static void kinetic_deriv_prim_general(
    const int a_ang[3], const int b_ang[3],
    double alpha_a, const double A[3],
    double alpha_b, const double B[3],
    double dT_dA[3])
{
    for (int d = 0; d < 3; ++d) {
        int a_plus[3] = {a_ang[0], a_ang[1], a_ang[2]};
        a_plus[d] += 1;
        double T_plus = native::kinetic_prim_3d(a_plus, b_ang, alpha_a, A, alpha_b, B);

        double result = 2.0 * alpha_a * T_plus;

        if (a_ang[d] > 0) {
            int a_minus[3] = {a_ang[0], a_ang[1], a_ang[2]};
            a_minus[d] -= 1;
            double T_minus = native::kinetic_prim_3d(a_minus, b_ang, alpha_a, A, alpha_b, B);
            result -= a_ang[d] * T_minus;
        }

        dT_dA[d] = result;
    }
}

/// @brief Derivative of primitive nuclear attraction w.r.t. center A (Gaussian center).
///
/// d/dA_d V(a,b) = 2*alpha_a * V(a+1_d, b) - a_d * V(a-1_d, b)
///
/// V is computed via native::nuclear_prim for a single nucleus.
/// The -Z_C prefactor is NOT included; caller must multiply.
///
/// @param a_ang Cartesian angular momenta on A [3].
/// @param b_ang Cartesian angular momenta on B [3].
/// @param gamma Sum of exponents.
/// @param alpha_a Exponent on center A.
/// @param A Center A position [3].
/// @param alpha_b Exponent on center B.
/// @param B Center B position [3].
/// @param C Nuclear position [3].
/// @param dV_dA Output: dV/dA in each direction [3].
static void nuclear_deriv_prim_center_general(
    const int a_ang[3], const int b_ang[3],
    double alpha_a, const double A[3],
    double alpha_b, const double B[3],
    const double C[3],
    double dV_dA[3])
{
    const double gamma = alpha_a + alpha_b;
    const double mu = alpha_a * alpha_b / gamma;
    const double R2_AB = (A[0]-B[0])*(A[0]-B[0]) + (A[1]-B[1])*(A[1]-B[1]) + (A[2]-B[2])*(A[2]-B[2]);
    const double K_AB = std::exp(-mu * R2_AB);

    double P[3], PA[3], PB[3], PC[3];
    for (int dd = 0; dd < 3; ++dd) {
        P[dd] = (alpha_a * A[dd] + alpha_b * B[dd]) / gamma;
        PA[dd] = P[dd] - A[dd];
        PB[dd] = P[dd] - B[dd];
        PC[dd] = P[dd] - C[dd];
    }

    for (int d = 0; d < 3; ++d) {
        int a_plus[3] = {a_ang[0], a_ang[1], a_ang[2]};
        a_plus[d] += 1;
        double V_plus = native::nuclear_prim(a_plus, b_ang, gamma, PA, PB, PC, K_AB);

        double result = 2.0 * alpha_a * V_plus;

        if (a_ang[d] > 0) {
            int a_minus[3] = {a_ang[0], a_ang[1], a_ang[2]};
            a_minus[d] -= 1;
            double V_minus = native::nuclear_prim(a_minus, b_ang, gamma, PA, PB, PC, K_AB);
            result -= a_ang[d] * V_minus;
        }

        dV_dA[d] = result;
    }
}

/// @brief Derivative of nuclear attraction w.r.t. nuclear position C.
///
/// The nucleus C appears in the Coulomb potential |r-C|^{-1}.
/// d/dC_d V(a,b;C) is computed by finite difference of the McMurchie-Davidson
/// R integrals w.r.t. C. Equivalently:
///   d/dC_d V = -d/dP_d V (for the P-C dependence only, keeping K_AB constant)
///
/// For the McMurchie-Davidson scheme:
///   V = (2pi/gamma) K_AB sum_{tuv} E[t]E[u]E[v] R^0_{tuv}(P-C)
///
/// d/dC_d only affects R^0_{tuv} via its PC dependence:
///   dR^0_{tuv}/dC_d = -dR^0_{tuv}/dPC_d = -(R^0_{t+delta_dx, u+delta_dy, v+delta_dz})
///     at next order = -R^0 with shifted index
///
/// Actually, the simplest correct approach: dV/dC_d = -dV/dP_d (for R dependence).
/// But dV/dP_d acts on both E and R. For the nuclear derivative w.r.t. C:
///   dV/dC_d = -(2pi/gamma) K_AB sum_{tuv} E[t]E[u]E[v] * dR^0_{tuv}/dPC_d
/// where dR^0_{tuv}/dPC_d = R^0_{t+delta_dx, u+delta_dy, v+delta_dz} (at the bumped index).
///
/// More precisely: R^0_{tuv} depends on PC, and the derivative is:
///   dR^0_{0,0,0}/dPC_d = R^1_{1_d, 0, 0} (via the downward recurrence)
///   But since we build R^0 from R^m, a cleaner way is:
///     dR^0_{tuv}/dC_d = -R^0_{t+delta_d(x), u+delta_d(y), v+delta_d(z)}  -- NO, this is wrong.
///
/// The correct way: from dV_pq^C/dC_d, the result is:
///   dV/dC_d = (2pi/gamma) K_AB sum_{tuv} E[t]E[u]E[v] * R^0_{t+delta_dx, u+delta_dy, v+delta_dz}
/// where delta_dx = (d==0), etc. This comes from differentiating the Hermite Coulomb integral
/// R^0_{tuv}(T) w.r.t. C_d.
///
/// Actually, the simplest and most correct formula uses translational invariance:
/// For a two-center integral V^C_ab summed over C:
///   d(total V)/dC_d = (contribution from C as a Gaussian center if C hosts a shell)
///                   + (contribution from C as a nuclear charge center)
///
/// The nuclear contribution dV_pq^C/dC_d has already been derived. We compute it
/// via a direct formula using Boys functions.
static void nuclear_deriv_prim_nuc_general(
    const int a_ang[3], const int b_ang[3],
    double alpha_a, const double A[3],
    double alpha_b, const double B[3],
    const double C[3],
    double dV_dC[3])
{
    const double gamma = alpha_a + alpha_b;
    const double mu = alpha_a * alpha_b / gamma;
    const double R2_AB = (A[0]-B[0])*(A[0]-B[0]) + (A[1]-B[1])*(A[1]-B[1]) + (A[2]-B[2])*(A[2]-B[2]);
    const double K_AB = std::exp(-mu * R2_AB);

    double P[3], PA[3], PB[3], PC[3];
    for (int dd = 0; dd < 3; ++dd) {
        P[dd] = (alpha_a * A[dd] + alpha_b * B[dd]) / gamma;
        PA[dd] = P[dd] - A[dd];
        PB[dd] = P[dd] - B[dd];
        PC[dd] = P[dd] - C[dd];
    }

    // The derivative dV/dC_d of the McMurchie-Davidson nuclear integral
    // w.r.t. nuclear position C is obtained by differentiating R^0_{tuv}:
    //   dR^0_{tuv}/dC_d = -dR^0_{tuv}/dPC_d
    // where dR^0_{tuv}/dPC_d involves the next-order R integrals.
    //
    // For the complete integral:
    //   V = (2pi/gamma) K_AB sum_{tuv} E[t]E[u]E[v] R^0_{tuv}
    //
    //   dV/dC_d = -(2pi/gamma) K_AB sum_{tuv} E[t]E[u]E[v] R^0_{t+d_x, u+d_y, v+d_z}
    //   where d_x = (d==0), d_y = (d==1), d_z = (d==2).
    //
    // This requires R^0 built to one higher order in direction d.
    // We implement this by evaluating V with a bumped angular momentum
    // on a fictitious direction, using the identity:
    //   dV/dC_d = -(2pi/gamma) K_AB sum_{tuv} E[t]E[u]E[v] R^0_{t+d_x, u+d_y, v+d_z}
    //
    // This equals: -V_prim(a_ang, b_ang + delta_d) evaluated with
    //   modified E coefficients for the b side... Actually this is not quite right.
    //
    // Simpler: just use the translational invariance of the total integral.
    // For a given nucleus C:
    //   d/dA_d V_C + d/dB_d V_C + d/dC_d V_C = 0  (translational invariance)
    // so:
    //   d/dC_d V_C = -(d/dA_d V_C + d/dB_d V_C)
    //
    // d/dA_d V_C uses angular momentum raising on center A.
    // d/dB_d V_C uses angular momentum raising on center B.
    // This is the cleanest approach.

    // Compute d/dA_d V_C (raising on center A)
    double dV_dA[3];
    for (int d = 0; d < 3; ++d) {
        int a_plus[3] = {a_ang[0], a_ang[1], a_ang[2]};
        a_plus[d] += 1;
        double V_plus = native::nuclear_prim(a_plus, b_ang, gamma, PA, PB, PC, K_AB);
        double result = 2.0 * alpha_a * V_plus;
        if (a_ang[d] > 0) {
            int a_minus[3] = {a_ang[0], a_ang[1], a_ang[2]};
            a_minus[d] -= 1;
            double V_minus = native::nuclear_prim(a_minus, b_ang, gamma, PA, PB, PC, K_AB);
            result -= a_ang[d] * V_minus;
        }
        dV_dA[d] = result;
    }

    // Compute d/dB_d V_C (raising on center B)
    double dV_dB[3];
    for (int d = 0; d < 3; ++d) {
        int b_plus[3] = {b_ang[0], b_ang[1], b_ang[2]};
        b_plus[d] += 1;
        double V_plus = native::nuclear_prim(a_ang, b_plus, gamma, PA, PB, PC, K_AB);
        double result = 2.0 * alpha_b * V_plus;
        if (b_ang[d] > 0) {
            int b_minus[3] = {b_ang[0], b_ang[1], b_ang[2]};
            b_minus[d] -= 1;
            double V_minus = native::nuclear_prim(a_ang, b_minus, gamma, PA, PB, PC, K_AB);
            result -= b_ang[d] * V_minus;
        }
        dV_dB[d] = result;
    }

    // Translational invariance: dV/dC = -(dV/dA + dV/dB)
    for (int d = 0; d < 3; ++d) {
        dV_dC[d] = -(dV_dA[d] + dV_dB[d]);
    }
}

/// @brief Derivative of primitive ERI w.r.t. center A in all 3 directions.
///
/// d/dA_d (ab|cd) = 2*alpha_a * (a+1_d b|cd) - a_d * (a-1_d b|cd)
///
/// @param a_ang, b_ang, c_ang, d_ang  Cartesian angular momenta [3] each.
/// @param alpha_a..alpha_d  Primitive exponents.
/// @param A..D  Center positions [3] each.
/// @param dERI_dA  Output: d(ERI)/dA in each direction [3].
static void eri_deriv_prim_general(
    const int a_ang[3], const int b_ang[3],
    const int c_ang[3], const int d_ang[3],
    double alpha_a, const double A[3],
    double alpha_b, const double B[3],
    double alpha_c, const double C[3],
    double alpha_d, const double D[3],
    double dERI_dA[3])
{
    for (int d = 0; d < 3; ++d) {
        int a_plus[3] = {a_ang[0], a_ang[1], a_ang[2]};
        a_plus[d] += 1;
        double ERI_plus = native::eri_prim(
            a_plus, b_ang, c_ang, d_ang,
            alpha_a, A, alpha_b, B, alpha_c, C, alpha_d, D);

        double result = 2.0 * alpha_a * ERI_plus;

        if (a_ang[d] > 0) {
            int a_minus[3] = {a_ang[0], a_ang[1], a_ang[2]};
            a_minus[d] -= 1;
            double ERI_minus = native::eri_prim(
                a_minus, b_ang, c_ang, d_ang,
                alpha_a, A, alpha_b, B, alpha_c, C, alpha_d, D);
            result -= a_ang[d] * ERI_minus;
        }

        dERI_dA[d] = result;
    }
}

// ── Main entry point ─────────────────────────────────────────────────

IntegralDerivPack build_integral_derivatives_sto3g(const Frame& frame, int charge) {
    const int natm = static_cast<int>(frame.atoms.size());

    // Count total expanded basis functions using sto3g_data.h
    Int n_basis = 0;
    for (Int i = 0; i < natm; ++i) {
        int Z = static_cast<int>(frame.atoms.type[i]);
        n_basis += sto3g_n_basis(Z);
    }

    IntegralDerivPack pack;
    pack.n_atoms = natm;
    pack.n_basis = n_basis;
    pack.dS     = xt::zeros<Real>({(Int)natm, (Int)3, n_basis, n_basis});
    pack.dT     = xt::zeros<Real>({(Int)natm, (Int)3, n_basis, n_basis});
    pack.dV     = xt::zeros<Real>({(Int)natm, (Int)3, n_basis, n_basis});
    pack.dERI   = xt::zeros<Real>({(Int)natm, (Int)3, n_basis, n_basis, n_basis, n_basis});
    pack.dV_nuc = xt::zeros<Real>({(Int)natm, (Int)3});

    // Build shells and atoms (same structure as integrals.cpp)
    std::vector<native::Shell> shells(n_basis);
    std::vector<native::Atom> atoms(natm);
    std::vector<int> shell_to_atom(n_basis);
    std::vector<std::array<double, 6>> norm_store(n_basis); // alpha[3] + norm_coeff[3]

    int shell_idx = 0;
    for (int a = 0; a < natm; ++a) {
        int Z = static_cast<int>(frame.atoms.type[a]);
        double cx = static_cast<double>(frame.atoms.x[a]) * kAngstromToBohr;
        double cy = static_cast<double>(frame.atoms.y[a]) * kAngstromToBohr;
        double cz = static_cast<double>(frame.atoms.z[a]) * kAngstromToBohr;

        atoms[a].Z = Z;
        atoms[a].pos[0] = cx;
        atoms[a].pos[1] = cy;
        atoms[a].pos[2] = cz;

        const auto& basis = sto3g_data[Z];
        for (int s = 0; s < basis.n_raw_shells; ++s) {
            const auto& raw = basis.raw[s];

            if (raw.type == 0) {
                // Pure s shell: 1 basis function
                for (int p = 0; p < 3; ++p) {
                    norm_store[shell_idx][p] = raw.alpha[p];
                    norm_store[shell_idx][p + 3] = raw.coeff_s[p] * native::phys_norm(raw.alpha[p]);
                }
                shells[shell_idx].n_prim = 3;
                shells[shell_idx].alpha = &norm_store[shell_idx][0];
                shells[shell_idx].coeff = &norm_store[shell_idx][3];
                shells[shell_idx].center[0] = cx;
                shells[shell_idx].center[1] = cy;
                shells[shell_idx].center[2] = cz;
                shells[shell_idx].l = 0;
                shells[shell_idx].cart = 0;
                shell_to_atom[shell_idx] = a;
                ++shell_idx;
            } else {
                // SP shell: expand to 4 basis functions (2s + 2px + 2py + 2pz)

                // s-component
                for (int p = 0; p < 3; ++p) {
                    norm_store[shell_idx][p] = raw.alpha[p];
                    norm_store[shell_idx][p + 3] = raw.coeff_s[p] * native::phys_norm(raw.alpha[p]);
                }
                shells[shell_idx].n_prim = 3;
                shells[shell_idx].alpha = &norm_store[shell_idx][0];
                shells[shell_idx].coeff = &norm_store[shell_idx][3];
                shells[shell_idx].center[0] = cx;
                shells[shell_idx].center[1] = cy;
                shells[shell_idx].center[2] = cz;
                shells[shell_idx].l = 0;
                shells[shell_idx].cart = 0;
                shell_to_atom[shell_idx] = a;
                ++shell_idx;

                // p-components: px (cart=0), py (cart=1), pz (cart=2)
                for (int c = 0; c < 3; ++c) {
                    for (int p = 0; p < 3; ++p) {
                        norm_store[shell_idx][p] = raw.alpha[p];
                        norm_store[shell_idx][p + 3] = raw.coeff_p[p] * native::phys_norm_p(raw.alpha[p]);
                    }
                    shells[shell_idx].n_prim = 3;
                    shells[shell_idx].alpha = &norm_store[shell_idx][0];
                    shells[shell_idx].coeff = &norm_store[shell_idx][3];
                    shells[shell_idx].center[0] = cx;
                    shells[shell_idx].center[1] = cy;
                    shells[shell_idx].center[2] = cz;
                    shells[shell_idx].l = 1;
                    shells[shell_idx].cart = c;
                    shell_to_atom[shell_idx] = a;
                    ++shell_idx;
                }
            }
        }
    }

    // ── 1e integral derivatives ──────────────────────────────────────

    for (int i = 0; i < n_basis; ++i) {
        int atom_i = shell_to_atom[i];
        int a_ang[3];
        cart_exponents(shells[i], a_ang);

        for (int j = 0; j <= i; ++j) {
            int atom_j = shell_to_atom[j];
            int b_ang[3];
            cart_exponents(shells[j], b_ang);

            // Overlap and kinetic derivatives
            double dS_total[3] = {0, 0, 0};
            double dT_total[3] = {0, 0, 0};

            for (int p = 0; p < shells[i].n_prim; ++p) {
                for (int q = 0; q < shells[j].n_prim; ++q) {
                    double ci = shells[i].coeff[p];
                    double cj = shells[j].coeff[q];

                    double dS_dA[3], dT_dA[3];
                    overlap_deriv_prim_general(
                        a_ang, b_ang,
                        shells[i].alpha[p], shells[i].center,
                        shells[j].alpha[q], shells[j].center,
                        dS_dA);
                    kinetic_deriv_prim_general(
                        a_ang, b_ang,
                        shells[i].alpha[p], shells[i].center,
                        shells[j].alpha[q], shells[j].center,
                        dT_dA);

                    for (int d = 0; d < 3; ++d) {
                        dS_total[d] += ci * cj * dS_dA[d];
                        dT_total[d] += ci * cj * dT_dA[d];
                    }
                }
            }

            // Assign to atom_i (derivative w.r.t. center of shell i)
            for (int d = 0; d < 3; ++d) {
                pack.dS(atom_i, d, i, j) += static_cast<Real>(dS_total[d]);
                pack.dT(atom_i, d, i, j) += static_cast<Real>(dT_total[d]);
                // Translational invariance: d/dB = -d/dA for S and T
                pack.dS(atom_j, d, i, j) += static_cast<Real>(-dS_total[d]);
                pack.dT(atom_j, d, i, j) += static_cast<Real>(-dT_total[d]);
                // Symmetry: S(i,j) = S(j,i)
                if (i != j) {
                    pack.dS(atom_i, d, j, i) += static_cast<Real>(dS_total[d]);
                    pack.dT(atom_i, d, j, i) += static_cast<Real>(dT_total[d]);
                    pack.dS(atom_j, d, j, i) += static_cast<Real>(-dS_total[d]);
                    pack.dT(atom_j, d, j, i) += static_cast<Real>(-dT_total[d]);
                }
            }

            // Nuclear attraction derivatives: contributions from each nucleus
            for (int c = 0; c < natm; ++c) {
                double dV_center_i[3] = {0, 0, 0};
                double dV_nuc_c[3] = {0, 0, 0};

                for (int p = 0; p < shells[i].n_prim; ++p) {
                    for (int q = 0; q < shells[j].n_prim; ++q) {
                        double ci = shells[i].coeff[p];
                        double cj = shells[j].coeff[q];

                        // d/dA (center of shell i) -- angular momentum raising on A
                        double dV_dA[3];
                        nuclear_deriv_prim_center_general(
                            a_ang, b_ang,
                            shells[i].alpha[p], shells[i].center,
                            shells[j].alpha[q], shells[j].center,
                            atoms[c].pos,
                            dV_dA);

                        // d/dC (nuclear position) -- via translational invariance
                        double dV_dC[3];
                        nuclear_deriv_prim_nuc_general(
                            a_ang, b_ang,
                            shells[i].alpha[p], shells[i].center,
                            shells[j].alpha[q], shells[j].center,
                            atoms[c].pos,
                            dV_dC);

                        for (int d = 0; d < 3; ++d) {
                            dV_center_i[d] += ci * cj * (-atoms[c].Z) * dV_dA[d];
                            dV_nuc_c[d]    += ci * cj * (-atoms[c].Z) * dV_dC[d];
                        }
                    }
                }

                // For V, the Gaussian center derivative d/dA_d applies to atom_i.
                // The derivative d/dB_d (center of shell j) is obtained via
                // translational invariance: dV/dA + dV/dB + dV/dC = 0
                // so dV/dB = -(dV/dA + dV/dC).
                // But we already computed dV/dC via translational invariance too.
                // Actually: dV_dC above already uses d/dC = -(d/dA + d/dB).
                // So we need dV/dB separately.
                //
                // Compute d/dB (center of shell j) by angular momentum raising on B.
                double dV_center_j[3] = {0, 0, 0};
                for (int p = 0; p < shells[i].n_prim; ++p) {
                    for (int q = 0; q < shells[j].n_prim; ++q) {
                        double ci = shells[i].coeff[p];
                        double cj = shells[j].coeff[q];

                        // d/dB uses raising on center B -- swap roles of A and B
                        double dV_dB[3];
                        nuclear_deriv_prim_center_general(
                            b_ang, a_ang,
                            shells[j].alpha[q], shells[j].center,
                            shells[i].alpha[p], shells[i].center,
                            atoms[c].pos,
                            dV_dB);

                        for (int d = 0; d < 3; ++d) {
                            dV_center_j[d] += ci * cj * (-atoms[c].Z) * dV_dB[d];
                        }
                    }
                }

                for (int d = 0; d < 3; ++d) {
                    pack.dV(atom_i, d, i, j) += static_cast<Real>(dV_center_i[d]);
                    pack.dV(atom_j, d, i, j) += static_cast<Real>(dV_center_j[d]);
                    pack.dV(c,      d, i, j) += static_cast<Real>(dV_nuc_c[d]);
                    if (i != j) {
                        pack.dV(atom_i, d, j, i) += static_cast<Real>(dV_center_i[d]);
                        pack.dV(atom_j, d, j, i) += static_cast<Real>(dV_center_j[d]);
                        pack.dV(c,      d, j, i) += static_cast<Real>(dV_nuc_c[d]);
                    }
                }
            }
        }
    }

    // ── 2e integral derivatives ──────────────────────────────────────

    for (int i = 0; i < n_basis; ++i) {
        int ai = shell_to_atom[i];
        int a_ang[3];
        cart_exponents(shells[i], a_ang);

        for (int j = 0; j < n_basis; ++j) {
            int aj = shell_to_atom[j];
            int b_ang[3];
            cart_exponents(shells[j], b_ang);

            for (int k = 0; k < n_basis; ++k) {
                int ak = shell_to_atom[k];
                int c_ang[3];
                cart_exponents(shells[k], c_ang);

                for (int l = 0; l < n_basis; ++l) {
                    int al = shell_to_atom[l];
                    int d_ang[3];
                    cart_exponents(shells[l], d_ang);

                    for (int p = 0; p < shells[i].n_prim; ++p) {
                      for (int q = 0; q < shells[j].n_prim; ++q) {
                        for (int r = 0; r < shells[k].n_prim; ++r) {
                          for (int s = 0; s < shells[l].n_prim; ++s) {
                            double coeff = shells[i].coeff[p] * shells[j].coeff[q]
                                         * shells[k].coeff[r] * shells[l].coeff[s];

                            // d/dA (center of shell i)
                            double dA[3];
                            eri_deriv_prim_general(
                                a_ang, b_ang, c_ang, d_ang,
                                shells[i].alpha[p], shells[i].center,
                                shells[j].alpha[q], shells[j].center,
                                shells[k].alpha[r], shells[k].center,
                                shells[l].alpha[s], shells[l].center,
                                dA);

                            // d/dB (center of shell j): raise angular momentum on B
                            double dB[3];
                            eri_deriv_prim_general(
                                b_ang, a_ang, c_ang, d_ang,
                                shells[j].alpha[q], shells[j].center,
                                shells[i].alpha[p], shells[i].center,
                                shells[k].alpha[r], shells[k].center,
                                shells[l].alpha[s], shells[l].center,
                                dB);

                            // d/dC (center of shell k): raise angular momentum on C
                            double dC[3];
                            eri_deriv_prim_general(
                                c_ang, d_ang, a_ang, b_ang,
                                shells[k].alpha[r], shells[k].center,
                                shells[l].alpha[s], shells[l].center,
                                shells[i].alpha[p], shells[i].center,
                                shells[j].alpha[q], shells[j].center,
                                dC);

                            // d/dD (center of shell l): raise angular momentum on D
                            double dD[3];
                            eri_deriv_prim_general(
                                d_ang, c_ang, a_ang, b_ang,
                                shells[l].alpha[s], shells[l].center,
                                shells[k].alpha[r], shells[k].center,
                                shells[i].alpha[p], shells[i].center,
                                shells[j].alpha[q], shells[j].center,
                                dD);

                            for (int d = 0; d < 3; ++d) {
                                pack.dERI(ai, d, i, j, k, l) += static_cast<Real>(coeff * dA[d]);
                                pack.dERI(aj, d, i, j, k, l) += static_cast<Real>(coeff * dB[d]);
                                pack.dERI(ak, d, i, j, k, l) += static_cast<Real>(coeff * dC[d]);
                                pack.dERI(al, d, i, j, k, l) += static_cast<Real>(coeff * dD[d]);
                            }
                          }
                        }
                      }
                    }
                }
            }
        }
    }

    // ── Nuclear repulsion derivative ─────────────────────────────────
    // V_nuc = sum_{A<B} Z_A * Z_B / R_AB
    // dV_nuc/dA_d = sum_{B!=A} Z_A * Z_B * (B_d - A_d) / R_AB^3
    for (int a = 0; a < natm; ++a) {
        for (int b = 0; b < natm; ++b) {
            if (a == b) continue;
            double dx = atoms[b].pos[0] - atoms[a].pos[0];
            double dy = atoms[b].pos[1] - atoms[a].pos[1];
            double dz = atoms[b].pos[2] - atoms[a].pos[2];
            double r = std::sqrt(dx*dx + dy*dy + dz*dz);
            double r3_inv = 1.0 / (r * r * r);
            double ZZ = atoms[a].Z * atoms[b].Z;
            pack.dV_nuc(a, 0) += static_cast<Real>(ZZ * dx * r3_inv);
            pack.dV_nuc(a, 1) += static_cast<Real>(ZZ * dy * r3_inv);
            pack.dV_nuc(a, 2) += static_cast<Real>(ZZ * dz * r3_inv);
        }
    }

    return pack;
}

} // namespace scf
} // namespace cpu
} // namespace ATV
