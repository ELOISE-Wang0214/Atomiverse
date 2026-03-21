/// @file cpu/scf/integral_derivatives.cpp
/// @brief Analytic derivatives of STO-3G molecular integrals w.r.t. nuclear positions.
///
/// For s-type Gaussians, the derivative w.r.t. center A acts as:
///   d/dA_x <a|O|b> = -2*alpha_a * <a+1_x|O|b>  (shift angular momentum)
///
/// where <a+1_x| denotes a p_x Gaussian on center A. Translational invariance
/// gives d/dB = -d/dA for two-center integrals.
///
/// References:
///   Obara & Saika, JCP 84, 3963 (1986) -- recurrence relations
///   Szabo & Ostlund, Appendix C -- RHF gradient formulas
///   Helgaker, Jorgensen, Olsen, Ch. 9 -- molecular gradients

#include "atomiverse/cpu/scf/integral_derivatives.h"
#include "atomiverse/cpu/scf/native_integrals.h"
#include "atomiverse/core/constants.h"
#include "atomiverse/core/gto.h"

#include <cmath>
#include <vector>
#include <array>

namespace ATV {
namespace cpu {
namespace scf {

static constexpr double PI = 3.14159265358979323846;

/// Squared distance between two 3D points.
static inline double dist2(const double a[3], const double b[3]) {
    double dx = a[0]-b[0], dy = a[1]-b[1], dz = a[2]-b[2];
    return dx*dx + dy*dy + dz*dz;
}

// ── Primitive-level derivative integrals ──────────────────────────────

/// @brief Derivative of primitive overlap w.r.t. center A.
///
/// d/dA_d <a|b> = -2*alpha_a * [(A_d - P_d) * S_pq + (d==x,y,z component)]
/// For s-type: d<s_a|s_b>/dA_d = 2*alpha_a * (B_d - A_d) * alpha_b/gamma * S_pq
///
/// Equivalent form using product center P = (alpha_a*A + alpha_b*B)/gamma:
///   d/dA_d S = 2*alpha_a*(P_d - A_d)*S/... = alpha_b/gamma * 2*alpha_a*(B_d-A_d)/gamma * ...
///
/// Simplest form for (ss) overlap:
///   dS_pq/dA_d = -2*mu*(A_d - B_d) * S_pq
/// where mu = alpha_a * alpha_b / gamma, S_pq = (pi/gamma)^{3/2} exp(-mu*R^2)
static inline void overlap_deriv_prim(
    double alpha_a, const double A[3],
    double alpha_b, const double B[3],
    double dS_dA[3])  // output: dS/dA [3]
{
    double gamma = alpha_a + alpha_b;
    double mu = alpha_a * alpha_b / gamma;
    double R2 = dist2(A, B);
    double S_pq = std::pow(PI / gamma, 1.5) * std::exp(-mu * R2);

    for (int d = 0; d < 3; ++d)
        dS_dA[d] = -2.0 * mu * (A[d] - B[d]) * S_pq;
}

/// @brief Derivative of primitive kinetic integral w.r.t. center A.
///
/// T_pq = mu * (3 - 2*mu*R^2) * S_pq
/// dT/dA_d = mu * [(3 - 2*mu*R^2) * dS/dA_d + S_pq * (-2*mu*dR^2/dA_d)]
/// where dR^2/dA_d = 2*(A_d - B_d)
static inline void kinetic_deriv_prim(
    double alpha_a, const double A[3],
    double alpha_b, const double B[3],
    double dT_dA[3])
{
    double gamma = alpha_a + alpha_b;
    double mu = alpha_a * alpha_b / gamma;
    double R2 = dist2(A, B);
    double S_pq = std::pow(PI / gamma, 1.5) * std::exp(-mu * R2);
    double T_factor = mu * (3.0 - 2.0 * mu * R2);

    double dS_dA[3];
    overlap_deriv_prim(alpha_a, A, alpha_b, B, dS_dA);

    for (int d = 0; d < 3; ++d) {
        double dR2_dA = 2.0 * (A[d] - B[d]);
        dT_dA[d] = T_factor * dS_dA[d] + mu * (-2.0 * mu * dR2_dA) * S_pq;
    }
}

/// @brief Derivative of primitive nuclear attraction w.r.t. center A (Gaussian center).
///
/// V_pq^C = (-Z_C) * (2*pi/gamma) * exp(-mu*R_AB^2) * F_0(gamma*|P-C|^2)
/// where P = (alpha_a*A + alpha_b*B)/gamma
///
/// dV_pq^C/dA_d = (-Z_C) * (2*pi/gamma) * [
///    exp(-mu*R_AB^2) * {-2*mu*(A_d-B_d)} * F_0(gamma*|P-C|^2)
///  + exp(-mu*R_AB^2) * dF_0/d(arg) * gamma * 2*(P_d-C_d) * alpha_a/gamma
/// ]
///
/// dF_0/d(arg) = -F_1(arg)  where F_1 is the Boys function of order 1
static inline void nuclear_deriv_prim_center(
    double alpha_a, const double A[3],
    double alpha_b, const double B[3],
    int Z_C, const double C[3],
    double dV_dA[3])
{
    double gamma = alpha_a + alpha_b;
    double mu = alpha_a * alpha_b / gamma;
    double R2_AB = dist2(A, B);
    double K = std::exp(-mu * R2_AB);

    double P[3];
    for (int d = 0; d < 3; ++d)
        P[d] = (alpha_a * A[d] + alpha_b * B[d]) / gamma;

    double PC2 = dist2(P, C);
    double arg = gamma * PC2;
    double F0 = boys_f0(arg);
    double F1 = boys_fn(1, arg);
    double prefactor = (-Z_C) * (2.0 * PI / gamma);

    for (int d = 0; d < 3; ++d) {
        double dK_dA = -2.0 * mu * (A[d] - B[d]) * K;
        double dP_dA = alpha_a / gamma;  // dP_d/dA_d
        double darg_dA = gamma * 2.0 * (P[d] - C[d]) * dP_dA;

        dV_dA[d] = prefactor * (dK_dA * F0 + K * (-F1) * darg_dA);
    }
}

/// @brief Derivative of nuclear attraction w.r.t. nuclear position C.
///
/// The nucleus C appears in the potential |r-C|^{-1}.
/// dV_pq^C/dC_d = (-Z_C) * (2*pi/gamma) * exp(-mu*R_AB^2)
///   * (-F_1(gamma*|P-C|^2)) * gamma * 2*(P_d-C_d) * (-1)
/// = (-Z_C) * (2*pi/gamma) * K * F_1(arg) * gamma * 2*(P_d-C_d)
static inline void nuclear_deriv_prim_nuc(
    double alpha_a, const double A[3],
    double alpha_b, const double B[3],
    int Z_C, const double C[3],
    double dV_dC[3])
{
    double gamma = alpha_a + alpha_b;
    double mu = alpha_a * alpha_b / gamma;
    double R2_AB = dist2(A, B);
    double K = std::exp(-mu * R2_AB);

    double P[3];
    for (int d = 0; d < 3; ++d)
        P[d] = (alpha_a * A[d] + alpha_b * B[d]) / gamma;

    double PC2 = dist2(P, C);
    double arg = gamma * PC2;
    double F1 = boys_fn(1, arg);
    double prefactor = (-Z_C) * (2.0 * PI / gamma) * K;

    // darg/dC_d = gamma * 2 * (P_d - C_d) * (-1) = -2*gamma*(P_d-C_d)
    // dF_0/darg = -F_1
    // dV/dC_d = prefactor * (-F_1) * (-2*gamma*(P_d-C_d))
    //         = prefactor * 2*gamma*(P_d-C_d)*F_1
    for (int d = 0; d < 3; ++d)
        dV_dC[d] = prefactor * 2.0 * gamma * (P[d] - C[d]) * F1;
}

/// @brief Derivative of primitive ERI w.r.t. center A.
///
/// (ab|cd) = K_AB * K_CD * prefactor * F_0(rho*|PQ|^2)
/// where K_AB = exp(-mu1*R_AB^2), P = (alpha_a*A+alpha_b*B)/gamma1
///
/// d/dA_d = d(K_AB)/dA_d * rest + K_AB * d(F_0(rho*|PQ|^2))/dA_d
///        = [-2*mu1*(A_d-B_d)*K_AB] * K_CD * pref * F_0
///          + K_AB * K_CD * pref * (-F_1) * rho * 2*(P_d-Q_d) * alpha_a/gamma1
static inline void eri_deriv_prim(
    double alpha_a, const double A[3],
    double alpha_b, const double B[3],
    double alpha_c, const double C[3],
    double alpha_d, const double D[3],
    double dERI_dA[3])
{
    double gamma1 = alpha_a + alpha_b;
    double gamma2 = alpha_c + alpha_d;
    double mu1 = alpha_a * alpha_b / gamma1;
    double mu2 = alpha_c * alpha_d / gamma2;
    double R2_AB = dist2(A, B);
    double R2_CD = dist2(C, D);
    double K_AB = std::exp(-mu1 * R2_AB);
    double K_CD = std::exp(-mu2 * R2_CD);

    double P[3], Q[3];
    for (int d = 0; d < 3; ++d) {
        P[d] = (alpha_a * A[d] + alpha_b * B[d]) / gamma1;
        Q[d] = (alpha_c * C[d] + alpha_d * D[d]) / gamma2;
    }

    double rho = gamma1 * gamma2 / (gamma1 + gamma2);
    double PQ2 = dist2(P, Q);
    double arg = rho * PQ2;
    double F0 = boys_f0(arg);
    double F1 = boys_fn(1, arg);
    double pref = 2.0 * std::pow(PI, 2.5)
                / (gamma1 * gamma2 * std::sqrt(gamma1 + gamma2));

    for (int d = 0; d < 3; ++d) {
        double dK_AB_dA = -2.0 * mu1 * (A[d] - B[d]) * K_AB;
        double dP_dA = alpha_a / gamma1;
        double darg_dA = rho * 2.0 * (P[d] - Q[d]) * dP_dA;

        dERI_dA[d] = pref * K_CD * (
            dK_AB_dA * F0 + K_AB * (-F1) * darg_dA
        );
    }
}

// ── Shell-level derivative integrals ──────────────────────────────────

// Reuse the STO-3G basis data from integrals.cpp

struct STO3GShell {
    int l;
    int n_prim;
    double alpha[3];
    double coeff[3];
};

struct STO3GBasis {
    int n_shells;
    STO3GShell shells[2];
};

static const STO3GBasis sto3g_basis[] = {
    {0, {}},
    {1, {{0, 3, {3.42525091, 0.62391373, 0.16885540}, {0.15432897, 0.53532814, 0.44463454}}}},
    {1, {{0, 3, {6.36242139, 1.15892300, 0.31364979}, {0.15432897, 0.53532814, 0.44463454}}}},
    {2, {
        {0, 3, {16.1195750, 2.93620070, 0.79465370}, { 0.15432897,  0.53532814,  0.44463454}},
        {0, 3, { 0.63628970, 0.14786010, 0.04808870}, {-0.09996723,  0.39951283,  0.70011547}},
    }},
};

IntegralDerivPack build_integral_derivatives_sto3g(const Frame& frame, int charge) {
    int natm = static_cast<int>(frame.n_atoms());
    Int n_basis = 0;
    for (Int i = 0; i < natm; ++i) {
        int Z = static_cast<int>(frame.type_id[i]);
        n_basis += sto3g_basis[Z].n_shells;
    }

    IntegralDerivPack pack;
    pack.n_atoms = natm;
    pack.n_basis = n_basis;
    pack.dS     = xt::zeros<Real>({(Int)natm, (Int)3, n_basis, n_basis});
    pack.dT     = xt::zeros<Real>({(Int)natm, (Int)3, n_basis, n_basis});
    pack.dV     = xt::zeros<Real>({(Int)natm, (Int)3, n_basis, n_basis});
    pack.dERI   = xt::zeros<Real>({(Int)natm, (Int)3, n_basis, n_basis, n_basis, n_basis});
    pack.dV_nuc = xt::zeros<Real>({(Int)natm, (Int)3});

    // Build shells and atoms (same as integrals.cpp)
    std::vector<std::array<double, 3>> norm_coeffs(n_basis);
    std::vector<native::Shell> shells(n_basis);
    std::vector<native::Atom> atoms(natm);
    std::vector<int> shell_to_atom(n_basis);  // which atom owns each shell

    int shell_idx = 0;
    for (int a = 0; a < natm; ++a) {
        int Z = static_cast<int>(frame.type_id[a]);
        double cx = static_cast<double>(frame.x[a]) * kAngstromToBohr;
        double cy = static_cast<double>(frame.y[a]) * kAngstromToBohr;
        double cz = static_cast<double>(frame.z[a]) * kAngstromToBohr;

        atoms[a].Z = Z;
        atoms[a].pos[0] = cx; atoms[a].pos[1] = cy; atoms[a].pos[2] = cz;

        const auto& basis = sto3g_basis[Z];
        for (int s = 0; s < basis.n_shells; ++s) {
            const auto& sh = basis.shells[s];
            for (int p = 0; p < 3; ++p)
                norm_coeffs[shell_idx][p] = sh.coeff[p] * native::phys_norm(sh.alpha[p]);
            shells[shell_idx].n_prim = 3;
            shells[shell_idx].alpha = sh.alpha;
            shells[shell_idx].coeff = norm_coeffs[shell_idx].data();
            shells[shell_idx].center[0] = cx;
            shells[shell_idx].center[1] = cy;
            shells[shell_idx].center[2] = cz;
            shell_to_atom[shell_idx] = a;
            ++shell_idx;
        }
    }

    // ── 1e integral derivatives ──────────────────────────────────────

    for (int i = 0; i < n_basis; ++i) {
        int atom_i = shell_to_atom[i];
        for (int j = 0; j <= i; ++j) {
            int atom_j = shell_to_atom[j];

            // Overlap and kinetic derivatives
            double dS_total[3] = {0, 0, 0};
            double dT_total[3] = {0, 0, 0};

            for (int p = 0; p < shells[i].n_prim; ++p)
                for (int q = 0; q < shells[j].n_prim; ++q) {
                    double ci = shells[i].coeff[p];
                    double cj = shells[j].coeff[q];
                    double dS_dA[3], dT_dA[3];

                    overlap_deriv_prim(
                        shells[i].alpha[p], shells[i].center,
                        shells[j].alpha[q], shells[j].center,
                        dS_dA);
                    kinetic_deriv_prim(
                        shells[i].alpha[p], shells[i].center,
                        shells[j].alpha[q], shells[j].center,
                        dT_dA);

                    for (int d = 0; d < 3; ++d) {
                        dS_total[d] += ci * cj * dS_dA[d];
                        dT_total[d] += ci * cj * dT_dA[d];
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

            // Nuclear attraction derivatives: two contributions per nucleus
            // (1) derivative w.r.t. Gaussian center (shell i or j)
            // (2) derivative w.r.t. nuclear position C
            for (int c = 0; c < natm; ++c) {
                double dV_center_i[3] = {0, 0, 0};
                double dV_center_j[3] = {0, 0, 0};
                double dV_nuc_c[3] = {0, 0, 0};

                for (int p = 0; p < shells[i].n_prim; ++p)
                    for (int q = 0; q < shells[j].n_prim; ++q) {
                        double ci = shells[i].coeff[p];
                        double cj = shells[j].coeff[q];

                        // d/dA (center of shell i)
                        double dV_dA[3];
                        nuclear_deriv_prim_center(
                            shells[i].alpha[p], shells[i].center,
                            shells[j].alpha[q], shells[j].center,
                            atoms[c].Z, atoms[c].pos,
                            dV_dA);

                        // d/dB (center of shell j) by swapping a<->b
                        double dV_dB[3];
                        nuclear_deriv_prim_center(
                            shells[j].alpha[q], shells[j].center,
                            shells[i].alpha[p], shells[i].center,
                            atoms[c].Z, atoms[c].pos,
                            dV_dB);

                        // d/dC (nuclear position)
                        double dV_dC[3];
                        nuclear_deriv_prim_nuc(
                            shells[i].alpha[p], shells[i].center,
                            shells[j].alpha[q], shells[j].center,
                            atoms[c].Z, atoms[c].pos,
                            dV_dC);

                        for (int d = 0; d < 3; ++d) {
                            dV_center_i[d] += ci * cj * dV_dA[d];
                            dV_center_j[d] += ci * cj * dV_dB[d];
                            dV_nuc_c[d]    += ci * cj * dV_dC[d];
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

    for (int i = 0; i < n_basis; ++i)
      for (int j = 0; j < n_basis; ++j)
        for (int k = 0; k < n_basis; ++k)
          for (int l = 0; l < n_basis; ++l) {
            int ai = shell_to_atom[i];
            int aj = shell_to_atom[j];
            int ak = shell_to_atom[k];
            int al = shell_to_atom[l];

            // Derivative w.r.t. each of the 4 shell centers
            for (int p = 0; p < shells[i].n_prim; ++p)
              for (int q = 0; q < shells[j].n_prim; ++q)
                for (int r = 0; r < shells[k].n_prim; ++r)
                  for (int s = 0; s < shells[l].n_prim; ++s) {
                    double coeff = shells[i].coeff[p] * shells[j].coeff[q]
                                 * shells[k].coeff[r] * shells[l].coeff[s];

                    // d/dA (center of shell i)
                    double dA[3];
                    eri_deriv_prim(
                        shells[i].alpha[p], shells[i].center,
                        shells[j].alpha[q], shells[j].center,
                        shells[k].alpha[r], shells[k].center,
                        shells[l].alpha[s], shells[l].center,
                        dA);

                    // d/dB (center of shell j): swap a<->b in bra
                    double dB[3];
                    eri_deriv_prim(
                        shells[j].alpha[q], shells[j].center,
                        shells[i].alpha[p], shells[i].center,
                        shells[k].alpha[r], shells[k].center,
                        shells[l].alpha[s], shells[l].center,
                        dB);

                    // d/dC (center of shell k): swap so k is "first" in ket
                    double dC[3];
                    eri_deriv_prim(
                        shells[k].alpha[r], shells[k].center,
                        shells[l].alpha[s], shells[l].center,
                        shells[i].alpha[p], shells[i].center,
                        shells[j].alpha[q], shells[j].center,
                        dC);

                    // d/dD (center of shell l): swap so l is "first" in ket
                    double dD[3];
                    eri_deriv_prim(
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

    // ── Nuclear repulsion derivative ─────────────────────────────────
    // V_nuc = sum_{A<B} Z_A * Z_B / R_AB
    // dV_nuc/dA_d = sum_{B!=A} Z_A * Z_B * (B_d - A_d) / R_AB^3
    for (int a = 0; a < natm; ++a)
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

    return pack;
}

} // namespace scf
} // namespace cpu
} // namespace ATV
