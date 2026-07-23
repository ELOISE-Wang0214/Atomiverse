/// @file cpu/scf/internal/mp2_gradient.cpp
/// @brief Closed-shell RMP2 analytic nuclear gradient (PySCF-compatible).
///
/// Algorithm mirrors pyscf.grad.mp2.grad_elec (no frozen orbitals):
///   T₂ → unrelaxed 1-RDM (doo,dvv) + non-separable 2-RDM
///   I-matrix + L_vo Lagrangian → Z-vector on ¹(A+B)
///   contract relaxed densities with skeleton integral derivatives.

#include "atomiverse/cpu/scf/internal/mp2_gradient.h"

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/containers/xtensor.hpp>
#include <xtensor/views/xview.hpp>

#include "atomiverse/core/units.h"
#include "atomiverse/cpu/scf/internal/integral_derivatives.h"
#include "atomiverse/cpu/scf/internal/orbital_hessian.h"
#include "atomiverse/cpu/scf/internal/shell_pair.h"
#include "atomiverse/cpu/scf/mp2.h"
#include "atomiverse/cpu/scf/rys_eri.h"

namespace ATV {
namespace cpu {
namespace scf {
namespace internal {

namespace {

/// Closed-shell G[X]_{μν} = Σ_{λσ} X_{λσ} [(μν|λσ) − ½(μλ|νσ)].
RealMatrix g_contraction(const xt::xtensor<double, 4>& eri,
                         const RealMatrix& X) {
    const int n = static_cast<int>(X.shape()[0]);
    RealMatrix G = xt::zeros<double>({n, n});
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            double acc = 0.0;
            for (int k = 0; k < n; ++k)
                for (int l = 0; l < n; ++l)
                    acc += X(k, l)
                           * (eri(i, j, k, l) - 0.5 * eri(i, k, j, l));
            G(i, j) = acc;
        }
    return G;
}

GradientResult make_force_result(const xt::xtensor<double, 2>& grad_ha_bohr,
                                 double energy_ha, int natm) {
    GradientResult out;
    out.forces = xt::zeros<double>({natm, 3});
    for (int A = 0; A < natm; ++A)
        for (int d = 0; d < 3; ++d)
            out.forces(A, d) = -units::convert(
                grad_ha_bohr(A, d), units::Force::HartreePerBohr,
                units::Force::eVPerAngstrom);
    out.energy =
        units::convert(energy_ha, units::Energy::Hartree, units::Energy::eV);
    return out;
}

} // namespace

GradientResult mp2_analytic_forces(const Frame& frame, const IntegralPack& pack,
                                   const SlaterDeterminant& det,
                                   double e_scf_ha, double e_corr_ha) {
    if (!det.restricted())
        throw std::invalid_argument(
            "mp2_analytic_forces: closed-shell RHF reference required");
    if (pack.ERI.size() == 0)
        throw std::invalid_argument(
            "mp2_analytic_forces: dense ERI required (compute_eri=true)");
    if (pack.spherical())
        throw std::invalid_argument(
            "mp2_analytic_forces: spherical-convention packs are not "
            "supported in this link — use a Cartesian basis (STO-3G / Pople)");

    using Spin = SlaterDeterminant::Spin;
    const RealMatrix C = det.coefficients(Spin::alpha);
    const RealVector eps = det.orbital_energies(Spin::alpha);
    const int n = static_cast<int>(C.shape()[0]);
    const int n_mo = static_cast<int>(C.shape()[1]);
    const int nocc = det.n_electrons(Spin::alpha);
    const int nvir = n_mo - nocc;
    const int natm = static_cast<int>(frame.nlocal());
    if (nvir <= 0)
        throw std::invalid_argument(
            "mp2_analytic_forces: no virtual orbitals (n_virt=0)");
    if (n != pack.n_basis)
        throw std::invalid_argument(
            "mp2_analytic_forces: AO dimension mismatch between det and pack");

    const xt::xtensor<double, 4>& eri = pack.ERI;
    const RealMatrix Co = xt::view(C, xt::all(), xt::range(0, nocc));
    const RealMatrix Cv = xt::view(C, xt::all(), xt::range(nocc, n_mo));

    // t2[i,j,a,b] = (ia|jb) / (ε_i + ε_j − ε_a − ε_b)
    const xt::xtensor<double, 4> ovov = ao2mo_ovov(eri, C, nocc, n_mo);
    xt::xtensor<double, 4> t2 = xt::zeros<double>({nocc, nocc, nvir, nvir});
    for (int i = 0; i < nocc; ++i)
        for (int j = 0; j < nocc; ++j)
            for (int a = 0; a < nvir; ++a)
                for (int b = 0; b < nvir; ++b) {
                    const double denom =
                        eps(i) + eps(j) - eps(nocc + a) - eps(nocc + b);
                    if (denom >= 0.0)
                        throw std::invalid_argument(
                            "mp2_analytic_forces: non-negative orbital-energy "
                            "denominator; reference is degenerate or "
                            "non-aufbau");
                    t2(i, j, a, b) = ovov(i, a, j, b) / denom;
                }

    // Unrelaxed 1-RDM (PySCF _gamma1_intermediates; returns −dm1occ, +dm1vir)
    RealMatrix doo = xt::zeros<double>({nocc, nocc});
    RealMatrix dvv = xt::zeros<double>({nvir, nvir});
    for (int i = 0; i < nocc; ++i) {
        // dm1occ(p,q) += Σ_ab [2 t[i,p,a,b] t[i,q,a,b] − t[i,p,a,b] t[i,q,b,a]]
        for (int p = 0; p < nocc; ++p)
            for (int q = 0; q < nocc; ++q) {
                double s2 = 0.0, sx = 0.0;
                for (int a = 0; a < nvir; ++a)
                    for (int b = 0; b < nvir; ++b) {
                        s2 += t2(i, p, a, b) * t2(i, q, a, b);
                        sx += t2(i, p, a, b) * t2(i, q, b, a);
                    }
                doo(p, q) += 2.0 * s2 - sx;
            }
        // dm1vir: einsum('jca,jcb->ba') → out(b,a) += t[i,j,c,a] t[i,j,c,b]
        for (int a = 0; a < nvir; ++a)
            for (int b = 0; b < nvir; ++b) {
                double s2 = 0.0, sx = 0.0;
                for (int j = 0; j < nocc; ++j)
                    for (int c = 0; c < nvir; ++c) {
                        s2 += t2(i, j, c, a) * t2(i, j, c, b);
                        sx += t2(i, j, c, a) * t2(i, j, b, c);
                    }
                dvv(b, a) += 2.0 * s2 - sx;
            }
    }
    doo = -doo;

    // part[i,j,u,v] = Σ_ab t[i,j,a,b] Cv[u,a] Cv[v,b]
    xt::xtensor<double, 4> part = xt::zeros<double>({nocc, nocc, n, n});
    for (int i = 0; i < nocc; ++i)
        for (int j = 0; j < nocc; ++j)
            for (int u = 0; u < n; ++u)
                for (int v = 0; v < n; ++v) {
                    double s = 0.0;
                    for (int a = 0; a < nvir; ++a)
                        for (int b = 0; b < nvir; ++b)
                            s += t2(i, j, a, b) * Cv(u, a) * Cv(v, b);
                    part(i, j, u, v) = s;
                }
    // part_dm2[i,u,v,j] = 4 part[i,j,u,v] − 2 part[i,j,v,u]
    xt::xtensor<double, 4> part_dm2 =
        xt::zeros<double>({nocc, n, n, nocc});
    for (int i = 0; i < nocc; ++i)
        for (int j = 0; j < nocc; ++j)
            for (int u = 0; u < n; ++u)
                for (int v = 0; v < n; ++v)
                    part_dm2(i, u, v, j) =
                        4.0 * part(i, j, u, v) - 2.0 * part(i, j, v, u);

    // Non-separable AO 2-RDM G[p,q,r,s]
    xt::xtensor<double, 4> Gamma = xt::zeros<double>({n, n, n, n});
    for (int p = 0; p < n; ++p)
        for (int q = 0; q < n; ++q)
            for (int r = 0; r < n; ++r)
                for (int s = 0; s < n; ++s) {
                    double acc = 0.0;
                    for (int i = 0; i < nocc; ++i)
                        for (int j = 0; j < nocc; ++j) {
                            acc += Co(p, i) * part_dm2(i, q, r, j) * Co(s, j);
                            acc += Co(q, i) * part_dm2(i, p, r, j) * Co(s, j);
                        }
                    Gamma(p, q, r, s) = acc;
                }
    {
        // (r,s) symmetrization as in PySCF dm2buf; full-ERI contraction then
        // needs the factor ½ that the s2kl pack path applies (diag×½ / tril
        // pairing) — otherwise Imat and de_ns are exactly 2× the PySCF values.
        const xt::xtensor<double, 4> Gs = Gamma;
        for (int p = 0; p < n; ++p)
            for (int q = 0; q < n; ++q)
                for (int r = 0; r < n; ++r)
                    for (int s = 0; s < n; ++s)
                        Gamma(p, q, r, s) =
                            0.5 * (Gs(p, q, r, s) + Gs(p, q, s, r));
    }

    // Imat[α,β] = Σ_{p,r,s} (p α | r s) Gamma[p, β, r, s]
    RealMatrix Imat = xt::zeros<double>({n, n});
    for (int p = 0; p < n; ++p)
        for (int alpha = 0; alpha < n; ++alpha)
            for (int beta = 0; beta < n; ++beta) {
                double acc = 0.0;
                for (int r = 0; r < n; ++r)
                    for (int s = 0; s < n; ++s)
                        acc += eri(p, alpha, r, s) * Gamma(p, beta, r, s);
                Imat(alpha, beta) += acc;
            }
    xt::xtensor<double, 2> de = xt::zeros<double>({natm, 3});
    const RealMatrix P_hf = det.density_total();

    // First-index dERI loop (int2e_ip1 analogue): non-separable 2-RDM
    // contribution + PySCF-style vhf1 = G-like skeleton from P_HF.
    // vhf1[A,d] is n×n, built exactly as pyscf.grad.mp2 (both index blocks).
    xt::xtensor<double, 4> vhf1 =
        xt::zeros<double>({natm, 3, n, n});
    {
        const ShellBlockPairList bpl{pack, /*thr=*/1e-14};
        const int n_blocks = static_cast<int>(bpl.blocks.size());
        std::vector<int> block_to_atom(bpl.blocks.size(), -1);
        for (std::size_t b = 0; b < bpl.blocks.size(); ++b) {
            const auto& Bl = bpl.blocks[b];
            for (int a = 0; a < natm; ++a) {
                const double dx = Bl.center[0] - pack.atoms[a].pos[0];
                const double dy = Bl.center[1] - pack.atoms[a].pos[1];
                const double dz = Bl.center[2] - pack.atoms[a].pos[2];
                if (dx * dx + dy * dy + dz * dz < 1e-20) {
                    block_to_atom[b] = a;
                    break;
                }
            }
        }
        for (int bi = 0; bi < n_blocks; ++bi) {
            const int A = block_to_atom[static_cast<std::size_t>(bi)];
            if (A < 0)
                continue;
            const auto& BA = bpl.blocks[static_cast<std::size_t>(bi)];
            for (int bj = 0; bj < n_blocks; ++bj) {
                const auto& BB = bpl.blocks[static_cast<std::size_t>(bj)];
                for (int bk = 0; bk < n_blocks; ++bk) {
                    const auto& BC = bpl.blocks[static_cast<std::size_t>(bk)];
                    for (int bl = 0; bl < n_blocks; ++bl) {
                        const auto& BD =
                            bpl.blocks[static_cast<std::size_t>(bl)];
                        const int na = BA.n_cart, nb = BB.n_cart,
                                  nc = BC.n_cart, nd = BD.n_cart;
                        const int block_n = na * nb * nc * nd;
                        std::vector<double> deriv(
                            static_cast<std::size_t>(3 * block_n), 0.0);
                        native::rys_eri_block_deriv(BA, BB, BC, BD,
                                                    /*deriv_center=*/0,
                                                    deriv.data());
                        for (int ia = 0; ia < na; ++ia) {
                            const int p = BA.ao_start + ia;
                            for (int ib = 0; ib < nb; ++ib) {
                                const int q = BB.ao_start + ib;
                                for (int ic = 0; ic < nc; ++ic) {
                                    const int r = BC.ao_start + ic;
                                    for (int id = 0; id < nd; ++id) {
                                        const int s = BD.ao_start + id;
                                        const int ao =
                                            ((ia * nb + ib) * nc + ic) * nd
                                            + id;
                                        const double g = Gamma(p, q, r, s);
                                        for (int ax = 0; ax < 3; ++ax) {
                                            // rys_eri_block_deriv ∂/∂A matches
                                            // the HF 2e-gradient engine; PySCF
                                            // int2e_ip1 is −∇_bra. Flip once so
                                            // the MP2 assembly matches PySCF's
                                            // nabla-on-bra sign convention.
                                            const double der =
                                                -deriv[static_cast<std::size_t>(
                                                    ax * block_n + ao)];
                                            // non-separable 2-RDM · dERI
                                            de(A, ax) -= 2.0 * der * g;
                                            // vhf1 from P_HF (PySCF mp2.py)
                                            vhf1(A, ax, r, s) +=
                                                der * P_hf(p, q);
                                            vhf1(A, ax, r, q) -=
                                                0.5 * der * P_hf(p, s);
                                            vhf1(A, ax, p, q) +=
                                                der * P_hf(r, s);
                                            vhf1(A, ax, p, s) -=
                                                0.5 * der * P_hf(q, r);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    IntegralDerivPack dpack{pack};

    // Unrelaxed dm1_mo
    RealMatrix dm1mo = xt::zeros<double>({n_mo, n_mo});
    for (int i = 0; i < nocc; ++i)
        for (int j = 0; j < nocc; ++j)
            dm1mo(i, j) = doo(i, j) + doo(j, i);
    for (int a = 0; a < nvir; ++a)
        for (int b = 0; b < nvir; ++b)
            dm1mo(nocc + a, nocc + b) = dvv(a, b) + dvv(b, a);

    // Imat_mo = − Cᵀ Imat S C
    const RealMatrix Imat_S = xt::linalg::dot(Imat, pack.S);
    RealMatrix Imat_mo =
        -xt::linalg::dot(xt::transpose(C), xt::linalg::dot(Imat_S, C));

    // Xvo = C_vᵀ (2 G[dm1_ao]) C_o + I_ovᵀ − I_vo
    RealMatrix dm1_ao =
        xt::linalg::dot(C, xt::linalg::dot(dm1mo, xt::transpose(C)));
    const RealMatrix vhf = g_contraction(eri, dm1_ao) * 2.0;
    RealMatrix Xvo =
        xt::linalg::dot(xt::transpose(Cv), xt::linalg::dot(vhf, Co));
    for (int a = 0; a < nvir; ++a)
        for (int i = 0; i < nocc; ++i)
            Xvo(a, i) += Imat_mo(i, nocc + a) - Imat_mo(nocc + a, i);

    // Z-vector: ¹(A+B) Z = −Xvo  (PySCF cphf.solve sign for h1=Xvo)
    OrbitalHessian oh{det, eri};
    const int nov = nocc * nvir;
    const RealMatrix Aplus = oh.dense_plus(Channel::singlet);
    RealVector rhs = xt::zeros<double>({static_cast<std::size_t>(nov)});
    for (int a = 0; a < nvir; ++a)
        for (int i = 0; i < nocc; ++i)
            rhs(a * nocc + i) = -Xvo(a, i);
    const RealVector zvec = xt::linalg::solve(Aplus, rhs);
    for (int a = 0; a < nvir; ++a)
        for (int i = 0; i < nocc; ++i) {
            dm1mo(nocc + a, i) = zvec(a * nocc + i);
            dm1mo(i, nocc + a) = zvec(a * nocc + i);
        }

    for (int i = 0; i < nocc; ++i)
        for (int a = 0; a < nvir; ++a)
            Imat_mo(nocc + a, i) = Imat_mo(i, nocc + a);
    const RealMatrix im1 =
        xt::linalg::dot(C, xt::linalg::dot(Imat_mo, xt::transpose(C)));

    // zeta_mo = ζ ⊙ dm1mo
    RealMatrix zeta_mo = xt::zeros<double>({n_mo, n_mo});
    for (int p = 0; p < n_mo; ++p)
        for (int q = 0; q < n_mo; ++q)
            zeta_mo(p, q) = 0.5 * (eps(p) + eps(q));
    for (int a = 0; a < nvir; ++a)
        for (int i = 0; i < nocc; ++i) {
            zeta_mo(nocc + a, i) = eps(i);
            zeta_mo(i, nocc + a) = eps(i);
        }
    for (int p = 0; p < n_mo; ++p)
        for (int q = 0; q < n_mo; ++q)
            zeta_mo(p, q) *= dm1mo(p, q);
    RealMatrix zeta =
        xt::linalg::dot(C, xt::linalg::dot(zeta_mo, xt::transpose(C)));

    dm1_ao = xt::linalg::dot(C, xt::linalg::dot(dm1mo, xt::transpose(C)));
    const RealMatrix P_occ = xt::linalg::dot(Co, xt::transpose(Co));
    const RealMatrix dm1_sym = dm1_ao + xt::transpose(dm1_ao);
    const RealMatrix g_dm1 = g_contraction(eri, dm1_sym);
    const RealMatrix vhf_s1occ =
        xt::linalg::dot(P_occ, xt::linalg::dot(g_dm1, P_occ));

    const RealMatrix dm1p = P_hf + dm1_ao * 2.0;
    const RealMatrix dm1_tot = P_hf + dm1_ao;

    // HF energy-weighted density (closed-shell factor-2 occupations)
    RealMatrix W_hf = xt::zeros<double>({n, n});
    for (int i = 0; i < nocc; ++i)
        for (int mu = 0; mu < n; ++mu)
            for (int nu = 0; nu < n; ++nu)
                W_hf(mu, nu) += 2.0 * eps(i) * C(mu, i) * C(nu, i);
    zeta += W_hf;

    // Final assembly. IntegralDerivPack dS/dH are FULL ∂/∂R_A (bra+ket).
    for (int A = 0; A < natm; ++A) {
        for (int d = 0; d < 3; ++d) {
            const RealMatrix dS =
                xt::view(dpack.dS, A, d, xt::all(), xt::all());
            const RealMatrix dH =
                xt::view(dpack.dT, A, d, xt::all(), xt::all())
                + xt::view(dpack.dV, A, d, xt::all(), xt::all());
            const RealMatrix v1 =
                xt::view(vhf1, A, d, xt::all(), xt::all());

            de(A, d) += xt::sum(dS * im1)();
            de(A, d) += xt::sum(dH * xt::transpose(dm1_tot))();
            de(A, d) -= xt::sum(dS * zeta)();
            de(A, d) -= xt::sum(dS * vhf_s1occ)();
            de(A, d) -= xt::sum(v1 * dm1p)();
            de(A, d) += dpack.dV_nuc(A, d);
        }
    }

    return make_force_result(de, e_scf_ha + e_corr_ha, natm);
}

} // namespace internal
} // namespace scf
} // namespace cpu
} // namespace ATV
