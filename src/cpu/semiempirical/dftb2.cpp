/// @file src/cpu/semiempirical/dftb2.cpp
/// SCC-DFTB2 multi-element: mio-1-1 H–H pin + sketch (spin + analytic F) + 3ob.
/// Closed-shell analytic nuclear forces AND closed-form pair-virial cell stress
/// on BOTH SK sets; MioSketch adds the collinear open-shell gradient (whose
/// stress is still multi-SCC affine strain-FD).
/// Γ-point PBC (MioSketch): lattice-summed H0/S, screened γ, and repulsive.
/// Multi-k (MioSketch restricted): Bloch-phased complex H0(k)/S(k), MP mesh.
/// Optional additive Grimme D3(BJ) dispersion (molecular sum, analytic forces).
#include "atomiverse/cpu/semiempirical/dftb2.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/core/mat3.h"
#include "atomiverse/core/units.h"
#include "atomiverse/cpu/gpw/internal/dftd3.h"
#include "atomiverse/cpu/gpw/internal/dftd3_params.h"
#include "atomiverse/cpu/linalg/eigensolve.h"
#include "atomiverse/cpu/linalg/matrix.h"
#include "atomiverse/cpu/semiempirical/basis_data/dftb_3ob_tables.h"
#include "atomiverse/cpu/semiempirical/basis_data/dftb_mio_tables.h"
#include "atomiverse/cpu/semiempirical/internal/mio_hh_table.h"
#include "atomiverse/cpu/semiempirical/internal/scc_mixer.h"

namespace ATV {
namespace cpu {
namespace semiempirical {
namespace {

using basis_data::DftbElementRow;
using basis_data::DftbPairRepRow;
using basis_data::DftbPairSkRow;
using basis_data::kDftbElements;
using basis_data::kDftbNElements;
using basis_data::kDftbNPairRep;
using basis_data::kDftbNPairSk;
using basis_data::kDftbPairRep;
using basis_data::kDftbPairSk;
using basis_data::Dftb3obElement;
using basis_data::Dftb3obPair;
using basis_data::Dftb3obSkPoint;
using basis_data::kDftb3obElements;
using basis_data::kDftb3obNElements;
using basis_data::kDftb3obNPairs;
using basis_data::kDftb3obPairs;

// Hartree → eV; Bohr → Å (CODATA-aligned with common DFTB pins).
constexpr double kHa2eV = 27.211386245988;
constexpr double kBohr2A = 0.529177210903;
constexpr double kA2Bohr = 1.0 / kBohr2A;

const DftbElementRow& element_row(int z) {
    for (std::size_t i = 0; i < kDftbNElements; ++i)
        if (kDftbElements[i].z == z)
            return kDftbElements[i];
    throw std::runtime_error("DFTB2: unsupported element Z=" + std::to_string(z)
                             + " (static mio-style tables cover H/C/N/O only)");
}

const DftbPairSkRow& pair_sk(int za, int zb) {
    for (std::size_t i = 0; i < kDftbNPairSk; ++i)
        if (kDftbPairSk[i].z_a == za && kDftbPairSk[i].z_b == zb)
            return kDftbPairSk[i];
    throw std::runtime_error("DFTB2: no SK table for pair Z="
                             + std::to_string(za) + "-" + std::to_string(zb));
}

const DftbPairRepRow& pair_rep(int za, int zb) {
    for (std::size_t i = 0; i < kDftbNPairRep; ++i)
        if (kDftbPairRep[i].z_a == za && kDftbPairRep[i].z_b == zb)
            return kDftbPairRep[i];
    throw std::runtime_error("DFTB2: no repulsive table for pair Z="
                             + std::to_string(za) + "-" + std::to_string(zb));
}

double h_radial(double v0, double r, double r0, double d) {
    return v0 * std::exp(-(r - r0) / d);
}
double s_radial(double s0, double r, double r0, double w) {
    const double t = (r - r0) / w;
    return s0 * std::exp(-0.5 * t * t);
}

/// Retained iterates in the Anderson charge mixer (dimensionless).
constexpr int kSccMixDepth = 6;

/// Σ_μν A_μν B_μν — Tr(AB) for the symmetric density/Hamiltonian pair.
double trace_product(const linalg::RealMatrix& a, const linalg::RealMatrix& b) {
    const int n = static_cast<int>(a.shape(0));
    double t = 0;
    for (int mu = 0; mu < n; ++mu)
        for (int nu = 0; nu < n; ++nu)
            t += a(mu, nu) * b(mu, nu);
    return t;
}

/// A non-converged SCC is an error, never a silently returned iterate: the
/// undamped Picard loop is a limit cycle, so "last iterate" is selected by the
/// parity of max_scc rather than by physics.
void require_scc_converged(bool converged, double residual, double scc_tol,
                           int max_scc) {
    if (converged)
        return;
    throw std::runtime_error(
        "DFTB2: SCC did not converge in " + std::to_string(max_scc)
        + " iterations (residual " + std::to_string(residual)
        + " e > scc_tol " + std::to_string(scc_tol)
        + "); raise max_scc or lower mix_beta");
}

/// γ_AB(R) = 1 / sqrt(R^2 + (0.5(1/U_A + 1/U_B))^2)  (eV, Å; Ohno–Klopman-like)
double gamma_ab(double r, double ua, double ub) {
    const double tau = 0.5 * (1.0 / ua + 1.0 / ub);
    return 1.0 / std::sqrt(r * r + tau * tau);
}

/// H–H mio cutoff [Å]; sketch pairs use element r_cut.
double pair_rcut_A(const DftbElementRow& ea, const DftbElementRow& eb, int za,
                   int zb) {
    if (za == 1 && zb == 1)
        return internal::MioHHTable{}.cutoff();
    return std::max(ea.r_cut, eb.r_cut);
}

double pair_repulsive_eV(int za, int zb, double r_A) {
    if (za == 1 && zb == 1)
        return internal::MioHHTable{}.repulsive(r_A);
    const auto& rep = pair_rep(za, zb);
    return rep.a * std::exp(-rep.b * r_A);
}

double pair_repulsive_dE_dr(int za, int zb, double r_A) {
    if (za == 1 && zb == 1)
        return internal::MioHHTable{}.d_repulsive_dr(r_A);
    const auto& rep = pair_rep(za, zb);
    return -rep.a * rep.b * std::exp(-rep.b * r_A);
}

/// @brief Write directed A→B block of H0/S for bond vector R = R_B − R_A (accumulate).
///
/// Does **not** write the transpose — lattice-sum callers walk all ordered pairs (a,b,T).
/// Min-image molecular forces keep a separate assign-both-triangles helper below.
void fill_sk_block_directed(const DftbElementRow& ea, const DftbElementRow& eb,
                            const DftbPairSkRow& sk, double dx, double dy,
                            double dz, double r, linalg::RealMatrix& H0,
                            linalg::RealMatrix& S, int ia0, int ib0) {
    if (r < 1e-12)
        return;

    // mio-1-1 H–H (s only): shared transcribed table, not the sketch radial.
    if (ea.z == 1 && eb.z == 1) {
        const internal::MioHHTable mio;
        H0(ia0, ib0) += mio.hopping_ss(r);
        S(ia0, ib0) += mio.overlap_ss(r);
        (void)sk;
        return;
    }

    const double inv_r = 1.0 / r;
    const double l = dx * inv_r;
    const double m = dy * inv_r;
    const double n = dz * inv_r;

    const double Hss = h_radial(sk.h_ss0, r, sk.h_r0, sk.h_d);
    const double Hsp = h_radial(sk.h_sp0, r, sk.h_r0, sk.h_d);
    const double Hpps = h_radial(sk.h_pps0, r, sk.h_r0, sk.h_d);
    const double Hppp = h_radial(sk.h_ppp0, r, sk.h_r0, sk.h_d);
    const double Sss = s_radial(sk.s_ss0, r, sk.s_r0, sk.s_w);
    const double Ssp = s_radial(sk.s_sp0, r, sk.s_r0, sk.s_w);
    const double Spps = s_radial(sk.s_pps0, r, sk.s_r0, sk.s_w);
    const double Sppp = s_radial(sk.s_ppp0, r, sk.s_r0, sk.s_w);

    auto add_hs = [&](int i, int j, double hv, double sv) {
        H0(ia0 + i, ib0 + j) += hv;
        S(ia0 + i, ib0 + j) += sv;
    };

    add_hs(0, 0, Hss, Sss);

    if (ea.n_ao >= 1 && eb.n_ao == 4) {
        add_hs(0, 1, l * Hsp, l * Ssp);
        add_hs(0, 2, m * Hsp, m * Ssp);
        add_hs(0, 3, n * Hsp, n * Ssp);
    }
    if (ea.n_ao == 4 && eb.n_ao >= 1) {
        add_hs(1, 0, -l * Hsp, -l * Ssp);
        add_hs(2, 0, -m * Hsp, -m * Ssp);
        add_hs(3, 0, -n * Hsp, -n * Ssp);
    }
    if (ea.n_ao == 4 && eb.n_ao == 4) {
        const double dir[3] = {l, m, n};
        for (int a = 0; a < 3; ++a) {
            for (int b = 0; b < 3; ++b) {
                const double d_ab = (a == b) ? 1.0 : 0.0;
                const double hv = dir[a] * dir[b] * (Hpps - Hppp) + d_ab * Hppp;
                const double sv = dir[a] * dir[b] * (Spps - Sppp) + d_ab * Sppp;
                add_hs(1 + a, 1 + b, hv, sv);
            }
        }
    }
}


struct LatVec {
    double tx, ty, tz;
    int n0 = 0, n1 = 0, n2 = 0; ///< integer lattice indices (for Bloch phase)
};

/// @brief Enumerate lattice translations with |T| ≲ max_cut + cell edge.
/// Columns of `cell` are lattice vectors (matches Box). Includes T = 0.
std::vector<LatVec> enumerate_lattice(const Mat3& cell, double max_cut) {
    double col_len[3];
    for (int d = 0; d < 3; ++d) {
        const double cx = static_cast<double>(cell(0, d));
        const double cy = static_cast<double>(cell(1, d));
        const double cz = static_cast<double>(cell(2, d));
        col_len[d] = std::sqrt(cx * cx + cy * cy + cz * cz);
        if (col_len[d] < 1e-9)
            throw std::runtime_error("DFTB2: degenerate cell edge");
    }
    int nmax[3];
    for (int d = 0; d < 3; ++d)
        nmax[d] = std::max(
            1, static_cast<int>(std::ceil(max_cut / col_len[d])) + 1);

    std::vector<LatVec> out;
    out.reserve(static_cast<std::size_t>((2 * nmax[0] + 1) * (2 * nmax[1] + 1)
                                         * (2 * nmax[2] + 1)));
    for (int n0 = -nmax[0]; n0 <= nmax[0]; ++n0) {
        for (int n1 = -nmax[1]; n1 <= nmax[1]; ++n1) {
            for (int n2 = -nmax[2]; n2 <= nmax[2]; ++n2) {
                const double tx = n0 * static_cast<double>(cell(0, 0))
                                  + n1 * static_cast<double>(cell(0, 1))
                                  + n2 * static_cast<double>(cell(0, 2));
                const double ty = n0 * static_cast<double>(cell(1, 0))
                                  + n1 * static_cast<double>(cell(1, 1))
                                  + n2 * static_cast<double>(cell(1, 2));
                const double tz = n0 * static_cast<double>(cell(2, 0))
                                  + n1 * static_cast<double>(cell(2, 1))
                                  + n2 * static_cast<double>(cell(2, 2));
                out.push_back({tx, ty, tz, n0, n1, n2});
            }
        }
    }
    return out;
}

std::vector<LatVec>
lattice_images_for(const cpu::Box& box,
                   const std::vector<const DftbElementRow*>& el) {
    double max_cut = 0;
    for (const auto* e : el)
        max_cut = std::max(max_cut, e->r_cut);
    // mio-1-1 H–H table may extend beyond sketch r_cut.
    max_cut = std::max(max_cut, internal::MioHHTable{}.cutoff());
    if (max_cut < 1e-12)
        max_cut = 1.0;
    return enumerate_lattice(box.cell(), max_cut);
}

/// dγ/dR = −R γ³.
double dgamma_dr(double r, double ua, double ub) {
    const double g = gamma_ab(r, ua, ub);
    return -r * g * g * g;
}

/// SK block derivatives w.r.t. bond vector components (R_B − R_A).
void fill_sk_block_deriv(const DftbElementRow& ea, const DftbElementRow& eb,
                         const DftbPairSkRow& sk, double dx, double dy,
                         double dz, double r, linalg::RealMatrix dH0[3],
                         linalg::RealMatrix dS[3], int ia0, int ib0) {
    if (r < 1e-12)
        return;
    const double inv_r = 1.0 / r;
    const double dir[3] = {dx * inv_r, dy * inv_r, dz * inv_r};

    // mio-1-1 H–H: radial derivative of the shared table, chained to bond vec.
    if (ea.z == 1 && eb.z == 1) {
        const internal::MioHHTable mio;
        const double dHss = mio.d_hopping_ss_dr(r);
        const double dSss = mio.d_overlap_ss_dr(r);
        for (int k = 0; k < 3; ++k) {
            dH0[k](ia0, ib0) = dHss * dir[k];
            dS[k](ia0, ib0) = dSss * dir[k];
        }
        (void)sk;
        return;
    }

    const double Hss = h_radial(sk.h_ss0, r, sk.h_r0, sk.h_d);
    const double Hsp = h_radial(sk.h_sp0, r, sk.h_r0, sk.h_d);
    const double Hpps = h_radial(sk.h_pps0, r, sk.h_r0, sk.h_d);
    const double Hppp = h_radial(sk.h_ppp0, r, sk.h_r0, sk.h_d);
    const double Sss = s_radial(sk.s_ss0, r, sk.s_r0, sk.s_w);
    const double Ssp = s_radial(sk.s_sp0, r, sk.s_r0, sk.s_w);
    const double Spps = s_radial(sk.s_pps0, r, sk.s_r0, sk.s_w);
    const double Sppp = s_radial(sk.s_ppp0, r, sk.s_r0, sk.s_w);

    const double dHss = -Hss / sk.h_d;
    const double dHsp = -Hsp / sk.h_d;
    const double dHpps = -Hpps / sk.h_d;
    const double dHppp = -Hppp / sk.h_d;
    const double dSss = Sss * (-(r - sk.s_r0) / (sk.s_w * sk.s_w));
    const double dSsp = Ssp * (-(r - sk.s_r0) / (sk.s_w * sk.s_w));
    const double dSpps = Spps * (-(r - sk.s_r0) / (sk.s_w * sk.s_w));
    const double dSppp = Sppp * (-(r - sk.s_r0) / (sk.s_w * sk.s_w));

    auto ddir = [&](int alpha, int k) {
        const double dab = (alpha == k) ? 1.0 : 0.0;
        return (dab - dir[alpha] * dir[k]) * inv_r;
    };

    for (int k = 0; k < 3; ++k) {
        const double dr_dk = dir[k];
        dH0[k](ia0, ib0) = dHss * dr_dk;
        dS[k](ia0, ib0) = dSss * dr_dk;

        if (ea.n_ao >= 1 && eb.n_ao == 4) {
            for (int a = 0; a < 3; ++a) {
                dH0[k](ia0, ib0 + 1 + a) =
                    ddir(a, k) * Hsp + dir[a] * dHsp * dr_dk;
                dS[k](ia0, ib0 + 1 + a) =
                    ddir(a, k) * Ssp + dir[a] * dSsp * dr_dk;
            }
        }
        if (ea.n_ao == 4 && eb.n_ao >= 1) {
            for (int a = 0; a < 3; ++a) {
                dH0[k](ia0 + 1 + a, ib0) =
                    -(ddir(a, k) * Hsp + dir[a] * dHsp * dr_dk);
                dS[k](ia0 + 1 + a, ib0) =
                    -(ddir(a, k) * Ssp + dir[a] * dSsp * dr_dk);
            }
        }
        if (ea.n_ao == 4 && eb.n_ao == 4) {
            const double dHsig = dHpps - dHppp;
            const double dSsig = dSpps - dSppp;
            const double Hsig = Hpps - Hppp;
            const double Ssig = Spps - Sppp;
            for (int a = 0; a < 3; ++a) {
                for (int b = 0; b < 3; ++b) {
                    const double dab = (a == b) ? 1.0 : 0.0;
                    dH0[k](ia0 + 1 + a, ib0 + 1 + b) =
                        (ddir(a, k) * dir[b] + dir[a] * ddir(b, k)) * Hsig
                        + dir[a] * dir[b] * dHsig * dr_dk + dab * dHppp * dr_dk;
                    dS[k](ia0 + 1 + a, ib0 + 1 + b) =
                        (ddir(a, k) * dir[b] + dir[a] * ddir(b, k)) * Ssig
                        + dir[a] * dir[b] * dSsig * dr_dk + dab * dSppp * dr_dk;
                }
            }
        }
    }
}

struct Occupancy {
    int n_alpha = 0;
    int n_beta = 0;
    bool unrestricted = false;
};

Occupancy resolve_occupancy(double n_elec, int multiplicity) {
    if (multiplicity < 1)
        throw std::runtime_error("DFTB2: multiplicity must be >= 1");
    if (std::fabs(n_elec - std::round(n_elec)) > 1e-9)
        throw std::runtime_error("DFTB2: non-integer valence electron count");
    const int n_e = static_cast<int>(std::lround(n_elec));
    if (n_e < 0)
        throw std::runtime_error("DFTB2: negative electron count after charge");
    const int n_unpaired = multiplicity - 1;
    if ((n_e - n_unpaired) % 2 != 0 || n_e < n_unpaired)
        throw std::runtime_error(
            "DFTB2: charge/multiplicity inconsistent with valence electron count");
    Occupancy occ;
    occ.n_alpha = (n_e + n_unpaired) / 2;
    occ.n_beta = (n_e - n_unpaired) / 2;
    occ.unrestricted = (multiplicity > 1) || (n_e % 2 != 0);
    return occ;
}

struct SccResult {
    double energy = 0;
    double residual = 0;
    int iterations = 0;
    std::vector<double> dq;
    std::vector<double> spin; // m_A = Nα − Nβ
    // Restricted-path orbitals for closed-form forces.
    linalg::RealMatrix C;
    linalg::RealVector eps;
    int n_occ = 0;
    // Unrestricted-path orbitals (α/β) for collinear open-shell analytic F.
    linalg::RealMatrix Ca;
    linalg::RealMatrix Cb;
    linalg::RealVector eps_a;
    linalg::RealVector eps_b;
    int n_alpha = 0;
    int n_beta = 0;
    bool unrestricted = false;
};

/// Build lattice-summed H0/S at the Γ-point (molecular limit when L ≫ r_cut).
void build_h0_s(int n_at, int n_ao, const std::vector<const DftbElementRow*>& el,
                const std::vector<int>& ao0, const int* type,
                const StorageReal* x, const StorageReal* y, const StorageReal* z,
                const cpu::Box& box, linalg::RealMatrix& H0,
                linalg::RealMatrix& S) {
    H0 = xt::zeros<double>({n_ao, n_ao});
    S = xt::zeros<double>({n_ao, n_ao});
    const auto images = lattice_images_for(box, el);

    // Onsite (T = 0, a = b): diagonal e_s / e_p and S = 1.
    for (int a = 0; a < n_at; ++a) {
        const auto& ea = *el[static_cast<std::size_t>(a)];
        const int ia = ao0[static_cast<std::size_t>(a)];
        H0(ia, ia) = ea.e_s;
        S(ia, ia) = 1.0;
        if (ea.n_ao == 4) {
            for (int p = 1; p < 4; ++p) {
                H0(ia + p, ia + p) = ea.e_p;
                S(ia + p, ia + p) = 1.0;
            }
        }
    }

    // Lattice sum over ordered pairs (a,b) and images T of B.
    for (int a = 0; a < n_at; ++a) {
        const auto& ea = *el[static_cast<std::size_t>(a)];
        const int ia = ao0[static_cast<std::size_t>(a)];
        const double xa = static_cast<double>(x[a]);
        const double ya = static_cast<double>(y[a]);
        const double za = static_cast<double>(z[a]);
        for (int b = 0; b < n_at; ++b) {
            const auto& eb = *el[static_cast<std::size_t>(b)];
            const int ib = ao0[static_cast<std::size_t>(b)];
            const double xb = static_cast<double>(x[b]);
            const double yb = static_cast<double>(y[b]);
            const double zb = static_cast<double>(z[b]);
            const double rcut = pair_rcut_A(ea, eb, type[a], type[b]);
            const auto& sk = pair_sk(type[a], type[b]);
            for (const auto& T : images) {
                const double dx = xb + T.tx - xa;
                const double dy = yb + T.ty - ya;
                const double dz = zb + T.tz - za;
                const double r = std::sqrt(dx * dx + dy * dy + dz * dz);
                if (r < 1e-12 || r > rcut)
                    continue;
                fill_sk_block_directed(ea, eb, sk, dx, dy, dz, r, H0, S, ia,
                                       ib);
            }
        }
    }
}

double repulsive_energy(int n_at, const int* type, const StorageReal* x,
                        const StorageReal* y, const StorageReal* z,
                        const cpu::Box& box,
                        const std::vector<const DftbElementRow*>& el) {
    const auto images = lattice_images_for(box, el);
    double e_rep = 0;
    for (int a = 0; a < n_at; ++a) {
        const auto& ea = *el[static_cast<std::size_t>(a)];
        const double xa = static_cast<double>(x[a]);
        const double ya = static_cast<double>(y[a]);
        const double za = static_cast<double>(z[a]);
        for (int b = 0; b < n_at; ++b) {
            const auto& eb = *el[static_cast<std::size_t>(b)];
            const double xb = static_cast<double>(x[b]);
            const double yb = static_cast<double>(y[b]);
            const double zb = static_cast<double>(z[b]);
            const double rcut = pair_rcut_A(ea, eb, type[a], type[b]);
            for (const auto& T : images) {
                const double dx = xb + T.tx - xa;
                const double dy = yb + T.ty - ya;
                const double dz = zb + T.tz - za;
                const double r = std::sqrt(dx * dx + dy * dy + dz * dz);
                if (r < 1e-12 || r > rcut)
                    continue;
                // 1/2 for double-counting over ordered (a,b,T).
                e_rep += 0.5 * pair_repulsive_eV(type[a], type[b], r);
            }
        }
    }
    return e_rep;
}

void apply_shift(int n_at, const std::vector<const DftbElementRow*>& el,
                 const std::vector<int>& ao0, const linalg::RealMatrix& H0,
                 const linalg::RealMatrix& S, const std::vector<double>& eps_at,
                 linalg::RealMatrix& H) {
    H = H0;
    for (int a = 0; a < n_at; ++a) {
        const int ia = ao0[static_cast<std::size_t>(a)];
        const int na = el[static_cast<std::size_t>(a)]->n_ao;
        for (int b = 0; b < n_at; ++b) {
            const int ib = ao0[static_cast<std::size_t>(b)];
            const int nb = el[static_cast<std::size_t>(b)]->n_ao;
            const double half = 0.5
                                * (eps_at[static_cast<std::size_t>(a)]
                                   + eps_at[static_cast<std::size_t>(b)]);
            for (int i = 0; i < na; ++i)
                for (int j = 0; j < nb; ++j)
                    H(ia + i, ib + j) += half * S(ia + i, ib + j);
        }
    }
}

double gamma_pair(int a, int b, const std::vector<const DftbElementRow*>& el,
                  const StorageReal* x, const StorageReal* y,
                  const StorageReal* z, const cpu::Box& box) {
    const double ua = el[static_cast<std::size_t>(a)]->hubbard_u;
    const double ub = el[static_cast<std::size_t>(b)]->hubbard_u;
    const double rcut = pair_rcut_A(*el[static_cast<std::size_t>(a)],
                                    *el[static_cast<std::size_t>(b)],
                                    el[static_cast<std::size_t>(a)]->z,
                                    el[static_cast<std::size_t>(b)]->z);
    const auto images = lattice_images_for(box, el);
    const double xa = static_cast<double>(x[a]);
    const double ya = static_cast<double>(y[a]);
    const double za = static_cast<double>(z[a]);
    const double xb = static_cast<double>(x[b]);
    const double yb = static_cast<double>(y[b]);
    const double zb = static_cast<double>(z[b]);
    double g = 0;
    for (const auto& T : images) {
        const double dx = xb + T.tx - xa;
        const double dy = yb + T.ty - ya;
        const double dz = zb + T.tz - za;
        const double r = std::sqrt(dx * dx + dy * dy + dz * dz);
        if (r < 1e-12) {
            if (a == b)
                g += ua; // onsite Hubbard
            continue;
        }
        if (r > rcut)
            continue;
        g += gamma_ab(r, ua, ub);
    }
    return g;
}

SccResult scc_energy_restricted(
    int n_at, int n_ao, int n_occ, const std::vector<const DftbElementRow*>& el,
    const std::vector<int>& ao0, const int* type, const StorageReal* x,
    const StorageReal* y, const StorageReal* z, const cpu::Box& box,
    double scc_tol, int max_scc, double mix_beta) {
    linalg::RealMatrix H0, S;
    build_h0_s(n_at, n_ao, el, ao0, type, x, y, z, box, H0, S);

    const std::size_t n = static_cast<std::size_t>(n_at);
    std::vector<double> dq(n, 0.0);
    std::vector<double> dq_out(n, 0.0);
    linalg::RealMatrix P = xt::zeros<double>({n_ao, n_ao});
    linalg::RealMatrix C = xt::zeros<double>({n_ao, n_ao});
    linalg::RealVector eps = xt::zeros<double>({n_ao});
    double residual = 0;
    int iterations = 0;
    bool converged = false;
    internal::ChargeMixer mixer(mix_beta, kSccMixDepth);

    for (int it = 0; it < max_scc; ++it) {
        iterations = it + 1;
        std::vector<double> dV(n, 0.0);
        for (int a = 0; a < n_at; ++a) {
            double shift = 0;
            for (int b = 0; b < n_at; ++b)
                shift += gamma_pair(a, b, el, x, y, z, box)
                         * dq[static_cast<std::size_t>(b)];
            dV[static_cast<std::size_t>(a)] = shift;
        }
        linalg::RealMatrix H;
        apply_shift(n_at, el, ao0, H0, S, dV, H);
        auto sol = linalg::diagonalize_fock_generalized(H, S);
        eps = sol.first;
        C = sol.second;

        P = xt::zeros<double>({n_ao, n_ao});
        for (int o = 0; o < n_occ; ++o)
            for (int mu = 0; mu < n_ao; ++mu)
                for (int nu = 0; nu < n_ao; ++nu)
                    P(mu, nu) += 2.0 * C(mu, o) * C(nu, o);

        linalg::RealMatrix PS = xt::linalg::dot(P, S);
        residual = 0;
        for (int a = 0; a < n_at; ++a) {
            double qa = 0;
            const int ia = ao0[static_cast<std::size_t>(a)];
            const int na = el[static_cast<std::size_t>(a)]->n_ao;
            for (int i = 0; i < na; ++i)
                qa += PS(ia + i, ia + i);
            dq_out[static_cast<std::size_t>(a)] =
                qa - el[static_cast<std::size_t>(a)]->n_val;
            residual = std::max(residual,
                                std::abs(dq_out[static_cast<std::size_t>(a)]
                                         - dq[static_cast<std::size_t>(a)]));
        }

        if (residual < scc_tol) {
            dq = dq_out;
            converged = true;
            break;
        }
        dq = mixer.next(dq, dq_out);
    }

    require_scc_converged(converged, residual, scc_tol, max_scc);

    double e_scc = 0;
    for (int a = 0; a < n_at; ++a)
        for (int b = 0; b < n_at; ++b)
            e_scc += 0.5 * gamma_pair(a, b, el, x, y, z, box)
                     * dq[static_cast<std::size_t>(a)]
                     * dq[static_cast<std::size_t>(b)];

    SccResult out;
    // E = Tr(ρH⁰) + E² + E_rep. NOT E_band − E²: the band energy is
    // Tr(ρH⁰) + Σ_A V_A q_A with the GROSS Mulliken population, so that route
    // leaves a spurious Σ_A V_A q⁰_A which is linear in Δq.
    out.energy = trace_product(P, H0) + e_scc
                 + repulsive_energy(n_at, type, x, y, z, box, el);
    out.residual = residual;
    out.iterations = iterations;
    out.dq = std::move(dq);
    out.spin.assign(static_cast<std::size_t>(n_at), 0.0);
    out.C = std::move(C);
    out.eps = std::move(eps);
    out.n_occ = n_occ;
    return out;
}

// ── Multi-k (Bloch) restricted SCC ──────────────────────────────────────

using Cplx = std::complex<double>;
using CplxMatrix = linalg::CplxMatrix;

struct KPoint {
    double frac[3] = {0, 0, 0};
    double weight = 1.0;
};

/// CP2K / DFTB+ Monkhorst–Pack: κ_d = (2r − N − 1)/(2N), r = 1..N.
std::vector<KPoint> monkhorst_pack_mesh(const int mesh[3], bool time_reversal) {
    for (int d = 0; d < 3; ++d)
        if (mesh[d] < 1)
            throw std::invalid_argument(
                "DFTB2: k_mesh subdivision must be >= 1 (axis "
                + std::to_string(d) + " got " + std::to_string(mesh[d]) + ")");

    std::vector<KPoint> full;
    const double w0 =
        1.0 / (static_cast<double>(mesh[0]) * mesh[1] * mesh[2]);
    for (int r0 = 1; r0 <= mesh[0]; ++r0)
        for (int r1 = 1; r1 <= mesh[1]; ++r1)
            for (int r2 = 1; r2 <= mesh[2]; ++r2) {
                KPoint kp;
                const int r[3] = {r0, r1, r2};
                for (int d = 0; d < 3; ++d) {
                    double f = (2.0 * r[d] - mesh[d] - 1.0)
                               / (2.0 * static_cast<double>(mesh[d]));
                    f -= std::floor(f + 0.5);
                    if (f <= -0.5)
                        f += 1.0;
                    kp.frac[d] = f;
                }
                kp.weight = w0;
                full.push_back(kp);
            }

    if (!time_reversal)
        return full;

    std::vector<KPoint> folded;
    std::vector<bool> used(full.size(), false);
    constexpr double tol = 1e-12;
    auto same = [&](const double a[3], const double b[3]) {
        for (int d = 0; d < 3; ++d) {
            double diff = a[d] - b[d];
            diff -= std::nearbyint(diff);
            if (std::abs(diff) > tol)
                return false;
        }
        return true;
    };
    for (std::size_t i = 0; i < full.size(); ++i) {
        if (used[i])
            continue;
        used[i] = true;
        KPoint kp = full[i];
        const double neg[3] = {-kp.frac[0], -kp.frac[1], -kp.frac[2]};
        for (std::size_t j = i + 1; j < full.size(); ++j) {
            if (used[j])
                continue;
            if (same(full[j].frac, neg)) {
                used[j] = true;
                kp.weight += full[j].weight;
                break;
            }
        }
        folded.push_back(kp);
    }
    return folded;
}

Cplx bloch_phase(const double kappa[3], int n0, int n1, int n2) {
    constexpr double kTwoPi = 6.28318530717958647692;
    const double ang =
        kTwoPi * (kappa[0] * n0 + kappa[1] * n1 + kappa[2] * n2);
    return Cplx(std::cos(ang), std::sin(ang));
}

/// Row-major complex product C = A B.
CplxMatrix cplx_mul(const CplxMatrix& A, const CplxMatrix& B) {
    const int n = static_cast<int>(A.shape(0));
    const int m = static_cast<int>(A.shape(1));
    const int p = static_cast<int>(B.shape(1));
    CplxMatrix C = xt::zeros<Cplx>({n, p});
    for (int i = 0; i < n; ++i)
        for (int k = 0; k < m; ++k) {
            const Cplx a = A(i, k);
            if (a == Cplx(0, 0))
                continue;
            for (int j = 0; j < p; ++j)
                C(i, j) += a * B(k, j);
        }
    return C;
}

/// Write directed real SK block scaled by complex phase into H0/S.
void fill_sk_block_phased(const DftbElementRow& ea, const DftbElementRow& eb,
                          const DftbPairSkRow& sk, double dx, double dy,
                          double dz, double r, const Cplx& phase,
                          CplxMatrix& H0, CplxMatrix& S, int ia0, int ib0) {
    // Stage into a tiny real scratch via the directed filler, then phase.
    // n_ao ≤ 4 per center; use full matrices only for the block region.
    if (r < 1e-12)
        return;

    auto add = [&](int i, int j, double hv, double sv) {
        H0(ia0 + i, ib0 + j) += phase * hv;
        S(ia0 + i, ib0 + j) += phase * sv;
    };

    if (ea.z == 1 && eb.z == 1) {
        const internal::MioHHTable mio;
        add(0, 0, mio.hopping_ss(r), mio.overlap_ss(r));
        (void)sk;
        return;
    }

    const double inv_r = 1.0 / r;
    const double l = dx * inv_r;
    const double m = dy * inv_r;
    const double n = dz * inv_r;

    const double Hss = h_radial(sk.h_ss0, r, sk.h_r0, sk.h_d);
    const double Hsp = h_radial(sk.h_sp0, r, sk.h_r0, sk.h_d);
    const double Hpps = h_radial(sk.h_pps0, r, sk.h_r0, sk.h_d);
    const double Hppp = h_radial(sk.h_ppp0, r, sk.h_r0, sk.h_d);
    const double Sss = s_radial(sk.s_ss0, r, sk.s_r0, sk.s_w);
    const double Ssp = s_radial(sk.s_sp0, r, sk.s_r0, sk.s_w);
    const double Spps = s_radial(sk.s_pps0, r, sk.s_r0, sk.s_w);
    const double Sppp = s_radial(sk.s_ppp0, r, sk.s_r0, sk.s_w);

    add(0, 0, Hss, Sss);
    if (ea.n_ao >= 1 && eb.n_ao == 4) {
        add(0, 1, l * Hsp, l * Ssp);
        add(0, 2, m * Hsp, m * Ssp);
        add(0, 3, n * Hsp, n * Ssp);
    }
    if (ea.n_ao == 4 && eb.n_ao >= 1) {
        add(1, 0, -l * Hsp, -l * Ssp);
        add(2, 0, -m * Hsp, -m * Ssp);
        add(3, 0, -n * Hsp, -n * Ssp);
    }
    if (ea.n_ao == 4 && eb.n_ao == 4) {
        const double dir[3] = {l, m, n};
        for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b) {
                const double d_ab = (a == b) ? 1.0 : 0.0;
                const double hv = dir[a] * dir[b] * (Hpps - Hppp) + d_ab * Hppp;
                const double sv = dir[a] * dir[b] * (Spps - Sppp) + d_ab * Sppp;
                add(1 + a, 1 + b, hv, sv);
            }
    }
}

void hermitize_in_place(CplxMatrix& M) {
    const int n = static_cast<int>(M.shape(0));
    for (int i = 0; i < n; ++i) {
        M(i, i) = Cplx(M(i, i).real(), 0.0);
        for (int j = i + 1; j < n; ++j) {
            const Cplx avg = 0.5 * (M(i, j) + std::conj(M(j, i)));
            M(i, j) = avg;
            M(j, i) = std::conj(avg);
        }
    }
}

/// Build lattice-summed complex H0(k)/S(k) for one fractional k-point.
void build_h0_s_k(int n_at, int n_ao,
                  const std::vector<const DftbElementRow*>& el,
                  const std::vector<int>& ao0, const int* type,
                  const StorageReal* x, const StorageReal* y,
                  const StorageReal* z, const cpu::Box& box,
                  const double kappa[3], CplxMatrix& H0, CplxMatrix& S) {
    H0 = xt::zeros<Cplx>({n_ao, n_ao});
    S = xt::zeros<Cplx>({n_ao, n_ao});
    const auto images = lattice_images_for(box, el);

    for (int a = 0; a < n_at; ++a) {
        const auto& ea = *el[static_cast<std::size_t>(a)];
        const int ia = ao0[static_cast<std::size_t>(a)];
        H0(ia, ia) = Cplx(ea.e_s, 0.0);
        S(ia, ia) = Cplx(1.0, 0.0);
        if (ea.n_ao == 4) {
            for (int p = 1; p < 4; ++p) {
                H0(ia + p, ia + p) = Cplx(ea.e_p, 0.0);
                S(ia + p, ia + p) = Cplx(1.0, 0.0);
            }
        }
    }

    for (int a = 0; a < n_at; ++a) {
        const auto& ea = *el[static_cast<std::size_t>(a)];
        const int ia = ao0[static_cast<std::size_t>(a)];
        const double xa = static_cast<double>(x[a]);
        const double ya = static_cast<double>(y[a]);
        const double za = static_cast<double>(z[a]);
        for (int b = 0; b < n_at; ++b) {
            const auto& eb = *el[static_cast<std::size_t>(b)];
            const int ib = ao0[static_cast<std::size_t>(b)];
            const double xb = static_cast<double>(x[b]);
            const double yb = static_cast<double>(y[b]);
            const double zb = static_cast<double>(z[b]);
            const double rcut = pair_rcut_A(ea, eb, type[a], type[b]);
            const auto& sk = pair_sk(type[a], type[b]);
            for (const auto& T : images) {
                const double dx = xb + T.tx - xa;
                const double dy = yb + T.ty - ya;
                const double dz = zb + T.tz - za;
                const double r = std::sqrt(dx * dx + dy * dy + dz * dz);
                if (r < 1e-12 || r > rcut)
                    continue;
                const Cplx phase = bloch_phase(kappa, T.n0, T.n1, T.n2);
                fill_sk_block_phased(ea, eb, sk, dx, dy, dz, r, phase, H0, S,
                                     ia, ib);
            }
        }
    }
    hermitize_in_place(H0);
    hermitize_in_place(S);
}

void apply_shift_k(int n_at, const std::vector<const DftbElementRow*>& el,
                   const std::vector<int>& ao0, const CplxMatrix& H0,
                   const CplxMatrix& S, const std::vector<double>& eps_at,
                   CplxMatrix& H) {
    H = H0;
    for (int a = 0; a < n_at; ++a) {
        const int ia = ao0[static_cast<std::size_t>(a)];
        const int na = el[static_cast<std::size_t>(a)]->n_ao;
        for (int b = 0; b < n_at; ++b) {
            const int ib = ao0[static_cast<std::size_t>(b)];
            const int nb = el[static_cast<std::size_t>(b)]->n_ao;
            const double half = 0.5
                                * (eps_at[static_cast<std::size_t>(a)]
                                   + eps_at[static_cast<std::size_t>(b)]);
            for (int i = 0; i < na; ++i)
                for (int j = 0; j < nb; ++j)
                    H(ia + i, ib + j) += half * S(ia + i, ib + j);
        }
    }
}

/// Generalized Hermitian eigen: H C = S C ε via Löwdin S^{-1/2} + zheevd.
std::pair<std::vector<double>, CplxMatrix>
diagonalize_fock_generalized_k(const CplxMatrix& H, const CplxMatrix& S) {
    const int n = static_cast<int>(H.shape(0));
    std::vector<double> s_eig;
    CplxMatrix Us;
    linalg::eigh_hermitian(S, s_eig, Us);
    // X = U s^{-1/2}
    CplxMatrix X = xt::zeros<Cplx>({n, n});
    for (int j = 0; j < n; ++j) {
        if (s_eig[static_cast<std::size_t>(j)] <= 1e-12)
            throw std::runtime_error(
                "DFTB2 multi-k: non-positive overlap eigenvalue at k");
        const double inv_sqrt = 1.0 / std::sqrt(s_eig[static_cast<std::size_t>(j)]);
        for (int i = 0; i < n; ++i)
            X(i, j) = Us(i, j) * inv_sqrt;
    }
    // H' = X† H X
    CplxMatrix HX = cplx_mul(H, X);
    CplxMatrix Hp = xt::zeros<Cplx>({n, n});
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            Cplx s(0, 0);
            for (int k = 0; k < n; ++k)
                s += std::conj(X(k, i)) * HX(k, j);
            Hp(i, j) = s;
        }
    hermitize_in_place(Hp);
    std::vector<double> eps;
    CplxMatrix Cp;
    linalg::eigh_hermitian(Hp, eps, Cp);
    return {std::move(eps), cplx_mul(X, Cp)};
}

double re_trace_product_k(const CplxMatrix& P, const CplxMatrix& H0) {
    const int n = static_cast<int>(P.shape(0));
    double t = 0;
    for (int mu = 0; mu < n; ++mu)
        for (int nu = 0; nu < n; ++nu)
            t += (P(mu, nu) * H0(nu, mu)).real();
    return t;
}

SccResult scc_energy_restricted_kmesh(
    int n_at, int n_ao, int n_occ, const std::vector<const DftbElementRow*>& el,
    const std::vector<int>& ao0, const int* type, const StorageReal* x,
    const StorageReal* y, const StorageReal* z, const cpu::Box& box,
    double scc_tol, int max_scc, double mix_beta, const int mesh[3],
    bool time_reversal) {
    if (n_occ > n_ao)
        throw std::runtime_error(
            "DFTB2 multi-k: n_occ exceeds basis size (no fractional fill)");

    const std::vector<KPoint> kpts = monkhorst_pack_mesh(mesh, time_reversal);
    const std::size_t nk = kpts.size();
    std::vector<CplxMatrix> H0k(nk), Sk(nk);
    for (std::size_t ik = 0; ik < nk; ++ik)
        build_h0_s_k(n_at, n_ao, el, ao0, type, x, y, z, box, kpts[ik].frac,
                     H0k[ik], Sk[ik]);

    const std::size_t n = static_cast<std::size_t>(n_at);
    std::vector<double> dq(n, 0.0);
    std::vector<double> dq_out(n, 0.0);
    double residual = 0;
    int iterations = 0;
    bool converged = false;
    internal::ChargeMixer mixer(mix_beta, kSccMixDepth);

    // Last-iterate densities for the energy assembly after SCC.
    std::vector<CplxMatrix> Pk(nk);

    for (int it = 0; it < max_scc; ++it) {
        iterations = it + 1;
        std::vector<double> dV(n, 0.0);
        for (int a = 0; a < n_at; ++a) {
            double shift = 0;
            for (int b = 0; b < n_at; ++b)
                shift += gamma_pair(a, b, el, x, y, z, box)
                         * dq[static_cast<std::size_t>(b)];
            dV[static_cast<std::size_t>(a)] = shift;
        }

        // Weighted Mulliken populations across the mesh.
        std::vector<double> q_ao(static_cast<std::size_t>(n_ao), 0.0);
        for (std::size_t ik = 0; ik < nk; ++ik) {
            CplxMatrix Hk;
            apply_shift_k(n_at, el, ao0, H0k[ik], Sk[ik], dV, Hk);
            auto sol = diagonalize_fock_generalized_k(Hk, Sk[ik]);
            const auto& eps = sol.first;
            const CplxMatrix& C = sol.second;
            (void)eps;

            CplxMatrix P = xt::zeros<Cplx>({n_ao, n_ao});
            for (int o = 0; o < n_occ; ++o)
                for (int mu = 0; mu < n_ao; ++mu)
                    for (int nu = 0; nu < n_ao; ++nu)
                        P(mu, nu) +=
                            2.0 * C(mu, o) * std::conj(C(nu, o));
            Pk[ik] = P;

            // PS = P S; Mulliken AO pop = Re((PS)_μμ)
            const CplxMatrix PS = cplx_mul(P, Sk[ik]);
            const double w = kpts[ik].weight;
            for (int mu = 0; mu < n_ao; ++mu)
                q_ao[static_cast<std::size_t>(mu)] += w * PS(mu, mu).real();
        }

        residual = 0;
        for (int a = 0; a < n_at; ++a) {
            double qa = 0;
            const int ia = ao0[static_cast<std::size_t>(a)];
            const int na = el[static_cast<std::size_t>(a)]->n_ao;
            for (int i = 0; i < na; ++i)
                qa += q_ao[static_cast<std::size_t>(ia + i)];
            dq_out[static_cast<std::size_t>(a)] =
                qa - el[static_cast<std::size_t>(a)]->n_val;
            residual = std::max(residual,
                                std::abs(dq_out[static_cast<std::size_t>(a)]
                                         - dq[static_cast<std::size_t>(a)]));
        }

        if (residual < scc_tol) {
            dq = dq_out;
            converged = true;
            break;
        }
        dq = mixer.next(dq, dq_out);
    }

    require_scc_converged(converged, residual, scc_tol, max_scc);

    // Rebuild final densities from the converged dq (last Pk is from the
    // residual check with dq_out; re-evaluate once for energy consistency).
    {
        std::vector<double> dV(n, 0.0);
        for (int a = 0; a < n_at; ++a) {
            double shift = 0;
            for (int b = 0; b < n_at; ++b)
                shift += gamma_pair(a, b, el, x, y, z, box)
                         * dq[static_cast<std::size_t>(b)];
            dV[static_cast<std::size_t>(a)] = shift;
        }
        for (std::size_t ik = 0; ik < nk; ++ik) {
            CplxMatrix Hk;
            apply_shift_k(n_at, el, ao0, H0k[ik], Sk[ik], dV, Hk);
            auto sol = diagonalize_fock_generalized_k(Hk, Sk[ik]);
            const CplxMatrix& C = sol.second;
            CplxMatrix P = xt::zeros<Cplx>({n_ao, n_ao});
            for (int o = 0; o < n_occ; ++o)
                for (int mu = 0; mu < n_ao; ++mu)
                    for (int nu = 0; nu < n_ao; ++nu)
                        P(mu, nu) +=
                            2.0 * C(mu, o) * std::conj(C(nu, o));
            Pk[ik] = std::move(P);
        }
    }

    double e_band = 0;
    for (std::size_t ik = 0; ik < nk; ++ik)
        e_band +=
            kpts[ik].weight * re_trace_product_k(Pk[ik], H0k[ik]);

    double e_scc = 0;
    for (int a = 0; a < n_at; ++a)
        for (int b = 0; b < n_at; ++b)
            e_scc += 0.5 * gamma_pair(a, b, el, x, y, z, box)
                     * dq[static_cast<std::size_t>(a)]
                     * dq[static_cast<std::size_t>(b)];

    SccResult out;
    out.energy = e_band + e_scc
                 + repulsive_energy(n_at, type, x, y, z, box, el);
    out.residual = residual;
    out.iterations = iterations;
    out.dq = std::move(dq);
    out.spin.assign(static_cast<std::size_t>(n_at), 0.0);
    out.n_occ = n_occ;
    // No real Γ orbitals for multi-k analytic forces (forces unsupported).
    return out;
}

SccResult scc_energy_unrestricted(
    int n_at, int n_ao, int n_alpha, int n_beta,
    const std::vector<const DftbElementRow*>& el, const std::vector<int>& ao0,
    const int* type, const StorageReal* x, const StorageReal* y,
    const StorageReal* z, const cpu::Box& box, double scc_tol, int max_scc,
    double mix_beta) {
    if (n_alpha > n_ao || n_beta > n_ao)
        throw std::runtime_error(
            "DFTB2: n_alpha/n_beta exceeds basis size");

    linalg::RealMatrix H0, S;
    build_h0_s(n_at, n_ao, el, ao0, type, x, y, z, box, H0, S);

    const std::size_t n = static_cast<std::size_t>(n_at);
    std::vector<double> dq(n, 0.0);
    std::vector<double> mspin(n, 0.0);
    std::vector<double> dq_out(n, 0.0);
    std::vector<double> mspin_out(n, 0.0);
    // Seed high-spin magnetization so the first iteration is not spin-null.
    if (n_alpha != n_beta) {
        const double m_tot = static_cast<double>(n_alpha - n_beta);
        for (int a = 0; a < n_at; ++a)
            mspin[static_cast<std::size_t>(a)] =
                m_tot / static_cast<double>(n_at);
    }

    double e_scc = 0, e_spin = 0, residual = 0;
    int iterations = 0;
    bool converged = false;
    // dq and m are coupled by the same density: mix them as one vector.
    internal::ChargeMixer mixer(mix_beta, kSccMixDepth);
    linalg::RealMatrix Pa = xt::zeros<double>({n_ao, n_ao});
    linalg::RealMatrix Pb = xt::zeros<double>({n_ao, n_ao});
    linalg::RealMatrix CA = xt::zeros<double>({n_ao, n_ao});
    linalg::RealMatrix CB = xt::zeros<double>({n_ao, n_ao});
    linalg::RealVector epsA = xt::zeros<double>({n_ao});
    linalg::RealVector epsB = xt::zeros<double>({n_ao});

    for (int it = 0; it < max_scc; ++it) {
        iterations = it + 1;
        std::vector<double> dVc(static_cast<std::size_t>(n_at), 0.0);
        std::vector<double> dVs(static_cast<std::size_t>(n_at), 0.0);
        for (int a = 0; a < n_at; ++a) {
            double shift = 0;
            for (int b = 0; b < n_at; ++b)
                shift += gamma_pair(a, b, el, x, y, z, box)
                         * dq[static_cast<std::size_t>(b)];
            dVc[static_cast<std::size_t>(a)] = shift;
            dVs[static_cast<std::size_t>(a)] =
                el[static_cast<std::size_t>(a)]->spin_w
                * mspin[static_cast<std::size_t>(a)];
        }

        // ε^σ_A = dVc + σ dVs  (σ = +1 α, −1 β)
        std::vector<double> eps_shift_a(static_cast<std::size_t>(n_at));
        std::vector<double> eps_shift_b(static_cast<std::size_t>(n_at));
        for (int a = 0; a < n_at; ++a) {
            eps_shift_a[static_cast<std::size_t>(a)] =
                dVc[static_cast<std::size_t>(a)]
                + dVs[static_cast<std::size_t>(a)];
            eps_shift_b[static_cast<std::size_t>(a)] =
                dVc[static_cast<std::size_t>(a)]
                - dVs[static_cast<std::size_t>(a)];
        }

        linalg::RealMatrix Ha, Hb;
        apply_shift(n_at, el, ao0, H0, S, eps_shift_a, Ha);
        apply_shift(n_at, el, ao0, H0, S, eps_shift_b, Hb);

        auto sol_a = linalg::diagonalize_fock_generalized(Ha, S);
        auto sol_b = linalg::diagonalize_fock_generalized(Hb, S);
        epsA = sol_a.first;
        CA = sol_a.second;
        epsB = sol_b.first;
        CB = sol_b.second;

        Pa = xt::zeros<double>({n_ao, n_ao});
        Pb = xt::zeros<double>({n_ao, n_ao});
        for (int o = 0; o < n_alpha; ++o)
            for (int mu = 0; mu < n_ao; ++mu)
                for (int nu = 0; nu < n_ao; ++nu)
                    Pa(mu, nu) += CA(mu, o) * CA(nu, o);
        for (int o = 0; o < n_beta; ++o)
            for (int mu = 0; mu < n_ao; ++mu)
                for (int nu = 0; nu < n_ao; ++nu)
                    Pb(mu, nu) += CB(mu, o) * CB(nu, o);

        linalg::RealMatrix PSa = xt::linalg::dot(Pa, S);
        linalg::RealMatrix PSb = xt::linalg::dot(Pb, S);

        residual = 0;
        for (int a = 0; a < n_at; ++a) {
            double qa = 0, qb = 0;
            const int ia = ao0[static_cast<std::size_t>(a)];
            const int na = el[static_cast<std::size_t>(a)]->n_ao;
            for (int i = 0; i < na; ++i) {
                qa += PSa(ia + i, ia + i);
                qb += PSb(ia + i, ia + i);
            }
            dq_out[static_cast<std::size_t>(a)] =
                (qa + qb) - el[static_cast<std::size_t>(a)]->n_val;
            mspin_out[static_cast<std::size_t>(a)] = qa - qb;
            residual = std::max(residual,
                                std::abs(dq_out[static_cast<std::size_t>(a)]
                                         - dq[static_cast<std::size_t>(a)]));
            residual = std::max(residual,
                                std::abs(mspin_out[static_cast<std::size_t>(a)]
                                         - mspin[static_cast<std::size_t>(a)]));
        }

        if (residual < scc_tol) {
            dq = dq_out;
            mspin = mspin_out;
            converged = true;
            break;
        }
        std::vector<double> in(dq);
        in.insert(in.end(), mspin.begin(), mspin.end());
        std::vector<double> out_vec(dq_out);
        out_vec.insert(out_vec.end(), mspin_out.begin(), mspin_out.end());
        const std::vector<double> mixed = mixer.next(in, out_vec);
        for (std::size_t a = 0; a < n; ++a) {
            dq[a] = mixed[a];
            mspin[a] = mixed[n + a];
        }
    }

    require_scc_converged(converged, residual, scc_tol, max_scc);

    e_scc = 0;
    for (int a = 0; a < n_at; ++a)
        for (int b = 0; b < n_at; ++b)
            e_scc += 0.5 * gamma_pair(a, b, el, x, y, z, box)
                     * dq[static_cast<std::size_t>(a)]
                     * dq[static_cast<std::size_t>(b)];
    e_spin = 0;
    for (int a = 0; a < n_at; ++a)
        e_spin += 0.5 * el[static_cast<std::size_t>(a)]->spin_w
                  * mspin[static_cast<std::size_t>(a)]
                  * mspin[static_cast<std::size_t>(a)];

    SccResult out;
    // E = Tr(ρH⁰) + E² + E_spin + E_rep (same convention as the restricted
    // path): Tr(ρ^σ H^σ) carries Σ_A V_A q_A on the GROSS population.
    out.energy = trace_product(Pa, H0) + trace_product(Pb, H0) + e_scc + e_spin
                 + repulsive_energy(n_at, type, x, y, z, box, el);
    out.residual = residual;
    out.iterations = iterations;
    out.dq = std::move(dq);
    out.spin = std::move(mspin);
    out.Ca = std::move(CA);
    out.Cb = std::move(CB);
    out.eps_a = std::move(epsA);
    out.eps_b = std::move(epsB);
    out.n_alpha = n_alpha;
    out.n_beta = n_beta;
    out.unrestricted = true;
    return out;
}


/// Closed-form nuclear gradient for **collinear unrestricted** SCC-DFTB2
/// (fixed Δq, m; dual α/β densities). Spin W is geometry-independent.
/// Lattice-image walk matches Γ-point energy (H⁰/S/γ/rep); self-images (a=a,
/// T≠0) have r = |T| only and contribute zero nuclear force.
void analytic_forces_unrestricted(
    int n_at, int n_ao, int n_alpha, int n_beta,
    const std::vector<const DftbElementRow*>& el, const std::vector<int>& ao0,
    const int* type, const StorageReal* x, const StorageReal* y,
    const StorageReal* z, const cpu::Box& box, const std::vector<double>& dq,
    const std::vector<double>& mspin, const linalg::RealMatrix& Ca,
    const linalg::RealMatrix& Cb, const linalg::RealVector& eps_a_orb,
    const linalg::RealVector& eps_b_orb, AccReal* fx, AccReal* fy,
    AccReal* fz) {
    // Atom shifts ε^σ_A = Σ_B γ_AB Δq_B ± W_A m_A (γ_AB already lattice-summed).
    std::vector<double> eps_at_a(static_cast<std::size_t>(n_at), 0.0);
    std::vector<double> eps_at_b(static_cast<std::size_t>(n_at), 0.0);
    for (int a = 0; a < n_at; ++a) {
        double shift = 0;
        for (int b = 0; b < n_at; ++b)
            shift += gamma_pair(a, b, el, x, y, z, box)
                     * dq[static_cast<std::size_t>(b)];
        const double vs = el[static_cast<std::size_t>(a)]->spin_w
                          * mspin[static_cast<std::size_t>(a)];
        eps_at_a[static_cast<std::size_t>(a)] = shift + vs;
        eps_at_b[static_cast<std::size_t>(a)] = shift - vs;
    }

    linalg::RealMatrix Pa = xt::zeros<double>({n_ao, n_ao});
    linalg::RealMatrix Pb = xt::zeros<double>({n_ao, n_ao});
    linalg::RealMatrix Wa = xt::zeros<double>({n_ao, n_ao});
    linalg::RealMatrix Wb = xt::zeros<double>({n_ao, n_ao});
    for (int o = 0; o < n_alpha; ++o)
        for (int mu = 0; mu < n_ao; ++mu)
            for (int nu = 0; nu < n_ao; ++nu) {
                Pa(mu, nu) += Ca(mu, o) * Ca(nu, o);
                Wa(mu, nu) += eps_a_orb(o) * Ca(mu, o) * Ca(nu, o);
            }
    for (int o = 0; o < n_beta; ++o)
        for (int mu = 0; mu < n_ao; ++mu)
            for (int nu = 0; nu < n_ao; ++nu) {
                Pb(mu, nu) += Cb(mu, o) * Cb(nu, o);
                Wb(mu, nu) += eps_b_orb(o) * Cb(mu, o) * Cb(nu, o);
            }

    linalg::RealMatrix dH0[3] = {xt::zeros<double>({n_ao, n_ao}),
                                 xt::zeros<double>({n_ao, n_ao}),
                                 xt::zeros<double>({n_ao, n_ao})};
    linalg::RealMatrix dS[3] = {xt::zeros<double>({n_ao, n_ao}),
                                xt::zeros<double>({n_ao, n_ao}),
                                xt::zeros<double>({n_ao, n_ao})};

    const auto images = lattice_images_for(box, el);

    for (int a = 0; a < n_at; ++a) {
        const auto& ea = *el[static_cast<std::size_t>(a)];
        const int ia = ao0[static_cast<std::size_t>(a)];
        const double xa = static_cast<double>(x[a]);
        const double ya = static_cast<double>(y[a]);
        const double za = static_cast<double>(z[a]);
        for (int b = a + 1; b < n_at; ++b) {
            const auto& eb = *el[static_cast<std::size_t>(b)];
            const int ib = ao0[static_cast<std::size_t>(b)];
            const double xb = static_cast<double>(x[b]);
            const double yb = static_cast<double>(y[b]);
            const double zb = static_cast<double>(z[b]);
            const double rcut = pair_rcut_A(ea, eb, type[a], type[b]);
            const auto& sk = pair_sk(type[a], type[b]);
            const double dqa = dq[static_cast<std::size_t>(a)];
            const double dqb = dq[static_cast<std::size_t>(b)];
            const double shift_a =
                0.5 * (eps_at_a[static_cast<std::size_t>(a)]
                       + eps_at_a[static_cast<std::size_t>(b)]);
            const double shift_b =
                0.5 * (eps_at_b[static_cast<std::size_t>(a)]
                       + eps_at_b[static_cast<std::size_t>(b)]);
            const int na = ea.n_ao;
            const int nb = eb.n_ao;

            for (const auto& T : images) {
                const double dx = xb + T.tx - xa;
                const double dy = yb + T.ty - ya;
                const double dz = zb + T.tz - za;
                const double r = std::sqrt(dx * dx + dy * dy + dz * dz);
                if (r < 1e-12)
                    continue;
                const double dir[3] = {dx / r, dy / r, dz / r};
                double dE_dbond[3] = {0, 0, 0};

                if (r <= rcut) {
                    for (int k = 0; k < 3; ++k) {
                        dH0[k] = xt::zeros<double>({n_ao, n_ao});
                        dS[k] = xt::zeros<double>({n_ao, n_ao});
                    }
                    fill_sk_block_deriv(ea, eb, sk, dx, dy, dz, r, dH0, dS, ia,
                                        ib);
                    for (int k = 0; k < 3; ++k) {
                        double acc = 0;
                        for (int i = 0; i < na; ++i)
                            for (int j = 0; j < nb; ++j) {
                                const double dh = dH0[k](ia + i, ib + j);
                                const double ds = dS[k](ia + i, ib + j);
                                acc += Pa(ia + i, ib + j) * (dh + shift_a * ds)
                                       - Wa(ia + i, ib + j) * ds;
                                acc += Pb(ia + i, ib + j) * (dh + shift_b * ds)
                                       - Wb(ia + i, ib + j) * ds;
                            }
                        // Factor 2: only the A–B upper block is contracted.
                        dE_dbond[k] += 2.0 * acc;
                    }
                    const double dErep_dr =
                        pair_repulsive_dE_dr(type[a], type[b], r);
                    for (int k = 0; k < 3; ++k)
                        dE_dbond[k] += dErep_dr * dir[k];
                }

                // Lattice γ: each image of γ(|R+T|) contributes independently.
                // E² = ½Σ_AB γ_AB Δq_A Δq_B with γ_AB = Σ_T γ(|R+T|); a<b
                // pair term is γ_AB Δq_A Δq_B so ∂E²/∂bond_T = Δq_A Δq_B dγ/dr.
                if (r <= rcut) {
                    const double dgdR =
                        dgamma_dr(r, ea.hubbard_u, eb.hubbard_u);
                    for (int k = 0; k < 3; ++k)
                        dE_dbond[k] += dqa * dqb * dgdR * dir[k];
                }

                // Spin energy −½ Σ W m²: W constant, m fixed ⇒ no geometry term.

                fx[a] += static_cast<AccReal>(dE_dbond[0]);
                fy[a] += static_cast<AccReal>(dE_dbond[1]);
                fz[a] += static_cast<AccReal>(dE_dbond[2]);
                fx[b] -= static_cast<AccReal>(dE_dbond[0]);
                fy[b] -= static_cast<AccReal>(dE_dbond[1]);
                fz[b] -= static_cast<AccReal>(dE_dbond[2]);
            }
        }
    }
}

/// Closed-form nuclear gradient (+ optional intensive stress) for **restricted**
/// SCC-DFTB2 (fixed Δq). Forces: F=−∇E [eV/Å]. Stress: σ=−(1/Ω)∂E/∂ε [eV/Å³]
/// Voigt (xx,yy,zz,xy,xz,yz) via pair virial of the same bond terms.
///
/// **Γ-PBC:** walks the same lattice images as `build_h0_s` / `gamma_pair` /
/// `repulsive_energy` (not min-image only). Self-images contribute zero F.
void analytic_forces_restricted(
    int n_at, int n_ao, int n_occ, const std::vector<const DftbElementRow*>& el,
    const std::vector<int>& ao0, const int* type, const StorageReal* x,
    const StorageReal* y, const StorageReal* z, const cpu::Box& box,
    const std::vector<double>& dq, const linalg::RealMatrix& C,
    const linalg::RealVector& eps, AccReal* fx, AccReal* fy, AccReal* fz,
    bool want_forces, double* stress_voigt /* nullable, length 6 */) {
    std::vector<double> dV(static_cast<std::size_t>(n_at), 0.0);
    for (int a = 0; a < n_at; ++a) {
        double shift = 0;
        for (int b = 0; b < n_at; ++b)
            shift += gamma_pair(a, b, el, x, y, z, box)
                     * dq[static_cast<std::size_t>(b)];
        dV[static_cast<std::size_t>(a)] = shift;
    }

    linalg::RealMatrix P = xt::zeros<double>({n_ao, n_ao});
    linalg::RealMatrix W = xt::zeros<double>({n_ao, n_ao});
    for (int o = 0; o < n_occ; ++o)
        for (int mu = 0; mu < n_ao; ++mu)
            for (int nu = 0; nu < n_ao; ++nu) {
                P(mu, nu) += 2.0 * C(mu, o) * C(nu, o);
                W(mu, nu) += 2.0 * eps(o) * C(mu, o) * C(nu, o);
            }

    linalg::RealMatrix dH0[3] = {xt::zeros<double>({n_ao, n_ao}),
                                 xt::zeros<double>({n_ao, n_ao}),
                                 xt::zeros<double>({n_ao, n_ao})};
    linalg::RealMatrix dS[3] = {xt::zeros<double>({n_ao, n_ao}),
                                xt::zeros<double>({n_ao, n_ao}),
                                xt::zeros<double>({n_ao, n_ao})};

    // Accumulate ∂E/∂ε (before −1/Ω) in Voigt order.
    double dE_de[6] = {0, 0, 0, 0, 0, 0};

    const auto images = lattice_images_for(box, el);

    for (int a = 0; a < n_at; ++a) {
        const auto& ea = *el[static_cast<std::size_t>(a)];
        const int ia = ao0[static_cast<std::size_t>(a)];
        const double xa = static_cast<double>(x[a]);
        const double ya = static_cast<double>(y[a]);
        const double za = static_cast<double>(z[a]);
        for (int b = a + 1; b < n_at; ++b) {
            const auto& eb = *el[static_cast<std::size_t>(b)];
            const int ib = ao0[static_cast<std::size_t>(b)];
            const double xb = static_cast<double>(x[b]);
            const double yb = static_cast<double>(y[b]);
            const double zb = static_cast<double>(z[b]);
            const double rcut = pair_rcut_A(ea, eb, type[a], type[b]);
            const auto& sk = pair_sk(type[a], type[b]);
            const double dqa = dq[static_cast<std::size_t>(a)];
            const double dqb = dq[static_cast<std::size_t>(b)];
            const double shift = 0.5 * (dV[static_cast<std::size_t>(a)]
                                        + dV[static_cast<std::size_t>(b)]);
            const int na = ea.n_ao;
            const int nb = eb.n_ao;

            for (const auto& T : images) {
                const double dx = xb + T.tx - xa;
                const double dy = yb + T.ty - ya;
                const double dz = zb + T.tz - za;
                const double r = std::sqrt(dx * dx + dy * dy + dz * dz);
                if (r < 1e-12)
                    continue;
                const double dir[3] = {dx / r, dy / r, dz / r};
                // Pair force contribution on atom a (F_b = −F_a).
                double F_a[3] = {0, 0, 0};

                if (r <= rcut) {
                    for (int k = 0; k < 3; ++k) {
                        dH0[k] = xt::zeros<double>({n_ao, n_ao});
                        dS[k] = xt::zeros<double>({n_ao, n_ao});
                    }
                    fill_sk_block_deriv(ea, eb, sk, dx, dy, dz, r, dH0, dS, ia,
                                        ib);
                    for (int k = 0; k < 3; ++k) {
                        double acc = 0;
                        for (int i = 0; i < na; ++i)
                            for (int j = 0; j < nb; ++j) {
                                const double p = P(ia + i, ib + j);
                                const double w = W(ia + i, ib + j);
                                const double dh = dH0[k](ia + i, ib + j);
                                const double ds = dS[k](ia + i, ib + j);
                                acc += p * (dh + shift * ds) - w * ds;
                            }
                        F_a[k] += 2.0 * acc;
                    }

                    const double dgdR =
                        dgamma_dr(r, ea.hubbard_u, eb.hubbard_u);
                    for (int k = 0; k < 3; ++k)
                        F_a[k] += dqa * dqb * dgdR * dir[k];

                    const double dErep_dr =
                        pair_repulsive_dE_dr(type[a], type[b], r);
                    for (int k = 0; k < 3; ++k)
                        F_a[k] += dErep_dr * dir[k];
                }

                if (want_forces && fx && fy && fz) {
                    fx[a] += static_cast<AccReal>(F_a[0]);
                    fy[a] += static_cast<AccReal>(F_a[1]);
                    fz[a] += static_cast<AccReal>(F_a[2]);
                    fx[b] -= static_cast<AccReal>(F_a[0]);
                    fy[b] -= static_cast<AccReal>(F_a[1]);
                    fz[b] -= static_cast<AccReal>(F_a[2]);
                }

                if (stress_voigt) {
                    const double R[3] = {dx, dy, dz};
                    dE_de[0] += F_a[0] * R[0];
                    dE_de[1] += F_a[1] * R[1];
                    dE_de[2] += F_a[2] * R[2];
                    dE_de[3] += F_a[0] * R[1] + F_a[1] * R[0]; // xy
                    dE_de[4] += F_a[0] * R[2] + F_a[2] * R[0]; // xz
                    dE_de[5] += F_a[1] * R[2] + F_a[2] * R[1]; // yz
                }
            }
        }
    }

    if (stress_voigt) {
        const double vol = box.volume();
        if (vol <= 0.0)
            throw std::runtime_error("DFTB2: non-positive box volume for stress");
        for (int v = 0; v < 6; ++v)
            stress_voigt[v] = -dE_de[v] / vol;
    }
}


/// @param mesh  Monkhorst–Pack subdivisions; {1,1,1} uses the real Γ path.
/// @param time_reversal Fold ±k when multi-k.
SccResult scc_total(int n_at, int n_ao, const Occupancy& occ,
                    const std::vector<const DftbElementRow*>& el,
                    const std::vector<int>& ao0, const int* type,
                    const StorageReal* x, const StorageReal* y,
                    const StorageReal* z, const cpu::Box& box, double scc_tol,
                    int max_scc, double mix_beta, const int mesh[3],
                    bool time_reversal) {
    const bool multi_k =
        mesh[0] > 1 || mesh[1] > 1 || mesh[2] > 1;
    if (multi_k) {
        if (occ.unrestricted)
            throw std::runtime_error(
                "DFTB2 multi-k: open shell not supported (restricted only)");
        return scc_energy_restricted_kmesh(
            n_at, n_ao, occ.n_alpha, el, ao0, type, x, y, z, box, scc_tol,
            max_scc, mix_beta, mesh, time_reversal);
    }
    if (!occ.unrestricted)
        return scc_energy_restricted(n_at, n_ao, occ.n_alpha, el, ao0, type, x,
                                     y, z, box, scc_tol, max_scc, mix_beta);
    return scc_energy_unrestricted(n_at, n_ao, occ.n_alpha, occ.n_beta, el, ao0,
                                   type, x, y, z, box, scc_tol, max_scc,
                                   mix_beta);
}


/// Affine strain of box + coordinates (mechanical ε; matches PW/GPW stress FD).
void apply_affine_strain(cpu::Frame& frame, int voigt, double h) {
    // Voigt: 0=xx,1=yy,2=zz,3=xy,4=xz,5=yz. Off-diag: ε_αβ = ε_βα = h.
    Mat3 F{};
    for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b)
            F(a, b) = (a == b) ? 1.0 : 0.0;
    if (voigt == 0)
        F(0, 0) += h;
    else if (voigt == 1)
        F(1, 1) += h;
    else if (voigt == 2)
        F(2, 2) += h;
    else if (voigt == 3) {
        F(0, 1) += h;
        F(1, 0) += h;
    } else if (voigt == 4) {
        F(0, 2) += h;
        F(2, 0) += h;
    } else if (voigt == 5) {
        F(1, 2) += h;
        F(2, 1) += h;
    }

    Mat3 cell = frame.box.cell();
    Mat3 cell_new;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            double s = 0;
            for (int k = 0; k < 3; ++k)
                s += F(i, k) * static_cast<double>(cell(k, j));
            cell_new(i, j) = static_cast<AccReal>(s);
        }
    frame.box.set_cell(cell_new);

    auto x = frame.span(kAtom_x);
    auto y = frame.span(kAtom_y);
    auto z = frame.span(kAtom_z);
    const int n = frame.nlocal();
    for (int a = 0; a < n; ++a) {
        const double rx = static_cast<double>(x[a]);
        const double ry = static_cast<double>(y[a]);
        const double rz = static_cast<double>(z[a]);
        x[a] = static_cast<StorageReal>(F(0, 0) * rx + F(0, 1) * ry
                                        + F(0, 2) * rz);
        y[a] = static_cast<StorageReal>(F(1, 0) * rx + F(1, 1) * ry
                                        + F(1, 2) * rz);
        z[a] = static_cast<StorageReal>(F(2, 0) * rx + F(2, 1) * ry
                                        + F(2, 2) * rz);
    }
}

double launch_mio(cpu::Frame& f, double scc_tol, int max_scc, double mix_beta,
                  bool with_forces, bool with_stress, int charge,
                  int multiplicity, const int k_mesh[3], bool k_time_reversal,
                  std::vector<double>& out_dq, std::vector<double>& out_spin,
                  std::array<double, 6>& out_stress, double& out_residual,
                  int& out_iterations) {
    const int n_at = f.nlocal();
    const int* type = f.span(kAtom_type).data();
    StorageReal* x = f.span(kAtom_x);
    StorageReal* y = f.span(kAtom_y);
    StorageReal* z = f.span(kAtom_z);
    AccReal* fx = f.span(kAtom_fx);
    AccReal* fy = f.span(kAtom_fy);
    AccReal* fz = f.span(kAtom_fz);
    double& energy_meta = f.meta.get(kEnergy);
    out_stress.fill(0.0);

    const bool multi_k =
        k_mesh[0] > 1 || k_mesh[1] > 1 || k_mesh[2] > 1;
    if (multi_k && (with_forces || with_stress))
        throw std::runtime_error(
            "DFTB2 multi-k: forces/stress not supported (energy-only path; "
            "set with_forces=false, with_stress=false)");

    std::vector<const DftbElementRow*> el(static_cast<std::size_t>(n_at));
    std::vector<int> ao0(static_cast<std::size_t>(n_at) + 1, 0);
    double n_elec = -static_cast<double>(charge);
    for (int a = 0; a < n_at; ++a) {
        el[static_cast<std::size_t>(a)] = &element_row(type[a]);
        n_elec += el[static_cast<std::size_t>(a)]->n_val;
        ao0[static_cast<std::size_t>(a) + 1] =
            ao0[static_cast<std::size_t>(a)]
            + el[static_cast<std::size_t>(a)]->n_ao;
    }
    const int n_ao = ao0[static_cast<std::size_t>(n_at)];
    if (n_ao < 1)
        throw std::runtime_error("DFTB2: no basis functions");

    const Occupancy occ = resolve_occupancy(n_elec, multiplicity);
    auto res = scc_total(n_at, n_ao, occ, el, ao0, type, x, y, z, f.box,
                         scc_tol, max_scc, mix_beta, k_mesh, k_time_reversal);
    out_dq = res.dq;
    out_spin = res.spin;
    out_residual = res.residual;
    out_iterations = res.iterations;
    const double energy = res.energy;
    energy_meta += energy;
    // The electronic per-term share of the total, same operator, same value
    // (esmm-contract-05): DFTB2 already accumulated into kEnergy.
    f.meta.get(kEnergyElectronic) += energy;

    // Multi-k is energy-only (no real Γ orbitals for analytic forces).
    if (multi_k)
        return energy;

    // Restricted closed shell has a closed-form gradient AND a closed-form pair
    // virial; collinear open shell has the closed-form gradient only, so its
    // stress still falls back to multi-SCC affine strain-FD.
    const bool analytic =
        !occ.unrestricted && res.n_occ > 0 && res.C.shape(0) > 0;
    const bool analytic_open = occ.unrestricted && res.Ca.shape(0) > 0;

    if ((with_forces || with_stress) && analytic) {
        double stress_buf[6] = {0, 0, 0, 0, 0, 0};
        analytic_forces_restricted(
            n_at, n_ao, res.n_occ, el, ao0, type, x, y, z, f.box, res.dq,
            res.C, res.eps, fx, fy, fz, with_forces,
            with_stress ? stress_buf : nullptr);
        if (with_stress) {
            for (int v = 0; v < 6; ++v)
                out_stress[static_cast<std::size_t>(v)] = stress_buf[v];
            auto& st = f.meta.get(kStress);
            for (int v = 0; v < 6; ++v)
                st[static_cast<std::size_t>(v)] += stress_buf[v];
        }
    } else {
        if (with_forces && analytic_open) {
            // Collinear open shell: closed-form dual α/β gradient.
            analytic_forces_unrestricted(
                n_at, n_ao, res.n_alpha, res.n_beta, el, ao0, type, x, y, z,
                f.box, res.dq, res.spin, res.Ca, res.Cb, res.eps_a, res.eps_b,
                fx, fy, fz);
        }
        // else: no orbitals (empty system) — leave forces unchanged
        if (with_stress) {
            // Multi-SCC affine strain-FD of E (open-shell fallback).
            constexpr double h_s = 2.0e-4;
            const double vol0 = f.box.volume();
            if (vol0 <= 0.0)
                throw std::runtime_error(
                    "DFTB2: non-positive box volume for stress");
            // Snapshot geometry
            std::vector<StorageReal> x0(static_cast<std::size_t>(n_at));
            std::vector<StorageReal> y0(static_cast<std::size_t>(n_at));
            std::vector<StorageReal> z0(static_cast<std::size_t>(n_at));
            for (int a = 0; a < n_at; ++a) {
                x0[static_cast<std::size_t>(a)] = x[a];
                y0[static_cast<std::size_t>(a)] = y[a];
                z0[static_cast<std::size_t>(a)] = z[a];
            }
            const Mat3 cell0 = f.box.cell();
            auto& st = f.meta.get(kStress);
            for (int v = 0; v < 6; ++v) {
                auto restore = [&]() {
                    for (int a = 0; a < n_at; ++a) {
                        x[a] = x0[static_cast<std::size_t>(a)];
                        y[a] = y0[static_cast<std::size_t>(a)];
                        z[a] = z0[static_cast<std::size_t>(a)];
                    }
                    f.box.set_cell(cell0);
                };
                restore();
                apply_affine_strain(f, v, h_s);
                // Stress FD is Γ-only (multi-k returns earlier).
                const int gamma_mesh[3] = {1, 1, 1};
                const double ep =
                    scc_total(n_at, n_ao, occ, el, ao0, type, x, y, z, f.box,
                              scc_tol, max_scc, mix_beta, gamma_mesh,
                              /*time_reversal=*/true)
                        .energy;
                restore();
                apply_affine_strain(f, v, -h_s);
                const double em =
                    scc_total(n_at, n_ao, occ, el, ao0, type, x, y, z, f.box,
                              scc_tol, max_scc, mix_beta, gamma_mesh,
                              /*time_reversal=*/true)
                        .energy;
                restore();
                const double dE_de = (ep - em) / (2.0 * h_s);
                const double sig = -dE_de / vol0;
                out_stress[static_cast<std::size_t>(v)] = sig;
                st[static_cast<std::size_t>(v)] += sig;
            }
        }
    }
    return energy;
}

// ── 3ob-3-1 path (internal Ha / Bohr; public eV / Å) ───────────────────

const Dftb3obElement& element_3ob(int z) {
    for (std::size_t i = 0; i < kDftb3obNElements; ++i)
        if (kDftb3obElements[i].z == z)
            return kDftb3obElements[i];
    throw std::runtime_error("DFTB2(ThreeOb): unsupported element Z="
                             + std::to_string(z) + " (H/C/N/O only)");
}

const Dftb3obPair& pair_3ob(int za, int zb) {
    for (std::size_t i = 0; i < kDftb3obNPairs; ++i)
        if (kDftb3obPairs[i].z_a == za && kDftb3obPairs[i].z_b == zb)
            return kDftb3obPairs[i];
    throw std::runtime_error("DFTB2(ThreeOb): no SK table for pair Z="
                             + std::to_string(za) + "-" + std::to_string(zb));
}

/// Linear interp of SK channels at r [Bohr]. Returns Hss,Hsp,Hpps,Hppp,Sss,...
void interp_3ob_sk(const Dftb3obPair& p, double r_bohr, double out[8]) {
    if (r_bohr <= p.grid) {
        const Dftb3obSkPoint& pt = p.sk[0];
        out[0] = pt.h_ss;
        out[1] = pt.h_sp;
        out[2] = pt.h_pps;
        out[3] = pt.h_ppp;
        out[4] = pt.s_ss;
        out[5] = pt.s_sp;
        out[6] = pt.s_pps;
        out[7] = pt.s_ppp;
        return;
    }
    const double x = r_bohr / p.grid; // 1-based continuous index
    if (x >= static_cast<double>(p.n_grid)) {
        for (int i = 0; i < 8; ++i)
            out[i] = 0;
        return;
    }
    int i0 = static_cast<int>(std::floor(x)) - 1;
    if (i0 < 0)
        i0 = 0;
    if (i0 > p.n_grid - 2)
        i0 = p.n_grid - 2;
    const int i1 = i0 + 1;
    const double t = x - static_cast<double>(i0 + 1);
    const Dftb3obSkPoint& a = p.sk[i0];
    const Dftb3obSkPoint& b = p.sk[i1];
    const double u = 1.0 - t;
    out[0] = u * a.h_ss + t * b.h_ss;
    out[1] = u * a.h_sp + t * b.h_sp;
    out[2] = u * a.h_pps + t * b.h_pps;
    out[3] = u * a.h_ppp + t * b.h_ppp;
    out[4] = u * a.s_ss + t * b.s_ss;
    out[5] = u * a.s_sp + t * b.s_sp;
    out[6] = u * a.s_pps + t * b.s_pps;
    out[7] = u * a.s_ppp + t * b.s_ppp;
}

/// d(channel)/dR [1/Bohr] of linear table interp (flat outside the grid).
void interp_3ob_sk_dr(const Dftb3obPair& p, double r_bohr, double dout[8]) {
    if (r_bohr <= p.grid) {
        for (int i = 0; i < 8; ++i)
            dout[i] = 0;
        return;
    }
    const double x = r_bohr / p.grid;
    if (x >= static_cast<double>(p.n_grid)) {
        for (int i = 0; i < 8; ++i)
            dout[i] = 0;
        return;
    }
    int i0 = static_cast<int>(std::floor(x)) - 1;
    if (i0 < 0)
        i0 = 0;
    if (i0 > p.n_grid - 2)
        i0 = p.n_grid - 2;
    const int i1 = i0 + 1;
    const Dftb3obSkPoint& a = p.sk[i0];
    const Dftb3obSkPoint& b = p.sk[i1];
    const double inv_g = 1.0 / p.grid;
    dout[0] = (static_cast<double>(b.h_ss) - static_cast<double>(a.h_ss)) * inv_g;
    dout[1] = (static_cast<double>(b.h_sp) - static_cast<double>(a.h_sp)) * inv_g;
    dout[2] = (static_cast<double>(b.h_pps) - static_cast<double>(a.h_pps)) * inv_g;
    dout[3] = (static_cast<double>(b.h_ppp) - static_cast<double>(a.h_ppp)) * inv_g;
    dout[4] = (static_cast<double>(b.s_ss) - static_cast<double>(a.s_ss)) * inv_g;
    dout[5] = (static_cast<double>(b.s_sp) - static_cast<double>(a.s_sp)) * inv_g;
    dout[6] = (static_cast<double>(b.s_pps) - static_cast<double>(a.s_pps)) * inv_g;
    dout[7] = (static_cast<double>(b.s_ppp) - static_cast<double>(a.s_ppp)) * inv_g;
}

double rep_3ob(const Dftb3obPair& p, double r_bohr) {
    if (r_bohr >= p.rep_cutoff)
        return 0.0;
    // First interval start
    const double r0 = p.rep_int[0];
    if (r_bohr < r0)
        return std::exp(-p.rep_a1 * r_bohr + p.rep_a2) + p.rep_a3;
    for (int i = 0; i < p.n_rep_int; ++i) {
        const double* iv = p.rep_int + static_cast<std::size_t>(i) * 8;
        const double ra = iv[0];
        const double rb = iv[1];
        if (r_bohr < rb || std::abs(r_bohr - rb) < 1e-14) {
            const double dr = r_bohr - ra;
            double e = 0;
            double pow = 1.0;
            for (int k = 0; k < 6; ++k) {
                e += iv[2 + k] * pow;
                pow *= dr;
            }
            return e;
        }
    }
    return 0.0;
}

/// dE_rep/dR [Ha/Bohr] of the 3ob repulsive (exp + spline).
double rep_3ob_dr(const Dftb3obPair& p, double r_bohr) {
    if (r_bohr >= p.rep_cutoff)
        return 0.0;
    const double r0 = p.rep_int[0];
    if (r_bohr < r0)
        return -p.rep_a1 * std::exp(-p.rep_a1 * r_bohr + p.rep_a2);
    for (int i = 0; i < p.n_rep_int; ++i) {
        const double* iv = p.rep_int + static_cast<std::size_t>(i) * 8;
        const double ra = iv[0];
        const double rb = iv[1];
        if (r_bohr < rb || std::abs(r_bohr - rb) < 1e-14) {
            const double dr = r_bohr - ra;
            // e = c0 + c1 dr + c2 dr^2 + ... + c5 dr^5
            // de/dr = c1 + 2 c2 dr + ... + 5 c5 dr^4
            double de = 0;
            double pow = 1.0;
            for (int k = 1; k < 6; ++k) {
                de += static_cast<double>(k) * iv[2 + k] * pow;
                pow *= dr;
            }
            return de;
        }
    }
    return 0.0;
}

void fill_sk_block_3ob(const Dftb3obElement& ea, const Dftb3obElement& eb,
                       int za, int zb, double dx, double dy, double dz, double r,
                       linalg::RealMatrix& H0, linalg::RealMatrix& S, int ia0,
                       int ib0) {
    if (r < 1e-12)
        return;
    const double inv_r = 1.0 / r;
    const double l = dx * inv_r;
    const double m = dy * inv_r;
    const double n = dz * inv_r;

    const Dftb3obPair& p_ab = pair_3ob(za, zb);
    double ch[8];
    interp_3ob_sk(p_ab, r, ch);
    double Hss = ch[0], Hsp = ch[1], Hpps = ch[2], Hppp = ch[3];
    double Sss = ch[4], Ssp = ch[5], Spps = ch[6], Sppp = ch[7];
    // X–H tables store only ss; s–p lives in H–X.
    if (za != 1 && zb == 1) {
        const Dftb3obPair& p_ba = pair_3ob(zb, za);
        double rev[8];
        interp_3ob_sk(p_ba, r, rev);
        Hsp = rev[1];
        Ssp = rev[5];
    }

    auto set_hs = [&](int i, int j, double hv, double sv) {
        H0(ia0 + i, ib0 + j) = hv;
        H0(ib0 + j, ia0 + i) = hv;
        S(ia0 + i, ib0 + j) = sv;
        S(ib0 + j, ia0 + i) = sv;
    };

    set_hs(0, 0, Hss, Sss);
    if (ea.n_ao >= 1 && eb.n_ao == 4) {
        set_hs(0, 1, l * Hsp, l * Ssp);
        set_hs(0, 2, m * Hsp, m * Ssp);
        set_hs(0, 3, n * Hsp, n * Ssp);
    }
    if (ea.n_ao == 4 && eb.n_ao >= 1) {
        set_hs(1, 0, -l * Hsp, -l * Ssp);
        set_hs(2, 0, -m * Hsp, -m * Ssp);
        set_hs(3, 0, -n * Hsp, -n * Ssp);
    }
    if (ea.n_ao == 4 && eb.n_ao == 4) {
        const double dir[3] = {l, m, n};
        for (int a = 0; a < 3; ++a) {
            for (int b = 0; b < 3; ++b) {
                const double d_ab = (a == b) ? 1.0 : 0.0;
                const double hv = dir[a] * dir[b] * (Hpps - Hppp) + d_ab * Hppp;
                const double sv = dir[a] * dir[b] * (Spps - Sppp) + d_ab * Sppp;
                set_hs(1 + a, 1 + b, hv, sv);
            }
        }
    }
}

/// SK block derivatives w.r.t. bond vector R_B − R_A [Bohr] (3ob tables).
void fill_sk_block_deriv_3ob(const Dftb3obElement& ea, const Dftb3obElement& eb,
                             int za, int zb, double dx, double dy, double dz,
                             double r, linalg::RealMatrix dH0[3],
                             linalg::RealMatrix dS[3], int ia0, int ib0) {
    if (r < 1e-12)
        return;
    const double inv_r = 1.0 / r;
    const double dir[3] = {dx * inv_r, dy * inv_r, dz * inv_r};

    const Dftb3obPair& p_ab = pair_3ob(za, zb);
    double ch[8], dch[8];
    interp_3ob_sk(p_ab, r, ch);
    interp_3ob_sk_dr(p_ab, r, dch);
    double Hss = ch[0], Hsp = ch[1], Hpps = ch[2], Hppp = ch[3];
    double Sss = ch[4], Ssp = ch[5], Spps = ch[6], Sppp = ch[7];
    double dHss = dch[0], dHsp = dch[1], dHpps = dch[2], dHppp = dch[3];
    double dSss = dch[4], dSsp = dch[5], dSpps = dch[6], dSppp = dch[7];
    // X–H tables store only ss; s–p lives in H–X.
    if (za != 1 && zb == 1) {
        const Dftb3obPair& p_ba = pair_3ob(zb, za);
        double rev[8], drev[8];
        interp_3ob_sk(p_ba, r, rev);
        interp_3ob_sk_dr(p_ba, r, drev);
        Hsp = rev[1];
        Ssp = rev[5];
        dHsp = drev[1];
        dSsp = drev[5];
    }

    auto ddir = [&](int alpha, int k) {
        const double dab = (alpha == k) ? 1.0 : 0.0;
        return (dab - dir[alpha] * dir[k]) * inv_r;
    };

    for (int k = 0; k < 3; ++k) {
        const double dr_dk = dir[k];
        dH0[k](ia0, ib0) = dHss * dr_dk;
        dS[k](ia0, ib0) = dSss * dr_dk;

        if (ea.n_ao >= 1 && eb.n_ao == 4) {
            for (int a = 0; a < 3; ++a) {
                dH0[k](ia0, ib0 + 1 + a) =
                    ddir(a, k) * Hsp + dir[a] * dHsp * dr_dk;
                dS[k](ia0, ib0 + 1 + a) =
                    ddir(a, k) * Ssp + dir[a] * dSsp * dr_dk;
            }
        }
        if (ea.n_ao == 4 && eb.n_ao >= 1) {
            for (int a = 0; a < 3; ++a) {
                dH0[k](ia0 + 1 + a, ib0) =
                    -(ddir(a, k) * Hsp + dir[a] * dHsp * dr_dk);
                dS[k](ia0 + 1 + a, ib0) =
                    -(ddir(a, k) * Ssp + dir[a] * dSsp * dr_dk);
            }
        }
        if (ea.n_ao == 4 && eb.n_ao == 4) {
            const double dHsig = dHpps - dHppp;
            const double dSsig = dSpps - dSppp;
            const double Hsig = Hpps - Hppp;
            const double Ssig = Spps - Sppp;
            for (int a = 0; a < 3; ++a) {
                for (int b = 0; b < 3; ++b) {
                    const double dab = (a == b) ? 1.0 : 0.0;
                    dH0[k](ia0 + 1 + a, ib0 + 1 + b) =
                        (ddir(a, k) * dir[b] + dir[a] * ddir(b, k)) * Hsig
                        + dir[a] * dir[b] * dHsig * dr_dk + dab * dHppp * dr_dk;
                    dS[k](ia0 + 1 + a, ib0 + 1 + b) =
                        (ddir(a, k) * dir[b] + dir[a] * ddir(b, k)) * Ssig
                        + dir[a] * dir[b] * dSsig * dr_dk + dab * dSppp * dr_dk;
                }
            }
        }
    }
}

double gamma_ab_ha(double r_bohr, double ua_ha, double ub_ha) {
    if (r_bohr < 1e-14)
        return ua_ha;
    const double tau = 0.5 * (1.0 / ua_ha + 1.0 / ub_ha);
    return 1.0 / std::sqrt(r_bohr * r_bohr + tau * tau);
}

/// dγ/dR [Ha/Bohr²] for Ohno–Klopman γ in Ha/Bohr.
double dgamma_dr_ha(double r_bohr, double ua_ha, double ub_ha) {
    if (r_bohr < 1e-14)
        return 0.0;
    const double g = gamma_ab_ha(r_bohr, ua_ha, ub_ha);
    return -r_bohr * g * g * g;
}

struct Scc3obResult {
    double energy_ha = 0;
    double residual = 0;
    int iterations = 0;
    std::vector<double> dq;
    linalg::RealMatrix C;
    linalg::RealVector eps;
    int n_occ = 0;
};

Scc3obResult scc_3ob(const std::vector<const Dftb3obElement*>& el,
                     const std::vector<int>& ao0, const std::vector<int>& Z,
                     const std::vector<double>& x_b,
                     const std::vector<double>& y_b,
                     const std::vector<double>& z_b, int n_at, int n_ao,
                     int n_occ, double scc_tol, int max_scc, double mix_beta) {
    auto dist = [&](int i, int j, double& dx, double& dy, double& dz) {
        dx = x_b[static_cast<std::size_t>(j)] - x_b[static_cast<std::size_t>(i)];
        dy = y_b[static_cast<std::size_t>(j)] - y_b[static_cast<std::size_t>(i)];
        dz = z_b[static_cast<std::size_t>(j)] - z_b[static_cast<std::size_t>(i)];
        return std::sqrt(dx * dx + dy * dy + dz * dz);
    };

    linalg::RealMatrix H0 = xt::zeros<double>({n_ao, n_ao});
    linalg::RealMatrix S = xt::zeros<double>({n_ao, n_ao});

    for (int a = 0; a < n_at; ++a) {
        const auto& ea = *el[static_cast<std::size_t>(a)];
        const int ia = ao0[static_cast<std::size_t>(a)];
        H0(ia, ia) = ea.e_s;
        S(ia, ia) = 1.0;
        if (ea.n_ao == 4) {
            for (int p = 1; p < 4; ++p) {
                H0(ia + p, ia + p) = ea.e_p;
                S(ia + p, ia + p) = 1.0;
            }
        }
        for (int b = a + 1; b < n_at; ++b) {
            const auto& eb = *el[static_cast<std::size_t>(b)];
            double dx, dy, dz;
            const double r = dist(a, b, dx, dy, dz);
            if (r > std::max(ea.r_cut, eb.r_cut))
                continue;
            fill_sk_block_3ob(ea, eb, Z[static_cast<std::size_t>(a)],
                              Z[static_cast<std::size_t>(b)], dx, dy, dz, r, H0,
                              S, ia, ao0[static_cast<std::size_t>(b)]);
        }
    }

    const std::size_t n = static_cast<std::size_t>(n_at);
    std::vector<double> dq(n, 0.0);
    std::vector<double> dq_out(n, 0.0);
    linalg::RealMatrix P = xt::zeros<double>({n_ao, n_ao});
    linalg::RealMatrix C = xt::zeros<double>({n_ao, n_ao});
    linalg::RealVector eps = xt::zeros<double>({n_ao});
    double e_scc = 0, residual = 0;
    int iterations = 0;
    bool converged = false;
    internal::ChargeMixer mixer(mix_beta, kSccMixDepth);

    for (int it = 0; it < max_scc; ++it) {
        iterations = it + 1;
        std::vector<double> dV(static_cast<std::size_t>(n_at), 0.0);
        for (int a = 0; a < n_at; ++a) {
            double shift = 0;
            for (int b = 0; b < n_at; ++b) {
                double dx, dy, dz;
                const double r = (a == b) ? 0.0 : dist(a, b, dx, dy, dz);
                const double g =
                    (a == b)
                        ? el[static_cast<std::size_t>(a)]->hubbard_u
                        : gamma_ab_ha(r,
                                      el[static_cast<std::size_t>(a)]->hubbard_u,
                                      el[static_cast<std::size_t>(b)]->hubbard_u);
                shift += g * dq[static_cast<std::size_t>(b)];
            }
            dV[static_cast<std::size_t>(a)] = shift;
        }

        linalg::RealMatrix H = H0;
        for (int a = 0; a < n_at; ++a) {
            const int ia = ao0[static_cast<std::size_t>(a)];
            const int na = el[static_cast<std::size_t>(a)]->n_ao;
            for (int b = 0; b < n_at; ++b) {
                const int ib = ao0[static_cast<std::size_t>(b)];
                const int nb = el[static_cast<std::size_t>(b)]->n_ao;
                const double half = 0.5
                                    * (dV[static_cast<std::size_t>(a)]
                                       + dV[static_cast<std::size_t>(b)]);
                for (int i = 0; i < na; ++i)
                    for (int j = 0; j < nb; ++j)
                        H(ia + i, ib + j) += half * S(ia + i, ib + j);
            }
        }

        auto sol = linalg::diagonalize_fock_generalized(H, S);
        eps = sol.first;
        C = sol.second;

        P = xt::zeros<double>({n_ao, n_ao});
        for (int o = 0; o < n_occ; ++o)
            for (int mu = 0; mu < n_ao; ++mu)
                for (int nu = 0; nu < n_ao; ++nu)
                    P(mu, nu) += 2.0 * C(mu, o) * C(nu, o);

        linalg::RealMatrix PS = xt::linalg::dot(P, S);
        residual = 0;
        for (int a = 0; a < n_at; ++a) {
            double qa = 0;
            const int ia = ao0[static_cast<std::size_t>(a)];
            const int na = el[static_cast<std::size_t>(a)]->n_ao;
            for (int i = 0; i < na; ++i)
                qa += PS(ia + i, ia + i);
            dq_out[static_cast<std::size_t>(a)] =
                qa - el[static_cast<std::size_t>(a)]->n_val;
            residual = std::max(residual,
                                std::abs(dq_out[static_cast<std::size_t>(a)]
                                         - dq[static_cast<std::size_t>(a)]));
        }

        if (residual < scc_tol) {
            dq = dq_out;
            converged = true;
            break;
        }
        dq = mixer.next(dq, dq_out);
    }

    require_scc_converged(converged, residual, scc_tol, max_scc);

    e_scc = 0;
    for (int a = 0; a < n_at; ++a)
        for (int b = 0; b < n_at; ++b) {
            double dx, dy, dz;
            const double r = (a == b) ? 0.0 : dist(a, b, dx, dy, dz);
            const double g =
                (a == b) ? el[static_cast<std::size_t>(a)]->hubbard_u
                         : gamma_ab_ha(
                               r, el[static_cast<std::size_t>(a)]->hubbard_u,
                               el[static_cast<std::size_t>(b)]->hubbard_u);
            e_scc += 0.5 * g * dq[static_cast<std::size_t>(a)]
                     * dq[static_cast<std::size_t>(b)];
        }

    // E_elec = Tr(P H0) + E_scc  (DFTB2; matches DFTB+ double-counting)
    const double e_h0 = trace_product(P, H0);

    double e_repulsive = 0;
    for (int a = 0; a < n_at; ++a)
        for (int b = a + 1; b < n_at; ++b) {
            double dx, dy, dz;
            const double r = dist(a, b, dx, dy, dz);
            const auto& pr = pair_3ob(Z[static_cast<std::size_t>(a)],
                                      Z[static_cast<std::size_t>(b)]);
            e_repulsive += rep_3ob(pr, r);
        }

    Scc3obResult out;
    out.energy_ha = e_h0 + e_scc + e_repulsive;
    out.residual = residual;
    out.iterations = iterations;
    out.dq = std::move(dq);
    out.C = std::move(C);
    out.eps = std::move(eps);
    out.n_occ = n_occ;
    return out;
}

/// Closed-form nuclear gradient (+ optional intensive stress) for **restricted**
/// SCC-DFTB2 on 3ob tables. Internally Ha/Bohr; writes Frame forces in eV/Å and
/// stress in eV/Å³ Voigt (xx,yy,zz,xy,xz,yz), mechanical convention
/// σ = −(1/Ω)∂E/∂ε — same assembly and same sign as the restricted MioSketch
/// path (`analytic_forces_restricted`).
///
/// `scc_3ob` is a raw-distance **cluster** sum (no lattice images, no
/// min-image), so E depends on geometry only through the pair bond vectors and
/// an affine strain of nuclei+cell gives ∂E/∂ε_αβ = Σ_pairs (∂E/∂R_α) R_β
/// exactly — the pair virial below is FD-exact, not an approximation.
///
/// @param vol_a3 Cell volume Ω [Å³]; only read when `stress_voigt` is non-null.
void analytic_forces_3ob(int n_at, int n_ao, int n_occ,
                         const std::vector<const Dftb3obElement*>& el,
                         const std::vector<int>& ao0, const std::vector<int>& Z,
                         const std::vector<double>& x_b,
                         const std::vector<double>& y_b,
                         const std::vector<double>& z_b,
                         const std::vector<double>& dq,
                         const linalg::RealMatrix& C,
                         const linalg::RealVector& eps, AccReal* fx, AccReal* fy,
                         AccReal* fz, bool want_forces,
                         double* stress_voigt /* nullable, length 6 */,
                         double vol_a3) {
    auto dist = [&](int i, int j, double& dx, double& dy, double& dz) {
        dx = x_b[static_cast<std::size_t>(j)] - x_b[static_cast<std::size_t>(i)];
        dy = y_b[static_cast<std::size_t>(j)] - y_b[static_cast<std::size_t>(i)];
        dz = z_b[static_cast<std::size_t>(j)] - z_b[static_cast<std::size_t>(i)];
        return std::sqrt(dx * dx + dy * dy + dz * dz);
    };

    std::vector<double> dV(static_cast<std::size_t>(n_at), 0.0);
    for (int a = 0; a < n_at; ++a) {
        double shift = 0;
        for (int b = 0; b < n_at; ++b) {
            double dx, dy, dz;
            const double r = (a == b) ? 0.0 : dist(a, b, dx, dy, dz);
            const double g =
                (a == b)
                    ? el[static_cast<std::size_t>(a)]->hubbard_u
                    : gamma_ab_ha(r, el[static_cast<std::size_t>(a)]->hubbard_u,
                                  el[static_cast<std::size_t>(b)]->hubbard_u);
            shift += g * dq[static_cast<std::size_t>(b)];
        }
        dV[static_cast<std::size_t>(a)] = shift;
    }

    linalg::RealMatrix P = xt::zeros<double>({n_ao, n_ao});
    linalg::RealMatrix W = xt::zeros<double>({n_ao, n_ao});
    for (int o = 0; o < n_occ; ++o)
        for (int mu = 0; mu < n_ao; ++mu)
            for (int nu = 0; nu < n_ao; ++nu) {
                P(mu, nu) += 2.0 * C(mu, o) * C(nu, o);
                W(mu, nu) += 2.0 * eps(o) * C(mu, o) * C(nu, o);
            }

    // Ha/Bohr → eV/Å: F_eV_A = F_Ha_Bohr * kHa2eV / kBohr2A
    const double fscale = kHa2eV / kBohr2A;

    linalg::RealMatrix dH0[3] = {xt::zeros<double>({n_ao, n_ao}),
                                 xt::zeros<double>({n_ao, n_ao}),
                                 xt::zeros<double>({n_ao, n_ao})};
    linalg::RealMatrix dS[3] = {xt::zeros<double>({n_ao, n_ao}),
                                xt::zeros<double>({n_ao, n_ao}),
                                xt::zeros<double>({n_ao, n_ao})};

    // ∂E/∂ε in Voigt order [Ha] (strain is dimensionless), before −1/Ω.
    double dE_de[6] = {0, 0, 0, 0, 0, 0};

    for (int a = 0; a < n_at; ++a) {
        const auto& ea = *el[static_cast<std::size_t>(a)];
        const int ia = ao0[static_cast<std::size_t>(a)];
        for (int b = a + 1; b < n_at; ++b) {
            const auto& eb = *el[static_cast<std::size_t>(b)];
            double dx, dy, dz;
            const double r = dist(a, b, dx, dy, dz);
            if (r < 1e-12)
                continue;
            const int ib = ao0[static_cast<std::size_t>(b)];
            const double dir[3] = {dx / r, dy / r, dz / r};
            double dE_dbond[3] = {0, 0, 0}; // ∂E/∂bond [Ha/Bohr]

            if (r <= std::max(ea.r_cut, eb.r_cut)) {
                for (int k = 0; k < 3; ++k) {
                    dH0[k] = xt::zeros<double>({n_ao, n_ao});
                    dS[k] = xt::zeros<double>({n_ao, n_ao});
                }
                fill_sk_block_deriv_3ob(ea, eb, Z[static_cast<std::size_t>(a)],
                                        Z[static_cast<std::size_t>(b)], dx, dy,
                                        dz, r, dH0, dS, ia, ib);
                const double shift = 0.5
                                     * (dV[static_cast<std::size_t>(a)]
                                        + dV[static_cast<std::size_t>(b)]);
                const int na = ea.n_ao;
                const int nb = eb.n_ao;
                for (int k = 0; k < 3; ++k) {
                    double acc = 0;
                    for (int i = 0; i < na; ++i)
                        for (int j = 0; j < nb; ++j) {
                            const double p = P(ia + i, ib + j);
                            const double w = W(ia + i, ib + j);
                            const double dh = dH0[k](ia + i, ib + j);
                            const double ds = dS[k](ia + i, ib + j);
                            acc += p * (dh + shift * ds) - w * ds;
                        }
                    dE_dbond[k] += 2.0 * acc;
                }
            }

            // E = Tr(PH0) + ½ Σ_AB γ_AB Δq_A Δq_B + E_rep; with
            // bond = R_B − R_A, ∂E²/∂bond = +Δq_A Δq_B dγ/dr.
            const double dgdR =
                dgamma_dr_ha(r, ea.hubbard_u, eb.hubbard_u);
            const double dqa = dq[static_cast<std::size_t>(a)];
            const double dqb = dq[static_cast<std::size_t>(b)];
            for (int k = 0; k < 3; ++k)
                dE_dbond[k] += dqa * dqb * dgdR * dir[k];

            const auto& pr = pair_3ob(Z[static_cast<std::size_t>(a)],
                                      Z[static_cast<std::size_t>(b)]);
            const double dErep_dr = rep_3ob_dr(pr, r);
            for (int k = 0; k < 3; ++k)
                dE_dbond[k] += dErep_dr * dir[k];

            // F = −∂E/∂R; bond = R_B − R_A ⇒ F_A += ∂E/∂bond, F_B −= ∂E/∂bond
            // after converting Ha/Bohr → eV/Å.
            if (want_forces && fx && fy && fz) {
                const double fax = dE_dbond[0] * fscale;
                const double fay = dE_dbond[1] * fscale;
                const double faz = dE_dbond[2] * fscale;
                fx[a] += static_cast<AccReal>(fax);
                fy[a] += static_cast<AccReal>(fay);
                fz[a] += static_cast<AccReal>(faz);
                fx[b] -= static_cast<AccReal>(fax);
                fy[b] -= static_cast<AccReal>(fay);
                fz[b] -= static_cast<AccReal>(faz);
            }

            if (stress_voigt) {
                // Pair virial, Voigt order and sign taken from the restricted
                // MioSketch path: off-diagonals are doubled because the affine
                // strain is symmetric (ε_αβ = ε_βα = h). [Ha], bond in Bohr.
                const double R[3] = {dx, dy, dz};
                dE_de[0] += dE_dbond[0] * R[0];
                dE_de[1] += dE_dbond[1] * R[1];
                dE_de[2] += dE_dbond[2] * R[2];
                dE_de[3] += dE_dbond[0] * R[1] + dE_dbond[1] * R[0]; // xy
                dE_de[4] += dE_dbond[0] * R[2] + dE_dbond[2] * R[0]; // xz
                dE_de[5] += dE_dbond[1] * R[2] + dE_dbond[2] * R[1]; // yz
            }
        }
    }

    if (stress_voigt) {
        if (vol_a3 <= 0.0)
            throw std::runtime_error(
                "DFTB2(ThreeOb): non-positive box volume for stress");
        // ∂E/∂ε [Ha] → [eV]; Ω [Å³] ⇒ σ [eV/Å³].
        for (int v = 0; v < 6; ++v)
            stress_voigt[v] = -dE_de[v] * kHa2eV / vol_a3;
    }
}

double launch_3ob(cpu::Frame& f, double scc_tol, int max_scc, double mix_beta,
                  bool with_forces, bool with_stress, std::vector<double>& out_dq,
                  std::array<double, 6>& out_stress, double& out_residual,
                  int& out_iterations) {
    const int n_at = f.nlocal();
    const auto type = f.span(kAtom_type);
    StorageReal* x = f.span(kAtom_x);
    StorageReal* y = f.span(kAtom_y);
    StorageReal* z = f.span(kAtom_z);
    AccReal* fx = f.span(kAtom_fx);
    AccReal* fy = f.span(kAtom_fy);
    AccReal* fz = f.span(kAtom_fz);
    double& energy_meta = f.meta.get(kEnergy);
    out_stress.fill(0.0);

    std::vector<const Dftb3obElement*> el(static_cast<std::size_t>(n_at));
    std::vector<int> ao0(static_cast<std::size_t>(n_at) + 1, 0);
    std::vector<int> Z(static_cast<std::size_t>(n_at));
    std::vector<double> xb(static_cast<std::size_t>(n_at));
    std::vector<double> yb(static_cast<std::size_t>(n_at));
    std::vector<double> zb(static_cast<std::size_t>(n_at));
    double n_elec = 0;
    for (int a = 0; a < n_at; ++a) {
        Z[static_cast<std::size_t>(a)] = type[a];
        el[static_cast<std::size_t>(a)] = &element_3ob(type[a]);
        n_elec += el[static_cast<std::size_t>(a)]->n_val;
        ao0[static_cast<std::size_t>(a) + 1] =
            ao0[static_cast<std::size_t>(a)]
            + el[static_cast<std::size_t>(a)]->n_ao;
        // Å → Bohr for SK/rep (molecular; no min-image in 3ob pin path)
        xb[static_cast<std::size_t>(a)] = static_cast<double>(x[a]) * kA2Bohr;
        yb[static_cast<std::size_t>(a)] = static_cast<double>(y[a]) * kA2Bohr;
        zb[static_cast<std::size_t>(a)] = static_cast<double>(z[a]) * kA2Bohr;
    }
    const int n_ao = ao0[static_cast<std::size_t>(n_at)];
    if (n_ao < 1)
        throw std::runtime_error("DFTB2(ThreeOb): no basis functions");
    if (std::fabs(n_elec - std::round(n_elec)) > 1e-9
        || static_cast<int>(std::lround(n_elec)) % 2 != 0)
        throw std::runtime_error(
            "DFTB2(ThreeOb): odd/non-integer valence electron count");
    const int n_occ = static_cast<int>(std::lround(n_elec)) / 2;

    auto res = scc_3ob(el, ao0, Z, xb, yb, zb, n_at, n_ao, n_occ, scc_tol,
                       max_scc, mix_beta);
    const double energy = res.energy_ha * kHa2eV;
    out_residual = res.residual;
    out_iterations = res.iterations;
    energy_meta += energy;
    // The electronic per-term share of the total, same operator, same value
    // (esmm-contract-05): DFTB2 already accumulated into kEnergy.
    f.meta.get(kEnergyElectronic) += energy;
    out_dq = std::move(res.dq);

    // One post-SCC pair contraction serves both F and σ: the virial is the
    // same bond terms as the gradient (no multi-SCC strain-FD product path).
    if (with_forces || with_stress) {
        double stress_buf[6] = {0, 0, 0, 0, 0, 0};
        analytic_forces_3ob(n_at, n_ao, res.n_occ, el, ao0, Z, xb, yb, zb,
                            out_dq, res.C, res.eps, fx, fy, fz, with_forces,
                            with_stress ? stress_buf : nullptr,
                            f.box.volume());
        if (with_stress) {
            auto& st = f.meta.get(kStress);
            for (int v = 0; v < 6; ++v) {
                out_stress[static_cast<std::size_t>(v)] = stress_buf[v];
                st[static_cast<std::size_t>(v)] += stress_buf[v];
            }
        }
    }
    return energy;
}


} // namespace

void DFTB2::build(const BuildContext& /*ctx*/) {}

void DFTB2::launch(const LaunchContext& ctx) {
    auto& f = ctx.frame;
    if (f.nlocal() < 1)
        throw std::invalid_argument("DFTB2: empty frame");

    if (!(mix_beta > 0.0) || mix_beta > 1.0)
        throw std::invalid_argument("DFTB2: mix_beta must be in (0, 1]");

    last_stress_.fill(0.0);
    for (int d = 0; d < 3; ++d)
        if (k_mesh[static_cast<std::size_t>(d)] < 1)
            throw std::invalid_argument(
                "DFTB2: k_mesh subdivision must be >= 1 (axis "
                + std::to_string(d) + " got "
                + std::to_string(k_mesh[static_cast<std::size_t>(d)]) + ")");
    const bool multi_k =
        k_mesh[0] > 1 || k_mesh[1] > 1 || k_mesh[2] > 1;
    if (sk_set == DftbSkSet::ThreeOb) {
        if (multi_k)
            throw std::runtime_error(
                "DFTB2(ThreeOb): multi-k not supported (MioSketch only)");
        if (multiplicity != 1 || charge != 0)
            throw std::runtime_error(
                "DFTB2(ThreeOb): open-shell/charge not supported on this path");
        energy_ = launch_3ob(f, scc_tol, max_scc, mix_beta, with_forces,
                             with_stress, last_dq_, last_stress_,
                             scc_residual_, scc_iterations_);
        last_spin_.assign(static_cast<std::size_t>(f.nlocal()), 0.0);
    } else {
        const int mesh[3] = {k_mesh[0], k_mesh[1], k_mesh[2]};
        energy_ = launch_mio(f, scc_tol, max_scc, mix_beta, with_forces,
                             with_stress, charge, multiplicity, mesh,
                             k_time_reversal, last_dq_, last_spin_,
                             last_stress_, scc_residual_, scc_iterations_);
    }

    // Optional Grimme D3(BJ): a purely geometric additive correction, so it
    // rides on top of whichever SK path ran above and needs no SCC re-entry.
    // Its own analytic gradient keeps F = -grad(E) for the total, independent of
    // whether the electronic forces above were analytic (closed-shell mio) or
    // finite-difference (open shell / 3ob).
    e_disp_ = 0.0;
    if (d3) {
        const int n_at = f.nlocal();
        const int* type = f.span(kAtom_type).data();
        const StorageReal* x = f.span(kAtom_x);
        const StorageReal* y = f.span(kAtom_y);
        const StorageReal* z = f.span(kAtom_z);

        std::vector<int> Z(static_cast<std::size_t>(n_at));
        std::vector<std::array<double, 3>> pos(static_cast<std::size_t>(n_at));
        for (int a = 0; a < n_at; ++a) {
            const int zi = type[a];
            if (!gpw::D3Params::supports_element(zi))
                throw std::runtime_error(
                    "DFTB2 D3: element Z=" + std::to_string(zi)
                    + " is not in the in-tree D3 table (no fallback)");
            Z[static_cast<std::size_t>(a)] = zi;
            pos[static_cast<std::size_t>(a)] = {
                static_cast<double>(x[a]) * units::kAngstromInBohr,
                static_cast<double>(y[a]) * units::kAngstromInBohr,
                static_cast<double>(z[a]) * units::kAngstromInBohr};
        }

        // Isolated cluster: the D3 sum takes no lattice images even when the
        // electronic Gamma-point path above is lattice-summed (documented
        // limitation on DFTB2::d3 -- periodic DFTB2-D3 is not implemented).
        const double cell[3] = {0.0, 0.0, 0.0};
        const gpw::D3Result res = gpw::dftd3_pbe_bj(
            Z, pos, cell, /*periodic=*/false, gpw::D3Params{},
            /*rcut_bohr=*/94.2, d3_s6, d3_s8, d3_a1, d3_a2,
            /*cp2k_cn=*/false, /*calculate_c9=*/false);

        e_disp_ = res.energy * units::kHartreeInEv;
        energy_ += e_disp_;
        // Total only: dispersion is an empirical geometric term, not part of the
        // electronic share reported on kEnergyElectronic.
        f.meta.get(kEnergy) += e_disp_;

        if (with_forces) {
            // Ha/Bohr -> eV/Angstrom.
            const double f_conv =
                units::kHartreeInEv * units::kAngstromInBohr;
            AccReal* fx = f.span(kAtom_fx);
            AccReal* fy = f.span(kAtom_fy);
            AccReal* fz = f.span(kAtom_fz);
            for (int a = 0; a < n_at; ++a) {
                const auto& fa = res.forces[static_cast<std::size_t>(a)];
                fx[a] += static_cast<AccReal>(fa[0] * f_conv);
                fy[a] += static_cast<AccReal>(fa[1] * f_conv);
                fz[a] += static_cast<AccReal>(fa[2] * f_conv);
            }
        }
    }
}

} // namespace semiempirical
} // namespace cpu
} // namespace ATV
