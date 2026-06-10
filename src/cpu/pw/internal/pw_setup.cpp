/// @file cpu/pw/internal/pw_setup.cpp
/// @brief Host-side plane-wave SCF problem assembly (moved verbatim from
///        the pre-loop section of pw_scf.cpp) shared by both backends.

#include "atomiverse/cpu/pw/internal/pw_setup.h"
#include "atomiverse/core/box.h"
#include "atomiverse/core/mat3.h"
#include "atomiverse/cpu/pw/internal/ewald.h"
#include "atomiverse/cpu/pw/internal/upf_reader.h"
#include <algorithm>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <memory>

namespace ATV {
namespace cpu {
namespace pw {
namespace internal {

namespace {

/// Radial → G-space transform of the local pseudopotential, QE-style
/// (vloc_mod.f90): V_loc(r) → −Z/r at large r, so the bare integrand
/// r²·V·j0(Gr) does not decay and the quadrature diverges. Add and
/// subtract an erf-screened Coulomb: the difference V + Z·erf(r)/r
/// decays like erfc, and the screened Coulomb has the analytic
/// transform −(4πZ/ΩG²)·exp(−G²/4).
///
/// The G=0 element gets the finite part only (QE's "αZ" term):
/// V(0) = (4π/Ω)∫ r²(V(r) + Z/r) dr — the divergent jellium piece is
/// cancelled against the Ewald/Hartree G=0 conventions exactly as for
/// the toy potential. All outputs in Hartree.
std::vector<Real> upf_vloc_to_gspace(const internal::UPFData& upf,
                                     const internal::GVectorList& gv,
                                     Real omega) {
    std::vector<Real> vg(gv.count(), Real(0));
    const auto& r = upf.r;
    const auto& vl = upf.v_loc;
    Real inv_omega = Real(1) / omega;
    Real ry_to_ha = Real(0.5); // UPF stores V_loc in Ry
    Real z = static_cast<Real>(upf.z_val);

    for (int ig = 0; ig < gv.count(); ++ig) {
        Real G = std::sqrt(gv.g[ig].norm_squared());
        Real g2 = gv.g[ig].norm_squared();
        bool g_zero = g2 < Real(1e-16);
        Real sum = 0;
        for (size_t ir = 1; ir < r.size(); ++ir) {
            Real r_mid = Real(0.5) * (r[ir] + r[ir - 1]);
            Real dr = r[ir] - r[ir - 1];
            Real v_mid = Real(0.5) * (vl[ir] + vl[ir - 1]) * ry_to_ha;
            if (r_mid < Real(1e-12)) {
                continue; // r²·(V + Z/r) → 0 at the origin
            }
            if (g_zero) {
                // r²(V + Z/r): exactly zero beyond the core radius.
                sum += r_mid * (r_mid * v_mid + z) * dr;
            } else {
                Real Gr = G * r_mid;
                Real j0 = (Gr < Real(1e-8)) ? Real(1) - Gr * Gr / Real(6)
                                            : std::sin(Gr) / Gr;
                Real v_scr = r_mid * v_mid + z * std::erf(r_mid);
                sum += r_mid * v_scr * j0 * dr;
            }
        }
        if (g_zero)
            vg[ig] = Real(4.0 * M_PI) * sum * inv_omega;
        else
            vg[ig] = Real(4.0 * M_PI) * inv_omega
                     * (sum - z * std::exp(-g2 / Real(4)) / g2);
    }
    return vg;
}

void setup_upf_vloc(PWWorkspace& ws, const internal::GVectorList& gv,
                    const internal::UPFData& upf, Real omega) {
    auto vg = upf_vloc_to_gspace(upf, gv, omega);
    // Atom form factor × structure factor Σ_a e^{-iG·R_a}.
    // vg[0] already holds only the finite part (QE "αZ" term); the
    // divergent −4πZ/ΩG² jellium piece is cancelled against the
    // Ewald/Hartree G=0 conventions, exactly as for the toy potential.
    const auto& pos = ws.atom_pos_bohr;
    for (int ig = 0; ig < ws.n_g; ++ig) {
        Real s_re = 0, s_im = 0;
        for (size_t ia = 0; ia < pos.size(); ++ia) {
            Real phase = gv.g[ig][0] * pos[ia][0] + gv.g[ig][1] * pos[ia][1]
                         + gv.g[ig][2] * pos[ia][2];
            s_re += std::cos(phase);
            s_im -= std::sin(phase);
        }
        ws.V_loc_G[ig] = vg[ig] * std::complex<Real>(s_re, s_im);
    }
}

void setup_toy_vloc(PWWorkspace& ws, const internal::GVectorList& gv,
                    const std::vector<Vec3>& positions,
                    const std::vector<Real>& charges, Real omega) {
    for (int ig = 0; ig < ws.n_g; ++ig) {
        Real g2 = gv.g[ig].norm_squared();
        Vec3 gvec = gv.g[ig];
        if (g2 < Real(1e-16)) {
            ws.V_loc_G[ig] = Real(0);
        } else {
            Real s_re = 0, s_im = 0;
            for (size_t ia = 0; ia < charges.size(); ++ia) {
                Real phase = gvec[0] * positions[ia][0]
                             + gvec[1] * positions[ia][1]
                             + gvec[2] * positions[ia][2];
                s_re += charges[ia] * std::cos(phase);
                s_im -= charges[ia] * std::sin(phase);
            }
            ws.V_loc_G[ig] = Real(-4.0 * M_PI) / (omega * g2)
                             * std::complex<Real>(s_re, s_im);
        }
    }
}

::ATV::Box to_core_box(const cpu::Frame& frame) {
    const auto& c = frame.box.cell();
    Mat3 cell;
    cell(0, 0) = c(0, 0);
    cell(0, 1) = c(0, 1);
    cell(0, 2) = c(0, 2);
    cell(1, 0) = c(1, 0);
    cell(1, 1) = c(1, 1);
    cell(1, 2) = c(1, 2);
    cell(2, 0) = c(2, 0);
    cell(2, 1) = c(2, 1);
    cell(2, 2) = c(2, 2);
    ::ATV::Box box;
    box.cell = cell;
    return box;
}

::ATV::Box box_to_bohr(const ::ATV::Box& box_ang) {
    Real a2b = Real(1.8897259886);
    ::ATV::Box box = box_ang;
    box.cell = box.cell * a2b;
    return box;
}

std::vector<Vec3> pos_to_bohr(const cpu::Frame& frame) {
    Real a2b = Real(1.8897259886);
    std::vector<Vec3> pos;
    for (size_t i = 0; i < frame.atoms.size(); ++i)
        pos.push_back({static_cast<Real>(frame.atoms.x[i]) * a2b,
                       static_cast<Real>(frame.atoms.y[i]) * a2b,
                       static_cast<Real>(frame.atoms.z[i]) * a2b});
    return pos;
}

// Spherical Bessel functions j_l(x)
Real sph_bessel_j0(Real x) {
    if (std::abs(x) < Real(1e-8))
        return Real(1) - x * x / Real(6);
    return std::sin(x) / x;
}
Real sph_bessel_j1(Real x) {
    if (std::abs(x) < Real(1e-8))
        return x / Real(3);
    return (std::sin(x) / (x * x)) - (std::cos(x) / x);
}
Real sph_bessel_j2(Real x) {
    if (std::abs(x) < Real(1e-8))
        return x * x / Real(15);
    Real x2 = x * x, x3 = x * x2;
    return (Real(3) / x3 - Real(1) / x) * std::sin(x)
           - (Real(3) / x2) * std::cos(x);
}

using SphBesselFn = Real (*)(Real);
SphBesselFn sph_bessel_for_l(int l) {
    switch (l) {
    case 0:
        return sph_bessel_j0;
    case 1:
        return sph_bessel_j1;
    case 2:
        return sph_bessel_j2;
    default:
        return sph_bessel_j0;
    }
}

/// Build Kleinman-Bylander non-local projectors in G-space for one
/// k-sphere. The radial table is evaluated at q = |k+G| — the projector
/// of a Bloch state lives on k+G, so the table is k-dependent and must
/// be rebuilt per k-point (size gv.count(), which also differs per k).
void setup_kb_projectors(std::vector<std::vector<PWWorkspace::KBProjG>>& out,
                         const internal::GVectorList& gv,
                         const internal::UPFData& upf,
                         const std::vector<Vec3>& pos_bohr) {
    int n_atoms = static_cast<int>(pos_bohr.size());
    int ng = gv.count();
    out.assign(n_atoms, {});

    // Precompute f_l(q) for each projector (radial integral, same for all atoms)
    int n_proj = upf.n_proj();
    std::vector<std::vector<Real>> fl_all(n_proj);
    std::vector<int> proj_l(n_proj);
    for (size_t ip = 0; ip < upf.kb.size(); ++ip) {
        auto& kb = upf.kb[ip];
        proj_l[ip] = kb.l;
        auto j_l = sph_bessel_for_l(kb.l);
        fl_all[ip].resize(ng, Real(0));
        const auto& r = upf.r;
        const auto& beta_r = kb.beta;
        for (int ig = 0; ig < ng; ++ig) {
            Vec3 q{gv.g[ig][0] + gv.kvec[0], gv.g[ig][1] + gv.kvec[1],
                   gv.g[ig][2] + gv.kvec[2]};
            Real Q = std::sqrt(q.norm_squared());
            Real sum = 0;
            for (size_t ir = 1; ir < r.size(); ++ir) {
                Real r_mid = Real(0.5) * (r[ir] + r[ir - 1]);
                Real dr = r[ir] - r[ir - 1];
                Real b_mid = Real(0.5) * (beta_r[ir] + beta_r[ir - 1]);
                sum += r_mid * b_mid * j_l(Q * r_mid) * dr;
            }
            fl_all[ip][ig] = sum;
        }
    }

    // Diagonalize D_ij per l-block, transform projectors to eigenbasis
    // This gives the correct effective projector strengths
    std::vector<Real> eff_D(n_proj, Real(0));
    std::vector<std::vector<Real>> fl_trans(n_proj); // transformed f_l(G)
    for (int ip = 0; ip < n_proj; ++ip)
        fl_trans[ip] = fl_all[ip];

    if (!upf.d_ij.empty()) {
        // Group projectors by l
        std::vector<int> l_groups[4]; // l=0,1,2,3
        for (int ip = 0; ip < n_proj; ++ip)
            if (proj_l[ip] < 4)
                l_groups[proj_l[ip]].push_back(ip);

        int ip_out = 0;
        for (int l = 0; l < 4; ++l) {
            int nl = static_cast<int>(l_groups[l].size());
            if (nl == 0)
                continue;

            // Extract D sub-block for this l
            std::vector<Real> D_sub(nl * nl);
            for (int i = 0; i < nl; ++i)
                for (int j = 0; j < nl; ++j) {
                    int ii = l_groups[l][i], jj = l_groups[l][j];
                    D_sub[i * nl + j] =
                        (ii * n_proj + jj < static_cast<int>(upf.d_ij.size()))
                            ? upf.d_ij[ii * n_proj + jj]
                            : Real(0);
                }

            // Diagonalize small matrix (nl ≤ 2, use analytic formula)
            std::vector<Real> eig(nl), vecs(nl * nl);
            if (nl == 1) {
                eig[0] = D_sub[0];
                vecs[0] = Real(1);
            } else {
                Real a = D_sub[0], b = D_sub[1], c = D_sub[2], d = D_sub[3];
                if (std::abs(b) > Real(1e-16)) {
                    Real tr = a + d, det = a * d - b * c,
                         disc = std::sqrt(
                             std::max(Real(0), tr * tr - Real(4) * det));
                    eig[0] = Real(0.5) * (tr - disc);
                    eig[1] = Real(0.5) * (tr + disc);
                    // Eigenvectors, column-major: vecs[i + nl*k]
                    Real n0 = std::sqrt(b * b + (eig[0] - a) * (eig[0] - a));
                    vecs[0 * nl + 0] = (eig[0] - a) / n0;
                    vecs[1 * nl + 0] = b / n0; // col 0
                    Real n1 = std::sqrt(b * b + (eig[1] - a) * (eig[1] - a));
                    vecs[0 * nl + 1] = (eig[1] - a) / n1;
                    vecs[1 * nl + 1] = b / n1; // col 1
                } else {
                    // Already diagonal: keep the ORIGINAL pairing. The old
                    // code sorted the eigenvalues here while leaving the
                    // eigenvectors as identity, swapping each projector's
                    // strength with its partner whenever D11 > D22 (true
                    // for the ONCV s and p channels — the l=2 block
                    // happened to be ascending and masked the bug).
                    eig[0] = a;
                    eig[1] = d;
                    vecs[0] = 1;
                    vecs[1] = 0;
                    vecs[2] = 0;
                    vecs[3] = 1;
                }
            }

            // Transform f_l(G): f'_k(G) = Σ_i V_{i,k} * f_i(G)
            // Store in fl_trans
            for (int k = 0; k < nl; ++k) {
                fl_trans[ip_out].resize(ng);
                for (int ig = 0; ig < ng; ++ig)
                    fl_trans[ip_out][ig] = Real(0);
                for (int i = 0; i < nl; ++i) {
                    Real v = vecs[i * nl + k];
                    for (int ig = 0; ig < ng; ++ig)
                        fl_trans[ip_out][ig] += v * fl_all[l_groups[l][i]][ig];
                }
                eff_D[ip_out] = eig[k] * Real(0.5); // Ry → Ha
                proj_l[ip_out] = l;
                ++ip_out;
            }
        }
        // Update number of effective projectors
        n_proj = ip_out;
    }

    auto ng_sz = static_cast<std::size_t>(ng);
    for (int ia = 0; ia < n_atoms; ++ia) {
        out[ia].resize(n_proj);
        for (int ip = 0; ip < n_proj; ++ip) {
            out[ia][ip].l = proj_l[ip];
            out[ia][ip].d = eff_D[ip];
            out[ia][ip].f_G.resize({ng_sz});
            for (int ig = 0; ig < ng; ++ig)
                out[ia][ip].f_G[ig] = fl_trans[ip][ig];
        }
    }
}

/// Build atomic wavefunctions in G-space for the initial guess on one
/// k-sphere; the radial transform is evaluated at q = |k+G|
/// (QE atwfc_mod.f90: chi_G = (4π/√Ω) ∫ χ(r)·j_l(qr)·r dr).
void setup_atomic_wfc(std::vector<std::vector<PWWorkspace::AtomWfc>>& out,
                      const internal::GVectorList& gv,
                      const internal::UPFData& upf, int n_atoms, Real omega) {
    int ng = gv.count();
    Real pref = Real(4.0 * M_PI) / std::sqrt(omega);

    std::vector<std::vector<Real>> chi_g_all(upf.nwfc);
    for (int iw = 0; iw < upf.nwfc; ++iw) {
        int l = upf.lchi[iw];
        auto j_l = sph_bessel_for_l(l);
        const auto& r = upf.r;
        const auto& chi_r = upf.chi[iw];

        chi_g_all[iw].resize(ng, Real(0));
        for (int ig = 0; ig < ng; ++ig) {
            Vec3 q{gv.g[ig][0] + gv.kvec[0], gv.g[ig][1] + gv.kvec[1],
                   gv.g[ig][2] + gv.kvec[2]};
            Real Q = std::sqrt(q.norm_squared());
            Real sum = 0;
            for (size_t ir = 1; ir < r.size(); ++ir) {
                Real r_mid = Real(0.5) * (r[ir] + r[ir - 1]);
                Real dr = r[ir] - r[ir - 1];
                Real chi_mid = Real(0.5) * (chi_r[ir] + chi_r[ir - 1]);
                sum += chi_mid * j_l(Q * r_mid) * r_mid * dr;
            }
            chi_g_all[iw][ig] = sum * pref;
        }
    }

    // Store for each atom (radial part shared; structure factor applied later)
    out.assign(n_atoms, {});
    auto ng_sz = static_cast<std::size_t>(ng);
    for (int ia = 0; ia < n_atoms; ++ia) {
        out[ia].resize(upf.nwfc);
        for (int iw = 0; iw < upf.nwfc; ++iw) {
            out[ia][iw].l = upf.lchi[iw];
            out[ia][iw].chi_G.resize({ng_sz});
            for (int ig = 0; ig < ng; ++ig)
                out[ia][iw].chi_G[ig] = chi_g_all[iw][ig];
        }
    }
}

} // namespace

int valence_from_z(int z) {
    switch (z) {
    case 1:
        return 1;
    case 2:
        return 2;
    case 3:
        return 1;
    case 4:
        return 2;
    case 5:
        return 3;
    case 6:
        return 4;
    case 7:
        return 5;
    case 8:
        return 6;
    case 9:
        return 7;
    case 10:
        return 8;
    case 11:
        return 1;
    case 12:
        return 2;
    case 13:
        return 3;
    case 14:
        return 4;
    case 15:
        return 5;
    case 16:
        return 6;
    case 17:
        return 7;
    case 18:
        return 8;
    default:
        return 0;
    }
}

void setup_pw_problem(const cpu::Frame& frame, const std::string& upf_path,
                      double ecut_ry, int k_mesh, PWProblem& out) {
    int n_atoms = static_cast<int>(frame.atoms.size());
    ::ATV::Box box_ang = to_core_box(frame);
    ::ATV::Box box_bohr = box_to_bohr(box_ang);

    Real ecut_ht = static_cast<Real>(ecut_ry) * Real(0.5);
    Real vol_bohr = box_bohr.volume();
    out.vol_bohr = vol_bohr;

    auto pos_bohr = pos_to_bohr(frame);

    int n_elec = 0;
    std::vector<Real> valence_charges, nuclear_Z;
    for (int i = 0; i < n_atoms; ++i) {
        int z = frame.atoms.type[i];
        int v = valence_from_z(z);
        n_elec += v;
        valence_charges.push_back(static_cast<Real>(v));
        nuclear_Z.push_back(static_cast<Real>(z));
    }
    out.n_elec = n_elec;
    out.n_bands = std::max(n_elec / 2, 1);

    bool has_upf = !upf_path.empty();

    std::cout << "\n========== PW SCF ==========" << std::endl;
    std::cout << "Atoms:" << n_atoms << " Elec:" << n_elec
              << " Bands:" << out.n_bands << " Ecut:" << ecut_ry << " Ry"
              << (has_upf ? " UPF" : " toy") << std::endl;

    // === Ewald ===
    std::vector<Vec3> pos_ang;
    for (int i = 0; i < n_atoms; ++i)
        pos_ang.push_back({static_cast<Real>(frame.atoms.x[i]),
                           static_cast<Real>(frame.atoms.y[i]),
                           static_cast<Real>(frame.atoms.z[i])});
    // Ewald: use same charges as V_loc (nuclear for toy, valence for UPF)
    auto& ewald_charges = has_upf ? valence_charges : nuclear_Z;
    internal::EwaldSummation ewald(pos_ang, ewald_charges, box_ang, ecut_ht);
    Real ew_total = ewald.total_energy();

    // Charged-system / jellium G=0 term: -π Q² / (2 α² Ω). NB: QE's
    // ewald.f90 writes -charge²/(4·alpha)·(2π/Ω), but QE's `alpha` is the
    // Gaussian exponent η² — our EwaldSummation::alpha() is η, so the
    // denominator needs α². (With α¹ the term acquires a spurious
    // cell-size dependence that tilted the Si E(V) curve by ~0.04 Ry
    // across the EOS scan — comparable to the entire physical span.)
    Real ew_charge = 0;
    for (auto& z : ewald_charges)
        ew_charge += z;
    Real ew_alpha = ewald.alpha();
    Real ew_g0 = -ew_charge * ew_charge / (Real(4) * ew_alpha * ew_alpha)
                 * Real(2.0 * M_PI) / vol_bohr;
    ew_total += ew_g0;
    out.ew_total = ew_total;

    std::cout << "\n--- Ewald (Ha) ---" << std::endl;
    std::cout << "  real=" << ewald.real_energy()
              << " recip=" << ewald.recip_energy()
              << " self=" << ewald.self_energy() << std::endl;
    std::cout << "  total=" << ew_total << " Ha (" << ew_total * 2.0 << " Ry)"
              << std::endl;

    // === G-vectors ===
    // gv0 is the DENSITY sphere: ρ = |ψ|² has spectral content out to
    // 2·G_max(wfc), i.e. ecutrho = 4·ecut (QE default). All density-side
    // quantities (ρ, V_H, V_xc, V_loc, V_eff, mixing, energy sums) live
    // on this sphere; the per-k wavefunction spheres stay at ecut.
    out.gv0 =
        internal::enumerate_g_vectors(box_bohr, Real(4) * ecut_ht,
                                      Vec3::zeros());
    auto& gv0 = out.gv0;
    std::cout << "\n--- G-vectors ---" << std::endl;
    std::cout << "  N_G=" << gv0.count()
              << " G_min=" << std::sqrt(gv0.g[1].norm_squared())
              << " G_max=" << std::sqrt(gv0.g[gv0.count() - 1].norm_squared())
              << std::endl;
    std::cout << "  FFT dims(hkl): " << gv0.nx << "x" << gv0.ny << "x" << gv0.nz
              << std::endl;

    // k-points
    out.kg = internal::generate_mp_mesh(k_mesh, k_mesh, k_mesh);

    // === Workspace ===
    PWWorkspace& ws = out.ws;
    ws.n_g = gv0.count();
    ws.gvec = &gv0;
    auto ng = static_cast<std::size_t>(ws.n_g);
    ws.rho_G.resize({ng});
    ws.rho_G_old.resize({ng});
    ws.V_H_G.resize({ng});
    ws.V_xc_G.resize({ng});
    ws.V_eff_G.resize({ng});
    ws.V_loc_G.resize({ng});
    ws.T_G.resize({ng});
    ws.mix_beta = Real(0.7); // QE default: mixing_beta = 0.7
    ws.omega = vol_bohr;     // cell volume for real-space integration
    ws.atom_pos_bohr =
        pos_bohr; // atom positions for augmentation structure factor

    // Initial density: read UPF first if needed
    internal::UPFData upf;
    if (has_upf)
        upf = internal::read_upf_v1(upf_path);

    for (int i = 0; i < ws.n_g; ++i)
        ws.rho_G[i] = std::complex<Real>(0, 0);

    if (has_upf && !upf.rho_atom.empty()
        && upf.rho_atom.size() == upf.r.size()) {
        // G-space transform of the atomic pseudocharge. UPF PP_RHOATOM
        // stores 4πr²ρ(r) (integrates to Z_val with plain dr), so
        // ρ_atom(G) = (1/Ω) ∫ PP_RHOATOM(r) j₀(Gr) dr.
        for (int ig = 0; ig < ws.n_g; ++ig) {
            Real G = std::sqrt(gv0.g[ig].norm_squared());
            Real sum = 0;
            for (size_t ir = 1; ir < upf.r.size(); ++ir) {
                Real r_mid = Real(0.5) * (upf.r[ir] + upf.r[ir - 1]);
                Real dr = upf.r[ir] - upf.r[ir - 1];
                Real rho_mid =
                    Real(0.5) * (upf.rho_atom[ir] + upf.rho_atom[ir - 1]);
                Real Gr = G * r_mid;
                Real j0 = (Gr < Real(1e-8)) ? Real(1) - Gr * Gr / Real(6)
                                            : std::sin(Gr) / Gr;
                sum += rho_mid * j0 * dr;
            }
            Real rho_g_atom = sum / vol_bohr;
            // Structure factor: Σ_a exp(-i G·R_a)
            Real s_re = 0, s_im = 0;
            for (int ia = 0; ia < n_atoms; ++ia) {
                Real phase = gv0.g[ig][0] * pos_bohr[ia][0]
                             + gv0.g[ig][1] * pos_bohr[ia][1]
                             + gv0.g[ig][2] * pos_bohr[ia][2];
                s_re += std::cos(phase);
                s_im -= std::sin(phase);
            }
            ws.rho_G[ig] = rho_g_atom * std::complex<Real>(s_re, s_im);
        }
        // Renormalize: total charge must be n_elec
        Real total_charge = ws.rho_G[0].real() * vol_bohr;
        if (std::abs(total_charge) > Real(1e-16)) {
            Real rescale = Real(n_elec) / total_charge;
            for (int ig = 0; ig < ws.n_g; ++ig)
                ws.rho_G[ig] *= rescale;
        }
    } else {
        // Fallback: uniform density
        ws.rho_G[0] = std::complex<Real>(Real(n_elec) / vol_bohr, 0);
    }

    // Density FFT grid: gv0 is already the ecutrho = 4·ecut sphere, so
    // 2× its index extent covers it (equals the old 4× wfc-sphere grid).
    ws.fft_nx = 2 * gv0.nx;
    ws.fft_ny = 2 * gv0.ny;
    ws.fft_nz = 2 * gv0.nz;
    if (ws.fft_nx % 2)
        ++ws.fft_nx;
    if (ws.fft_ny % 2)
        ++ws.fft_ny;
    if (ws.fft_nz % 2)
        ++ws.fft_nz;
    ws.nr = ws.fft_nx * ws.fft_ny * ws.fft_nz;
    ws.nc = (ws.fft_nx / 2 + 1) * ws.fft_ny * ws.fft_nz;
    ws.rho_r.resize({static_cast<std::size_t>(ws.nr)});
    ws.V_xc_r.resize({static_cast<std::size_t>(ws.nr)});
    ws.fft =
        std::make_unique<internal::FFTPlan3D>(ws.fft_nx, ws.fft_ny, ws.fft_nz);

    std::cout << "  Density FFT grid: " << ws.fft_nx << "x" << ws.fft_ny << "x"
              << ws.fft_nz << " (nr=" << ws.nr << ")" << std::endl;

    // === V_loc ===
    if (has_upf) {
        std::cout << "\n--- UPF: Z_val=" << upf.z_val << " l_max=" << upf.l_max
                  << " n_proj=" << upf.kb.size() << " ---" << std::endl;
        setup_upf_vloc(ws, gv0, upf, vol_bohr);

        // Setup augmentation form factor: rho_atom(G) / Q_atom
        ws.q_ij = upf.q_ij;
        ws.n_proj_total = upf.n_proj();
        ws.aug_form_factor_G.resize({static_cast<std::size_t>(ws.n_g)});
        if (!upf.rho_atom.empty() && upf.rho_atom.size() == upf.r.size()) {
            // PP_RHOATOM stores 4πr²ρ(r): Q = ∫ PP_RHOATOM dr,
            // form factor (G) = (1/ΩQ) ∫ PP_RHOATOM j₀(Gr) dr.
            Real Q_atom = 0;
            for (size_t ir = 1; ir < upf.r.size(); ++ir) {
                Real dr = upf.r[ir] - upf.r[ir - 1];
                Real rho_mid =
                    Real(0.5) * (upf.rho_atom[ir] + upf.rho_atom[ir - 1]);
                Q_atom += rho_mid * dr;
            }
            Real inv_Q =
                (std::abs(Q_atom) > Real(1e-16)) ? Real(1) / Q_atom : Real(0);
            for (int ig = 0; ig < ws.n_g; ++ig) {
                Real G = std::sqrt(gv0.g[ig].norm_squared());
                Real sum = 0;
                for (size_t ir = 1; ir < upf.r.size(); ++ir) {
                    Real r_mid = Real(0.5) * (upf.r[ir] + upf.r[ir - 1]);
                    Real dr = upf.r[ir] - upf.r[ir - 1];
                    Real rho_mid =
                        Real(0.5) * (upf.rho_atom[ir] + upf.rho_atom[ir - 1]);
                    Real Gr = G * r_mid;
                    Real j0 = (Gr < Real(1e-8)) ? Real(1) - Gr * Gr / Real(6)
                                                : std::sin(Gr) / Gr;
                    sum += rho_mid * j0 * dr;
                }
                ws.aug_form_factor_G[ig] = sum * inv_Q / vol_bohr;
            }
        } else {
            for (int ig = 0; ig < ws.n_g; ++ig)
                ws.aug_form_factor_G[ig] = Real(0);
        }
    } else {
        setup_toy_vloc(ws, gv0, pos_bohr, valence_charges, vol_bohr);
    }
    std::cout << "  V_loc(G) [1]=" << ws.V_loc_G[1].real()
              << " [10]=" << ws.V_loc_G[10].real() << std::endl;

    // Pre-compute k-point G-vectors and the per-k nonlocal / wfc tables
    // (radial parts live on q = k+G, so size and values are k-dependent).
    out.gv_k.clear();
    out.gv_k.reserve(out.kg.kpoints.size());
    {
        // MP mesh k-points are fractional (crystal) coordinates;
        // enumerate_g_vectors needs Cartesian [1/Bohr].
        auto recip = box_bohr.recip_cell();
        for (auto& kf : out.kg.kpoints) {
            Vec3 kc{recip(0, 0) * kf[0] + recip(1, 0) * kf[1]
                        + recip(2, 0) * kf[2],
                    recip(0, 1) * kf[0] + recip(1, 1) * kf[1]
                        + recip(2, 1) * kf[2],
                    recip(0, 2) * kf[0] + recip(1, 2) * kf[1]
                        + recip(2, 2) * kf[2]};
            out.gv_k.push_back(
                internal::enumerate_g_vectors(box_bohr, ecut_ht, kc));
        }
    }

    out.kb_proj_k.clear();
    out.atom_wfc_k.clear();
    if (has_upf) {
        auto n_k = out.gv_k.size();
        out.kb_proj_k.resize(n_k);
        out.atom_wfc_k.resize(n_k);
        for (size_t ik = 0; ik < n_k; ++ik) {
            setup_kb_projectors(out.kb_proj_k[ik], out.gv_k[ik], upf,
                                pos_bohr);
            if (upf.nwfc > 0)
                setup_atomic_wfc(out.atom_wfc_k[ik], out.gv_k[ik], upf,
                                 n_atoms, vol_bohr);
        }
        std::cout << "  KB projectors: " << upf.kb.size() << " proj, "
                  << n_atoms << " atoms, " << n_k << " k-spheres"
                  << " (atomic wfc: " << upf.nwfc << ")" << std::endl;
    }
}

void pulay_mix(PWWorkspace& ws) {
    // Pulay/DIIS mixing (QE's Anderson/Broyden simplified)
    // Store residual F = ρ_out - ρ_in
    xt::xarray<std::complex<Real>> F = ws.rho_G - ws.rho_G_old;
    ws.pulay_res.push_back(F);
    ws.pulay_rho.push_back(ws.rho_G_old);

    int depth = std::min(ws.pulay_depth, static_cast<int>(ws.pulay_res.size()));
    if (depth > 8) {
        ws.pulay_res.erase(ws.pulay_res.begin());
        ws.pulay_rho.erase(ws.pulay_rho.begin());
        depth = 8;
    }

    if (depth >= 2) {
        // Build B_{ij} = ⟨F_i|F_j⟩
        std::vector<Real> B(depth * depth, Real(0));
        for (int i = 0; i < depth; ++i) {
            for (int j = i; j < depth; ++j) {
                Real dot = 0;
                for (int ig = 0; ig < ws.n_g; ++ig) {
                    auto& fi = ws.pulay_res[i];
                    auto& fj = ws.pulay_res[j];
                    dot += fi[ig].real() * fj[ig].real()
                           + fi[ig].imag() * fj[ig].imag();
                }
                B[i * depth + j] = B[j * depth + i] = dot;
            }
        }
        // Solve B·c = 1 (with constraint Σc=1 via Lagrangian)
        // Augmented system: [B  1] [c] = [0]
        //                    [1  0] [λ]   [1]
        int n = depth + 1;
        std::vector<Real> A(n * n, Real(0));
        for (int i = 0; i < depth; ++i)
            for (int j = 0; j < depth; ++j)
                A[i * n + j] = B[i * depth + j];
        for (int i = 0; i < depth; ++i)
            A[i * n + depth] = A[depth * n + i] = Real(1);

        std::vector<Real> rhs(n, Real(0));
        rhs[depth] = Real(1);

        // Solve via Gaussian elimination
        std::vector<Real> c(n, Real(0));
        auto aug = A;
        auto b = rhs;
        for (int k = 0; k < n; ++k) {
            // Pivot
            int imax = k;
            for (int i = k + 1; i < n; ++i)
                if (std::abs(aug[i * n + k]) > std::abs(aug[imax * n + k]))
                    imax = i;
            if (imax != k) {
                for (int j = 0; j < n; ++j)
                    std::swap(aug[k * n + j], aug[imax * n + j]);
                std::swap(b[k], b[imax]);
            }
            if (std::abs(aug[k * n + k]) < Real(1e-14))
                continue;
            for (int i = k + 1; i < n; ++i) {
                Real factor = aug[i * n + k] / aug[k * n + k];
                for (int j = k; j < n; ++j)
                    aug[i * n + j] -= factor * aug[k * n + j];
                b[i] -= factor * b[k];
            }
        }
        for (int k = n - 1; k >= 0; --k) {
            Real sum = b[k];
            for (int j = k + 1; j < n; ++j)
                sum -= aug[k * n + j] * c[j];
            if (std::abs(aug[k * n + k]) > Real(1e-14))
                c[k] = sum / aug[k * n + k];
        }
        // Build mixed density: ρ_new = Σ c_i (ρ_in^(i) + α·F^(i))
        Real alpha = ws.mix_beta;
        for (int ig = 0; ig < ws.n_g; ++ig)
            ws.rho_G[ig] = std::complex<Real>(0, 0);
        for (int i = 0; i < depth; ++i) {
            for (int ig = 0; ig < ws.n_g; ++ig)
                ws.rho_G[ig] +=
                    c[i] * (ws.pulay_rho[i][ig] + alpha * ws.pulay_res[i][ig]);
        }
    } else {
        // Simple linear mixing for first iteration
        for (int i = 0; i < ws.n_g; ++i)
            ws.rho_G[i] = ws.rho_G_old[i] + ws.mix_beta * F[i];
    }
}

Real report_iteration_energy(const PWWorkspace& ws, Real e_band, Real ew_total,
                             Real vol_bohr, int iter, int n_bands) {
    // Energy decomposition
    // Parseval: ∫ρ(r)v(r)d³r = Ω·Σ_G Re[ρ*(G)·v(G)] — both real AND
    // imaginary parts contribute (ρ_G is complex for any cell without
    // inversion symmetry about the origin; taking Re·Re only silently
    // dropped the Im·Im half and tilted E(V)).
    Real e_h = 0, e_xc_int = 0, e_loc = 0;
    for (int i = 0; i < ws.n_g; ++i) {
        std::complex<Real> rh = ws.rho_G[i];
        e_h += Real(0.5) * std::real(std::conj(rh) * ws.V_H_G[i]);
        e_xc_int += std::real(std::conj(rh) * ws.V_xc_G[i]);
        e_loc += std::real(std::conj(rh) * ws.V_loc_G[i]);
    }
    e_loc *= vol_bohr;
    e_h *= vol_bohr;
    e_xc_int *= vol_bohr;

    // one-electron = Σ ε_i - ∫ρ V_H - ∫ρ V_xc
    Real e_one = e_band - Real(2.0) * e_h - e_xc_int;
    // total = one-electron + Hartree + XC + Ewald + augmentation
    Real e_total_ha = e_one + e_h + ws.E_xc + ew_total + ws.E_aug;

    std::cout << "\n-- iter " << (iter + 1) << " --" << std::endl;
    std::cout << "  one-electron = " << std::setw(12) << std::setprecision(6)
              << e_one << " Ha (" << e_one * 2.0 << " Ry)" << std::endl;
    std::cout << "  Hartree      = " << std::setw(12) << e_h << " Ha ("
              << e_h * 2.0 << " Ry)" << std::endl;
    std::cout << "  E_loc(rho)   = " << std::setw(12) << e_loc << " Ha ("
              << e_loc * 2.0 << " Ry)  [T+E_NL = "
              << (e_one - e_loc) << " Ha]" << std::endl;
    std::cout << "  XC           = " << std::setw(12) << ws.E_xc << " Ha ("
              << ws.E_xc * 2.0 << " Ry)" << std::endl;
    std::cout << "  Ewald        = " << std::setw(12) << ew_total << " Ha ("
              << ew_total * 2.0 << " Ry)" << std::endl;
    std::cout << "  Aug          = " << std::setw(12) << ws.E_aug << " Ha ("
              << ws.E_aug * 2.0 << " Ry)" << std::endl;
    std::cout << "  Total        = " << std::setw(12) << e_total_ha << " Ha ("
              << e_total_ha * 2.0 << " Ry)" << std::endl;
    std::cout << "  E_band=" << e_band << " Ha  eig:";
    for (int b = 0; b < std::min(4, n_bands); ++b)
        std::cout << " " << ws.eigenvalues[b];
    std::cout << ")" << std::endl;

    return e_total_ha;
}

} // namespace internal
} // namespace pw
} // namespace cpu
} // namespace ATV
