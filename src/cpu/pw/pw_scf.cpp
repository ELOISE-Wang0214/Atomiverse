/// @file cpu/pw/pw_scf.cpp
/// @brief pw::pw_scf() — generic plane-wave DFT total energy.

#include "atomiverse/cpu/pw/pw_scf.h"
#include "atomiverse/core/box.h"
#include "atomiverse/core/constants.h"
#include "atomiverse/core/mat3.h"
#include "atomiverse/cpu/context.h"
#include "atomiverse/cpu/pw/internal/ewald.h"
#include "atomiverse/cpu/pw/internal/g_vector.h"
#include "atomiverse/cpu/pw/internal/kpoint.h"
#include "atomiverse/cpu/pw/internal/upf_reader.h"
#include "atomiverse/cpu/pw/pw_convergence.h"
#include "atomiverse/cpu/pw/pw_davidson.h"
#include "atomiverse/cpu/pw/pw_density.h"
#include "atomiverse/cpu/pw/pw_hamiltonian.h"
#include "atomiverse/cpu/pw/pw_lda.h"
#include "atomiverse/cpu/pw/pw_poisson.h"
#include "atomiverse/cpu/pw/pw_veff.h"
#include "atomiverse/cpu/pw/pw_workspace.h"
#include <algorithm>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <memory>

namespace ATV {
namespace cpu {
namespace pw {

namespace {

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

std::vector<Real> upf_vloc_to_gspace(const internal::UPFData& upf,
                                     const internal::GVectorList& gv,
                                     Real omega) {
    std::vector<Real> vg(gv.count(), Real(0));
    const auto& r = upf.r;
    const auto& vl = upf.v_loc;
    Real inv_omega = Real(1) / omega;
    Real ry_to_ha = Real(0.5); // UPF stores V_loc in Ry

    for (int ig = 0; ig < gv.count(); ++ig) {
        Real G = std::sqrt(gv.g[ig].norm_squared());
        Real sum = 0;
        for (size_t ir = 1; ir < r.size(); ++ir) {
            Real r_mid = Real(0.5) * (r[ir] + r[ir - 1]);
            Real dr = r[ir] - r[ir - 1];
            Real v_mid = Real(0.5) * (vl[ir] + vl[ir - 1]) * ry_to_ha;
            Real Gr = G * r_mid;
            Real j0 = (Gr < Real(1e-8)) ? Real(1) - Gr * Gr / Real(6)
                                        : std::sin(Gr) / Gr;
            sum += r_mid * r_mid * v_mid * j0 * dr;
        }
        vg[ig] = Real(4.0 * M_PI) * sum * inv_omega;
    }
    return vg;
}

void setup_upf_vloc(PWWorkspace& ws, const internal::GVectorList& gv,
                    const internal::UPFData& upf, Real omega) {
    auto vg = upf_vloc_to_gspace(upf, gv, omega);
    for (int i = 0; i < ws.n_g; ++i)
        ws.V_loc_G_real[i] = vg[i];
    ws.V_loc_G_real[0] = Real(0); // G=0 cancelled against Ewald
}

void setup_toy_vloc(PWWorkspace& ws, const internal::GVectorList& gv,
                    const std::vector<Vec3>& positions,
                    const std::vector<Real>& charges, Real omega) {
    for (int ig = 0; ig < ws.n_g; ++ig) {
        Real g2 = gv.g[ig].norm_squared();
        Vec3 gvec = gv.g[ig];
        if (g2 < Real(1e-16)) {
            ws.V_loc_G_real[ig] = Real(0);
        } else {
            Real s_real = 0;
            for (size_t ia = 0; ia < charges.size(); ++ia) {
                Real phase = gvec[0] * positions[ia][0]
                             + gvec[1] * positions[ia][1]
                             + gvec[2] * positions[ia][2];
                s_real += charges[ia] * std::cos(phase);
            }
            ws.V_loc_G_real[ig] = Real(-4.0 * M_PI) * s_real / (omega * g2);
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

/// Setup Kleinman-Bylander non-local projectors in G-space.
void setup_kb_projectors(PWWorkspace& ws, const internal::GVectorList& gv,
                         const internal::UPFData& upf,
                         const std::vector<Vec3>& pos_bohr, Real omega) {
    int n_atoms = static_cast<int>(pos_bohr.size());
    ws.kb_proj.resize(n_atoms);

    Real inv_sqrt_omega = Real(1) / std::sqrt(omega);
    Real four_pi = Real(4.0 * M_PI);
    Real sqrt_four_pi = std::sqrt(four_pi);

    std::cout << "  KB projectors: " << upf.kb.size() << " proj, " << n_atoms
              << " atoms" << std::endl;

    // Precompute f_l(G) for each projector (radial integral, same for all atoms)
    int n_proj = upf.n_proj();
    std::vector<std::vector<Real>> fl_all(n_proj);
    std::vector<int> proj_l(n_proj);
    for (size_t ip = 0; ip < upf.kb.size(); ++ip) {
        auto& kb = upf.kb[ip];
        proj_l[ip] = kb.l;
        auto j_l = sph_bessel_for_l(kb.l);
        fl_all[ip].resize(ws.n_g, Real(0));
        const auto& r = upf.r;
        const auto& beta_r = kb.beta;
        for (int ig = 0; ig < ws.n_g; ++ig) {
            Real G = std::sqrt(gv.g[ig].norm_squared());
            Real sum = 0;
            for (size_t ir = 1; ir < r.size(); ++ir) {
                Real r_mid = Real(0.5) * (r[ir] + r[ir - 1]);
                Real dr = r[ir] - r[ir - 1];
                Real b_mid = Real(0.5) * (beta_r[ir] + beta_r[ir - 1]);
                Real Gr = G * r_mid;
                sum += r_mid * b_mid * j_l(Gr) * dr;
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
                Real tr = a + d, det = a * d - b * c,
                     disc =
                         std::sqrt(std::max(Real(0), tr * tr - Real(4) * det));
                eig[0] = Real(0.5) * (tr - disc);
                eig[1] = Real(0.5) * (tr + disc);
                // Eigenvectors for symmetric 2×2, column-major: vecs[i + nl*k]
                if (std::abs(b) > Real(1e-16)) {
                    Real n0 = std::sqrt(b * b + (eig[0] - a) * (eig[0] - a));
                    vecs[0 * nl + 0] = (eig[0] - a) / n0;
                    vecs[1 * nl + 0] = b / n0; // col 0
                    Real n1 = std::sqrt(b * b + (eig[1] - a) * (eig[1] - a));
                    vecs[0 * nl + 1] = (eig[1] - a) / n1;
                    vecs[1 * nl + 1] = b / n1; // col 1
                } else {
                    vecs[0] = 1;
                    vecs[1] = 0;
                    vecs[2] = 0;
                    vecs[3] = 1;
                }
            }

            // Transform f_l(G): f'_k(G) = Σ_i V_{i,k} * f_i(G)
            // Store in fl_trans
            for (int k = 0; k < nl; ++k) {
                int ip_orig = l_groups[l][k];
                fl_trans[ip_out].resize(ws.n_g);
                for (int ig = 0; ig < ws.n_g; ++ig)
                    fl_trans[ip_out][ig] = Real(0);
                for (int i = 0; i < nl; ++i) {
                    Real v = vecs[i * nl + k];
                    for (int ig = 0; ig < ws.n_g; ++ig)
                        fl_trans[ip_out][ig] += v * fl_all[l_groups[l][i]][ig];
                }
                eff_D[ip_out] = eig[k] * Real(0.5); // Ry → Ha
                proj_l[ip_out] = l;
                ++ip_out;
            }
        }
        // Update number of effective projectors
        n_proj = ip_out;
        std::cout << "  D_ij diagonalized: " << n_proj << " eff proj (was "
                  << upf.n_proj() << ")" << std::endl;
    }

    auto ng_sz = static_cast<std::size_t>(ws.n_g);
    for (int ia = 0; ia < n_atoms; ++ia) {
        ws.kb_proj[ia].resize(n_proj);
        for (int ip = 0; ip < n_proj; ++ip) {
            ws.kb_proj[ia][ip].l = proj_l[ip];
            ws.kb_proj[ia][ip].d = eff_D[ip];
            ws.kb_proj[ia][ip].f_G.resize({ng_sz});
            for (int ig = 0; ig < ws.n_g; ++ig)
                ws.kb_proj[ia][ip].f_G[ig] = fl_trans[ip][ig];
        }
    }
}

/// Setup atomic wavefunctions in G-space for initial guess.
/// Follows QE's atwfc_mod.f90: chi_G(G) = (4π/√Ω) * ∫ χ(r) * j_l(Gr) * r * dr
void setup_atomic_wfc(PWWorkspace& ws, const internal::GVectorList& gv,
                      const internal::UPFData& upf,
                      const std::vector<Vec3>& pos_bohr, Real omega) {
    int n_atoms = static_cast<int>(pos_bohr.size());
    Real pref = Real(4.0 * M_PI) / std::sqrt(omega);

    std::cout << "  Atomic wfc: nwfc=" << upf.nwfc << " atoms=" << n_atoms
              << std::endl;

    // Precompute chi_G(G) for each atomic state (radial part, same for all atoms)
    // chi_G(G) = (4π/√Ω) * ∫ χ(r) * j_l(Gr) * r * dr  (QE atwfc_mod.f90:113-116)
    std::vector<std::vector<Real>> chi_g_all(upf.nwfc);
    for (int iw = 0; iw < upf.nwfc; ++iw) {
        int l = upf.lchi[iw];
        auto j_l = sph_bessel_for_l(l);
        const auto& r = upf.r;
        const auto& chi_r = upf.chi[iw];

        chi_g_all[iw].resize(ws.n_g, Real(0));
        for (int ig = 0; ig < ws.n_g; ++ig) {
            Real G = std::sqrt(gv.g[ig].norm_squared());
            Real sum = 0;
            for (size_t ir = 1; ir < r.size(); ++ir) {
                Real r_mid = Real(0.5) * (r[ir] + r[ir - 1]);
                Real dr = r[ir] - r[ir - 1];
                Real chi_mid = Real(0.5) * (chi_r[ir] + chi_r[ir - 1]);
                Real Gr = G * r_mid;
                sum += chi_mid * j_l(Gr) * r_mid * dr;
            }
            chi_g_all[iw][ig] = sum * pref;
        }
        std::cout << "    wfc[" << iw << "] l=" << l
                  << " chi_g[0]=" << chi_g_all[iw][0]
                  << " chi_g[1]=" << chi_g_all[iw][1] << std::endl;
    }

    // Store for each atom (radial part is the same, structure factor applied later)
    ws.atom_wfc.resize(n_atoms);
    auto ng_sz = static_cast<std::size_t>(ws.n_g);
    for (int ia = 0; ia < n_atoms; ++ia) {
        ws.atom_wfc[ia].resize(upf.nwfc);
        for (int iw = 0; iw < upf.nwfc; ++iw) {
            ws.atom_wfc[ia][iw].l = upf.lchi[iw];
            ws.atom_wfc[ia][iw].chi_G.resize({ng_sz});
            for (int ig = 0; ig < ws.n_g; ++ig)
                ws.atom_wfc[ia][iw].chi_G[ig] = chi_g_all[iw][ig];
        }
    }
}

} // namespace

PWResult pw_scf(const cpu::Frame& frame, const std::string& upf_path,
                double ecut_ry, int k_mesh, double tol, int max_iter) {
    PWResult result;

    int n_atoms = static_cast<int>(frame.atoms.size());
    ::ATV::Box box_ang = to_core_box(frame);
    ::ATV::Box box_bohr = box_to_bohr(box_ang);

    Real ecut_ht = static_cast<Real>(ecut_ry) * Real(0.5);
    Real vol_bohr = box_bohr.volume();

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
    int n_bands = std::max(n_elec / 2, 1);

    bool has_upf = !upf_path.empty();

    std::cout << "\n========== PW SCF ==========" << std::endl;
    std::cout << "Atoms:" << n_atoms << " Elec:" << n_elec
              << " Bands:" << n_bands << " Ecut:" << ecut_ry << " Ry"
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

    // QE Ewald G=0 term: ewaldg = -charge²/(4*alpha) * (2π/Ω)  (ewald.f90:112,129)
    Real ew_charge = 0;
    for (auto& z : ewald_charges)
        ew_charge += z;
    Real ew_alpha = ewald.alpha();
    Real ew_g0 = -ew_charge * ew_charge / (Real(4) * ew_alpha)
                 * Real(2.0 * M_PI) / vol_bohr;
    ew_total += ew_g0;

    // Alpha correction note: UPF stores V_loc in Ry, so V_loc(G=0) is in Ry.
    // The proper Ewald alpha correction requires careful handling of the
    // Coulomb singularity. Deferred for now — the raw Ewald (-15.52 Ry) is
    // close to QE's (-16.80 Ry). The 1.28 Ry difference will be addressed
    // after the SCF density and eigenvalues are fixed.

    std::cout << "\n--- Ewald (Ha) ---" << std::endl;
    std::cout << "  real=" << ewald.real_energy()
              << " recip=" << ewald.recip_energy()
              << " self=" << ewald.self_energy() << std::endl;
    std::cout << "  total=" << ew_total << " Ha (" << ew_total * 2.0 << " Ry)"
              << std::endl;

    // === G-vectors ===
    auto gv0 = internal::enumerate_g_vectors(box_bohr, ecut_ht, Vec3::zeros());
    std::cout << "\n--- G-vectors ---" << std::endl;
    std::cout << "  N_G=" << gv0.count()
              << " G_min=" << std::sqrt(gv0.g[1].norm_squared())
              << " G_max=" << std::sqrt(gv0.g[gv0.count() - 1].norm_squared())
              << std::endl;
    std::cout << "  FFT dims(hkl): " << gv0.nx << "x" << gv0.ny << "x" << gv0.nz
              << std::endl;

    // k-points
    internal::KPointGrid kg =
        internal::generate_mp_mesh(k_mesh, k_mesh, k_mesh);
    int n_kpts = static_cast<int>(kg.kpoints.size());

    // === Workspace ===
    PWWorkspace ws;
    ws.n_g = gv0.count();
    ws.gvec = &gv0;
    auto ng = static_cast<std::size_t>(ws.n_g);
    ws.rho_G.resize({ng});
    ws.rho_G_old.resize({ng});
    ws.V_H_G.resize({ng});
    ws.V_xc_G.resize({ng});
    ws.V_eff_G.resize({ng});
    ws.V_loc_G_real.resize({ng});
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
        // G-space transform of atomic pseudocharge: rho_atom(G) = 4π ∫ r² ρ(r) j₀(Gr) dr
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
                sum += r_mid * r_mid * rho_mid * j0 * dr;
            }
            Real rho_g_atom = Real(4.0 * M_PI) * sum / vol_bohr;
            // Structure factor: sum over atoms exp(-i G·R_a)
            Real s_real = 0;
            for (int ia = 0; ia < n_atoms; ++ia) {
                Real phase = gv0.g[ig][0] * pos_bohr[ia][0]
                             + gv0.g[ig][1] * pos_bohr[ia][1]
                             + gv0.g[ig][2] * pos_bohr[ia][2];
                s_real += std::cos(phase);
            }
            ws.rho_G[ig] = std::complex<Real>(rho_g_atom * s_real, Real(0));
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

    // Density FFT grid: QE default ecutrho = 4*ecutwfc
    ws.fft_nx = 4 * gv0.nx;
    ws.fft_ny = 4 * gv0.ny;
    ws.fft_nz = 4 * gv0.nz;
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
        setup_kb_projectors(ws, gv0, upf, pos_bohr, vol_bohr);
        if (upf.nwfc > 0)
            setup_atomic_wfc(ws, gv0, upf, pos_bohr, vol_bohr);

        // Setup augmentation form factor: rho_atom(G) / Q_atom
        ws.q_ij = upf.q_ij;
        ws.n_proj_total = upf.n_proj();
        ws.aug_form_factor_G.resize({static_cast<std::size_t>(ws.n_g)});
        if (!upf.rho_atom.empty() && upf.rho_atom.size() == upf.r.size()) {
            Real Q_atom = 0;
            for (size_t ir = 1; ir < upf.r.size(); ++ir) {
                Real r_mid = Real(0.5) * (upf.r[ir] + upf.r[ir - 1]);
                Real dr = upf.r[ir] - upf.r[ir - 1];
                Real rho_mid =
                    Real(0.5) * (upf.rho_atom[ir] + upf.rho_atom[ir - 1]);
                Q_atom += r_mid * r_mid * rho_mid * dr;
            }
            Q_atom *= Real(4.0 * M_PI);
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
                    sum += r_mid * r_mid * rho_mid * j0 * dr;
                }
                ws.aug_form_factor_G[ig] =
                    Real(4.0 * M_PI) * sum * inv_Q / vol_bohr;
            }
        } else {
            for (int ig = 0; ig < ws.n_g; ++ig)
                ws.aug_form_factor_G[ig] = Real(0);
        }
    } else {
        setup_toy_vloc(ws, gv0, pos_bohr, valence_charges, vol_bohr);
    }
    std::cout << "  V_loc(G) [1]=" << ws.V_loc_G_real[1]
              << " [10]=" << ws.V_loc_G_real[10] << std::endl;

    // === Context ===
    cpu::Context ctx;
    cpu::Frame mutable_frame = frame;
    cpu::BuildContext bctx{ctx, mutable_frame, &ws};

    // Pre-compute k-point G-vectors
    std::vector<internal::GVectorList> gv_k;
    gv_k.reserve(n_kpts);
    for (auto& k : kg.kpoints)
        gv_k.push_back(internal::enumerate_g_vectors(box_bohr, ecut_ht, k));

    // === SCF loop ===
    for (int iter = 0; iter < max_iter; ++iter) {
        PW_Poisson po;
        po.build(bctx);
        po.launch({ctx, mutable_frame, 0, iter});
        PW_LDA lda;
        lda.build(bctx);
        lda.launch({ctx, mutable_frame, 0, iter});
        PW_Veff v;
        v.build(bctx);
        v.launch({ctx, mutable_frame, 0, iter});

        // Accumulate density over k-points in real space
        auto* rho_r_kpt = static_cast<Real*>(
            fftw_malloc(sizeof(Real) * static_cast<size_t>(ws.nr)));
        for (int i = 0; i < ws.nr; ++i)
            rho_r_kpt[i] = Real(0);

        Real e_band = 0;
        for (int ik = 0; ik < n_kpts; ++ik) {
            int n_g_k = gv_k[ik].count();
            auto* saved_gvec = ws.gvec;
            ws.gvec = &gv_k[ik];
            ws.n_g = n_g_k;
            for (int ig = 0; ig < ws.n_g; ++ig)
                ws.T_G[ig] = (ig < n_g_k) ? gv_k[ik].ekin[ig] : Real(1e10);

            PW_Hamiltonian h;
            h.build(bctx);
            h.launch({ctx, mutable_frame, 1, iter});
            PW_Davidson d{n_bands, Real(tol) * Real(0.1), 40};
            d.build(bctx);
            d.launch({ctx, mutable_frame, 1, iter});

            Real weight = kg.weights[ik];
            for (int b = 0; b < n_bands; ++b)
                e_band += weight * ws.eigenvalues[b] * Real(2);

            // Accumulate density from this k-point in real space
            if (ws.nr > 0 && ws.fft) {
                int nb =
                    std::min(n_bands, static_cast<int>(ws.eigenvalues.size()));
                int ng_k = ws.n_g;
                int nx = ws.fft_nx, ny = ws.fft_ny, nz = ws.fft_nz;
                auto& gvk = *ws.gvec;
                auto* psi_cc = static_cast<fftw_complex*>(fftw_malloc(
                    sizeof(fftw_complex) * static_cast<size_t>(ws.nr)));
                for (int b = 0; b < nb; ++b) {
                    for (int i = 0; i < ws.nr; ++i)
                        psi_cc[i][0] = psi_cc[i][1] = 0;
                    for (int ig = 0; ig < ng_k; ++ig) {
                        int ix = ((gvk.h[ig] % nx) + nx) % nx;
                        int iy = ((gvk.k[ig] % ny) + ny) % ny;
                        int iz = ((gvk.l[ig] % nz) + nz) % nz;
                        int idx = ix + nx * (iy + ny * iz);
                        if (idx >= 0 && idx < ws.nr) {
                            auto val = ws.wfc(b, ig);
                            psi_cc[idx][0] = static_cast<double>(val.real());
                            psi_cc[idx][1] = static_cast<double>(val.imag());
                        }
                    }
                    ws.fft->c2c_bwd(psi_cc);
                    for (int i = 0; i < ws.nr; ++i) {
                        Real psi_sq =
                            static_cast<Real>(psi_cc[i][0] * psi_cc[i][0]
                                              + psi_cc[i][1] * psi_cc[i][1]);
                        // ψ_phys(r) = ψ_fft(r)/√Ω, so |ψ|²_phys = |ψ_fft|²/Ω
                        rho_r_kpt[i] += weight * Real(2) * psi_sq / ws.omega;
                    }
                }
                fftw_free(psi_cc);
            }

            ws.gvec = saved_gvec;
            ws.n_g = gv0.count();
        }

        // FFT total density to G-space; slow mixing preserves atomic density
        if (ws.nr > 0 && ws.fft) {
            int nc = ws.nc;
            auto* rho_g_full = static_cast<fftw_complex*>(
                fftw_malloc(sizeof(fftw_complex) * static_cast<size_t>(nc)));
            ws.fft->r2c(rho_r_kpt, rho_g_full);

            for (int i = 0; i < ws.n_g; ++i)
                ws.rho_G_old[i] = ws.rho_G[i];

            int nx = ws.fft_nx, ny = ws.fft_ny, nz = ws.fft_nz,
                nxc = nx / 2 + 1;
            Real inv_nr = Real(1) / Real(ws.nr);
            for (int ig = 0; ig < ws.n_g; ++ig) {
                int h = gv0.h[ig], k = gv0.k[ig], l = gv0.l[ig];
                int kx = ((h % nx) + nx) % nx, ky = ((k % ny) + ny) % ny,
                    kz = ((l % nz) + nz) % nz;
                bool conj = false;
                if (kx > nx / 2) {
                    kx = nx - kx;
                    ky = (ny - ky) % ny;
                    kz = (nz - kz) % nz;
                    conj = true;
                }
                int idx = kx + nxc * (ky + ny * kz);
                if (idx >= 0 && idx < nc) {
                    ws.rho_G[ig] = std::complex<Real>(
                        static_cast<Real>(rho_g_full[idx][0]) * inv_nr,
                        static_cast<Real>(conj ? -rho_g_full[idx][1]
                                               : rho_g_full[idx][1])
                            * inv_nr);
                } else
                    ws.rho_G[ig] = std::complex<Real>(0, 0);
            }
            fftw_free(rho_g_full);

            // Pulay/DIIS mixing (QE's Anderson/Broyden simplified)
            // Store residual F = ρ_out - ρ_in
            xt::xarray<std::complex<Real>> F = ws.rho_G - ws.rho_G_old;
            ws.pulay_res.push_back(F);
            ws.pulay_rho.push_back(ws.rho_G_old);

            int depth =
                std::min(ws.pulay_depth, static_cast<int>(ws.pulay_res.size()));
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
                        if (std::abs(aug[i * n + k])
                            > std::abs(aug[imax * n + k]))
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
                        ws.rho_G[ig] += c[i]
                                        * (ws.pulay_rho[i][ig]
                                           + alpha * ws.pulay_res[i][ig]);
                }
            } else {
                // Simple linear mixing for first iteration
                for (int i = 0; i < ws.n_g; ++i)
                    ws.rho_G[i] = ws.rho_G_old[i] + ws.mix_beta * F[i];
            }
        }
        fftw_free(rho_r_kpt);

        PW_Convergence conv{Real(tol), 2};
        conv.build(bctx);
        conv.launch({ctx, mutable_frame, 3, iter});

        result.iterations = iter + 1;

        // Energy decomposition
        Real e_h = 0, e_xc_int = 0;
        for (int i = 0; i < ws.n_g; ++i) {
            Real rh = ws.rho_G[i].real();
            e_h += Real(0.5) * rh * ws.V_H_G[i].real();
            e_xc_int += rh * ws.V_xc_G[i].real();
        }
        e_h *= vol_bohr;
        e_xc_int *= vol_bohr;

        // one-electron = Σ ε_i - ∫ρ V_H - ∫ρ V_xc
        Real e_one = e_band - Real(2.0) * e_h - e_xc_int;
        // total = one-electron + Hartree + XC + Ewald + augmentation
        Real e_total_ha = e_one + e_h + ws.E_xc + ew_total + ws.E_aug;
        result.energy = static_cast<double>(e_total_ha) * 2.0;

        std::cout << "\n-- iter " << (iter + 1) << " --" << std::endl;
        std::cout << "  one-electron = " << std::setw(12)
                  << std::setprecision(6) << e_one << " Ha (" << e_one * 2.0
                  << " Ry)" << std::endl;
        std::cout << "  Hartree      = " << std::setw(12) << e_h << " Ha ("
                  << e_h * 2.0 << " Ry)" << std::endl;
        std::cout << "  XC           = " << std::setw(12) << ws.E_xc << " Ha ("
                  << ws.E_xc * 2.0 << " Ry)" << std::endl;
        std::cout << "  Ewald        = " << std::setw(12) << ew_total << " Ha ("
                  << ew_total * 2.0 << " Ry)" << std::endl;
        std::cout << "  Aug          = " << std::setw(12) << ws.E_aug << " Ha ("
                  << ws.E_aug * 2.0 << " Ry)" << std::endl;
        std::cout << "  Total        = " << std::setw(12) << e_total_ha
                  << " Ha (" << result.energy << " Ry)" << std::endl;
        std::cout << "  E_band=" << e_band << " Ha  eig:";
        for (int b = 0; b < std::min(4, n_bands); ++b)
            std::cout << " " << ws.eigenvalues[b];
        std::cout << ")" << std::endl;

        if (ws.converged) {
            result.converged = true;
            std::cout << "\n>>> Converged in " << (iter + 1) << " iters"
                      << std::endl;
            break;
        }
    }

    if (!result.converged)
        std::cout << "\n>>> NOT converged after " << max_iter << " iters"
                  << std::endl;

    return result;
}

} // namespace pw
} // namespace cpu
} // namespace ATV
