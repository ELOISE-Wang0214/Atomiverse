/// @file cpu/pw/pw_solver.cpp
/// @brief PWSolver implementation — orchestrates the internal step helpers
///        through the plane-wave SCF self-consistency cycle.
///
/// `pw_scf()` is a thin functional wrapper that constructs a transient
/// `PWSolver` and runs it; both produce identical results.
///
/// MPI: when built with `ATV_HAS_MPI` and MPI is initialized, k-points
/// are distributed round-robin across `MPI_COMM_WORLD`; the band-energy
/// scalar and the real-space density accumulation are Allreduced each
/// SCF iteration, so all ranks hold the identical mixed density and
/// take identical convergence decisions.

#include "atomiverse/cpu/pw/pw_solver.h"

#include "atomiverse/core/box.h"
#include "atomiverse/core/constants.h"
#include "atomiverse/core/log.h"
#include "atomiverse/core/mat3.h"
#include "atomiverse/core/units.h"
#include "atomiverse/cpu/pw/internal/convergence.h"
#include "atomiverse/cpu/pw/internal/davidson.h"
#include "atomiverse/cpu/pw/internal/density_builder.h"
#include "atomiverse/cpu/pw/internal/ewald.h"
#include "atomiverse/cpu/pw/internal/g_vector.h"
#include "atomiverse/cpu/pw/internal/kpoint.h"
#include "atomiverse/cpu/pw/internal/poisson.h"
#include "atomiverse/cpu/pw/internal/pw_workspace.h"
#include "atomiverse/cpu/pw/internal/upf_reader.h"
#include "atomiverse/cpu/pw/internal/veff_builder.h"
#include "atomiverse/cpu/pw/internal/xc_builder.h"

#ifdef ATV_HAS_MPI
#include <mpi.h>
#endif

#include <algorithm>
#include <cmath>
#include <complex>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace ATV {
namespace cpu {
namespace pw {

namespace {

/// @brief Valence electron count for the toy (UPF-less) pseudopotential.
/// @throws std::invalid_argument for elements outside H..Ar — no silent
///         fallback (a 0-electron system would be silently wrong).
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
        throw std::invalid_argument(
            "pw: no valence table entry for Z=" + std::to_string(z)
            + " (toy V_loc supports H..Ar; provide a UPF pseudopotential)");
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

void setup_upf_vloc(internal::PWWorkspace& ws,
                    const internal::GVectorList& gv,
                    const internal::UPFData& upf, Real omega) {
    auto vg = upf_vloc_to_gspace(upf, gv, omega);
    for (int i = 0; i < ws.n_g; ++i)
        ws.V_loc_G_real[i] = vg[i];
    ws.V_loc_G_real[0] = Real(0); // G=0 cancelled against Ewald
}

void setup_toy_vloc(internal::PWWorkspace& ws,
                    const internal::GVectorList& gv,
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
    ::ATV::Box box = box_ang;
    box.cell = box.cell * Real(kAngstrom2Bohr);
    return box;
}

std::vector<Vec3> pos_to_bohr(const cpu::Frame& frame) {
    Real a2b = Real(kAngstrom2Bohr);
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
/// @throws std::invalid_argument for l > 2 — no silent j0 fallback.
SphBesselFn sph_bessel_for_l(int l) {
    switch (l) {
    case 0:
        return sph_bessel_j0;
    case 1:
        return sph_bessel_j1;
    case 2:
        return sph_bessel_j2;
    default:
        throw std::invalid_argument(
            "pw: spherical Bessel j_l implemented for l<=2, got l="
            + std::to_string(l));
    }
}

/// Setup Kleinman-Bylander non-local projectors in G-space.
void setup_kb_projectors(internal::PWWorkspace& ws,
                         const internal::GVectorList& gv,
                         const internal::UPFData& upf,
                         const std::vector<Vec3>& pos_bohr, Real omega) {
    int n_atoms = static_cast<int>(pos_bohr.size());
    ws.kb_proj.resize(n_atoms);

    ATV_DEBUG("[pw] KB projectors: {} proj, {} atoms", upf.kb.size(), n_atoms);

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
                (void)ip_orig;
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
        ATV_DEBUG("[pw] D_ij diagonalized: {} eff proj (was {})", n_proj,
                  upf.n_proj());
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
void setup_atomic_wfc(internal::PWWorkspace& ws,
                      const internal::GVectorList& gv,
                      const internal::UPFData& upf,
                      const std::vector<Vec3>& pos_bohr, Real omega) {
    int n_atoms = static_cast<int>(pos_bohr.size());
    Real pref = Real(4.0 * M_PI) / std::sqrt(omega);

    ATV_DEBUG("[pw] atomic wfc: nwfc={} atoms={}", upf.nwfc, n_atoms);

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

/// @brief MPI layout for the k-point round-robin distribution.
struct MpiLayout {
    int rank = 0;
    int nprocs = 1;
    bool active = false; ///< true only when MPI is built AND initialized.
};

MpiLayout mpi_layout() {
    MpiLayout layout;
#ifdef ATV_HAS_MPI
    int initialized = 0;
    MPI_Initialized(&initialized);
    if (initialized) {
        MPI_Comm_rank(MPI_COMM_WORLD, &layout.rank);
        MPI_Comm_size(MPI_COMM_WORLD, &layout.nprocs);
        layout.active = (layout.nprocs > 1);
    }
#endif
    return layout;
}

#ifdef ATV_HAS_MPI
MPI_Datatype mpi_real_type() {
    return sizeof(Real) == sizeof(double) ? MPI_DOUBLE : MPI_FLOAT;
}
#endif

} // namespace

// ── PWSolver ────────────────────────────────────────────────────────────

PWSolver::PWSolver(PWConfig config)
    : config_(std::move(config)),
      upf_(std::make_unique<internal::UPFData>()),
      ws_(std::make_unique<internal::PWWorkspace>()),
      poisson_(std::make_unique<internal::PoissonSolver>()),
      xc_(std::make_unique<internal::XCBuilder>()),
      veff_(std::make_unique<internal::VeffBuilder>()),
      density_(std::make_unique<internal::DensityBuilder>(
          Real(config_.mix_beta), config_.pulay_depth)) {}

// Out-of-line where internal::* are complete types (required for the
// unique_ptr members over forward-declared helpers).
PWSolver::~PWSolver() = default;
PWSolver::PWSolver(PWSolver&&) noexcept = default;
PWSolver& PWSolver::operator=(PWSolver&&) noexcept = default;

void PWSolver::build(const cpu::Frame& reference) {
    if (!config_.upf_path.empty()) {
        *upf_ = internal::read_upf_v1(config_.upf_path);
    } else {
        // Validate the toy valence table covers every atom type up front
        // (throws on unsupported Z — no silent 0-electron fallback).
        for (size_t i = 0; i < reference.atoms.size(); ++i)
            (void)valence_from_z(reference.atoms.type[i]);
    }
    built_ = true;
}

PWResult& PWSolver::run(const cpu::Frame& frame) {
    // Lazy build — convenience for callers that skip explicit build().
    if (!built_)
        build(frame);

    result_ = PWResult{};
    auto& ws = *ws_;

    const auto layout = mpi_layout();
    const bool log_rank0 = (layout.rank == 0);

    int n_atoms = static_cast<int>(frame.atoms.size());
    ::ATV::Box box_ang = to_core_box(frame);
    ::ATV::Box box_bohr = box_to_bohr(box_ang);

    Real ecut_ht = static_cast<Real>(config_.ecut_ry) * Real(kRy2Hartree);
    Real vol_bohr = box_bohr.volume();

    auto pos_bohr = pos_to_bohr(frame);

    const bool has_upf = !config_.upf_path.empty();

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

    if (log_rank0)
        ATV_INFO("[pw] start atoms={} elec={} bands={} ecut_ry={} {} "
                 "k_mesh={} ranks={}",
                 n_atoms, n_elec, n_bands, config_.ecut_ry,
                 has_upf ? "upf" : "toy", config_.k_mesh, layout.nprocs);

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

    if (log_rank0)
        ATV_DEBUG("[pw] ewald_ha real={} recip={} self={} total={}",
                  ewald.real_energy(), ewald.recip_energy(),
                  ewald.self_energy(), ew_total);

    // === G-vectors (density sphere at Γ) ===
    ws.gvec0 = internal::enumerate_g_vectors(box_bohr, ecut_ht, Vec3::zeros());
    auto& gv0 = ws.gvec0;
    if (log_rank0)
        ATV_DEBUG("[pw] n_g={} fft_hkl={}x{}x{}", gv0.count(), gv0.nx, gv0.ny,
                  gv0.nz);

    // k-points: full mesh on every rank, round-robin ownership.
    ws.kgrid =
        internal::generate_mp_mesh(config_.k_mesh, config_.k_mesh,
                                   config_.k_mesh);
    int n_kpts = static_cast<int>(ws.kgrid.kpoints.size());
    ws.owned_k.clear();
    ws.gvec_k.clear();
    for (int ik = layout.rank; ik < n_kpts; ik += layout.nprocs) {
        ws.owned_k.push_back(ik);
        ws.gvec_k.push_back(internal::enumerate_g_vectors(
            box_bohr, ecut_ht, ws.kgrid.kpoints[ik]));
    }

    // === Workspace sizing ===
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
    ws.mix_beta = Real(config_.mix_beta);
    ws.pulay_depth = config_.pulay_depth;
    ws.omega = vol_bohr; // cell volume for real-space integration
    ws.atom_pos_bohr = pos_bohr;
    ws.E_aug = Real(0);
    ws.E_xc = Real(0);

    const auto& upf = *upf_;

    // Initial density
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
        // Uniform initial density (no UPF atomic pseudocharge available)
        ws.rho_G[0] = std::complex<Real>(Real(n_elec) / vol_bohr, 0);
    }

    // Density FFT grid: QE default ecutrho = 4*ecutwfc
    int fft_nx = 4 * gv0.nx;
    int fft_ny = 4 * gv0.ny;
    int fft_nz = 4 * gv0.nz;
    if (fft_nx % 2)
        ++fft_nx;
    if (fft_ny % 2)
        ++fft_ny;
    if (fft_nz % 2)
        ++fft_nz;
    const bool grid_changed = (fft_nx != ws.fft_nx || fft_ny != ws.fft_ny
                               || fft_nz != ws.fft_nz || !ws.fft);
    ws.fft_nx = fft_nx;
    ws.fft_ny = fft_ny;
    ws.fft_nz = fft_nz;
    ws.nr = ws.fft_nx * ws.fft_ny * ws.fft_nz;
    ws.nc = (ws.fft_nx / 2 + 1) * ws.fft_ny * ws.fft_nz;
    ws.rho_r.resize({static_cast<std::size_t>(ws.nr)});
    ws.V_xc_r.resize({static_cast<std::size_t>(ws.nr)});
    if (grid_changed)
        ws.fft = std::make_unique<internal::FFTPlan3D>(ws.fft_nx, ws.fft_ny,
                                                       ws.fft_nz);

    if (log_rank0)
        ATV_DEBUG("[pw] density fft grid {}x{}x{} (nr={})", ws.fft_nx,
                  ws.fft_ny, ws.fft_nz, ws.nr);

    // === V_loc ===
    if (has_upf) {
        if (log_rank0)
            ATV_DEBUG("[pw] upf z_val={} l_max={} n_proj={}", upf.z_val,
                      upf.l_max, upf.kb.size());
        setup_upf_vloc(ws, gv0, upf, vol_bohr);
        setup_kb_projectors(ws, gv0, upf, pos_bohr, vol_bohr);
        if (upf.nwfc > 0)
            setup_atomic_wfc(ws, gv0, upf, pos_bohr, vol_bohr);

        // Augmentation form factor: rho_atom(G) / Q_atom. Parsed and
        // stored as fundamental pseudopotential data; the augmentation
        // density itself is not yet applied by this solver (E_aug = 0).
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

    // === Step-helper preparation (idempotent for unchanged grids) ===
    xc_->prepare(ws);
    veff_->prepare(ws);
    density_->prepare(ws);
    density_->reset(ws);
    if (!davidson_ || davidson_bands_ != n_bands) {
        davidson_ = std::make_unique<internal::DavidsonSolver>(
            n_bands, Real(config_.tol) * Real(0.1), 40);
        davidson_bands_ = n_bands;
    }
    davidson_->prepare(ws);
    internal::ConvergenceCheck convergence{Real(config_.tol), 2};

    // === SCF loop ===
    bool converged = false;
    for (int iter = 0; iter < config_.max_iter; ++iter) {
        poisson_->step(ws);
        xc_->step(ws);
        veff_->step(ws);

        // Accumulate density and band energy over this rank's k-points.
        density_->begin(ws);
        Real e_band = 0;
        for (size_t i_owned = 0; i_owned < ws.owned_k.size(); ++i_owned) {
            int ik = ws.owned_k[i_owned];
            auto& gvk = ws.gvec_k[i_owned];
            int n_g_k = gvk.count();

            // Swap the workspace onto this k-point's G-sphere.
            ws.gvec = &gvk;
            ws.n_g = n_g_k;
            ws.T_G.resize({static_cast<std::size_t>(n_g_k)});
            for (int ig = 0; ig < n_g_k; ++ig)
                ws.T_G[ig] = gvk.ekin[ig];

            davidson_->step(ws);

            Real weight = ws.kgrid.weights[ik];
            for (int b = 0; b < n_bands; ++b)
                e_band += weight * ws.eigenvalues[b] * Real(2);

            density_->accumulate_kpoint(ws, weight, n_bands);

            // Restore the density G-sphere.
            ws.gvec = &gv0;
            ws.n_g = gv0.count();
        }

#ifdef ATV_HAS_MPI
        if (layout.active) {
            MPI_Allreduce(MPI_IN_PLACE, density_->accumulation(), ws.nr,
                          mpi_real_type(), MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &e_band, 1, mpi_real_type(), MPI_SUM,
                          MPI_COMM_WORLD);
        }
#endif

        // FFT the (globally summed) density to G-space + Pulay mixing.
        density_->finish(ws);

        converged = convergence.step(ws);
        result_.iterations = iter + 1;

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
        result_.energy_ha = static_cast<double>(e_total_ha);
        result_.energy = result_.energy_ha / kRy2Hartree;

        if (log_rank0)
            ATV_DEBUG("[pw] iter={} one_e={:.6f} hartree={:.6f} xc={:.6f} "
                      "ewald={:.6f} total_ha={:.8f} ({:.8f} Ry)",
                      iter + 1, e_one, e_h, ws.E_xc, ew_total, e_total_ha,
                      result_.energy);

        if (converged) {
            result_.converged = true;
            if (log_rank0)
                ATV_INFO("[pw] converged iter={} energy_ry={:.8f}", iter + 1,
                         result_.energy);
            break;
        }
    }

    if (!result_.converged && log_rank0)
        ATV_WARN("[pw] not converged after {} iters energy_ry={:.8f}",
                 config_.max_iter, result_.energy);

    return result_;
}

} // namespace pw
} // namespace cpu
} // namespace ATV
