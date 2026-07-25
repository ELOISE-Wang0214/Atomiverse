/// @file cpu/gpw/internal/gpw_stress.cpp
/// @brief GPW cell stress: Hellmann–Feynman fixed-P product path + multi-SCF FD.
///
/// Analytic product path (gpw_stress / gpw_stress_analytic_impl):
///   1. One full SCF at the reference geometry → density matrix P.
///   2. Central FD of the *non-self-consistent* total energy E[P; ε] under
///      orthorhombic diagonal strain at fixed fractional coordinates. Geometry-
///      dependent T, V_pp, grid, Poisson, XC (and optional D3) are rebuilt;
///      the SCF is *not* re-converged. This is the GPW form of the Hellmann–
///      Feynman cell theorem used for plane-wave Nielsen–Martin in this tree
///      (Nielsen & Martin, PRB 32, 3780 (1985); Quickstep energy surface).
///
/// Oracle path (gpw_stress_fd / gpw_stress_fd_impl): full multi-SCF strain-FD
/// of E_tot (2 SCFs per diagonal axis).
///
/// Off-diagonal shear remains 0 (orthorhombic GPW cell only). Sign matches
/// CP2K STRESS_TENSOR: σ_aa = −(1/Ω) ∂E/∂ε_aa.

#include "atomiverse/cpu/gpw/internal/gpw_stress.h"

#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "atomiverse/cpu/gpw/internal/collocate.h"
#include "atomiverse/cpu/gpw/internal/dftd3.h"
#include "atomiverse/cpu/gpw/internal/gth_integrals.h"
#include "atomiverse/cpu/gpw/internal/local_pp.h"
#include "atomiverse/cpu/gpw/internal/mt_poisson.h"
#include "atomiverse/cpu/gpw/internal/one_electron.h"
#include "atomiverse/cpu/gpw/internal/periodic.h"
#include "atomiverse/cpu/gpw/internal/poisson_fft.h"
#include "atomiverse/cpu/gpw/internal/xc_grid.h"
#include "atomiverse/cpu/gpw/internal/xc_grid_gga.h"
#include "atomiverse/cpu/linalg/matrix.h"

namespace ATV {
namespace cpu {
namespace gpw {
namespace internal {
namespace {

using linalg::RealMatrix;

GpwDft copy_gpw(const GpwDft& m) {
    GpwDft p;
    for (int d = 0; d < 3; ++d)
        p.cell[d] = m.cell[d];
    p.cutoff_ry = m.cutoff_ry;
    p.ngrids = m.ngrids;
    p.rel_cutoff = m.rel_cutoff;
    p.n_occ = m.n_occ;
    p.n_alpha = m.n_alpha;
    p.n_beta = m.n_beta;
    p.max_iter = m.max_iter;
    p.tol = m.tol;
    p.eps_default = m.eps_default;
    p.mix = m.mix;
    p.diis_depth = m.diis_depth;
    p.lindep_threshold = m.lindep_threshold;
    p.level_shift = m.level_shift;
    p.xc = m.xc;
    p.nonperiodic = m.nonperiodic;
    p.local_pp_on_grid = m.local_pp_on_grid;
    p.d3 = m.d3;
    p.d3bj = m.d3bj;
    p.d4 = m.d4;
    p.use_ot = m.use_ot;
    p.plus_u = m.plus_u;
    p.atom_spin_guess = m.atom_spin_guess;
    for (int d = 0; d < 3; ++d)
        p.k_mesh[d] = m.k_mesh[d];
    p.k_time_reversal = m.k_time_reversal;
    p.k_points = m.k_points;
    p.k_point_weights = m.k_point_weights;
    return p;
}

/// Scale atom positions and all basis-function centers along axis `ax` by
/// factor `s` (fixed fractional coordinates under orthorhombic strain).
std::vector<GpwAtom> scale_axis(const std::vector<GpwAtom>& atoms, int ax,
                                double s) {
    std::vector<GpwAtom> out = atoms;
    for (auto& a : out) {
        a.pos[ax] *= s;
        for (auto& bf : a.basis)
            for (auto& p : bf)
                p.center[ax] *= s;
    }
    return out;
}

double scf_energy_ha(const std::vector<GpwAtom>& atoms, const GpwDft& method,
                     const linalg::RealMatrix* guess_C) {
    GpwDft solver = copy_gpw(method);
    solver.build(atoms);
    const GpwResult& r = solver.run(atoms, guess_C);
    if (!r.converged)
        throw std::runtime_error("gpw_stress_fd: SCF did not converge");
    return r.e_total;
}

Prim to_prim(const GridPrimitive& p) {
    Prim q;
    q.alpha = p.alpha;
    q.center[0] = p.center[0];
    q.center[1] = p.center[1];
    q.center[2] = p.center[2];
    q.l[0] = p.l[0];
    q.l[1] = p.l[1];
    q.l[2] = p.l[2];
    return q;
}

std::vector<GridBasisFunction> global_basis(const std::vector<GpwAtom>& atoms) {
    std::vector<GridBasisFunction> b;
    for (const auto& a : atoms)
        for (const auto& bf : a.basis)
            b.push_back(bf);
    return b;
}

// Γ-point periodic <a | Σ_C v_C | b> — same convention as gpw_forces / SCF.
double vpp_images(const GridBasisFunction& a, const GridBasisFunction& b,
                  const std::vector<const GpwAtom*>& centers, const double L[3],
                  const int m[3], double reach, bool with_local) {
    if (a.empty() || b.empty())
        return 0.0;
    const double* ac = a[0].center;
    const double* bc = b[0].center;
    const double reach2 = reach * reach;
    const bool per = m[0] || m[1] || m[2];
    double v = 0.0;
    for (const auto& p : a) {
        const Prim pp = to_prim(p);
        for (const auto& q0 : b) {
            const double cc = p.coeff * q0.coeff;
            if (cc == 0.0)
                continue;
            const double aa = p.alpha;
            for (int tx = -m[0]; tx <= m[0]; ++tx)
                for (int ty = -m[1]; ty <= m[1]; ++ty)
                    for (int tz = -m[2]; tz <= m[2]; ++tz) {
                        const double sx = ac[0] - bc[0] - tx * L[0];
                        const double sy = ac[1] - bc[1] - ty * L[1];
                        const double sz = ac[2] - bc[2] - tz * L[2];
                        if (sx * sx + sy * sy + sz * sz > reach2)
                            continue;
                        GridPrimitive q = q0;
                        q.center[0] += tx * L[0];
                        q.center[1] += ty * L[1];
                        q.center[2] += tz * L[2];
                        const double bb = q.alpha, g = aa + bb;
                        double P[3];
                        for (int d = 0; d < 3; ++d)
                            P[d] = (aa * p.center[d] + bb * q.center[d]) / g;
                        const Prim qq = to_prim(q);
                        for (const GpwAtom* Cp : centers) {
                            const GpwAtom& C = *Cp;
                            double apos[3];
                            for (int d = 0; d < 3; ++d) {
                                apos[d] = C.pos[d];
                                if (per && L[d] > 0.0)
                                    apos[d] -= L[d]
                                               * std::nearbyint(
                                                   (C.pos[d] - P[d]) / L[d]);
                            }
                            if (with_local && !C.cloc.empty())
                                v += cc
                                     * gth_local_sr(pp, qq, apos, C.rloc,
                                                    C.cloc);
                            for (const auto& ch : C.nonlocal)
                                v += cc * gth_nonlocal(pp, qq, apos, ch);
                        }
                    }
        }
    }
    return v;
}

double core_self_energy(const GpwAtom& a) {
    return a.z_val * a.z_val / (2.0 * a.rloc * std::sqrt(M_PI));
}

/// @brief Hellmann–Feynman Lagrangian energy at fixed (P, W) and geometry.
///
///   E_L = Tr(P T) + Tr(P V_pp) + E_es[ρ_P] + E_xc[ρ_P] + E_disp − Tr(W S)
///
/// At the SCF stationary point the total-energy derivative w.r.t. a geometry
/// parameter equals ∂E_L/∂λ at fixed P and energy-weighted density W
/// (Pulay −W·dS). Dropping −Tr(W S) leaves an incomplete fixed-P surface whose
/// strain FD does not match multi-SCF FD.
double energy_fixed_P(const std::vector<GpwAtom>& atoms, const GpwDft& cfg,
                      const RealMatrix& P, const RealMatrix& W) {
    if (cfg.is_uks())
        throw std::invalid_argument(
            "gpw_stress analytic: UKS not supported in this foothold");
    if (cfg.k_mesh[0] != 1 || cfg.k_mesh[1] != 1 || cfg.k_mesh[2] != 1
        || !cfg.k_points.empty())
        throw std::invalid_argument(
            "gpw_stress analytic: k-point GPW not supported in this foothold");
    if (!cfg.plus_u.empty())
        throw std::invalid_argument(
            "gpw_stress analytic: DFT+U not supported in this foothold");
    if (cfg.d4)
        throw std::invalid_argument(
            "gpw_stress analytic: D4 not supported in this foothold");
    const int n_disp = (cfg.d3 ? 1 : 0) + (cfg.d3bj ? 1 : 0) + (cfg.d4 ? 1 : 0);
    if (n_disp > 1)
        throw std::invalid_argument(
            "gpw_stress analytic: d3, d3bj, and d4 are mutually exclusive");

    const std::vector<GridBasisFunction> basis = global_basis(atoms);
    const int nb = static_cast<int>(basis.size());
    if (static_cast<int>(P.shape(0)) != nb || static_cast<int>(P.shape(1)) != nb
        || static_cast<int>(W.shape(0)) != nb
        || static_cast<int>(W.shape(1)) != nb)
        throw std::invalid_argument(
            "gpw_stress analytic: P/W shape does not match basis");

    const bool periodic = !cfg.nonperiodic;
    const double eps = cfg.eps_default;
    const bool analytic_loc = !cfg.local_pp_on_grid;

    std::vector<double> radii(static_cast<std::size_t>(nb));
    double rmax = 0.0;
    for (int u = 0; u < nb; ++u) {
        radii[static_cast<std::size_t>(u)] =
            basis_function_radius(basis[u], eps);
        rmax = std::max(rmax, radii[static_cast<std::size_t>(u)]);
    }
    int mimg[3];
    pair_image_counts(cfg.cell, periodic, 2.0 * rmax, mimg);

    std::vector<const GpwAtom*> all(atoms.size());
    for (std::size_t a = 0; a < atoms.size(); ++a)
        all[a] = &atoms[a];

    // Analytic one-electron: Tr(P T) + Tr(P V_pp) − Tr(W S) (Pulay).
    double e_kin = 0.0;
    double e_pp = 0.0;
    double e_pulay = 0.0;
    for (int mu = 0; mu < nb; ++mu) {
        if (basis[static_cast<std::size_t>(mu)].empty())
            continue;
        for (int nu = 0; nu < nb; ++nu) {
            if (basis[static_cast<std::size_t>(nu)].empty())
                continue;
            const double reach = radii[static_cast<std::size_t>(mu)]
                                 + radii[static_cast<std::size_t>(nu)];
            const double s =
                overlap(basis[mu], basis[nu], cfg.cell, mimg, reach);
            const double t =
                kinetic(basis[mu], basis[nu], cfg.cell, mimg, reach);
            const double vpp = vpp_images(basis[mu], basis[nu], all, cfg.cell,
                                          mimg, reach, analytic_loc);
            e_kin += P(mu, nu) * t;
            e_pp += P(mu, nu) * vpp;
            e_pulay += W(mu, nu) * s;
        }
    }

    UniformGrid grid = make_uniform_grid(cfg.cell, cfg.cutoff_ry);
    const long ng = grid.npts();
    const double wgt = grid.weight();
    const CompactBasisGrid cbg =
        build_compact_basis_grid(grid, basis, periodic, eps);

    std::vector<double> Pvec(static_cast<std::size_t>(nb) * nb);
    for (int i = 0; i < nb; ++i)
        for (int j = 0; j < nb; ++j)
            Pvec[static_cast<std::size_t>(i) * nb + j] = P(i, j);
    std::vector<double> rho = collocate_density(cbg, nb, Pvec);

    // Ionic core density n_c (min-image; matches SCF / forces).
    std::vector<double> n_ion(static_cast<std::size_t>(ng), 0.0);
    double e_self = 0.0;
    for (const auto& a : atoms)
        e_self += core_self_energy(a);
    for (long g = 0; g < ng; ++g) {
        double r[3];
        grid.point(g, r);
        double s = 0.0;
        for (const auto& a : atoms) {
            const double s2 = a.rloc * a.rloc;
            const double norm = a.z_val * std::pow(2.0 * M_PI * s2, -1.5);
            double dx = r[0] - a.pos[0], dy = r[1] - a.pos[1],
                   dz = r[2] - a.pos[2];
            if (periodic) {
                dx -= grid.L[0] * std::nearbyint(dx / grid.L[0]);
                dy -= grid.L[1] * std::nearbyint(dy / grid.L[1]);
                dz -= grid.L[2] * std::nearbyint(dz / grid.L[2]);
            }
            s += norm * std::exp(-(dx * dx + dy * dy + dz * dz) / (2.0 * s2));
        }
        n_ion[static_cast<std::size_t>(g)] = s;
    }

    std::vector<double> n_tot(rho.size());
    for (std::size_t g = 0; g < rho.size(); ++g)
        n_tot[g] = rho[g] - n_ion[g];

    const PoissonResult es = cfg.nonperiodic ? solve_poisson_mt(grid, n_tot)
                                             : solve_poisson(grid, n_tot);
    const XcOnGrid xc = cfg.xc == XcKind::PBE ? xc_grid_pbe(grid, rho)
                                              : (cfg.xc == XcKind::PBEsol
                                                     ? xc_grid_pbesol(grid, rho)
                                                     : xc_grid_pade(grid, rho));

    double e_locgrid = 0.0;
    if (!analytic_loc) {
        const std::vector<double> v_locsr =
            local_pp_grid(grid, atoms, periodic);
        for (std::size_t g = 0; g < rho.size(); ++g)
            e_locgrid += rho[g] * v_locsr[g];
        e_locgrid *= wgt;
    }

    const double e_es = es.e_hartree - e_self;
    const double e_xc = xc.e_xc;

    double e_disp = 0.0;
    if (cfg.d3 || cfg.d3bj) {
        const int na = static_cast<int>(atoms.size());
        std::vector<int> Z(static_cast<std::size_t>(na));
        std::vector<std::array<double, 3>> pos(static_cast<std::size_t>(na));
        for (int i = 0; i < na; ++i) {
            Z[static_cast<std::size_t>(i)] =
                atoms[static_cast<std::size_t>(i)].z;
            pos[static_cast<std::size_t>(i)] = {
                atoms[static_cast<std::size_t>(i)].pos[0],
                atoms[static_cast<std::size_t>(i)].pos[1],
                atoms[static_cast<std::size_t>(i)].pos[2]};
        }
        constexpr double kAng2Bohr = 1.0 / 0.52917720859;
        constexpr double kRcutBjBohr = 2.0 * 20.0 * kAng2Bohr;
        const D3Result d3 =
            cfg.d3 ? dftd3_pbe(Z, pos, cfg.cell, periodic, D3Params{})
                   : dftd3_pbe_bj(Z, pos, cfg.cell, periodic, D3Params{},
                                  kRcutBjBohr, 1.0, 0.7875, 0.4289, 4.4407,
                                  /*cp2k_cn=*/true);
        e_disp = d3.energy;
    }

    // −Tr(W S): Pulay contribution to the Lagrangian (constant −Tr(W S0) drops
    // out of the strain FD; only ∂S/∂ε matters).
    return e_kin + e_pp + e_locgrid + e_es + e_xc + e_disp - e_pulay;
}

void validate_stress_args(const std::vector<GpwAtom>& atoms, double h,
                          const std::array<bool, 3>& axes, const char* tag) {
    if (atoms.empty())
        throw std::invalid_argument(std::string(tag) + ": empty atom list");
    if (!(h > 0.0) || !(h < 0.1))
        throw std::invalid_argument(std::string(tag)
                                    + ": strain step h out of range");
    if (!axes[0] && !axes[1] && !axes[2])
        throw std::invalid_argument(std::string(tag) + ": all axes masked");
}

} // namespace

std::array<double, 6> gpw_stress_analytic_impl(
    const std::vector<GpwAtom>& atoms, const GpwDft& method, double h,
    const std::array<bool, 3>& axes, const linalg::RealMatrix* guess_C) {
    validate_stress_args(atoms, h, axes, "gpw_stress");

    const double vol = method.cell[0] * method.cell[1] * method.cell[2];
    if (!(vol > 0.0))
        throw std::invalid_argument("gpw_stress: non-positive cell volume");

    // Reference SCF → fixed density matrix P.
    GpwDft solver = copy_gpw(method);
    solver.build(atoms);
    const GpwResult scf = solver.run(atoms, guess_C);
    if (!scf.converged)
        throw std::runtime_error("gpw_stress: reference SCF did not converge");
    if (scf.P.size() == 0)
        throw std::runtime_error(
            "gpw_stress: reference SCF produced empty density matrix");

    // Closed-shell RKS density + energy-weighted density W = Σ_i 2 ε_i c_i c_iᵀ.
    // UKS rejected in energy_fixed_P.
    const RealMatrix& P = scf.P;
    const int nb = static_cast<int>(P.shape(0));
    RealMatrix W = xt::zeros<double>({nb, nb});
    const int n_occ = method.n_occ;
    if (n_occ <= 0)
        throw std::invalid_argument(
            "gpw_stress analytic: n_occ must be positive for RKS");
    if (static_cast<int>(scf.eps.size()) < n_occ
        || static_cast<int>(scf.C.shape(0)) != nb
        || static_cast<int>(scf.C.shape(1)) < n_occ)
        throw std::runtime_error(
            "gpw_stress analytic: converged C/eps missing for W build");
    for (int i = 0; i < n_occ; ++i) {
        const double w = 2.0 * scf.eps[static_cast<std::size_t>(i)];
        for (int mu = 0; mu < nb; ++mu) {
            const double cmu = w * scf.C(mu, i);
            for (int nu = 0; nu < nb; ++nu)
                W(mu, nu) += cmu * scf.C(nu, i);
        }
    }

    std::array<double, 6> stress{};
    for (int a = 0; a < 3; ++a) {
        if (!axes[static_cast<std::size_t>(a)])
            continue;
        const double sp = 1.0 + h;
        const double sm = 1.0 - h;

        std::vector<GpwAtom> plus = scale_axis(atoms, a, sp);
        std::vector<GpwAtom> minus = scale_axis(atoms, a, sm);

        GpwDft mp = copy_gpw(method);
        GpwDft mm = copy_gpw(method);
        mp.cell[a] = method.cell[a] * sp;
        mm.cell[a] = method.cell[a] * sm;

        const double ep = energy_fixed_P(plus, mp, P, W);
        const double em = energy_fixed_P(minus, mm, P, W);
        // σ_aa = −(1/Ω) ∂E/∂ε_aa  (CP2K STRESS_TENSOR convention).
        stress[static_cast<std::size_t>(a)] = -(ep - em) / (2.0 * h) / vol;
    }
    stress[3] = stress[4] = stress[5] = 0.0;
    return stress;
}

std::array<double, 6> gpw_stress_fd_impl(const std::vector<GpwAtom>& atoms,
                                         const GpwDft& method, double h,
                                         const std::array<bool, 3>& axes,
                                         const linalg::RealMatrix* guess_C) {
    validate_stress_args(atoms, h, axes, "gpw_stress_fd");

    const double vol = method.cell[0] * method.cell[1] * method.cell[2];
    if (!(vol > 0.0))
        throw std::invalid_argument("gpw_stress_fd: non-positive cell volume");

    std::array<double, 6> stress{};
    for (int a = 0; a < 3; ++a) {
        if (!axes[static_cast<std::size_t>(a)])
            continue;
        const double sp = 1.0 + h;
        const double sm = 1.0 - h;

        std::vector<GpwAtom> plus = scale_axis(atoms, a, sp);
        std::vector<GpwAtom> minus = scale_axis(atoms, a, sm);

        GpwDft mp = copy_gpw(method);
        GpwDft mm = copy_gpw(method);
        mp.cell[a] = method.cell[a] * sp;
        mm.cell[a] = method.cell[a] * sm;

        const double ep = scf_energy_ha(plus, mp, guess_C);
        const double em = scf_energy_ha(minus, mm, guess_C);
        stress[static_cast<std::size_t>(a)] = -(ep - em) / (2.0 * h) / vol;
    }
    stress[3] = stress[4] = stress[5] = 0.0;
    return stress;
}

} // namespace internal

std::array<double, 6> gpw_stress(const std::vector<GpwAtom>& atoms,
                                 const GpwDft& method, double h,
                                 const std::array<bool, 3>& axes,
                                 const linalg::RealMatrix* guess_C) {
    return internal::gpw_stress_analytic_impl(atoms, method, h, axes, guess_C);
}

std::array<double, 6> gpw_stress_fd(const std::vector<GpwAtom>& atoms,
                                    const GpwDft& method, double h,
                                    const std::array<bool, 3>& axes,
                                    const linalg::RealMatrix* guess_C) {
    return internal::gpw_stress_fd_impl(atoms, method, h, axes, guess_C);
}

} // namespace gpw
} // namespace cpu
} // namespace ATV
