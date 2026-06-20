/// @file src/cpu/gpw/internal/gpw_forces.cpp
/// @brief Analytic GPW forces (see header).

#include "atomiverse/cpu/gpw/internal/gpw_forces.h"

#include <cmath>
#include <stdexcept>
#include <vector>

#include "atomiverse/cpu/gpw/internal/collocate.h"
#include "atomiverse/cpu/gpw/internal/gth_integrals.h"
#include "atomiverse/cpu/gpw/internal/mt_poisson.h"
#include "atomiverse/cpu/gpw/internal/one_electron.h"
#include "atomiverse/cpu/gpw/internal/poisson_fft.h"
#include "atomiverse/cpu/gpw/internal/xc_grid.h"
#include "atomiverse/cpu/gpw/internal/xc_grid_gga.h"
#include "atomiverse/cpu/linalg/matrix.h"

namespace ATV {
namespace cpu {
namespace gpw {

namespace {

using linalg::RealMatrix;

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
        for (const auto& bf : a.basis) b.push_back(bf);
    return b;
}

// d φ / d(center_k) = 2α φ_{l_k+1} − l_k φ_{l_k−1}  (Gaussian center derivative).
GridBasisFunction deriv_bf(const GridBasisFunction& bf, int k) {
    GridBasisFunction d;
    d.reserve(bf.size() * 2);
    for (const auto& p : bf) {
        GridPrimitive up = p;
        up.l[k] += 1;
        up.coeff = p.coeff * 2.0 * p.alpha;
        d.push_back(up);
        if (p.l[k] > 0) {
            GridPrimitive dn = p;
            dn.coeff = p.coeff * (-static_cast<double>(p.l[k]));
            dn.l[k] -= 1;
            d.push_back(dn);
        }
    }
    return d;
}

// <a| v_C |b> for a single GTH center C (local short-range + nonlocal channels).
double vpp_center(const GridBasisFunction& a, const GridBasisFunction& b,
                  const GpwAtom& C) {
    double v = 0.0;
    for (const auto& p : a) {
        const Prim pp = to_prim(p);
        for (const auto& q : b) {
            const double cc = p.coeff * q.coeff;
            const Prim qq = to_prim(q);
            if (!C.cloc.empty()) v += cc * gth_local_sr(pp, qq, C.pos, C.rloc, C.cloc);
            for (const auto& ch : C.nonlocal) v += cc * gth_nonlocal(pp, qq, C.pos, ch);
        }
    }
    return v;
}

// <a| Σ_C v_C |b> over all GTH centers.
double vpp_all(const GridBasisFunction& a, const GridBasisFunction& b,
               const std::vector<GpwAtom>& atoms) {
    double v = 0.0;
    for (const auto& C : atoms) v += vpp_center(a, b, C);
    return v;
}

} // namespace

std::vector<std::array<double, 3>>
gpw_forces_analytic(const std::vector<GpwAtom>& atoms, const GpwConfig& cfg) {
    const GpwResult scf = gpw_scf(atoms, cfg);
    if (!scf.converged)
        throw std::runtime_error("gpw_forces_analytic: SCF did not converge");

    const std::vector<GridBasisFunction> basis = global_basis(atoms);
    const int nb = static_cast<int>(basis.size());
    const int na = static_cast<int>(atoms.size());

    // Map each atom to its contiguous block of global basis indices, and each
    // global basis index back to its owning atom.
    std::vector<int> atom_off(na + 1, 0);
    for (int a = 0; a < na; ++a)
        atom_off[a + 1] = atom_off[a] + static_cast<int>(atoms[a].basis.size());
    std::vector<int> owner(nb, 0);
    for (int a = 0; a < na; ++a)
        for (int u = atom_off[a]; u < atom_off[a + 1]; ++u) owner[u] = a;

    const RealMatrix& P = scf.P;
    // Energy-weighted density matrix W_μν = Σ_i 2 ε_i C_μi C_νi (occupied).
    RealMatrix W = xt::zeros<double>({nb, nb});
    for (int i = 0; i < cfg.n_occ; ++i) {
        const double w = 2.0 * scf.eps[static_cast<std::size_t>(i)];
        for (int mu = 0; mu < nb; ++mu) {
            const double cmu = w * scf.C(mu, i);
            for (int nu = 0; nu < nb; ++nu) W(mu, nu) += cmu * scf.C(nu, i);
        }
    }

    // Precompute derivative basis functions d φ_u / dR_u along each axis.
    std::vector<std::array<GridBasisFunction, 3>> dphi(nb);
    for (int u = 0; u < nb; ++u)
        for (int k = 0; k < 3; ++k) dphi[u][k] = deriv_bf(basis[u], k);

    // Recompute the converged grid quantities once (one pass, not 6N+1):
    // ρ → V_es (Poisson of the neutral total charge) and V_eff = V_es + V_xc.
    UniformGrid grid = make_uniform_grid(cfg.cell, cfg.cutoff_ry);
    const long ng = grid.npts();
    const double wgt = grid.weight();
    const std::vector<double> phi = eval_basis_on_grid(grid, basis);
    const std::vector<long> active = active_grid_points(phi, ng, nb);

    // Ionic core density n_c(r) = Σ_A z_A (2π s_A²)^{-3/2} exp(-|r-R_A|²/2 s_A²),
    // matching gpw_solver::ionic_density (no min-image; Gaussians well inside box).
    std::vector<double> n_ion(static_cast<std::size_t>(ng), 0.0);
    for (long g = 0; g < ng; ++g) {
        double r[3];
        grid.point(g, r);
        double s = 0.0;
        for (const auto& a : atoms) {
            const double s2 = a.rloc * a.rloc;
            const double norm = a.z_val * std::pow(2.0 * M_PI * s2, -1.5);
            const double dx = r[0] - a.pos[0], dy = r[1] - a.pos[1], dz = r[2] - a.pos[2];
            s += norm * std::exp(-(dx * dx + dy * dy + dz * dz) / (2.0 * s2));
        }
        n_ion[static_cast<std::size_t>(g)] = s;
    }

    std::vector<double> Pvec(static_cast<std::size_t>(nb) * nb);
    for (int i = 0; i < nb; ++i)
        for (int j = 0; j < nb; ++j)
            Pvec[static_cast<std::size_t>(i) * nb + j] = P(i, j);
    std::vector<double> rho = collocate_density(grid, phi, nb, Pvec, active);
    std::vector<double> n_tot(rho.size());
    for (std::size_t g = 0; g < rho.size(); ++g) n_tot[g] = rho[g] - n_ion[static_cast<std::size_t>(g)];

    const PoissonResult es =
        cfg.nonperiodic ? solve_poisson_mt(grid, n_tot) : solve_poisson(grid, n_tot);
    const XcOnGrid xc =
        cfg.xc == XcKind::PBE ? xc_grid_pbe(grid, rho) : xc_grid_pade(grid, rho);
    const std::vector<double>& v_es = es.v_r;
    std::vector<double> veff(rho.size());
    for (std::size_t g = 0; g < rho.size(); ++g) veff[g] = es.v_r[g] + xc.vxc_r[g];

    // Pφ_μ(g) = Σ_ν P_μν φ_ν(g), over the active grid points (for the grid term).
    // Stored as Pphi[active_index][mu].
    const std::size_t nact = active.size();
    std::vector<double> Pphi(nact * static_cast<std::size_t>(nb), 0.0);
    for (std::size_t ai = 0; ai < nact; ++ai) {
        const long g = active[ai];
        const double* phig = &phi[static_cast<std::size_t>(g) * nb];
        for (int mu = 0; mu < nb; ++mu) {
            double acc = 0.0;
            for (int nu = 0; nu < nb; ++nu) acc += P(mu, nu) * phig[nu];
            Pphi[ai * static_cast<std::size_t>(nb) + mu] = acc;
        }
    }

    std::vector<std::array<double, 3>> forces(na, {0.0, 0.0, 0.0});

    for (int A = 0; A < na; ++A) {
        for (int k = 0; k < 3; ++k) {
            double f = 0.0;

            // --- AO-basis Pulay terms (symmetric operators: 2× the μ∈A part) ---
            for (int mu = atom_off[A]; mu < atom_off[A + 1]; ++mu) {
                const GridBasisFunction& dmu = dphi[mu][k];
                for (int nu = 0; nu < nb; ++nu) {
                    const double s = overlap(dmu, basis[nu]);
                    const double t = kinetic(dmu, basis[nu]);
                    f += 2.0 * W(mu, nu) * s;        // + W dS
                    f -= 2.0 * P(mu, nu) * t;        // - P dT
                    f -= 2.0 * P(mu, nu) * vpp_all(dmu, basis[nu], atoms); // - P dV_pp (basis)
                }
            }

            // --- GTH operator-center derivative (∂_C = −(∂_μ+∂_ν), trans. inv.) ---
            // + P_μν ∂_A<μ|v_A|ν> summed; via symmetry = 2 Σ_{all μ,ν} P_μν <dμ|v_A|ν>.
            const GpwAtom& Aatom = atoms[A];
            for (int mu = 0; mu < nb; ++mu) {
                const GridBasisFunction& dmu = dphi[mu][k];
                for (int nu = 0; nu < nb; ++nu)
                    f += 2.0 * P(mu, nu) * vpp_center(dmu, basis[nu], Aatom);
            }

            // --- Grid term: -2 Σ_{μ∈A} Σ_g w dφ_μ(g) V_eff(g) (Σ_ν P_μν φ_ν(g)) ---
            for (int mu = atom_off[A]; mu < atom_off[A + 1]; ++mu) {
                const std::vector<double> dmu_g =
                    eval_basis_on_grid(grid, {dphi[mu][k]}); // length ng (nbasis=1)
                double acc = 0.0;
                for (std::size_t ai = 0; ai < nact; ++ai) {
                    const long g = active[ai];
                    acc += dmu_g[static_cast<std::size_t>(g)] * veff[static_cast<std::size_t>(g)] *
                           Pphi[ai * static_cast<std::size_t>(nb) + mu];
                }
                f -= 2.0 * wgt * acc;
            }

            // --- Ionic Hellmann–Feynman: + ∫ V_es ∂n_c^A/∂R_A,k ---
            // ∂n_c^A/∂R_A,k = n_c^A(r) (r_k - R_A,k)/s_A².
            {
                const double s2 = Aatom.rloc * Aatom.rloc;
                const double norm = Aatom.z_val * std::pow(2.0 * M_PI * s2, -1.5);
                double acc = 0.0;
                for (long g = 0; g < ng; ++g) {
                    double r[3];
                    grid.point(g, r);
                    const double dx = r[0] - Aatom.pos[0];
                    const double dy = r[1] - Aatom.pos[1];
                    const double dz = r[2] - Aatom.pos[2];
                    const double gauss = norm * std::exp(-(dx * dx + dy * dy + dz * dz) / (2.0 * s2));
                    const double dr_k = (k == 0 ? dx : (k == 1 ? dy : dz));
                    acc += v_es[static_cast<std::size_t>(g)] * gauss * dr_k / s2;
                }
                f += wgt * acc;
            }

            forces[A][k] = f;
        }
    }

    return forces;
}

} // namespace gpw
} // namespace cpu
} // namespace ATV
