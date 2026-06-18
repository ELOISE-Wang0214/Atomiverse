/// @file src/cpu/gpw/gpw_solver.cpp
/// @brief Implementation of the closed-shell GPW SCF (see header).
///
/// Electrostatics use the GTH Gaussian core charge: the long-range local PP
/// -Z erf(r/(sqrt2 rloc))/r is the potential of a negative Gaussian -Z g(rloc).
/// Solving one periodic Poisson of the NEUTRAL total charge n = rho_e - n_ion
/// yields V_es = v_H + v_ext^lr in a single FFT, and its energy
/// 1/2 integral n V_es = E_H + E_ext^lr + E_ii^smeared. The spurious Gaussian
/// self-energy of each core is removed analytically; the erf-vs-point ion-ion
/// correction is a documented TODO for the CP2K bit-parity milestone.

#include "atomiverse/cpu/gpw/gpw_solver.h"

#include <array>
#include <cmath>
#include <vector>

#include "atomiverse/cpu/gpw/internal/gth_integrals.h"
#include "atomiverse/cpu/gpw/internal/hartree_grid.h"
#include "atomiverse/cpu/gpw/internal/mt_poisson.h"
#include "atomiverse/cpu/gpw/internal/one_electron.h"
#include "atomiverse/cpu/gpw/internal/poisson_fft.h"
#include "atomiverse/cpu/gpw/internal/xc_grid.h"
#include "atomiverse/cpu/gpw/internal/xc_grid_gga.h"
#include "atomiverse/cpu/linalg/density.h"    // build_density_matrix
#include "atomiverse/cpu/linalg/eigensolve.h"  // diagonalize_fock_generalized
#include "atomiverse/cpu/linalg/matrix.h"      // RealMatrix

#ifdef ATV_GPW_FFTW_THREADS
#include <fftw3.h>
#include <omp.h>
#endif

namespace ATV {
namespace cpu {
namespace gpw {

namespace {

using ATV::cpu::linalg::RealMatrix;

// Flatten all atoms' basis functions into one global list.
std::vector<GridBasisFunction> global_basis(const std::vector<GpwAtom>& atoms) {
    std::vector<GridBasisFunction> b;
    for (const auto& a : atoms)
        for (const auto& bf : a.basis) b.push_back(bf);
    return b;
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

// <mu| V_loc^sr + V_nl |nu> summed over all atoms acting as PP centers.
double gth_matrix_element(const GridBasisFunction& mu, const GridBasisFunction& nu,
                          const std::vector<GpwAtom>& atoms) {
    double v = 0.0;
    for (const auto& a : atoms) {
        for (const auto& p : mu) {
            for (const auto& q : nu) {
                const double cc = p.coeff * q.coeff;
                const Prim pp = to_prim(p);
                const Prim qq = to_prim(q);
                if (!a.cloc.empty())
                    v += cc * gth_local_sr(pp, qq, a.pos, a.rloc, a.cloc);
                for (const auto& ch : a.nonlocal)
                    v += cc * gth_nonlocal(pp, qq, a.pos, ch);
            }
        }
    }
    return v;
}

// Ionic core density n_ion(r) = sum_A z_A (2 pi rloc^2)^{-3/2} exp(-|r-R|^2/2 rloc^2).
std::vector<double> ionic_density(const UniformGrid& grid,
                                  const std::vector<GpwAtom>& atoms) {
    const long ng = grid.npts();
    std::vector<double> n(static_cast<std::size_t>(ng), 0.0);
#pragma omp parallel for schedule(static)
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
        n[g] = s;
    }
    return n;
}

// Analytic Hartree self-energy of one normalized Gaussian core of charge z and
// width rloc: 1/2 integral rho V[rho] = z^2 / (2 rloc sqrt(pi)).
double core_self_energy(const GpwAtom& a) {
    return a.z_val * a.z_val / (2.0 * a.rloc * std::sqrt(M_PI));
}

RealMatrix vec_to_mat(const std::vector<double>& v, int n) {
    RealMatrix m = xt::zeros<double>({n, n});
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) m(i, j) = v[static_cast<std::size_t>(i) * n + j];
    return m;
}

std::vector<double> mat_to_vec(const RealMatrix& m, int n) {
    std::vector<double> v(static_cast<std::size_t>(n) * n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) v[static_cast<std::size_t>(i) * n + j] = m(i, j);
    return v;
}

double trace_product(const RealMatrix& P, const RealMatrix& M, int n) {
    double t = 0.0;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) t += P(i, j) * M(i, j);
    return t;
}

} // namespace

namespace {
// One Cartesian term of a real solid harmonic: powers + weight.
struct CartTerm {
    int l[3];
    double w;
};
// Real solid harmonics as Cartesian combinations, per angular momentum.
// s: 1; p: 3 (x,y,z); d: 5 spherical (xy, yz, xz, x^2-y^2, 2z^2-x^2-y^2) so the
// spanned space matches CP2K (Cartesian d would add a spurious s-like 6th).
// Internal relative weights matter (they define each function); the overall
// scale of each function is irrelevant to the energy (basis-norm invariance).
const std::vector<std::vector<std::vector<CartTerm>>>& harmonics() {
    static const std::vector<std::vector<std::vector<CartTerm>>> H = {
        {{{{0, 0, 0}, 1.0}}},
        {{{{1, 0, 0}, 1.0}}, {{{0, 1, 0}, 1.0}}, {{{0, 0, 1}, 1.0}}},
        {{{{1, 1, 0}, 1.0}},
         {{{0, 1, 1}, 1.0}},
         {{{1, 0, 1}, 1.0}},
         {{{2, 0, 0}, 1.0}, {{0, 2, 0}, -1.0}},
         {{{0, 0, 2}, 2.0}, {{2, 0, 0}, -1.0}, {{0, 2, 0}, -1.0}}},
    };
    return H;
}
} // namespace

GpwAtom make_gpw_atom(const double pos[3], const GaussianBasis& basis,
                      const GthPotential& pot) {
    GpwAtom atom;
    atom.z_val = static_cast<double>(pot.zion);
    atom.pos[0] = pos[0];
    atom.pos[1] = pos[1];
    atom.pos[2] = pos[2];
    atom.rloc = pot.rloc;
    atom.cloc = pot.cloc;
    atom.nonlocal = pot.channels;

    for (const auto& set : basis.sets) {
        int col = 0;
        for (int l = set.lmin; l <= set.lmax; ++l) {
            const int L = l;
            for (int sh = 0; sh < set.nshell_l(l); ++sh, ++col) {
                for (const auto& harm : harmonics()[static_cast<std::size_t>(L)]) {
                    GridBasisFunction bf;
                    for (int e = 0; e < set.nexp(); ++e) {
                        const double a = set.exps[static_cast<std::size_t>(e)];
                        const double c = set.coef_at(e, col);
                        if (c == 0.0) continue;
                        // N(a,L) = (2a/pi)^{3/4} (4a)^{L/2}: the L alpha-power fixes
                        // the spanned space to CP2K's; constant factors cancel.
                        const double norm =
                            std::pow(2.0 * a / M_PI, 0.75) * std::pow(4.0 * a, 0.5 * L);
                        for (const auto& term : harm) {
                            GridPrimitive p;
                            p.alpha = a;
                            p.center[0] = pos[0];
                            p.center[1] = pos[1];
                            p.center[2] = pos[2];
                            p.l[0] = term.l[0];
                            p.l[1] = term.l[1];
                            p.l[2] = term.l[2];
                            p.coeff = c * norm * term.w;
                            bf.push_back(p);
                        }
                    }
                    if (!bf.empty()) atom.basis.push_back(std::move(bf));
                }
            }
        }
    }
    return atom;
}

GpwResult gpw_scf(const std::vector<GpwAtom>& atoms, const GpwConfig& cfg) {
#ifdef ATV_GPW_FFTW_THREADS
    // One-time global FFTW threading init; fftw_plan_with_nthreads is a global
    // setting, so every FFT plan created downstream (Poisson, GGA grad/div, MT)
    // runs multithreaded. The density-grid FFTs dominate GPW's wall time.
    static const int kFftwThreadsOk = fftw_init_threads();
    if (kFftwThreadsOk) fftw_plan_with_nthreads(omp_get_max_threads());
#endif
    const std::vector<GridBasisFunction> basis = global_basis(atoms);
    const int nb = static_cast<int>(basis.size());
    UniformGrid grid = make_uniform_grid(cfg.cell, cfg.cutoff_ry);

    // Geometry-fixed matrices.
    RealMatrix S = xt::zeros<double>({nb, nb});
    RealMatrix T = xt::zeros<double>({nb, nb});
    RealMatrix Vpp = xt::zeros<double>({nb, nb}); // V_loc^sr + V_nl
    for (int i = 0; i < nb; ++i) {
        for (int j = 0; j < nb; ++j) {
            S(i, j) = overlap(basis[i], basis[j]);
            T(i, j) = kinetic(basis[i], basis[j]);
            Vpp(i, j) = gth_matrix_element(basis[i], basis[j], atoms);
        }
    }
    const RealMatrix H0 = T + Vpp;

    // Geometry-fixed ionic density and the constant self-energy to remove.
    const std::vector<double> n_ion = ionic_density(grid, atoms);
    double e_self = 0.0;
    for (const auto& a : atoms) e_self += core_self_energy(a);

    // Core (H0) initial guess.
    auto [eig0, C0] = linalg::diagonalize_fock_generalized(H0, S);
    RealMatrix P = linalg::build_density_matrix(C0, cfg.n_occ, nb, 2.0);

    // Basis-on-grid is geometry-fixed: evaluate the Gaussians once and reuse
    // across all SCF iterations. The molecule fills a small part of the box, so
    // also precompute the "active" grid points (where any basis function is
    // non-negligible) and collocate/integrate only over those.
    const std::vector<double> phi = eval_basis_on_grid(grid, basis);
    const std::vector<long> active = active_grid_points(phi, grid.npts(), nb);

    GpwResult res;
    double e_prev = 0.0;
    for (int iter = 0; iter < cfg.max_iter; ++iter) {
        // Density on the grid and the neutral total charge.
        const std::vector<double> Pvec = mat_to_vec(P, nb);
        std::vector<double> rho = collocate_density(grid, phi, nb, Pvec, active);
        std::vector<double> n_tot(rho.size());
        for (std::size_t g = 0; g < rho.size(); ++g) n_tot[g] = rho[g] - n_ion[g];

        // Electrostatics (one Poisson) and XC — periodic or isolated, LDA or GGA.
        const PoissonResult es =
            cfg.nonperiodic ? solve_poisson_mt(grid, n_tot) : solve_poisson(grid, n_tot);
        const XcOnGrid xc =
            cfg.xc == XcKind::PBE ? xc_grid_pbe(grid, rho) : xc_grid_pade(grid, rho);

        std::vector<double> veff(rho.size());
        for (std::size_t g = 0; g < rho.size(); ++g) veff[g] = es.v_r[g] + xc.vxc_r[g];
        const RealMatrix Vgrid = vec_to_mat(integrate_potential(grid, phi, nb, veff, active), nb);

        const RealMatrix H = H0 + Vgrid;

        // Energy (e_es already holds E_H + E_ext^lr + E_ii^smeared; remove the
        // spurious Gaussian self-energy of the cores).
        res.e_kinetic = trace_product(P, T, nb);
        res.e_nonlocal = 0.0; // folded into e_loc_sr below for the first cut
        res.e_loc_sr = trace_product(P, Vpp, nb);
        res.e_es = es.e_hartree - e_self;
        res.e_xc = xc.e_xc;
        const double e_total = res.e_kinetic + res.e_loc_sr + res.e_es + res.e_xc;

        // New density.
        auto [eig, C] = linalg::diagonalize_fock_generalized(H, S);
        RealMatrix P_new = linalg::build_density_matrix(C, cfg.n_occ, nb, 2.0);

        res.iters = iter + 1;
        res.e_total = e_total;
        if (iter > 0 && std::abs(e_total - e_prev) < cfg.tol) {
            res.converged = true;
            break;
        }
        e_prev = e_total;

        // Linear density mixing.
        P = (1.0 - cfg.mix) * P + cfg.mix * P_new;
    }
    return res;
}

namespace {
// Displace an atom (its position AND all its basis-function centers) by d.
GpwAtom displaced(GpwAtom a, int axis, double d) {
    a.pos[axis] += d;
    for (auto& bf : a.basis)
        for (auto& p : bf) p.center[axis] += d;
    return a;
}
} // namespace

std::vector<std::array<double, 3>>
gpw_forces_fd(const std::vector<GpwAtom>& atoms, const GpwConfig& cfg, double h) {
    const int na = static_cast<int>(atoms.size());
    std::vector<std::array<double, 3>> forces(na, {0.0, 0.0, 0.0});
    for (int a = 0; a < na; ++a) {
        for (int k = 0; k < 3; ++k) {
            std::vector<GpwAtom> plus = atoms, minus = atoms;
            plus[a] = displaced(atoms[a], k, +h);
            minus[a] = displaced(atoms[a], k, -h);
            const double e_plus = gpw_scf(plus, cfg).e_total;
            const double e_minus = gpw_scf(minus, cfg).e_total;
            forces[a][k] = -(e_plus - e_minus) / (2.0 * h); // F = -dE/dR
        }
    }
    return forces;
}

} // namespace gpw
} // namespace cpu
} // namespace ATV
