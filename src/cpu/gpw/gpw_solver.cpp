#include <cstdio>
#include <cstdlib>
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

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

#include "atomiverse/cpu/gpw/internal/diis.h"
#include "atomiverse/cpu/gpw/internal/gth_integrals.h"
#include "atomiverse/cpu/gpw/internal/hartree_grid.h"
#include "atomiverse/cpu/gpw/internal/mt_poisson.h"
#include "atomiverse/cpu/gpw/internal/one_electron.h"
#include "atomiverse/cpu/gpw/internal/ot.h"
#include "atomiverse/cpu/gpw/internal/poisson_fft.h"
#include "atomiverse/cpu/gpw/internal/real_harmonics.h"
#include "atomiverse/cpu/gpw/internal/xc_grid.h"
#include "atomiverse/cpu/gpw/internal/xc_grid_gga.h"
#include "atomiverse/cpu/linalg/density.h"    // build_density_matrix
#include "atomiverse/cpu/linalg/eigensolve.h"  // diagonalize_fock_generalized
#include "atomiverse/cpu/linalg/matrix.h"      // RealMatrix

#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/views/xview.hpp>

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

// Screening tolerance for the geometry-fixed GTH/one-electron build. A term is
// skipped only when a rigorous Gaussian-magnitude bound falls below this — the
// exponential decay dominates, so dropped terms are far below the SCF energy
// tolerance (validated bit-unchanged against the small-system CP2K-parity tests).
// This is a compute-time saving in the O(nb^2*natoms) build; the stored matrices
// remain dense (the overlap matrices are not sparse for this geometry).
constexpr double kBuildScreen = 1e-12;

// <mu| V_loc^sr + V_nl |nu> summed over all atoms acting as PP centers.
//
// Screened: a primitive pair (p,q) is skipped when its Gaussian product
// prefactor exp(-a*b/(a+b)|Rp-Rq|^2) is negligible; for a surviving pair, an
// atom A is skipped when the combined product-PP decay exp(-(g*k/(g+k))|P-A|^2)
// (g=a+b, k=1/(2 rmax^2), P the product center, rmax the atom's largest PP
// width) puts its contribution below the screen. The |P-A| local-style bound
// over-estimates the separable nonlocal magnitude when p,q overlap, so it is
// conservative for both the local and nonlocal terms. Turns the per-pair atom
// loop from O(natoms) into O(local neighbours) -> overall O(N) build.
double gth_matrix_element(const GridBasisFunction& mu, const GridBasisFunction& nu,
                          const std::vector<GpwAtom>& atoms) {
    double v = 0.0;
    for (const auto& p : mu) {
        const Prim pp = to_prim(p);
        for (const auto& q : nu) {
            const double cc = p.coeff * q.coeff;
            if (cc == 0.0) continue;
            const double acc = std::abs(cc);
            const double a = p.alpha, b = q.alpha, g = a + b;
            double dpq2 = 0.0, P[3];
            for (int d = 0; d < 3; ++d) {
                const double diff = p.center[d] - q.center[d];
                dpq2 += diff * diff;
                P[d] = (a * p.center[d] + b * q.center[d]) / g;
            }
            const double Kpq = std::exp(-a * b / g * dpq2);
            if (acc * Kpq < kBuildScreen) continue; // overlap-negligible primitive pair
            const Prim qq = to_prim(q);
            for (const auto& at : atoms) {
                double rmax = at.rloc;
                for (const auto& ch : at.nonlocal) rmax = std::max(rmax, ch.r);
                const double kappa = 1.0 / (2.0 * rmax * rmax);
                double dPA2 = 0.0;
                for (int d = 0; d < 3; ++d) {
                    const double e = P[d] - at.pos[d];
                    dPA2 += e * e;
                }
                if (acc * Kpq * std::exp(-(g * kappa / (g + kappa)) * dPA2) < kBuildScreen)
                    continue; // atom too far for any non-negligible PP contribution
                if (!at.cloc.empty())
                    v += cc * gth_local_sr(pp, qq, at.pos, at.rloc, at.cloc);
                for (const auto& ch : at.nonlocal)
                    v += cc * gth_nonlocal(pp, qq, at.pos, ch);
            }
        }
    }
    return v;
}

// Max Gaussian-product prefactor over all primitive pairs of (mu,nu): bounds the
// magnitude of S, T, and Vpp for this pair. Below kBuildScreen the whole (i,j)
// element is negligible and can be skipped (overlap pair screening).
double pair_overlap_bound(const GridBasisFunction& mu, const GridBasisFunction& nu) {
    double best = 0.0;
    for (const auto& p : mu)
        for (const auto& q : nu) {
            const double a = p.alpha, b = q.alpha;
            double d2 = 0.0;
            for (int k = 0; k < 3; ++k) {
                const double e = p.center[k] - q.center[k];
                d2 += e * e;
            }
            best = std::max(best, std::abs(p.coeff * q.coeff) * std::exp(-a * b / (a + b) * d2));
        }
    return best;
}

// Ionic core density n_ion(r) = sum_A z_A (2 pi rloc^2)^{-3/2} exp(-|r-R|^2/2 rloc^2).
// Periodic: minimum-image displacement so unwrapped slab cores collocate through
// their nearest image (matching the electron-density collocation in collocate.cpp).
std::vector<double> ionic_density(const UniformGrid& grid,
                                  const std::vector<GpwAtom>& atoms, bool periodic) {
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
            double dx = r[0] - a.pos[0], dy = r[1] - a.pos[1], dz = r[2] - a.pos[2];
            if (periodic) {
                dx -= grid.L[0] * std::nearbyint(dx / grid.L[0]);
                dy -= grid.L[1] * std::nearbyint(dy / grid.L[1]);
                dz -= grid.L[2] * std::nearbyint(dz / grid.L[2]);
            }
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

// On convergence, diagonalize the converged KS matrix once more to record a
// self-consistent (C, eps, P) triple in the result (P == density(C)). This is
// what the analytic-force assembly needs to build the energy-weighted density
// matrix W = Σ_i 2 ε_i c_i c_iᵀ.
void store_converged_state(GpwResult& res, const RealMatrix& H, const RealMatrix& S,
                           int n_occ, int nb) {
    auto [eig, C] = linalg::diagonalize_fock_generalized(H, S);
    res.C = C;
    res.eps.assign(eig.begin(), eig.end());
    res.P = linalg::build_density_matrix(C, n_occ, nb, 2.0);
}

} // namespace

namespace {
// One Cartesian term of a real solid harmonic: powers + weight.
struct CartTerm {
    int l[3];
    double w;
};
// Real solid harmonics as Cartesian combinations, per angular momentum, for
// arbitrary l. Generated by real_solid_harmonics() (Use A) and cached: the
// (2l+1) pure-spherical functions span the same space as CP2K's, and the overall
// per-function scale is irrelevant to the energy (basis-norm invariance). l<=2
// reproduces the canonical s; p:{x,y,z}; d:{xy,yz,xz,x^2-y^2,2z^2-x^2-y^2}.
constexpr int kHarmLMax = 6; // s..i; covers every GTH/MOLOPT basis shell.
const std::vector<std::vector<std::vector<CartTerm>>>& harmonics() {
    static const std::vector<std::vector<std::vector<CartTerm>>> H = [] {
        std::vector<std::vector<std::vector<CartTerm>>> table;
        table.reserve(kHarmLMax + 1);
        for (int l = 0; l <= kHarmLMax; ++l) {
            std::vector<std::vector<CartTerm>> shell;
            for (const auto& harm : real_solid_harmonics(l)) {
                std::vector<CartTerm> terms;
                for (const auto& t : harm)
                    terms.push_back({{t.lx, t.ly, t.lz}, t.coeff});
                shell.push_back(std::move(terms));
            }
            table.push_back(std::move(shell));
        }
        return table;
    }();
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
    // The geometry-fixed one-electron / GTH matrices dominate the GPW SCF wall
    // time at scale (profiled: ~75-80% for hundreds of atoms; the eigensolve is
    // <2%). Each (i,j) is an independent primitive contraction → parallelize the
    // outer index. NOTE: gth_matrix_element still loops ALL atoms per pair
    // (O(nb²·natoms)); spatial screening of that atom loop is the next, larger
    // algorithmic win (the PP contribution is local).
#pragma omp parallel for schedule(dynamic, 8)
    for (int i = 0; i < nb; ++i) {
        for (int j = 0; j < nb; ++j) {
            // Overlap pair screening: S, T, Vpp are all bounded by the Gaussian
            // product magnitude of (i,j). If even the best primitive pair is
            // below the screen, the whole element is negligible -> leave 0.
            if (pair_overlap_bound(basis[i], basis[j]) < kBuildScreen) continue;
            S(i, j) = overlap(basis[i], basis[j]);
            T(i, j) = kinetic(basis[i], basis[j]);
            Vpp(i, j) = gth_matrix_element(basis[i], basis[j], atoms);
        }
    }
    const RealMatrix H0 = T + Vpp;

    // Geometry-fixed ionic density and the constant self-energy to remove.
    const std::vector<double> n_ion = ionic_density(grid, atoms, !cfg.nonperiodic);
    double e_self = 0.0;
    for (const auto& a : atoms) e_self += core_self_energy(a);

    // Core (H0) initial guess.
    auto [eig0, C0] = linalg::diagonalize_fock_generalized(H0, S);
    RealMatrix P = linalg::build_density_matrix(C0, cfg.n_occ, nb, 2.0);

    // Basis-on-grid is geometry-fixed: evaluate the Gaussians once and reuse
    // across all SCF iterations. The molecule fills a small part of the box, so
    // also precompute the "active" grid points (where any basis function is
    // non-negligible) and collocate/integrate only over those.
    const std::vector<double> phi = eval_basis_on_grid(grid, basis, !cfg.nonperiodic);
    const std::vector<long> active = active_grid_points(phi, grid.npts(), nb);

    GpwResult res;
    double e_prev = 0.0;

    // Convergence acceleration. With diis_depth > 0 we run Pulay DIIS on the KS
    // matrix (X = S^{-1/2} is geometry-fixed → built once); with diis_depth == 0
    // we fall back to the legacy linear density mixing. DIIS converges
    // transition-metal / semicore systems (Zn-q12, Na-q9) that linear mixing
    // stalls on, and is the prerequisite for scaling to large cells.
    const bool use_diis = cfg.diis_depth > 0;
    const RealMatrix X =
        use_diis ? linalg::build_orthogonalization_matrix(S) : RealMatrix{};
    internal::GpwDiis diis(use_diis ? cfg.diis_depth : 1);
    // DIIS commutator-residual threshold: a converged step must have both a
    // stable energy and a near-zero [H,P]_S residual (guards against declaring
    // convergence on an energy plateau before self-consistency).
    constexpr double kDiisErrTol = 1e-5;
    // OT convergence: max element of the projected gradient (the SCF error
    // vector). Paired with the energy-stability test for a tight fixed point.
    constexpr double kOtGradTol = 1e-6;

    // KS matrix + total energy from a density matrix P (collocate -> Poisson +
    // XC -> integrate). Fills res's energy terms; shared by the diagonalization
    // and OT paths.
    auto build_H_and_energy = [&](const RealMatrix& Pin) -> std::pair<RealMatrix, double> {
        const std::vector<double> Pvec = mat_to_vec(Pin, nb);
        std::vector<double> rho = collocate_density(grid, phi, nb, Pvec, active);
        std::vector<double> n_tot(rho.size());
        for (std::size_t g = 0; g < rho.size(); ++g) n_tot[g] = rho[g] - n_ion[g];
        const PoissonResult es =
            cfg.nonperiodic ? solve_poisson_mt(grid, n_tot) : solve_poisson(grid, n_tot);
        const XcOnGrid xc =
            cfg.xc == XcKind::PBE ? xc_grid_pbe(grid, rho) : xc_grid_pade(grid, rho);
        std::vector<double> veff(rho.size());
        for (std::size_t g = 0; g < rho.size(); ++g) veff[g] = es.v_r[g] + xc.vxc_r[g];
        const RealMatrix Vgrid =
            vec_to_mat(integrate_potential(grid, phi, nb, veff, active), nb);
        const RealMatrix H = H0 + Vgrid;
        res.e_kinetic = trace_product(Pin, T, nb);
        res.e_nonlocal = 0.0; // folded into e_loc_sr below for the first cut
        res.e_loc_sr = trace_product(Pin, Vpp, nb);
        res.e_es = es.e_hartree - e_self;
        res.e_xc = xc.e_xc;
        return {H, res.e_kinetic + res.e_loc_sr + res.e_es + res.e_xc};
    };

    // ── OT (orbital transformation) path ────────────────────────────────────
    // Diagonalization-free SCF: each outer step rebuilds H, then takes one
    // preconditioned-CG step on the OT manifold (gradient O(nb²·n_occ), no
    // O(nb³) eigensolve). Converges to the same gapped ground state as the
    // diagonalization path. The line search uses the frozen-H band energy
    // 2 Tr(CᵀHC) (matrix-only); self-consistency comes from the H rebuild.
    if (cfg.use_ot) {
        RealMatrix Cocc = internal::ot_orthonormalize(
            RealMatrix(xt::view(C0, xt::all(), xt::range(0, cfg.n_occ))), S);
        auto frob = [](const RealMatrix& a, const RealMatrix& b) {
            return xt::sum(a * b)();
        };
        RealMatrix d, pz_prev;
        double denom_prev = 1.0;
        bool have_prev = false;
        for (int iter = 0; iter < cfg.max_iter; ++iter) {
            P = linalg::build_density_matrix(Cocc, cfg.n_occ, nb, 2.0); // 2 Cocc Coccᵀ
            auto [H, e_total] = build_H_and_energy(P);
            res.iters = iter + 1;
            res.e_total = e_total;

            const RealMatrix g = internal::ot_gradient(H, S, Cocc);
            double gnorm = 0.0;
            for (double v : g) gnorm = std::max(gnorm, std::abs(v));
            if (std::getenv("ATV_GPW_TRACE"))
                std::fprintf(stderr,
                             "[gpw OT %2d] E=%.8f  T=%.4f loc=%.4f es=%.4f xc=%.4f  "
                             "dE=%.2e gnorm=%.2e\n",
                             iter, e_total, res.e_kinetic, res.e_loc_sr, res.e_es,
                             res.e_xc, std::abs(e_total - e_prev), gnorm);
            if (iter > 0 && gnorm < kOtGradTol && std::abs(e_total - e_prev) < cfg.tol) {
                res.converged = true;
                store_converged_state(res, H, S, cfg.n_occ, nb);
                break;
            }
            e_prev = e_total;

            // FULL_SINGLE_INVERSE preconditioner: H-aware, rebuilt each step
            // (CP2K's default for gapped systems; FULL_KINETIC T+S was too weak
            // and OT crawled). M = H + a S with a below the spectrum -> SPD.
            const internal::OtPrecondSingleInverse precond(H, S, Cocc);
            const RealMatrix pz = internal::ot_tangent(Cocc, S, precond.apply(g));
            if (!have_prev) {
                d = -pz;
            } else {
                const double beta = std::max(0.0, frob(g, pz - pz_prev) / denom_prev);
                d = internal::ot_tangent(Cocc, S, -pz + beta * d);
            }
            pz_prev = pz;
            denom_prev = frob(g, pz);
            have_prev = true;

            // Monotonic line search on the SCF TOTAL energy E[rho(C)]. The
            // frozen-H band energy 2 Tr(CᵀHC) is only the local quadratic model;
            // minimizing *it* overshoots the self-consistent E and let OT wander
            // (descend then drift up). Armijo backtracking on the rebuilt-H total
            // energy enforces monotonic descent — every accepted step lowers E.
            auto etot = [&](const RealMatrix& C) {
                return build_H_and_energy(
                           linalg::build_density_matrix(C, cfg.n_occ, nb, 2.0))
                    .second;
            };
            double slope = frob(g, d);
            if (slope >= 0.0) { // not a descent direction -> restart to (precond) SD
                d = internal::ot_tangent(Cocc, S, -pz);
                slope = frob(g, d);
            }
            // Curvature-aware initial step from the frozen-H band energy (cheap,
            // no grid), then Armijo-backtrack on the rebuilt-H total energy. The
            // band parabola gives a near-optimal first alpha so backtracking
            // (each trial = a full grid build) is rare.
            auto band = [&](const RealMatrix& C) {
                return 2.0 * frob(C, xt::linalg::dot(H, C));
            };
            const double Eb0 = band(Cocc);
            const double Ebt = band(internal::ot_parameterize(Cocc, d, S));
            const double curv = 2.0 * (Ebt - Eb0 - slope); // parabola curvature
            double alpha = (curv > 1e-300) ? -slope / curv : 1.0; // vertex
            if (!(alpha > 0.0)) alpha = 1.0;
            alpha = std::min(alpha, 2.0);
            constexpr double c1 = 1e-4;
            RealMatrix Cnew = internal::ot_parameterize(Cocc, alpha * d, S);
            for (int bt = 0; bt < 50; ++bt) {
                if (etot(Cnew) <= e_total + c1 * alpha * slope) break;
                alpha *= 0.5;
                Cnew = internal::ot_parameterize(Cocc, alpha * d, S);
            }
            Cocc = internal::ot_orthonormalize(Cnew, S);
        }
        return res;
    }

    for (int iter = 0; iter < cfg.max_iter; ++iter) {
        const auto [H, e_total] = build_H_and_energy(P);

        res.iters = iter + 1;
        res.e_total = e_total;

        if (use_diis) {
            // Pulay DIIS: extrapolate the KS matrix from the (H, residual)
            // history and rebuild P directly from its eigenvectors (no density
            // mixing — DIIS does the averaging).
            const RealMatrix err = internal::gpw_diis_error(H, P, S, X);
            double err_norm = 0.0;
            for (const double v : err) err_norm = std::max(err_norm, std::abs(v));

            if (std::getenv("ATV_GPW_TRACE"))
                std::fprintf(stderr,
                             "[gpw it %2d] E=%.8f  T=%.4f loc=%.4f es=%.4f xc=%.4f  "
                             "dE=%.2e err=%.2e\n",
                             iter, e_total, res.e_kinetic, res.e_loc_sr, res.e_es,
                             res.e_xc, std::abs(e_total - e_prev), err_norm);

            const bool no_stop = std::getenv("ATV_GPW_NOSTOP") != nullptr;
            if (!no_stop && iter > 0 && std::abs(e_total - e_prev) < cfg.tol &&
                err_norm < kDiisErrTol) {
                res.converged = true;
                store_converged_state(res, H, S, cfg.n_occ, nb);
                break;
            }
            e_prev = e_total;

            const RealMatrix H_ext = diis.step(H, err);
            auto [eig, C] = linalg::diagonalize_fock_generalized(H_ext, S);
            P = linalg::build_density_matrix(C, cfg.n_occ, nb, 2.0);
        } else {
            // Legacy linear density mixing.
            auto [eig, C] = linalg::diagonalize_fock_generalized(H, S);
            RealMatrix P_new = linalg::build_density_matrix(C, cfg.n_occ, nb, 2.0);
            if (iter > 0 && std::abs(e_total - e_prev) < cfg.tol) {
                res.converged = true;
                store_converged_state(res, H, S, cfg.n_occ, nb);
                break;
            }
            e_prev = e_total;
            P = (1.0 - cfg.mix) * P + cfg.mix * P_new;
        }
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
