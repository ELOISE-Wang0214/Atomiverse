#include <chrono>
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
#include "atomiverse/cpu/gpw/internal/dftd3.h"
#include "atomiverse/cpu/gpw/internal/gth_integrals.h"
#include "atomiverse/cpu/gpw/internal/hartree_grid.h"
#include "atomiverse/cpu/gpw/internal/mt_poisson.h"
#include "atomiverse/cpu/gpw/internal/multigrid.h"
#include "atomiverse/cpu/gpw/internal/one_electron.h"
#include "atomiverse/cpu/gpw/internal/ot.h"
#include "atomiverse/cpu/gpw/internal/periodic.h"
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
// Γ-point periodic GTH pseudopotential matrix element. The basis is diffuse
// (reach > L/2 in a thin cell), so the bra φ_ν is summed over its lattice images
// T_ν (outer loop, screened by the Gaussian-product prefactor K_pq). The GTH
// potential itself is short-ranged (rloc / projector radii ~0.3-0.5 Bohr ≪ L/2),
// so for each atom only its single nearest periodic image to the product center P
// contributes — taken by min-image, no atom-image loop. `m` (from image_counts)
// must match the basis-collocation / overlap image set for consistency.
double gth_matrix_element(const GridBasisFunction& mu, const GridBasisFunction& nu,
                          const std::vector<GpwAtom>& atoms, const double L[3],
                          const int m[3], double reach) {
    if (mu.empty() || nu.empty()) return 0.0;
    const double* muc = mu[0].center;
    const double* nuc = nu[0].center;
    const double reach2 = reach * reach;
    double v = 0.0;
    for (const auto& p : mu) {
        const Prim pp = to_prim(p);
        for (const auto& q0 : nu) {
            const double cc = p.coeff * q0.coeff;
            if (cc == 0.0) continue;
            const double acc = std::abs(cc);
            const double a = p.alpha, b = q0.alpha, g = a + b;
            for (int tx = -m[0]; tx <= m[0]; ++tx)
                for (int ty = -m[1]; ty <= m[1]; ++ty)
                    for (int tz = -m[2]; tz <= m[2]; ++tz) {
                        // Skip ν-images whose shifted center cannot overlap μ:
                        // no μ-ν overlap ⇒ no 3-centre PP contribution.
                        const double sx = muc[0] - nuc[0] - tx * L[0];
                        const double sy = muc[1] - nuc[1] - ty * L[1];
                        const double sz = muc[2] - nuc[2] - tz * L[2];
                        if (sx * sx + sy * sy + sz * sz > reach2) continue;
                        GridPrimitive q = q0; // image T_ν of the bra primitive
                        q.center[0] += tx * L[0];
                        q.center[1] += ty * L[1];
                        q.center[2] += tz * L[2];
                        double dpq2 = 0.0, P[3];
                        for (int d = 0; d < 3; ++d) {
                            const double diff = p.center[d] - q.center[d];
                            dpq2 += diff * diff;
                            P[d] = (a * p.center[d] + b * q.center[d]) / g;
                        }
                        const double Kpq = std::exp(-a * b / g * dpq2);
                        if (acc * Kpq < kBuildScreen) continue; // image too far
                        const Prim qq = to_prim(q);
                        for (const auto& at : atoms) {
                            double rmax = at.rloc;
                            for (const auto& ch : at.nonlocal) rmax = std::max(rmax, ch.r);
                            const double kappa = 1.0 / (2.0 * rmax * rmax);
                            // nearest periodic image of this atom to the product center
                            double apos[3], dPA2 = 0.0;
                            for (int d = 0; d < 3; ++d) {
                                apos[d] = at.pos[d];
                                if (L[d] > 0.0 && (m[0] || m[1] || m[2]))
                                    apos[d] -= L[d] * std::nearbyint((at.pos[d] - P[d]) / L[d]);
                                const double e = P[d] - apos[d];
                                dPA2 += e * e;
                            }
                            if (acc * Kpq * std::exp(-(g * kappa / (g + kappa)) * dPA2) <
                                kBuildScreen)
                                continue; // atom too far for any PP contribution
                            if (!at.cloc.empty())
                                v += cc * gth_local_sr(pp, qq, apos, at.rloc, at.cloc);
                            for (const auto& ch : at.nonlocal)
                                v += cc * gth_nonlocal(pp, qq, apos, ch);
                        }
                    }
        }
    }
    return v;
}

// Max Gaussian-product prefactor over all primitive pairs of (mu,nu): bounds the
// magnitude of S, T, and Vpp for this pair. Below kBuildScreen the whole (i,j)
// element is negligible and can be skipped (overlap pair screening).
double pair_overlap_bound(const GridBasisFunction& mu, const GridBasisFunction& nu,
                          const double L[3], bool periodic) {
    double best = 0.0;
    for (const auto& p : mu)
        for (const auto& q : nu) {
            const double a = p.alpha, b = q.alpha;
            double d2 = 0.0;
            for (int k = 0; k < 3; ++k) {
                double e = p.center[k] - q.center[k];
                // Min-image: a pair straddling a periodic boundary is a true near
                // neighbour through its nearest image — keep it past the screen.
                if (periodic && L[k] > 0.0) e -= L[k] * std::nearbyint(e / L[k]);
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
                      const GthPotential& pot, int z) {
    GpwAtom atom;
    atom.z_val = static_cast<double>(pot.zion);
    atom.z = z;
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

    // Γ-point lattice-image extent, shared by the analytic S/T/Vpp lattice sums
    // and the grid collocation. A Gaussian *pair* reaches up to twice the largest
    // single-function radius; images within that distance contribute. EPS_DEFAULT
    // sets the radius (CP2K's Gaussian-tail tolerance).
    const bool periodic = !cfg.nonperiodic;
    const double eps_default = cfg.eps_default;
    std::vector<double> radii(static_cast<std::size_t>(nb));
    double rmax = 0.0;
    for (int i = 0; i < nb; ++i) {
        radii[static_cast<std::size_t>(i)] = basis_function_radius(basis[i], eps_default);
        rmax = std::max(rmax, radii[static_cast<std::size_t>(i)]);
    }
    int mimg[3];
    image_counts(cfg.cell, periodic, 2.0 * rmax, mimg);

    // Grimme DFT-D3 (zero damping, PBE) dispersion — a geometry-only additive
    // correction (independent of the SCF), folded into the final total energy.
    double e_d3 = 0.0;
    if (cfg.d3) {
        std::vector<int> Z(atoms.size());
        std::vector<std::array<double, 3>> pos(Z.size());
        for (std::size_t i = 0; i < Z.size(); ++i) {
            Z[i] = atoms[i].z;
            pos[i] = {atoms[i].pos[0], atoms[i].pos[1], atoms[i].pos[2]};
        }
        e_d3 = dftd3_pbe(Z, pos, cfg.cell, periodic, cfg.d3_param_file).energy;
    }

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
    const bool sprof = std::getenv("ATV_GPW_PROF") != nullptr;
    auto sclk = std::chrono::steady_clock::now();
    auto smark = [&](const char* what) {
        if (!sprof) return;
        auto now = std::chrono::steady_clock::now();
        std::fprintf(stderr, "[gpw setup] %-14s %.0f ms\n", what,
                     std::chrono::duration<double, std::milli>(now - sclk).count());
        sclk = now;
    };
#pragma omp parallel for schedule(dynamic, 8)
    for (int i = 0; i < nb; ++i) {
        for (int j = 0; j < nb; ++j) {
            // Overlap pair screening: S, T, Vpp are all bounded by the Gaussian
            // product magnitude of (i,j). If even the best primitive pair is
            // below the screen, the whole element is negligible -> leave 0.
            if (pair_overlap_bound(basis[i], basis[j], cfg.cell, periodic) < kBuildScreen)
                continue;
            const double reach = radii[static_cast<std::size_t>(i)] +
                                 radii[static_cast<std::size_t>(j)];
            S(i, j) = overlap(basis[i], basis[j], cfg.cell, mimg, reach);
            T(i, j) = kinetic(basis[i], basis[j], cfg.cell, mimg, reach);
            Vpp(i, j) = gth_matrix_element(basis[i], basis[j], atoms, cfg.cell, mimg, reach);
        }
    }
    const RealMatrix H0 = T + Vpp;
    smark("STVpp_build");

    // Geometry-fixed ionic density and the constant self-energy to remove.
    const std::vector<double> n_ion = ionic_density(grid, atoms, !cfg.nonperiodic);
    double e_self = 0.0;
    for (const auto& a : atoms) e_self += core_self_energy(a);
    smark("ionic_density");

    // Core (H0) initial guess.
    auto [eig0, C0] = linalg::diagonalize_fock_generalized(H0, S);
    RealMatrix P = linalg::build_density_matrix(C0, cfg.n_occ, nb, 2.0);
    smark("init_eigensolve");

    // Basis-on-grid is geometry-fixed: evaluate the Gaussians once and reuse
    // across all SCF iterations. The molecule fills a small part of the box, so
    // also precompute the "active" grid points (where any basis function is
    // non-negligible) and collocate/integrate only over those.
    // Compact-support basis-on-grid: each Gaussian touches only the grid points
    // within its radius (CSR-by-point), built once. Replaces the dense
    // npts x nb matrix — prohibitive in memory at scale (ngrid*nb*8 bytes) and
    // O(na*nb^2) per iter; the compact path is O(sum_g k_g^2), k_g << nb.
    const CompactBasisGrid cbg =
        build_compact_basis_grid(grid, basis, periodic, eps_default);
    smark("cbg_build");

    // CP2K-style REL_CUTOFF multigrid (gated on ngrids>1 && rel_cutoff>0). Each
    // primitive-pair product of the density is assigned (by its combined exponent)
    // to the coarsest grid level that resolves it; diffuse products collocate on a
    // coarse grid (far fewer points). The per-level densities are prolonged to the
    // finest grid for the (nonlinear) Hartree+XC, and the fine potential is
    // restricted back to each level to integrate the KS matrix. ngrids==1 leaves
    // the single-grid path above untouched (bit-identical). The product list is
    // geometry-fixed (excludes P); P is applied per SCF iteration.
    const bool use_mg = cfg.ngrids > 1 && cfg.rel_cutoff > 0.0;
    MultiGrid mg;
    ProductSet products;
    // Inter-level transfer is spectral (band-limited FFT interpolation): a diffuse
    // density collocated on a coarse level is embedded into the fine grid's Fourier
    // spectrum at the same physical frequencies and zero-padded above the coarse
    // Nyquist (shape-exact, unlike trilinear). prolong_fft conserves mass natively
    // (the G=0 mode is preserved) and restrict_fft is its exact adjoint, so the KS
    // integration stays the exact gradient of the collocated energy and the SCF is
    // variational — no per-level mass-conservation correction is needed.
    if (use_mg) {
        mg = make_multigrid(cfg.cell, cfg.cutoff_ry, cfg.ngrids, cfg.rel_cutoff);
        products = build_products(mg, basis, cfg.cell, cfg.rel_cutoff, eps_default,
                                  periodic);
    }

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
    const bool prof = std::getenv("ATV_GPW_PROF") != nullptr;
    using ProfClk = std::chrono::steady_clock;
    auto pms = [](ProfClk::time_point a, ProfClk::time_point b) {
        return std::chrono::duration<double, std::milli>(b - a).count();
    };
    auto build_H_and_energy = [&](const RealMatrix& Pin) -> std::pair<RealMatrix, double> {
        auto t0p = ProfClk::now();
        const std::vector<double> Pvec = mat_to_vec(Pin, nb);
        // Density on the finest grid. Single-grid: direct compact collocation.
        // Multigrid: collocate each level's products, prolong to the fine grid.
        std::vector<double> rho;
        if (use_mg) {
            rho.assign(static_cast<std::size_t>(mg.levels[0].npts()), 0.0);
            for (std::size_t l = 0; l < mg.levels.size(); ++l) {
                std::vector<double> rho_l = collocate_products(
                    mg.levels[l], products.by_level[l], nb, Pvec, eps_default);
                if (l == 0) {
                    rho = std::move(rho_l); // finest level is the reference grid
                } else {
                    // Spectral (band-limited) prolongation to the fine grid:
                    // shape-exact and mass-conserving (no D_l correction).
                    const std::vector<double> up = prolong_fft(
                        mg.levels[l], mg.levels[0], rho_l, periodic);
                    for (std::size_t g = 0; g < rho.size(); ++g) rho[g] += up[g];
                }
            }
        } else {
            rho = collocate_density(cbg, nb, Pvec);
        }
        auto t1p = ProfClk::now();
        if (use_mg && std::getenv("ATV_GPW_MG_TRACE")) {
            double qe = 0.0;
            for (double v : rho) qe += v;
            qe *= grid.weight();
            std::fprintf(stderr, "[gpw mg] int rho_fine = %.10f (Tr(P S) check)\n", qe);
        }
        std::vector<double> n_tot(rho.size());
        for (std::size_t g = 0; g < rho.size(); ++g) n_tot[g] = rho[g] - n_ion[g];
        const PoissonResult es =
            cfg.nonperiodic ? solve_poisson_mt(grid, n_tot) : solve_poisson(grid, n_tot);
        auto t2p = ProfClk::now();
        const XcOnGrid xc =
            cfg.xc == XcKind::PBE ? xc_grid_pbe(grid, rho) : xc_grid_pade(grid, rho);
        auto t3p = ProfClk::now();
        std::vector<double> veff(rho.size());
        for (std::size_t g = 0; g < rho.size(); ++g) veff[g] = es.v_r[g] + xc.vxc_r[g];
        // KS grid matrix. Single-grid: integrate veff directly. Multigrid:
        // restrict the fine potential to each level (the transpose of prolong) and
        // integrate each level's products — the adjoint of the density assembly,
        // so the KS matrix stays the exact gradient of the collocated energy.
        RealMatrix Vgrid;
        if (use_mg) {
            std::vector<double> V_ks(static_cast<std::size_t>(nb) * nb, 0.0);
            for (std::size_t l = 0; l < mg.levels.size(); ++l) {
                std::vector<double> V_l;
                if (l == 0) {
                    V_l = veff;
                } else {
                    // Exact adjoint of the spectral prolongation used for the
                    // density: restrict the fine potential by spectral truncation,
                    // so the KS matrix is the exact gradient of the collocated
                    // energy (variational SCF).
                    V_l = restrict_fft(mg.levels[0], mg.levels[l], veff, periodic);
                }
                const std::vector<double> Vl_ks = integrate_products(
                    mg.levels[l], products.by_level[l], nb, V_l, eps_default);
                for (std::size_t kk = 0; kk < V_ks.size(); ++kk) V_ks[kk] += Vl_ks[kk];
            }
            Vgrid = vec_to_mat(V_ks, nb);
        } else {
            Vgrid = vec_to_mat(integrate_potential(cbg, nb, veff), nb);
        }
        auto t4p = ProfClk::now();
        if (prof)
            std::fprintf(stderr,
                         "[gpw prof] collocate=%.0f poisson=%.0f xc=%.0f integrate=%.0f ms "
                         "(ngrid=%ld)\n",
                         pms(t0p, t1p), pms(t1p, t2p), pms(t2p, t3p), pms(t3p, t4p),
                         grid.npts());
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
        const bool prof = std::getenv("ATV_GPW_OT_PROF") != nullptr;
        using Clk = std::chrono::steady_clock;
        auto ms = [](Clk::time_point a, Clk::time_point b) {
            return std::chrono::duration<double, std::milli>(b - a).count();
        };
        for (int iter = 0; iter < cfg.max_iter; ++iter) {
            auto tA = Clk::now();
            P = linalg::build_density_matrix(Cocc, cfg.n_occ, nb, 2.0); // 2 Cocc Coccᵀ
            auto [H, e_total] = build_H_and_energy(P);
            auto tB = Clk::now();
            res.iters = iter + 1;
            res.e_total = e_total;
            int ls_rebuilds = 0;

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
            auto tC = Clk::now();
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
                ++ls_rebuilds;
                if (etot(Cnew) <= e_total + c1 * alpha * slope) break;
                alpha *= 0.5;
                Cnew = internal::ot_parameterize(Cocc, alpha * d, S);
            }
            Cocc = internal::ot_orthonormalize(Cnew, S);
            if (prof) {
                auto tD = Clk::now();
                std::fprintf(stderr,
                             "[OT prof %2d] H_grad=%.0fms precond=%.0fms "
                             "linesearch=%.0fms (%d rebuilds) | total=%.0fms\n",
                             iter, ms(tA, tB), ms(tB, tC), ms(tC, tD), ls_rebuilds,
                             ms(tA, tD));
            }
        }
        res.e_total += e_d3;
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
    res.e_total += e_d3;
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
