/// @file src/cpu/gpw/internal/multigrid.cpp
/// @brief Implementation of the CP2K-style CUTOFF/REL_CUTOFF multigrid helpers
///        (see multigrid.h). Level construction, exponent->level assignment, and
///        the adjoint trilinear prolong/restrict transfer operators.

#include "atomiverse/cpu/gpw/internal/multigrid.h"

#include <algorithm>
#include <cmath>

namespace ATV {
namespace cpu {
namespace gpw {

MultiGrid make_multigrid(const double cell_abc[3], double cutoff_ry, int ngrids,
                         double rel_cutoff) {
    MultiGrid mg;
    const int ng = std::max(1, ngrids);
    if (ng == 1 || rel_cutoff <= 0.0 || rel_cutoff >= cutoff_ry) {
        mg.levels.push_back(make_uniform_grid(cell_abc, cutoff_ry));
        mg.level_cutoff_ry.push_back(cutoff_ry);
        return mg;
    }
    // Geometric ladder cutoff_ry = c0 > c1 > ... > c_{ng-1} ~ rel_cutoff.
    const double ratio = std::pow(cutoff_ry / rel_cutoff, 1.0 / (ng - 1));
    double c = cutoff_ry;
    for (int l = 0; l < ng; ++l) {
        mg.levels.push_back(make_uniform_grid(cell_abc, c));
        mg.level_cutoff_ry.push_back(c);
        c /= ratio;
    }
    return mg;
}

std::vector<int> assign_levels(const MultiGrid& mg,
                               const std::vector<GridBasisFunction>& basis,
                               double rel_cutoff) {
    std::vector<int> level(basis.size(), 0);
    const int nlev = static_cast<int>(mg.levels.size());
    if (nlev <= 1) return level;
    // A primitive of exponent alpha has density exponent 2*alpha; the plane-wave
    // cutoff needed to resolve it scales linearly with alpha. Calibrate so the
    // reference exponent alpha_ref = rel_cutoff (in Ry, as an exponent proxy)
    // lands a Gaussian at the rel_cutoff level: need_cutoff = rel_cutoff*alpha.
    // A basis function goes on the COARSEST level whose cutoff still meets the
    // need of its SHARPEST primitive (so the function is never under-resolved).
    for (std::size_t u = 0; u < basis.size(); ++u) {
        double amax = 0.0;
        for (const auto& p : basis[u]) amax = std::max(amax, p.alpha);
        const double need = rel_cutoff * amax; // [Ry] resolution need
        // Coarsest level whose cutoff still resolves the Gaussian (most
        // economical). If none does (a very sharp Gaussian needs more than the
        // finest cutoff), clamp to the finest level (0).
        int lev = 0;
        for (int l = nlev - 1; l >= 0; --l) {
            if (mg.level_cutoff_ry[l] >= need) { lev = l; break; }
        }
        level[u] = lev;
    }
    return level;
}

namespace {
// Periodic wrap of an integer grid index into [0, n).
inline int wrap(int i, int n) {
    i %= n;
    return i < 0 ? i + n : i;
}
} // namespace

std::vector<double> prolong_trilinear(const UniformGrid& coarse,
                                      const UniformGrid& fine,
                                      const std::vector<double>& f_coarse,
                                      bool periodic) {
    std::vector<double> out(static_cast<std::size_t>(fine.npts()), 0.0);
    const int cn0 = coarse.n[0], cn1 = coarse.n[1], cn2 = coarse.n[2];
    const double cstep0 = coarse.L[0] / cn0, cstep1 = coarse.L[1] / cn1,
                 cstep2 = coarse.L[2] / cn2;
#pragma omp parallel for schedule(static)
    for (long g = 0; g < fine.npts(); ++g) {
        double r[3];
        fine.point(g, r);
        // Fractional position in coarse-grid index space.
        double fx = (r[0] - coarse.origin[0]) / cstep0;
        double fy = (r[1] - coarse.origin[1]) / cstep1;
        double fz = (r[2] - coarse.origin[2]) / cstep2;
        int i0 = static_cast<int>(std::floor(fx));
        int j0 = static_cast<int>(std::floor(fy));
        int k0 = static_cast<int>(std::floor(fz));
        double tx = fx - i0, ty = fy - j0, tz = fz - k0;
        double val = 0.0;
        for (int di = 0; di < 2; ++di)
            for (int dj = 0; dj < 2; ++dj)
                for (int dk = 0; dk < 2; ++dk) {
                    const double w = (di ? tx : 1 - tx) * (dj ? ty : 1 - ty) *
                                     (dk ? tz : 1 - tz);
                    int ci = periodic ? wrap(i0 + di, cn0)
                                      : std::clamp(i0 + di, 0, cn0 - 1);
                    int cj = periodic ? wrap(j0 + dj, cn1)
                                      : std::clamp(j0 + dj, 0, cn1 - 1);
                    int ck = periodic ? wrap(k0 + dk, cn2)
                                      : std::clamp(k0 + dk, 0, cn2 - 1);
                    const long ci_lin =
                        (static_cast<long>(ci) * cn1 + cj) * cn2 + ck;
                    val += w * f_coarse[static_cast<std::size_t>(ci_lin)];
                }
        out[static_cast<std::size_t>(g)] = val;
    }
    return out;
}

std::vector<double> restrict_transpose(const UniformGrid& fine,
                                       const UniformGrid& coarse,
                                       const std::vector<double>& f_fine,
                                       bool periodic) {
    // Exact transpose of prolong_trilinear, times the per-point volume ratio so
    // the operator is the adjoint w.r.t. the grid integration inner products:
    // sum_fine w_f V_f (P g)_f = sum_coarse w_c (R V)_c g_c.
    std::vector<double> out(static_cast<std::size_t>(coarse.npts()), 0.0);
    const int cn0 = coarse.n[0], cn1 = coarse.n[1], cn2 = coarse.n[2];
    const double cstep0 = coarse.L[0] / cn0, cstep1 = coarse.L[1] / cn1,
                 cstep2 = coarse.L[2] / cn2;
    const double vol_ratio = fine.weight() / coarse.weight();
    // Serial scatter (atomic adds would dominate); coarse grids are small.
    for (long g = 0; g < fine.npts(); ++g) {
        double r[3];
        fine.point(g, r);
        double fx = (r[0] - coarse.origin[0]) / cstep0;
        double fy = (r[1] - coarse.origin[1]) / cstep1;
        double fz = (r[2] - coarse.origin[2]) / cstep2;
        int i0 = static_cast<int>(std::floor(fx));
        int j0 = static_cast<int>(std::floor(fy));
        int k0 = static_cast<int>(std::floor(fz));
        double tx = fx - i0, ty = fy - j0, tz = fz - k0;
        const double vf = f_fine[static_cast<std::size_t>(g)] * vol_ratio;
        for (int di = 0; di < 2; ++di)
            for (int dj = 0; dj < 2; ++dj)
                for (int dk = 0; dk < 2; ++dk) {
                    const double w = (di ? tx : 1 - tx) * (dj ? ty : 1 - ty) *
                                     (dk ? tz : 1 - tz);
                    int ci = periodic ? wrap(i0 + di, cn0)
                                      : std::clamp(i0 + di, 0, cn0 - 1);
                    int cj = periodic ? wrap(j0 + dj, cn1)
                                      : std::clamp(j0 + dj, 0, cn1 - 1);
                    int ck = periodic ? wrap(k0 + dk, cn2)
                                      : std::clamp(k0 + dk, 0, cn2 - 1);
                    const long ci_lin =
                        (static_cast<long>(ci) * cn1 + cj) * cn2 + ck;
                    out[static_cast<std::size_t>(ci_lin)] += w * vf;
                }
    }
    return out;
}

} // namespace gpw
} // namespace cpu
} // namespace ATV
