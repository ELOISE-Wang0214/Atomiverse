/// @file src/cpu/gpw/internal/collocate.cpp
/// @brief Implementation of uniform-grid Gaussian collocation / integration.

#include "atomiverse/cpu/gpw/internal/collocate.h"

#include <cmath>

namespace ATV {
namespace cpu {
namespace gpw {

void UniformGrid::point(long idx, double r[3]) const {
    // k fastest, then j, then i (row-major i,j,k).
    const long nk = n[2];
    const long nj = n[1];
    const long k = idx % nk;
    const long j = (idx / nk) % nj;
    const long i = idx / (nk * nj);
    r[0] = origin[0] + (L[0] * static_cast<double>(i)) / n[0];
    r[1] = origin[1] + (L[1] * static_cast<double>(j)) / n[1];
    r[2] = origin[2] + (L[2] * static_cast<double>(k)) / n[2];
}

UniformGrid make_uniform_grid(const double cell_abc[3], double cutoff_ry) {
    UniformGrid g;
    const double gmax = std::sqrt(cutoff_ry); // |G_max| [1/Bohr], E_cut[Ry] = G^2
    for (int d = 0; d < 3; ++d) {
        g.L[d] = cell_abc[d];
        g.origin[d] = 0.0;
        // Nyquist: G_max = pi * n / L  ->  n = L * gmax / pi. Round up to even.
        int n = static_cast<int>(std::ceil(cell_abc[d] * gmax / M_PI));
        if (n < 2) n = 2;
        if (n % 2 != 0) ++n;
        g.n[d] = n;
    }
    return g;
}

namespace {

/// @brief Evaluate one primitive Gaussian at a point.
inline double eval_primitive(const GridPrimitive& p, const double r[3]) {
    const double dx = r[0] - p.center[0];
    const double dy = r[1] - p.center[1];
    const double dz = r[2] - p.center[2];
    const double r2 = dx * dx + dy * dy + dz * dz;
    double v = p.coeff * std::exp(-p.alpha * r2);
    for (int i = 0; i < p.l[0]; ++i) v *= dx;
    for (int i = 0; i < p.l[1]; ++i) v *= dy;
    for (int i = 0; i < p.l[2]; ++i) v *= dz;
    return v;
}

inline double eval_basis_function(const GridBasisFunction& bf, const double r[3]) {
    double v = 0.0;
    for (const auto& p : bf) v += eval_primitive(p, r);
    return v;
}

} // namespace

std::vector<double> eval_basis_on_grid(const UniformGrid& grid,
                                       const std::vector<GridBasisFunction>& basis) {
    const long ng = grid.npts();
    const int nb = static_cast<int>(basis.size());
    std::vector<double> phi(static_cast<std::size_t>(ng) * nb);
#pragma omp parallel for schedule(static)
    for (long g = 0; g < ng; ++g) {
        double r[3];
        grid.point(g, r);
        double* row = &phi[static_cast<std::size_t>(g) * nb];
        for (int u = 0; u < nb; ++u) row[u] = eval_basis_function(basis[u], r);
    }
    return phi;
}

std::vector<long> active_grid_points(const std::vector<double>& phi, long ng,
                                     int nb, double eps) {
    std::vector<long> active;
    for (long g = 0; g < ng; ++g) {
        const double* row = &phi[static_cast<std::size_t>(g) * nb];
        for (int u = 0; u < nb; ++u) {
            if (std::fabs(row[u]) > eps) {
                active.push_back(g);
                break;
            }
        }
    }
    return active;
}

std::vector<double> collocate_density(const UniformGrid& grid,
                                      const std::vector<double>& phi, int nb,
                                      const std::vector<double>& P,
                                      const std::vector<long>& active) {
    std::vector<double> rho(static_cast<std::size_t>(grid.npts()), 0.0);
    const long na = static_cast<long>(active.size());
#pragma omp parallel for schedule(static)
    for (long i = 0; i < na; ++i) {
        const long g = active[static_cast<std::size_t>(i)];
        const double* row = &phi[static_cast<std::size_t>(g) * nb];
        double acc = 0.0;
        for (int u = 0; u < nb; ++u) {
            double s = 0.0;
            const double* Pu = &P[static_cast<std::size_t>(u) * nb];
            for (int v = 0; v < nb; ++v) s += Pu[v] * row[v];
            acc += row[u] * s;
        }
        rho[g] = acc;
    }
    return rho;
}

std::vector<double> collocate_density(const UniformGrid& grid,
                                      const std::vector<double>& phi, int nb,
                                      const std::vector<double>& P) {
    return collocate_density(grid, phi, nb, P,
                             active_grid_points(phi, grid.npts(), nb));
}

std::vector<double> collocate_density(const UniformGrid& grid,
                                      const std::vector<GridBasisFunction>& basis,
                                      const std::vector<double>& P) {
    return collocate_density(grid, eval_basis_on_grid(grid, basis),
                             static_cast<int>(basis.size()), P);
}

std::vector<double> integrate_potential(const UniformGrid& grid,
                                        const std::vector<double>& phi, int nb,
                                        const std::vector<double>& V_r,
                                        const std::vector<long>& active) {
    const double w = grid.weight();
    const long na = static_cast<long>(active.size());
    const std::size_t nn = static_cast<std::size_t>(nb) * nb;
    std::vector<double> V(nn, 0.0);
    // Per-thread partial KS blocks, merged at the end (V is a shared reduction).
#pragma omp parallel
    {
        std::vector<double> Vloc(nn, 0.0);
#pragma omp for schedule(static) nowait
        for (long i = 0; i < na; ++i) {
            const long g = active[static_cast<std::size_t>(i)];
            const double* row = &phi[static_cast<std::size_t>(g) * nb];
            const double vw = V_r[g] * w;
            for (int u = 0; u < nb; ++u) {
                const double fu = row[u] * vw;
                double* Vu = &Vloc[static_cast<std::size_t>(u) * nb];
                for (int v = 0; v < nb; ++v) Vu[v] += fu * row[v];
            }
        }
#pragma omp critical
        for (std::size_t k = 0; k < nn; ++k) V[k] += Vloc[k];
    }
    return V;
}

std::vector<double> integrate_potential(const UniformGrid& grid,
                                        const std::vector<double>& phi, int nb,
                                        const std::vector<double>& V_r) {
    return integrate_potential(grid, phi, nb, V_r,
                               active_grid_points(phi, grid.npts(), nb));
}

std::vector<double> integrate_potential(const UniformGrid& grid,
                                        const std::vector<GridBasisFunction>& basis,
                                        const std::vector<double>& V_r) {
    return integrate_potential(grid, eval_basis_on_grid(grid, basis),
                               static_cast<int>(basis.size()), V_r);
}

} // namespace gpw
} // namespace cpu
} // namespace ATV
