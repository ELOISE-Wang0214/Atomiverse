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
    for (long g = 0; g < ng; ++g) {
        double r[3];
        grid.point(g, r);
        double* row = &phi[static_cast<std::size_t>(g) * nb];
        for (int u = 0; u < nb; ++u) row[u] = eval_basis_function(basis[u], r);
    }
    return phi;
}

std::vector<double> collocate_density(const UniformGrid& grid,
                                      const std::vector<GridBasisFunction>& basis,
                                      const std::vector<double>& P) {
    const long ng = grid.npts();
    const int nb = static_cast<int>(basis.size());
    const std::vector<double> phi = eval_basis_on_grid(grid, basis);
    std::vector<double> rho(static_cast<std::size_t>(ng), 0.0);
    std::vector<double> pv(nb); // P * phi_row
    for (long g = 0; g < ng; ++g) {
        const double* row = &phi[static_cast<std::size_t>(g) * nb];
        // pv_u = sum_v P_uv phi_v ; rho = sum_u phi_u pv_u
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

std::vector<double> integrate_potential(const UniformGrid& grid,
                                        const std::vector<GridBasisFunction>& basis,
                                        const std::vector<double>& V_r) {
    const long ng = grid.npts();
    const int nb = static_cast<int>(basis.size());
    const double w = grid.weight();
    const std::vector<double> phi = eval_basis_on_grid(grid, basis);
    std::vector<double> V(static_cast<std::size_t>(nb) * nb, 0.0);
    for (long g = 0; g < ng; ++g) {
        const double* row = &phi[static_cast<std::size_t>(g) * nb];
        const double vw = V_r[g] * w;
        for (int u = 0; u < nb; ++u) {
            const double fu = row[u] * vw;
            double* Vu = &V[static_cast<std::size_t>(u) * nb];
            for (int v = 0; v < nb; ++v) Vu[v] += fu * row[v];
        }
    }
    return V;
}

} // namespace gpw
} // namespace cpu
} // namespace ATV
