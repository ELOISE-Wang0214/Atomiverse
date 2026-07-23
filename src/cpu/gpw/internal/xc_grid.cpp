/// @file src/cpu/gpw/internal/xc_grid.cpp
/// @brief Implementation of XC energy/potential on the GPW uniform grid.

#include "atomiverse/cpu/gpw/internal/xc_grid.h"

#include "atomiverse/cpu/xc/pade.h"

namespace ATV {
namespace cpu {
namespace gpw {

XcOnGrid xc_grid_pade(const UniformGrid& grid,
                      const std::vector<double>& rho_r) {
    const long ng = grid.npts();
    const double w = grid.weight();
    XcOnGrid out;
    out.vxc_r.assign(static_cast<std::size_t>(ng), 0.0);
    double e = 0.0;
#pragma omp parallel for reduction(+ : e) schedule(static)
    for (long g = 0; g < ng; ++g) {
        const double rho = rho_r[g];
        // Spin-unpolarized: split evenly. pade_lsd applies its own low-density
        // cutoff, returning zeros below threshold.
        const ATV::cpu::xc::PadeXC xc =
            ATV::cpu::xc::pade_lsd(0.5 * rho, 0.5 * rho);
        e += w * rho * xc.exc;
        out.vxc_r[g] = xc.vrho_a; // = vrho_b in the unpolarized case
    }
    out.e_xc = e;
    return out;
}

XcOnGridSpin xc_grid_pade_lsd(const UniformGrid& grid,
                              const std::vector<double>& rho_a,
                              const std::vector<double>& rho_b) {
    const long ng = grid.npts();
    const double w = grid.weight();
    XcOnGridSpin out;
    out.vxc_a.assign(static_cast<std::size_t>(ng), 0.0);
    out.vxc_b.assign(static_cast<std::size_t>(ng), 0.0);
    double e = 0.0;
#pragma omp parallel for reduction(+ : e) schedule(static)
    for (long g = 0; g < ng; ++g) {
        const double ra = rho_a[static_cast<std::size_t>(g)];
        const double rb = rho_b[static_cast<std::size_t>(g)];
        const ATV::cpu::xc::PadeXC xc = ATV::cpu::xc::pade_lsd(ra, rb);
        e += w * (ra + rb) * xc.exc;
        out.vxc_a[static_cast<std::size_t>(g)] = xc.vrho_a;
        out.vxc_b[static_cast<std::size_t>(g)] = xc.vrho_b;
    }
    out.e_xc = e;
    return out;
}

} // namespace gpw
} // namespace cpu
} // namespace ATV
