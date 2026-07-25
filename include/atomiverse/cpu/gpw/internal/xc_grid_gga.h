/// @file cpu/gpw/internal/xc_grid_gga.h
/// @brief GGA (PBE / PBEsol / BLYP) exchange-correlation energy and potential
///        on the GPW grid.
///
/// The GGA counterpart of `xc_grid_pade` (LDA). Density gradients are computed
/// spectrally (`grid_fft`), the pointwise functional is reused from the
/// project's libxc-parity `cpu::scf::evaluate_xc`, and the gradient-dependent
/// part of the potential is assembled via the divergence term. Closed-shell
/// (spin-unpolarized): rho split evenly between spins.
#pragma once

#include <vector>

#include "atomiverse/cpu/gpw/internal/collocate.h" // UniformGrid
#include "atomiverse/cpu/gpw/internal/xc_grid.h"   // XcOnGrid
#include "atomiverse/cpu/scf/xc_functional.h"      // XCType

namespace ATV {
namespace cpu {
namespace gpw {

/// @brief Evaluate a pure GGA XC energy and potential on the grid.
///
/// E_xc = sum_g w rho eps_xc; the real-space potential is
///   v_xc = vrho_a - div[ (vsigma_aa + 1/2 vsigma_ab) grad rho ]
/// (the spin-unpolarized collapse of v_a = vrho_a - div[2 vsigma_aa grad rho_a
///  + vsigma_ab grad rho_b] with rho_a=rho_b=rho/2, grad rho_a=grad rho/2).
///
/// @param grid  Uniform grid.
/// @param rho_r Total electron density on the grid, length grid.npts().
/// @param xc    Pure GGA type (`PBE`, `PBEsol`, or `BLYP`).
/// @return XcOnGrid with the total energy and the real-space potential.
XcOnGrid xc_grid_gga(const UniformGrid& grid, const std::vector<double>& rho_r,
                     cpu::scf::XCType xc);

/// @brief PBE GGA on the grid (`xc_grid_gga` with `XCType::PBE`).
XcOnGrid xc_grid_pbe(const UniformGrid& grid, const std::vector<double>& rho_r);

/// @brief PBEsol GGA on the grid (`xc_grid_gga` with `XCType::PBEsol`).
XcOnGrid xc_grid_pbesol(const UniformGrid& grid,
                        const std::vector<double>& rho_r);

/// @brief BLYP GGA on the grid (`xc_grid_gga` with `XCType::BLYP`).
XcOnGrid xc_grid_blyp(const UniformGrid& grid, const std::vector<double>& rho_r);

} // namespace gpw
} // namespace cpu
} // namespace ATV
