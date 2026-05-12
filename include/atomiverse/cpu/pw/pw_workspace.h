/// @file cpu/pw/pw_workspace.h
/// @brief Shared workspace for plane-wave SCF components.
#pragma once

#include <xtensor/xarray.hpp>
#include <vector>
#include <complex>
#include "atomiverse/core/real.h"
#include "atomiverse/cpu/pw/internal/g_vector.h"

namespace ATV {
namespace cpu {
namespace pw {

struct PWWorkspace {
    int n_g = 0;
    int nr = 0, nc = 0; // real and complex FFT grid sizes
    const internal::GVectorList* gvec = nullptr;

    xt::xarray<std::complex<Real>> rho_G;      // charge density in G-space
    xt::xarray<std::complex<Real>> rho_G_old;  // previous iteration density
    xt::xarray<std::complex<Real>> V_H_G;      // Hartree potential
    xt::xarray<std::complex<Real>> V_xc_G;     // XC potential
    xt::xarray<std::complex<Real>> V_eff_G;    // total effective potential
    xt::xarray<Real> rho_r;                    // real-space density
    xt::xarray<Real> V_xc_r;                   // real-space XC potential
    xt::xarray<Real> T_G;                      // kinetic energy diagonal
    xt::xarray<Real> eigenvalues;              // per k-point eigenvalues
    xt::xarray<Real> V_loc_G_real;             // local PP in G-space (real part)

    // Pulay mixer
    std::vector<xt::xarray<std::complex<Real>>> pulay_res;
    std::vector<xt::xarray<std::complex<Real>>> pulay_rho;
    int pulay_depth = 5;
    Real mix_beta = Real(0.4);

    bool converged = false;
    Real total_energy = 0;
};

} // namespace pw
} // namespace cpu
} // namespace ATV
