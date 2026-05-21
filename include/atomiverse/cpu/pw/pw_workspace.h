/// @file cpu/pw/pw_workspace.h
/// @brief Shared workspace for plane-wave SCF components.
#pragma once

#include "atomiverse/core/real.h"
#include "atomiverse/cpu/pw/internal/fft.h"
#include "atomiverse/cpu/pw/internal/g_vector.h"
#include <complex>
#include <memory>
#include <vector>
#include <xtensor/xarray.hpp>

namespace ATV {
namespace cpu {
namespace pw {

struct PWWorkspace {
    int n_g = 0;
    int nr = 0, nc = 0; // real and complex FFT grid sizes
    int fft_nx = 0, fft_ny = 0, fft_nz = 0;
    const internal::GVectorList* gvec = nullptr;

    // FFT plan for r2c / c2r transforms
    std::unique_ptr<internal::FFTPlan3D> fft;

    xt::xarray<std::complex<Real>> rho_G;     // charge density in G-space
    xt::xarray<std::complex<Real>> rho_G_old; // previous iteration density
    xt::xarray<std::complex<Real>> V_H_G;     // Hartree potential
    xt::xarray<std::complex<Real>> V_xc_G;    // XC potential
    xt::xarray<std::complex<Real>> V_eff_G;   // total effective potential
    xt::xarray<Real> rho_r;                   // real-space density
    xt::xarray<Real> V_xc_r;                  // real-space XC potential
    xt::xarray<Real> T_G;                     // kinetic energy diagonal
    xt::xarray<Real> eigenvalues;             // per-band eigenvalues
    xt::xarray<Real> V_loc_G_real;            // local PP in G-space (real part)
    xt::xarray<Real> V_eff_r; // V_eff on real-space grid [Hartree]
    xt::xarray<Real> exc_r;   // real-space XC energy density ε_xc(r)

    // Wavefunctions: wfc[b][ig] for band b, G-vector ig
    xt::xarray<std::complex<Real>> wfc;

    // Cell volume [Bohr³] for real-space integration
    Real omega = Real(0);

    // Accumulated XC energy [Hartree]
    Real E_xc = Real(0);

    // Non-local KB projector data in G-space
    struct KBProjG {
        int l;
        Real d;               // D_ij coefficient (diagonal approximation)
        xt::xarray<Real> f_G; // f_l(|G|) = ∫ r·β_upf(r)·j_l(Gr) dr
    };
    // kb_proj[atom][proj] — same f_G for all atoms (radial part only)
    std::vector<std::vector<KBProjG>> kb_proj;

    // Atomic wavefunctions for initial guess (QE atomic_wfc_mod.f90)
    struct AtomWfc {
        int l;
        xt::xarray<Real> chi_G; // (4π/√Ω) * ∫χ(r)j_l(Gr)rdr  for all G
    };
    std::vector<std::vector<AtomWfc>> atom_wfc; // [atom][iwfc]

    // USPP augmentation: Q_ij matrix (n_proj × n_proj, row-major)
    std::vector<Real> q_ij;
    int n_proj_total = 0;

    // Normalized atomic form factor: rho_atom(G) / ∫ rho_atom(r) d³r
    xt::xarray<Real> aug_form_factor_G;

    // Atom positions [Bohr] for structure factors
    std::vector<Vec3> atom_pos_bohr;

    // Augmentation energy [Hartree] from current density
    Real E_aug = Real(0);

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
