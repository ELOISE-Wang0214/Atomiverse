/// @file cpu/pw/internal/pw_workspace.h
/// @brief PWWorkspace — facade-private scratch for the plane-wave SCF.
///
/// Owned by `cpu::pw::PWSolver` (the inner solver), never by `cpu::Frame`
/// and never threaded through Build/Launch contexts. Mirrors the gaussian
/// SCF precedent: the workspace is the solver's private state, and the
/// step helpers (`internal::PoissonSolver`, `XCBuilder`, `VeffBuilder`,
/// `DavidsonSolver`, `DensityBuilder`, `ConvergenceCheck`) receive it as
/// an explicit parameter.
#pragma once

#include "atomiverse/core/real.h"
#include "atomiverse/cpu/pw/internal/fft.h"
#include "atomiverse/cpu/pw/internal/g_vector.h"
#include "atomiverse/cpu/pw/internal/kpoint.h"
#include <complex>
#include <memory>
#include <vector>
#include <xtensor/containers/xarray.hpp>

namespace ATV {
namespace cpu {
namespace pw {
namespace internal {

struct PWWorkspace {
    int n_g = 0;
    int nr = 0, nc = 0; // real and complex FFT grid sizes
    int fft_nx = 0, fft_ny = 0, fft_nz = 0;
    const GVectorList* gvec = nullptr; // current G-sphere (density or k)

    // Owned G-vector storage the `gvec` pointer aliases into:
    // the density sphere (Γ) plus one sphere per locally-owned k-point.
    GVectorList gvec0;
    KPointGrid kgrid;             // full Monkhorst-Pack mesh
    std::vector<int> owned_k;     // global k indices this rank computes
    std::vector<GVectorList> gvec_k; // spheres for owned_k (same order)

    // FFT plan for r2c / c2r transforms
    std::unique_ptr<FFTPlan3D> fft;

    xt::xarray<std::complex<Real>> rho_G;     // charge density in G-space
    xt::xarray<std::complex<Real>> rho_G_old; // previous iteration density
    xt::xarray<std::complex<Real>> V_H_G;     // Hartree potential
    xt::xarray<std::complex<Real>> V_xc_G;    // XC potential
    xt::xarray<std::complex<Real>> V_eff_G;   // total effective potential
    xt::xarray<Real> rho_r;                   // real-space density
    xt::xarray<Real> V_xc_r;                  // real-space XC potential
    xt::xarray<Real> T_G;                     // kinetic energy diagonal
    xt::xarray<Real> eigenvalues;             // per-band eigenvalues
    // Local PP in G-space, incl. structure factor Σ_a e^{-iG·R_a}.
    xt::xarray<std::complex<Real>> V_loc_G;
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

    // USPP augmentation: Q_ij matrix (n_proj × n_proj, row-major).
    // Parsed from the UPF and kept as fundamental pseudopotential data;
    // the augmentation density is not yet applied by the solver, so
    // E_aug stays 0 in the current implementation.
    std::vector<Real> q_ij;
    int n_proj_total = 0;

    // Normalized atomic form factor: rho_atom(G) / ∫ rho_atom(r) d³r
    xt::xarray<Real> aug_form_factor_G;

    // Atom positions [Bohr] for structure factors
    std::vector<Vec3> atom_pos_bohr;

    // Augmentation energy [Hartree] from current density
    Real E_aug = Real(0);

    // Pulay mixer history (managed by internal::DensityBuilder)
    std::vector<xt::xarray<std::complex<Real>>> pulay_res;
    std::vector<xt::xarray<std::complex<Real>>> pulay_rho;
    int pulay_depth = 5;
    Real mix_beta = Real(0.4);
};

} // namespace internal
} // namespace pw
} // namespace cpu
} // namespace ATV
