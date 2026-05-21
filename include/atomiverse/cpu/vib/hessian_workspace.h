/// @file cpu/vib/hessian_workspace.h
/// @brief cpu::vib::HessianWorkspace — shared scratch for the four
///        numerical-Hessian structural components.
///
/// The four Response-category components (`ReferenceGeometry`,
/// `DisplaceStep`, `HessianAccumulator`, `Symmetrize`) cooperate on one
/// `[3N,3N]` matrix plus the equilibrium-geometry cache. The driver
/// (`cpu::vib::Hessian`) owns one `HessianWorkspace` and hands each
/// component a non-owning pointer at construction — exactly as the
/// driver owns its component set today. No state lives on `cpu::Frame`.
#pragma once

#include <vector>

#include <xtensor/containers/xtensor.hpp>

namespace ATV {
namespace cpu {
namespace vib {

/// @brief Finite-difference step size in Angstrom.
///
/// Applied directly in the frame's Angstrom coordinates. The equivalent
/// in Bohr is `1.889726e-3` (non-binding note); the binding value is
/// `1.0e-3 Angstrom`. Central differencing has `O(h^2)` truncation
/// error with leading coefficient `g'''/6`.
inline constexpr double kHessianStepAngstrom = 1.0e-3;

/// @brief Shared scratch for the numerical-Hessian components.
///
/// Sized once by the driver in `build()` from `N = frame.atoms.size()`.
/// Lifetime: owned by `cpu::vib::Hessian`; components hold raw pointers.
struct HessianWorkspace {
    /// @brief Number of atoms `N`.
    int n_atoms = 0;

    /// @brief Equilibrium Cartesian coordinates `R0`, flattened length `3N`
    ///        (atom-major: `[x0,y0,z0, x1,y1,z1, ...]`), in Angstrom.
    std::vector<double> r0;

    /// @brief True once `ReferenceGeometry` has cached `R0`.
    bool r0_ready = false;

    /// @brief Dense Hessian under construction, shape `[3N,3N]`, eV/Angstrom^2.
    xt::xtensor<double, 2> hessian;

    /// @brief Force vector from the most recent `+h` evaluation of the
    ///        current coordinate, flattened length `3N` (eV/Angstrom).
    ///        `HessianAccumulator` stores the `+h` forces here, then pairs
    ///        them with the next `-h` forces to write one Hessian row.
    std::vector<double> f_plus;

    // ── Vibrational-frequency link (phase4-vibrations-02) ──────────────
    //
    // The `MassWeighting` and `EigenSolver` components cooperate on the
    // fields below during the finalize pass, after `Symmetrize`.

    /// @brief Mass-weighted Hessian `H^mw` in atomic units
    ///        (Hartree/Bohr^2/amu), shape `[3N,3N]`. `MassWeighting`
    ///        derives it from `hessian` without mutating `hessian` — so
    ///        `HessianResult.hessian` keeps the raw link-1 eV/Angstrom^2
    ///        matrix. `EigenSolver` consumes this buffer.
    xt::xtensor<double, 2> hessian_mw;

    /// @brief True once `MassWeighting` has populated `hessian_mw`.
    ///        Guards `EigenSolver` against running on an empty buffer.
    bool mass_weighted = false;

    /// @brief Per-atom masses [amu] used for mass-weighting, length `N`.
    ///        Populated by `MassWeighting` from `mass_amu(frame type)`.
    std::vector<double> mass_amu;

    /// @brief Vibrational wavenumbers [cm^-1] produced by `EigenSolver`,
    ///        length `3N-6` / `3N-5`. Negative for imaginary modes.
    xt::xtensor<double, 1> frequencies;

    /// @brief Normal-mode eigenvectors `[3N, n_modes]` produced by
    ///        `EigenSolver` (columns lifted to the full Cartesian basis).
    xt::xtensor<double, 2> normal_modes;

    /// @brief Imaginary-mode count produced by `EigenSolver`.
    int n_imaginary = 0;

    /// @brief Linear-vs-nonlinear classification produced by `EigenSolver`.
    bool is_linear = false;

    /// @brief Allocate all buffers for an `N`-atom system. Called once
    ///        from the driver's `build()`; no allocation occurs in any
    ///        component `launch()`.
    /// @param n Number of atoms.
    void allocate(int n) {
        n_atoms = n;
        const int dim = 3 * n;
        r0.assign(dim, 0.0);
        r0_ready = false;
        hessian = xt::zeros<double>({dim, dim});
        f_plus.assign(dim, 0.0);
        hessian_mw = xt::zeros<double>({dim, dim});
        mass_weighted = false;
        mass_amu.assign(n, 0.0);
        frequencies = xt::zeros<double>({0});
        normal_modes = xt::zeros<double>({dim, 0});
        n_imaginary = 0;
        is_linear = false;
    }
};

} // namespace vib
} // namespace cpu
} // namespace ATV
