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
    }
};

} // namespace vib
} // namespace cpu
} // namespace ATV
