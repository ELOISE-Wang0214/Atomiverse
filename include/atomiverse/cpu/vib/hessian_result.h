/// @file cpu/vib/hessian_result.h
/// @brief cpu::vib::HessianResult — output of the numerical-Hessian driver.
///
/// Frame purity is preserved by NOT adding a Hessian field to
/// `cpu::Frame`. Instead `cpu::vib::Hessian::run()` returns this plain
/// struct by value (mirroring `OptimizeDriver::run()` returning `int`).
/// Vibrational-frequency and thermochemistry links extend this struct
/// with mass-weighted frequencies / normal modes / thermo quantities.
#pragma once

#include <string>

#include <xtensor/containers/xtensor.hpp>

namespace ATV {
namespace cpu {
namespace vib {

/// @brief Raw numerical (central-difference) molecular Hessian.
///
/// The Hessian is the matrix of second derivatives of the
/// Born-Oppenheimer energy with respect to nuclear Cartesian
/// coordinates, `H[i.alpha, j.beta] = d2E / dR_{i.alpha} dR_{j.beta}`.
/// It is assembled directly in the frame's working units, eV/Angstrom^2
/// (positions in Angstrom, forces in eV/Angstrom). Conversion to
/// Hartree/Bohr^2 is deferred to the mass-weighting / frequency link.
struct HessianResult {
    /// @brief Dense symmetric Hessian, shape `[3N, 3N]`, in eV/Angstrom^2.
    ///        Row/column index `c` maps to atom `c/3`, axis `c%3` (x,y,z).
    xt::xtensor<double, 2> hessian;

    /// @brief Number of atoms `N` (so the matrix is `3N x 3N`).
    int n_atoms = 0;

    /// @brief Physical unit tag of `hessian`. Always "eV/Angstrom^2"
    ///        for the numerical-Hessian link.
    std::string units = "eV/Angstrom^2";
};

} // namespace vib
} // namespace cpu
} // namespace ATV
