/// @file cpu/planewave/dfpt_phonon.h
/// @brief Γ-point DFPT phonons for the plane-wave method (CPU).
///
/// Analytic dynamical matrix from density-functional perturbation theory
/// (self-consistent Sternheimer response) instead of `6N` finite-displacement
/// SCF solves. Restricted LSDA with a UPF norm-conserving pseudopotential and
/// fixed occupations; q = 0 only (full dispersions are a follow-on).
/// Parameters are members — no public `*Config` bag (`Bands` precedent).
///
/// Reference: Baroni, de Gironcoli, Dal Corso, Giannozzi,
/// *Rev. Mod. Phys.* **73**, 515 (2001).
///
/// @code
///   cpu::pw::PW pw;
///   pw.upf_path = "Si.upf";
///   pw.ecut_ry = 15.0;
///   pw.k_mesh = 2;
///   pw.xc = cpu::scf::XCType::LSDA;
///
///   cpu::pw::DfptPhonon dfpt;
///   auto r = dfpt.compute(pw, frame);   // ω [cm⁻¹] at Γ
/// @endcode
#pragma once

#include <xtensor/containers/xtensor.hpp>

#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/planewave/pw_solver.h"

namespace ATV {
namespace cpu {
namespace pw {

/// @brief Result of a Γ-point DFPT phonon calculation.
///
/// Frame purity: spectra and force constants live here, never on Frame.
struct DfptPhononResult {
    /// @brief Wavenumbers [cm⁻¹], length `3N`, ascending; imaginary modes
    ///        (λ < 0) are reported as negative values.
    xt::xtensor<double, 1> frequencies;

    /// @brief Symmetrized Cartesian force-constant matrix [eV/Å²],
    ///        shape `[3N, 3N]`, index `3·I + α`.
    xt::xtensor<double, 2> force_constants;

    /// @brief Whether the ground-state SCF converged.
    bool scf_converged = false;

    /// @brief Whether every one of the `3N` response SCFs converged.
    bool response_converged = false;

    /// @brief Max |C − Cᵀ|/2 before symmetrization [eV/Å²] (diagnostic).
    double max_asymmetry = 0;

    /// @brief Atoms in the cell.
    int n_atoms = 0;
};

/// @brief Γ-point DFPT phonon calculator over a `PW` ground state.
///
/// `compute` runs one ground-state SCF, then `3N` self-consistent
/// Sternheimer responses, and diagonalizes the mass-weighted dynamical
/// matrix (masses from the IUPAC table by atomic number, the same
/// convention as `cpu::vib`). All workspace is private and rebuilt per
/// call; the input Frame is not mutated.
class DfptPhonon {
  public:
    /// @brief Response-density L2 convergence per perturbation
    ///        [electrons/Bohr³ per Bohr, volume-weighted grid norm].
    double response_tol = 1e-8;
    /// @brief Response SCF iteration cap per perturbation (dimensionless).
    int response_max_iter = 80;
    /// @brief Sternheimer CG relative residual tolerance (dimensionless).
    double cg_tol = 1e-10;
    /// @brief Sternheimer CG iteration cap (dimensionless).
    int cg_max_iter = 200;
    /// @brief Linear mixing factor for the response density in (0, 1].
    double mix_beta = 0.5;

    /// @brief DFPT phonons at Γ for the geometry in `frame`.
    ///
    /// Throws `std::invalid_argument` unless `method.xc == LSDA`,
    /// `method.upf_path` is set, smearing is off, and spin is restricted.
    ///
    /// @param method Plane-wave method (members only; not mutated).
    /// @param frame  Atom positions [Å], types (Z), simulation cell.
    /// @return Owned result with ω [cm⁻¹] and force constants [eV/Å²].
    DfptPhononResult compute(const PW& method, const cpu::Frame& frame) const;
};

} // namespace pw
} // namespace cpu
} // namespace ATV
