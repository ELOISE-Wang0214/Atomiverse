/// @file cpu/path/dimer.h
/// @brief Henkelman–Jónsson dimer method skeleton (CPU).
///
/// Finds a first-order saddle by maintaining a finite-difference dimer
/// (two images separated by δ along unit mode N), rotating N to the lowest
/// curvature direction, then translating the midpoint uphill along N and
/// downhill perpendicular to N.
///
/// Breadth-first skeleton: single-rank, Cartesian, user-supplied force via
/// ProgramBuilder (or built-in inverted-harmonic test potential). Not a new
/// Driver category — a method preset like Remd/Pimd.
///
/// Refs: Henkelman & Jónsson, J. Chem. Phys. 111, 7010 (1999),
///       DOI 10.1063/1.480097.
#pragma once

#include <cstdint>
#include <functional>
#include <string>
#include <vector>

#include "atomiverse/cpu/frame/frame.h"

namespace ATV {
namespace cpu {
namespace path {

/// @brief Configuration for one dimer search.
struct DimerConfig {
    double delta = 0.01;     ///< Half-separation δ [Å] (dimer length = 2δ).
    double dR = 0.05;        ///< Translation step [Å] along modified force.
    int max_rotations = 8;   ///< Rotation iterations per translation step.
    double rot_tol = 1e-3;   ///< Stop rotating when |F_perp| < rot_tol [eV/Å].
    double f_tol = 0.05;     ///< Converge when |F_mod| < f_tol [eV/Å].
    std::uint64_t seed = 1;  ///< RNG seed for initial mode N.
};

/// @brief Dimer saddle search over one Frame (midpoint = saddle estimate).
class Dimer {
  public:
    explicit Dimer(DimerConfig cfg = {});

    /// @brief Force evaluator: write forces [eV/Å] and energy [eV] on frame.
    using ForceFn = std::function<void(cpu::Frame&)>;

    /// @brief Install the force callback (required before run unless test pot).
    void set_force_fn(ForceFn fn);

    /// @brief Use built-in saddle V = ½ k_x x² − ½ k_y y² (z free free particle).
    /// @param kx Curvature along x [eV/Å²] (> 0, bound).
    /// @param ky Curvature along y [eV/Å²] (> 0, inverted).
    void use_inverted_harmonic(double kx = 1.0, double ky = 1.0);

    /// @brief Run dimer search; overwrites frame with the saddle estimate.
    /// @param frame Starting geometry (and masses/types/box).
    /// @param max_iter Maximum translation steps.
    /// @return Number of translation steps performed.
    int run(cpu::Frame& frame, int max_iter);

    /// @brief True if last run met f_tol.
    bool converged() const { return converged_; }

    /// @brief Final |F_mod| [eV/Å].
    double max_force() const { return f_mod_; }

  private:
    void eval_force_(cpu::Frame& f);
    void normalize3_(double& x, double& y, double& z) const;
    void random_mode_(std::uint64_t seed);

    DimerConfig cfg_;
    ForceFn force_fn_;
    bool use_test_pot_ = false;
    double kx_ = 1.0, ky_ = 1.0;
    double nx_ = 1, ny_ = 0, nz_ = 0; ///< Unit dimer mode.
    double f_mod_ = 0;
    bool converged_ = false;
};

} // namespace path
} // namespace cpu
} // namespace ATV
