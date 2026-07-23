/// @file cpu/path/string_method.h
/// @brief Simplified string method (E–Ren–Vanden-Eijnden) skeleton (CPU).
///
/// Fixed-N image chain; evolve perpendicular to the path then reparameterize
/// by equal arc length. Not a Driver category.
///
/// Ref: E, Ren, Vanden-Eijnden, Phys. Rev. B 66, 052301 (2002).
#pragma once

#include <functional>
#include <vector>

#include "atomiverse/cpu/frame/frame.h"

namespace ATV {
namespace cpu {
namespace path {

/// @brief Configuration for one simplified string run.
struct StringConfig {
    int n_images = 12;     ///< Interior + endpoints count (≥ 3).
    double step = 0.01;    ///< Steepest-descent step [coord units].
    int max_iter = 200;    ///< Outer iterations.
    double f_tol = 0.05;   ///< Stop when max |F_perp| < f_tol.
};

/// @brief Simplified string method over a chain of 1-particle Frames (2D/3D).
class StringMethod {
  public:
    explicit StringMethod(StringConfig cfg = {});

    using ForceFn = std::function<void(cpu::Frame&)>;

    void set_force_fn(ForceFn fn);

    /// @brief Built-in Müller–Brown force (same constants as NEB regression).
    void use_muller_brown();

    /// @brief Run string from endpoint geometries; returns max |F_⊥| at end.
    /// @param start Left endpoint Frame (copied).
    /// @param end Right endpoint Frame (copied).
    /// @return Number of iterations performed.
    int run(const cpu::Frame& start, const cpu::Frame& end);

    /// @brief Highest image energy after last run.
    double max_energy() const { return max_energy_; }

    /// @brief Energy of image `i` after last run.
    double energy(int i) const { return energies_.at(static_cast<size_t>(i)); }

    /// @brief Copy image `i` into `out`.
    void image(int i, cpu::Frame& out) const;

    int n_images() const { return cfg_.n_images; }
    bool converged() const { return converged_; }

  private:
    void reparameterize_();
    void evolve_once_();

    StringConfig cfg_;
    ForceFn force_;
    std::vector<cpu::Frame> images_;
    std::vector<double> energies_;
    double max_energy_ = 0;
    bool converged_ = false;
};

} // namespace path
} // namespace cpu
} // namespace ATV
