/// @file cpu/path/irc.h
/// @brief Intrinsic reaction coordinate (Fukui) skeleton — CPU.
///
/// Mass-weighted steepest descent from a saddle in both directions.
/// Ref: Fukui, Acc. Chem. Res. 14, 363 (1981).
#pragma once

#include <functional>
#include <vector>

#include "atomiverse/cpu/frame/frame.h"

namespace ATV {
namespace cpu {
namespace path {

/// @brief IRC configuration.
struct IrcConfig {
    double step = 0.01;     ///< SD step [coord units].
    int max_steps = 500;    ///< Max steps per branch.
    double f_tol = 1e-3;    ///< Stop when |F| < f_tol.
    double kick = 0.05;     ///< Initial displacement along mode [coord].
};

/// @brief Fukui IRC from a first-order saddle (1-particle / few-atom Frames).
class Irc {
  public:
    explicit Irc(IrcConfig cfg = {});

    using ForceFn = std::function<void(cpu::Frame&)>;
    void set_force_fn(ForceFn fn);
    void use_muller_brown();

    /// @brief Run both branches from `saddle`; fills forward/backward paths.
    /// @param saddle Starting geometry (near first-order saddle).
    /// @param mode_x Optional initial mode x-component (unit free); if 0, use F.
    /// @param mode_y Optional initial mode y-component.
    void run(const cpu::Frame& saddle, double mode_x = 0.0,
             double mode_y = 1.0);

    double forward_final_energy() const { return E_fwd_; }
    double backward_final_energy() const { return E_bwd_; }
    const std::vector<double>& forward_energies() const { return path_E_fwd_; }
    const std::vector<double>& backward_energies() const { return path_E_bwd_; }

  private:
    void branch_(cpu::Frame f, double mx, double my, bool forward);

    IrcConfig cfg_;
    ForceFn force_;
    std::vector<double> path_E_fwd_, path_E_bwd_;
    double E_fwd_ = 0, E_bwd_ = 0;
};

} // namespace path
} // namespace cpu
} // namespace ATV
