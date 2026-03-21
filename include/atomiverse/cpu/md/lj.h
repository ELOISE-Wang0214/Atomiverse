/// @file cpu/md/lj.h
/// @brief CPU Lennard-Jones 12-6 (shifted) force component.
#pragma once

#include <string>
#include <vector>

#include "atomiverse/cpu/component.h"

namespace ATV {
namespace cpu {
namespace md {

/// @brief Lennard-Jones 12-6 with shifted potential.
///
/// V(r) = 4*eps * [(sig/r)^12 - (sig/r)^6] - V_shift
/// F(r) = 24*eps/r * [2*(sig/r)^12 - (sig/r)^6] * r_hat
///
/// Thread-parallel with per-thread force reduction buffers (Newton's 3rd law).
class LJ : public cpu::ComponentBase {
public:
    /// @param epsilon Well depth [eV].
    /// @param sigma Size parameter [Angstrom].
    /// @param cutoff Cutoff distance [Angstrom].
    LJ(Real epsilon, Real sigma, Real cutoff);

    std::string name() const override { return "lj"; }
    void build(const cpu::BuildContext& ctx) override;
    void launch(const cpu::LaunchContext& ctx) override;
    void destroy(const cpu::DestroyContext& ctx) override;

private:
    Real epsilon_, sigma_, cutoff_;
    Real sigma6_, cutoff_sq_, v_shift_;

    // Per-thread force reduction buffers
    int n_threads_ = 1;
    std::vector<std::vector<Real>> fx_thr_, fy_thr_, fz_thr_;
};

} // namespace md
} // namespace cpu
} // namespace ATV
