/// @file cpu/potentials/angle_harmonic.h
/// @brief LAMMPS-style AngleHarmonic — V = K(θ − θ₀)² (CPU).
///
/// Was `ATV::HarmonicAngle` before the LAMMPS-naming refactor.
#pragma once

#include <string>
#include <vector>

#include "atomiverse/cpu/potentials/angle.h"

namespace ATV {
namespace cpu {

class AngleHarmonic : public Angle {
  public:
    /// @param K Force constants per angle type [eV/rad²].
    /// @param theta0 Equilibrium angles per angle type [radians].
    AngleHarmonic(std::vector<double> K, std::vector<double> theta0);

    std::string name() const override { return "harmonic_angle"; }
    Capability capabilities() const override {
        // Angles contribute to the virial tensor — narrow the base's
        // default (Energy|Force) to include Stress.
        return Capability::Energy | Capability::Force | Capability::Stress;
    }

    void build(const BuildContext& ctx) override;
    void destroy(const DestroyContext& ctx) override;

  protected:
    void compute(const LaunchContext& ctx) override;

  private:
    std::vector<double> K_, theta0_;

    int n_threads_ = 1;
    std::vector<std::vector<double>> fx_thr_, fy_thr_, fz_thr_;
};

} // namespace cpu
} // namespace ATV
