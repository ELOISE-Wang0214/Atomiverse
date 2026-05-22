/// @file cuda/md/init.cuh
/// @brief Initializer components for forces, energy, kinetic energy, and velocity.
///
/// Initializers prepare accumulator fields before the computation phase:
/// - ZeroForces / ZeroEnergy / ZeroKinetic: zero accumulators via cudaMemsetAsync.
/// - InitVelocity: Maxwell-Boltzmann velocity initialization (LAMMPS `velocity create`).
///
/// The Zero* components are **phase-agnostic**: they declare no
/// `default_bindings()`, so the concrete driver binds them to its own
/// phase enum via explicit `on(phase).add()`. `cuda::md::MD` binds them
/// at `md::Phase`; `cuda::minimizer` drivers reuse the same classes at
/// `minimizer::Phase`. This mirrors the CPU resolution in `cpu/md/init.h`.
/// CUDA-free header.
///
/// @see md::MD for the MD driver that explicitly binds the Zero* initializers
/// @see md::KineticEnergy for the component that accumulates kinetic energy
#pragma once

#include <string>
#include <vector>

#include "atomiverse/cuda/component.cuh"
#include "atomiverse/facade/component_role.h"
#include "atomiverse/facade/stage.h"

namespace ATV {
namespace cuda {
namespace md {

/// @brief Zeros force arrays (fx, fy, fz) on device using cudaMemsetAsync.
///
/// Phase-agnostic infrastructure component — the driver binds it to its
/// own phase via explicit `on(phase).add()` (md::MD binds it at
/// pre_force, order -1000). Uses the DMA engine (cudaMemsetAsync), not a
/// compute kernel.
class ZeroForces : public ComponentBase {
  public:
    /// @brief Component name.
    /// @return "zero_forces".
    std::string name() const override { return "zero_forces"; }

    /// @brief Structural role — none. Zeroing an accumulator is not a
    ///        loop slot; the driver decides where it binds.
    /// @return ComponentRole::None.
    ComponentRole role() const override { return ComponentRole::None; }

    /// @brief Zero fx, fy, fz arrays via cudaMemsetAsync.
    /// @param ctx Launch context.
    void launch(const LaunchContext& ctx) override;
};

/// @brief Zeros potential energy scalar on device using cudaMemsetAsync.
///
/// Phase-agnostic infrastructure component — the driver binds it to its
/// own phase via explicit `on(phase).add()` (md::MD binds it at
/// pre_force, order -999). Uses the DMA engine (cudaMemsetAsync), not a
/// compute kernel.
class ZeroEnergy : public ComponentBase {
  public:
    /// @brief Component name.
    /// @return "zero_energy".
    std::string name() const override { return "zero_energy"; }

    /// @brief Structural role — none. Zeroing an accumulator is not a
    ///        loop slot; the driver decides where it binds.
    /// @return ComponentRole::None.
    ComponentRole role() const override { return ComponentRole::None; }

    /// @brief Zero energy scalar via cudaMemsetAsync.
    /// @param ctx Launch context.
    void launch(const LaunchContext& ctx) override;
};

/// @brief Zeros kinetic energy scalar on device using cudaMemsetAsync.
///
/// Phase-agnostic infrastructure component — the driver binds it to its
/// own phase via explicit `on(phase).add()` (md::MD binds it at
/// end_of_step, order -1000). Uses the DMA engine (cudaMemsetAsync), not
/// a compute kernel.
class ZeroKinetic : public ComponentBase {
  public:
    /// @brief Component name.
    /// @return "zero_kinetic".
    std::string name() const override { return "zero_kinetic"; }

    /// @brief Structural role — none. Zeroing an accumulator is not a
    ///        loop slot; the driver decides where it binds.
    /// @return ComponentRole::None.
    ComponentRole role() const override { return ComponentRole::None; }

    /// @brief Zero kinetic energy scalar via cudaMemsetAsync.
    /// @param ctx Launch context.
    void launch(const LaunchContext& ctx) override;
};

/// @brief Maxwell-Boltzmann velocity initialization (LAMMPS `velocity create`).
///
/// Generates Gaussian-distributed velocities at the target temperature
/// and removes center-of-mass velocity (zero total momentum).
///
/// One-time initialization in build(). launch() is a no-op.
///
/// `role()` stays `None` and `capabilities()` stays `None` (the
/// `ComponentBase` defaults), mirroring `cpu::md::InitVelocity`: it is a
/// structural one-shot setup component, neither a loop slot (Integrator
/// / Thermostat / …) nor a per-step writer of a `Capability` field — it
/// touches velocities once in `build()`, not in `launch()`.
///
/// @par Units
/// Temperature in Kelvin. Velocities in Angstrom/fs.
/// v_i ~ N(0, sqrt(kB*T / (kMvv2e * m_i))) per component.
///
/// @par Usage
/// @code{.cpp}
/// sim.add(md::InitVelocity{300.0, 42});  // 300 K, seed 42
/// @endcode
///
/// @see LAMMPS velocity create: https://docs.lammps.org/velocity.html
class InitVelocity : public ComponentBase {
  public:
    /// @brief Construct velocity initializer.
    /// @param temperature Target temperature [K].
    /// @param seed Random number seed.
    InitVelocity(double temperature, unsigned long long seed)
        : temperature_(temperature), seed_(seed) {}

    /// @brief Component name.
    /// @return "init_velocity".
    std::string name() const override { return "init_velocity"; }

    /// @brief Default bindings: initial_integrate at order -2000 (before integrator).
    /// @return Single binding to initial_integrate.
    std::vector<DefaultBinding> default_bindings() const override {
        return {{static_cast<int>(ATV::md::Phase::initial_integrate), -2000}};
    }

    /// @brief Generate Maxwell-Boltzmann velocities on device, remove COM velocity.
    /// @param ctx Build context.
    void build(const BuildContext& ctx) override;

    /// @brief No-op (velocities initialized once in build).
    /// @param ctx Launch context.
    void launch(const LaunchContext& ctx) override;

    /// @brief Target temperature.
    /// @return Temperature [K].
    double temperature() const { return temperature_; }

    /// @brief Random seed.
    /// @return Seed value.
    unsigned long long seed() const { return seed_; }

  private:
    double temperature_;      ///< Target temperature [K].
    unsigned long long seed_; ///< RNG seed.
};

} // namespace md
} // namespace cuda
} // namespace ATV
