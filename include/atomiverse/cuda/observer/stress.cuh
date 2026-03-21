/// @file cuda/observer/stress.cuh
/// @brief VirialStress -- computes virial stress tensor from kinetic + pair contributions.
///
/// σ_αβ = (1/V) [Σ_i m_i v_iα v_iβ + Σ_{pairs} r_ij_α F_ij_β]
///
/// Cross-driver component: empty default_bindings(), bind at post_force or end_of_step.
/// CUDA-free header.
#pragma once

#include <memory>
#include <string>

#include "atomiverse/cuda/component.cuh"

namespace ATV {

/// @brief Virial stress tensor computation.
///
/// Computes 6-component Voigt stress (xx, yy, zz, xy, xz, yz) and writes to
/// DeviceFrame::stress. Includes kinetic and pair-wise contributions.
class VirialStress : public ComponentBase {
public:
    VirialStress();
    ~VirialStress() override;
    VirialStress(VirialStress&&) noexcept;
    VirialStress& operator=(VirialStress&&) noexcept;

    std::string name() const override { return "virial_stress"; }
    void build(const BuildContext& ctx) override;
    void launch(const LaunchContext& ctx) override;
    void destroy(const DestroyContext& ctx) override;

private:
    struct State;
    std::unique_ptr<State> state_;
};

} // namespace ATV
