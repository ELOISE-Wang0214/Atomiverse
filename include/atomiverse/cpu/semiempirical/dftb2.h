/// @file cpu/semiempirical/dftb2.h
/// @brief SCC-DFTB2 skeleton — H-only textbook path (CPU).
#pragma once

#include <string>
#include <vector>

#include "atomiverse/cpu/component/semiempirical.h"
#include "atomiverse/cpu/frame/frame.h"

namespace ATV {
namespace cpu {
namespace semiempirical {

/// @brief Minimal SCC-DFTB2 energy (+ optional forces) for H systems.
///
/// Textbook skeleton (Elstner et al., PRB 58, 7260 (1998) shape):
/// electronic TB + self-consistent charge Coulomb. In-tree H–H SK-like
/// tables only; other elements throw. Not a product DFTB suite.
class DFTB2 : public cpu::Semiempirical {
  public:
    /// @brief SCC residual tolerance on max |Δq| (dimensionless e).
    double scc_tol = 1e-8;
    /// @brief Maximum SCC iterations (dimensionless).
    int max_scc = 50;
    /// @brief When true, write analytic forces [eV/Å] on launch.
    bool with_forces = true;

    DFTB2() = default;

    std::string name() const override { return "dftb2"; }

    void build(const BuildContext& ctx) override;
    void launch(const LaunchContext& ctx) override;

    /// @brief Last electronic energy [eV].
    double last_energy() const { return energy_; }

  private:
    double energy_ = 0;
};

} // namespace semiempirical
} // namespace cpu
} // namespace ATV
