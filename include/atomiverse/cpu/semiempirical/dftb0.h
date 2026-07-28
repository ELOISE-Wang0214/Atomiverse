/// @file cpu/semiempirical/dftb0.h
/// @brief Non-SCC DFTB (DFTB0) method — static in-tree SK tables, CPU.
#pragma once

#include <string>

#include "atomiverse/cpu/component/semiempirical.h"
#include "atomiverse/cpu/frame/frame.h"

namespace ATV {
namespace cpu {
namespace semiempirical {

/// @brief Non-SCC DFTB / DFTB0 energy (+ optional forces) for H/C/N/O.
///
/// Electronic tight-binding without self-consistent charge: Slater–Koster
/// H0/S from static mio-style tables, generalized eigen, band energy + pair
/// repulsive. No Hubbard / γ / Δq iteration (contrast @ref DFTB2). Tables:
/// `basis_data/dftb_mio_tables.h` (not bit-identical to DFTB.org mio-1-1).
/// Supported Z: 1, 6, 7, 8; any other element or missing pair throws.
///
/// Family oracle for pins: **DFTB+** (non-SCC / DFTB0 path). Units on the
/// public surface: energy [eV], positions [Å], forces [eV/Å].
///
/// @code
///   cpu::semiempirical::DFTB0 dftb0;
///   dftb0.with_forces = true;
///   eval.on(evaluate::Phase::evaluate).add(std::move(dftb0));
/// @endcode
class DFTB0 : public cpu::Semiempirical {
  public:
    /// @brief When true, write forces F=−∇E [eV/Å] on launch (central FD of E).
    bool with_forces = true;

    DFTB0() = default;

    std::string name() const override { return "dftb0"; }

    void build(const BuildContext& ctx) override;
    void launch(const LaunchContext& ctx) override;

    /// @brief Last electronic+repulsive energy [eV].
    double last_energy() const { return energy_; }

  private:
    double energy_ = 0;
};

} // namespace semiempirical
} // namespace cpu
} // namespace ATV
