/// @file cpu/semiempirical/dftb3.h
/// @brief SCC-DFTB3 multi-element method (static in-tree SK tables, CPU).
#pragma once

#include <string>

#include "atomiverse/cpu/component/semiempirical.h"
#include "atomiverse/cpu/frame/frame.h"

namespace ATV {
namespace cpu {
namespace semiempirical {

/// @brief SCC-DFTB3 energy (+ optional FD forces) for H/C/N/O.
///
/// Third-order self-consistent-charge DFTB (Gaus et al., J. Chem. Theory
/// Comput. **7**, 931 (2011), DOI 10.1021/ct100684s) on top of the
/// Elstner DFTB2 skeleton: electronic tight-binding + second-order SCC +
/// third-order charge fluctuation (\(U^{\mathrm{d}}\), \(\Gamma_{AB}\)) +
/// pair repulsive. Integrals and onsite/Hubbard data come from static
/// mio-style tables in `basis_data/dftb_mio_tables.h` (no runtime `.skf`
/// parser, no public factory). Supported Z: 1, 6, 7, 8; any other element
/// or missing pair throws.
///
/// Units on the public surface: energy [eV], positions [Å], forces [eV/Å].
///
/// @code
///   cpu::semiempirical::DFTB3 dftb3;
///   dftb3.with_forces = false;
///   eval.on(evaluate::Phase::evaluate).add(std::move(dftb3));
/// @endcode
class DFTB3 : public cpu::Semiempirical {
  public:
    /// @brief SCC residual tolerance on max |Δq| (dimensionless e).
    double scc_tol = 1e-8;
    /// @brief Maximum SCC iterations (dimensionless).
    int max_scc = 50;
    /// @brief When true, write forces F=−∇E [eV/Å] on launch (central FD of E).
    bool with_forces = true;

    DFTB3() = default;

    std::string name() const override { return "dftb3"; }

    void build(const BuildContext& ctx) override;
    void launch(const LaunchContext& ctx) override;

    /// @brief Last electronic+repulsive energy [eV].
    double last_energy() const { return energy_; }

    /// @brief Last second-order SCC double-counting term \(E^{(2)}\) [eV].
    double last_e2() const { return e2_; }

    /// @brief Last third-order term \(E^{(3)}\) [eV].
    double last_e3() const { return e3_; }

  private:
    double energy_ = 0;
    double e2_ = 0;
    double e3_ = 0;
};

} // namespace semiempirical
} // namespace cpu
} // namespace ATV
