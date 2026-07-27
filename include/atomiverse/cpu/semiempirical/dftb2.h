/// @file cpu/semiempirical/dftb2.h
/// @brief SCC-DFTB2 multi-element method (static in-tree SK tables, CPU).
#pragma once

#include <string>
#include <vector>

#include "atomiverse/cpu/component/semiempirical.h"
#include "atomiverse/cpu/frame/frame.h"

namespace ATV {
namespace cpu {
namespace semiempirical {

/// @brief SCC-DFTB2 energy (+ optional analytic nuclear forces) for H/C/N/O.
///
/// Self-consistent-charge DFTB2 (Elstner et al., Phys. Rev. B **58**, 7260
/// (1998), DOI 10.1103/PhysRevB.58.7260): electronic tight-binding + second-
/// order SCC Coulomb + pair repulsive. Integrals and onsite/Hubbard data come
/// from static mio-style tables in `basis_data/dftb_mio_tables.h` (no runtime
/// `.skf` parser, no public factory). Supported Z: 1, 6, 7, 8; any other
/// element or missing pair throws.
///
/// Energy is \(E=\sum_i n_i\varepsilon_i - E_{\mathrm{SCC}} + E_{\mathrm{rep}}\)
/// with SCC Hamiltonian \(H=H_0+\tfrac12 S_{\mu\nu}(\Delta V_A+\Delta V_B)\).
/// When @c with_forces is true, nuclear forces are the closed-form gradient
/// (Hellmann–Feynman + Pulay + SCC shift on \(\partial S\) + fixed-\(\Delta q\)
/// \(\partial\gamma\) + repulsive) — the DFTB+ force module boundary — not a
/// re-SCC finite difference.
///
/// Units on the public surface: energy [eV], positions [Å], forces [eV/Å].
///
/// @code
///   cpu::semiempirical::DFTB2 dftb;
///   dftb.with_forces = true;
///   eval.on(evaluate::Phase::evaluate).add(std::move(dftb));
/// @endcode
class DFTB2 : public cpu::Semiempirical {
  public:
    /// @brief SCC residual tolerance on max |Δq| (dimensionless e).
    double scc_tol = 1e-8;
    /// @brief Maximum SCC iterations (dimensionless).
    int max_scc = 50;
    /// @brief When true, write analytic F=−∇E [eV/Å] on launch.
    bool with_forces = true;

    DFTB2() = default;

    std::string name() const override { return "dftb2"; }

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
