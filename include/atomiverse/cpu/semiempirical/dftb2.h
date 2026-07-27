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

/// @brief SCC-DFTB2 energy (+ optional D3 + forces) for H/C/N/O.
///
/// Self-consistent-charge DFTB2 (Elstner et al., Phys. Rev. B **58**, 7260
/// (1998), DOI 10.1103/PhysRevB.58.7260): electronic tight-binding + second-
/// order SCC Coulomb + pair repulsive. Integrals and onsite/Hubbard data come
/// from static mio-style tables in `basis_data/dftb_mio_tables.h` (no runtime
/// `.skf` parser, no public factory). Supported Z: 1, 6, 7, 8; any other
/// element or missing pair throws.
///
/// Optional **Grimme DFT-D3(BJ)** (`d3 = true`) is a geometry-only additive
/// correction using in-tree Grimme tables and the published DFTB3-D3(BJ)
/// damping set (s6, s8, a1, a2) as used by DFTB+ / simple-dftd3 `--bj dftb3`.
/// Three-body ATM is off; evaluation is non-periodic (molecular path).
///
/// Units on the public surface: energy [eV], positions [Å], forces [eV/Å].
///
/// @code
///   cpu::semiempirical::DFTB2 dftb;
///   dftb.d3 = true;
///   dftb.with_forces = true;
///   eval.on(evaluate::Phase::evaluate).add(std::move(dftb));
/// @endcode
class DFTB2 : public cpu::Semiempirical {
  public:
    /// @brief SCC residual tolerance on max |Δq| (dimensionless e).
    double scc_tol = 1e-8;
    /// @brief Maximum SCC iterations (dimensionless).
    int max_scc = 50;
    /// @brief When true, write forces F=−∇E [eV/Å] on launch (central FD of E).
    bool with_forces = true;
    /// @brief Opt-in Grimme DFT-D3 with Becke–Johnson damping (two-body).
    bool d3 = false;
    /// @brief D3 C6 scale (dimensionless). DFTB3-D3(BJ) default = 1.
    double d3_s6 = 1.0;
    /// @brief D3 C8 scale (dimensionless). DFTB3-D3(BJ) default = 0.5883.
    double d3_s8 = 0.5883;
    /// @brief BJ steepness (dimensionless). DFTB3-D3(BJ) default = 0.5719.
    double d3_a1 = 0.5719;
    /// @brief BJ offset [Bohr]. DFTB3-D3(BJ) default = 3.6017.
    double d3_a2 = 3.6017;

    DFTB2() = default;

    std::string name() const override { return "dftb2"; }

    void build(const BuildContext& ctx) override;
    void launch(const LaunchContext& ctx) override;

    /// @brief Last total energy E_elec+rep+disp [eV].
    double last_energy() const { return energy_; }
    /// @brief Last D3(BJ) dispersion energy [eV] (0 when `d3` is false).
    double last_e_disp() const { return e_disp_; }

  private:
    double energy_ = 0;
    double e_disp_ = 0;
};

} // namespace semiempirical
} // namespace cpu
} // namespace ATV
