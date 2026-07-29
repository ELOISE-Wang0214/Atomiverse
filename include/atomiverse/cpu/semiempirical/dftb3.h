/// @file cpu/semiempirical/dftb3.h
/// @brief SCC-DFTB3 multi-element method (static in-tree SK tables, CPU).
#pragma once

#include <string>

#include "atomiverse/cpu/component/semiempirical.h"
#include "atomiverse/cpu/frame/frame.h"

namespace ATV {
namespace cpu {
namespace semiempirical {

/// @brief SCC-DFTB3 energy (+ optional analytic nuclear forces) for H/C/N/O.
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
/// The SCC potential is the exact derivative
/// \f$V_A=\partial(E^{(2)}+E^{(3)})/\partial\Delta q_A\f$, so the band
/// energy double-counts by \f$E^{(2)}+2E^{(3)}\f$ (Euler:
/// \f$\sum_A\Delta q_A V_A = 2E^{(2)}+3E^{(3)}\f$). Both live in
/// `internal/dftb3_scc.h`.
///
/// Nuclear forces are the **closed-form** gradient at the converged SCC
/// density — Hellmann–Feynman \f$\partial H^0/\partial R\f$, Pulay
/// \f$-W\,\partial S/\partial R\f$ with the energy-weighted density
/// \f$W_{\mu\nu}=\sum_o n_o\varepsilon_o C_{\mu o}C_{\nu o}\f$, the SCC shift
/// on \f$\partial S/\partial R\f$, the fixed-Δq
/// \f$\partial(E^{(2)}+E^{(3)})/\partial R\f$ (both the \f$\gamma\f$ and the
/// non-symmetric \f$\Gamma\f$ radial derivatives), and the pair repulsive
/// gradient. Because \f$H=\partial E/\partial P\f$ exactly, the density is
/// stationary and no coupled-perturbed charge response is needed; gated
/// against re-SCC finite differences by
/// `regressions/semiempirical/dftb3-analytic-forces/`. This is **not** a
/// finite-difference path: one SCC solve per force evaluation, not \f$6N+1\f$.
///
/// **PBC:** Γ-point lattice-image walk for H⁰/S, γ/Γ, repulsive, and analytic
/// forces (same shape as DFTB2 MioSketch). Closed shell only; no multi-k,
/// no cell stress, no dispersion.
///
/// Units on the public surface: energy [eV], positions [Å], forces [eV/Å].
///
/// @code
///   cpu::semiempirical::DFTB3 dftb3;
///   dftb3.with_forces = true;
///   eval.on(evaluate::Phase::evaluate).add(std::move(dftb3));
/// @endcode
class DFTB3 : public cpu::Semiempirical {
  public:
    /// @brief SCC residual tolerance on max |Δq| (dimensionless e).
    double scc_tol = 1e-8;
    /// @brief Maximum SCC iterations (dimensionless). Exhausting it throws.
    int max_scc = 50;
    /// @brief Anderson charge-mixing relaxation β ∈ (0,1] (dimensionless).
    ///
    /// Undamped Picard (no mixing) is a period-2 limit cycle for polar
    /// systems: it never meets `scc_tol` and the reported energy is picked by
    /// the parity of `max_scc`. Charges go through `internal::ChargeMixer`.
    /// The default converges every in-tree case in ≤ 14 iterations; β ≳ 0.7
    /// stops converging for hydrogen-bonded systems.
    double mix_beta = 0.3;
    /// @brief When true, write forces F=−∇E [eV/Å] on launch (closed form).
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

    /// @brief SCC residual max|Δq_out − Δq_in| at exit [e]; below `scc_tol`
    ///        whenever `launch` returned (it throws otherwise).
    double last_scc_residual() const { return scc_residual_; }

    /// @brief SCC iterations the last solve took (dimensionless).
    int last_scc_iterations() const { return scc_iterations_; }

  private:
    double energy_ = 0;
    double e2_ = 0;
    double e3_ = 0;
    double scc_residual_ = 0;
    int scc_iterations_ = 0;
};

} // namespace semiempirical
} // namespace cpu
} // namespace ATV
