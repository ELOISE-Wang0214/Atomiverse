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

/// @brief Which static SK parameter set DFTB2 uses.
///
/// Parameters are a public member on the method (no `*Config` bag).
enum class DftbSkSet {
    /// Mio path: H–H uses transcribed mio-1-1 grids; other H/C/N/O pairs use
    /// exponential sketches (+ spin_W).
    MioSketch = 0,
    /// DFTB.org 3ob-3-1 organic CHNO radial tables (static header).
    ThreeOb = 1,
};

/// @brief SCC-DFTB2 energy (+ optional nuclear forces) for H/C/N/O.
///
/// Self-consistent-charge DFTB2 (Elstner et al., Phys. Rev. B **58**, 7260
/// (1998), DOI 10.1103/PhysRevB.58.7260) with optional **collinear** spin
/// polarization (Köhler / DFTB+ spin channel): dual α/β Fock matrices,
/// Mulliken Δq and spin population m = Nα−Nβ.
///
/// - `sk_set = MioSketch` (default): H–H from static mio-1-1
///   (`basis_data/dftb_mio_hh_sk.h`, DFTB+ pin
///   `regressions/semiempirical/dftb2-dftbplus-energy-pin/`); other pairs from
///   exponential sketch (`dftb_mio_tables.h` + spin_W). **Γ-point PBC:**
///   lattice-summed H⁰/S, screened γ, and repulsive within static SK cutoffs
///   (molecular limit when L ≫ r_cut). Closed shell (`multiplicity == 1`):
///   analytic nuclear forces (HF + Pulay + SCC γ + repulsive) under the
///   min-image molecular limit. Open shell: re-SCC central FD.
/// - `sk_set = ThreeOb`: `basis_data/dftb_3ob_tables.h` 3ob-3-1 radial H/S +
///   repulsive (Gaus et al., JCTC 2013; organic CHNO). Closed-shell energy;
///   forces are re-SCC FD. Spin channel is MioSketch-only.
///
/// Open shell uses `charge` and `multiplicity` to set nα, nβ. Units: energy
/// [eV], positions [Å], forces [eV/Å].
///
/// @code
///   cpu::semiempirical::DFTB2 dftb;
///   dftb.sk_set = cpu::semiempirical::DftbSkSet::ThreeOb;
///   dftb.multiplicity = 3;   // O2 triplet (MioSketch)
///   dftb.with_forces = true;
///   eval.on(evaluate::Phase::evaluate).add(std::move(dftb));
/// @endcode
class DFTB2 : public cpu::Semiempirical {
  public:
    /// @brief SCC residual tolerance on max |Δq| (and |Δm| when open shell).
    double scc_tol = 1e-8;
    /// @brief Maximum SCC iterations (dimensionless).
    int max_scc = 50;
    /// @brief When true, write F=−∇E [eV/Å] on launch.
    bool with_forces = true;
    /// @brief Static SK parameter set (member — not a Config bag).
    DftbSkSet sk_set = DftbSkSet::MioSketch;
    /// @brief Net molecular charge [e] (subtracted from valence electron count).
    int charge = 0;
    /// @brief Spin multiplicity 2S+1 (dimensionless). Default 1 = closed shell.
    int multiplicity = 1;

    DFTB2() = default;

    std::string name() const override { return "dftb2"; }

    void build(const BuildContext& ctx) override;
    void launch(const LaunchContext& ctx) override;

    /// @brief Last electronic+repulsive energy [eV].
    double last_energy() const { return energy_; }

    /// @brief Per-atom Mulliken charge fluctuation Δq_A [e] after last launch.
    const std::vector<double>& last_mulliken_dq() const { return last_dq_; }

    /// @brief Per-atom Mulliken spin population m_A = Nα−Nβ [e] after last launch.
    const std::vector<double>& last_mulliken_spin() const { return last_spin_; }

  private:
    double energy_ = 0;
    std::vector<double> last_dq_;
    std::vector<double> last_spin_;
};

} // namespace semiempirical
} // namespace cpu
} // namespace ATV
