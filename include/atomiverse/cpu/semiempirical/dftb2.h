/// @file cpu/semiempirical/dftb2.h
/// @brief SCC-DFTB2 multi-element method (static in-tree SK tables, CPU).
#pragma once

#include <array>
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

/// @brief SCC-DFTB2 energy (+ optional D3 dispersion, nuclear forces, cell
///        stress) for H/C/N/O.
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
///   (molecular limit when L ≫ r_cut). **k-mesh beyond Γ:** set
///   `k_mesh = {n1,n2,n3}` (Monkhorst–Pack, CP2K/DFTB+ convention) for
///   Bloch-phased complex \(H^0(\mathbf{k})\), \(S(\mathbf{k})\), and
///   k-weighted SCC energy (`reg_dftb2-kpoints`). Default `{1,1,1}` keeps the
///   real Γ solver (pin-safe). Multi-k is **energy-only** (restricted
///   MioSketch); `with_forces` / `with_stress` / open shell / ThreeOb throw.
///   **Analytic** nuclear forces for closed shell **and** collinear open shell
///   (HF + Pulay + SCC γ + repulsive; spin shifts on ∂S) under the min-image
///   molecular / Γ limit — fixed-Δq/m Hellmann–Feynman (see
///   `reg_dftb2-analytic-forces`, `reg_dftb2-spin-analytic-forces`). Closed
///   shell also offers an analytic intensive cell stress on `kStress` [eV/Å³]
///   — the pair virial of those same bond terms (`reg_dftb2-analytic-stress`);
///   open shell falls back to multi-SCC affine strain-FD.
/// - `sk_set = ThreeOb`: `basis_data/dftb_3ob_tables.h` 3ob-3-1 radial H/S +
///   repulsive (Gaus et al., JCTC 2013; organic CHNO). Closed shell: **analytic**
///   nuclear forces (same HF + Pulay + γ + repulsive assembly as MioSketch, on
///   the Ha/Bohr 3ob tables; `reg_dftb2-analytic-forces-threeob`). Stress is
///   multi-SCC strain-FD when requested. Spin / nonzero charge are not
///   supported on the ThreeOb path.
///
/// Optional **Grimme DFT-D3(BJ)** (`d3 = true`) is a geometry-only additive
/// correction using in-tree Grimme tables and the published DFTB3-D3(BJ)
/// damping set (s6, s8, a1, a2) as used by DFTB+ / simple-dftd3 `--bj dftb3`.
/// It is added on top of **either** `sk_set` path, with its own analytic
/// gradient, so `F = −∇E` holds for the total. Three-body ATM is off.
/// **Limitation:** the dispersion sum is evaluated as an *isolated cluster*
/// (no lattice images) even when the electronic Γ-point path is lattice-summed,
/// so `d3` is only physical for molecules / large boxes — a periodic D3 lattice
/// sum for DFTB2 is not implemented.
///
/// Open shell uses `charge` and `multiplicity` to set nα, nβ. Units: energy
/// [eV], positions [Å], forces [eV/Å], stress [eV/Å³] Voigt
/// `(xx,yy,zz,xy,xz,yz)` with mechanical convention
/// \(\sigma=-(1/\Omega)\partial E/\partial\varepsilon\).
///
/// @code
///   cpu::semiempirical::DFTB2 dftb;
///   dftb.sk_set = cpu::semiempirical::DftbSkSet::ThreeOb;
///   dftb.multiplicity = 3;   // O2 triplet (MioSketch)
///   dftb.d3 = true;
///   dftb.with_forces = true;
///   dftb.with_stress = true;
///   dftb.k_mesh = {4, 1, 1};   // multi-k energy (forces off)
///   eval.on(evaluate::Phase::evaluate).add(std::move(dftb));
/// @endcode
class DFTB2 : public cpu::Semiempirical {
  public:
    /// @brief SCC residual tolerance on max |Δq| (and |Δm| when open shell).
    double scc_tol = 1e-8;
    /// @brief Maximum SCC iterations (dimensionless). Exhausting it throws.
    int max_scc = 50;
    /// @brief Anderson charge-mixing relaxation β ∈ (0,1] (dimensionless).
    ///
    /// Undamped Picard (no mixing) is a period-2 limit cycle for polar
    /// systems: it never meets `scc_tol` and the reported energy is picked by
    /// the parity of `max_scc`. Charges (and the spin population m, mixed as
    /// one vector with them) go through `internal::ChargeMixer`. The default
    /// converges every in-tree case in ≤ 14 iterations; β ≳ 0.7 stops
    /// converging for hydrogen-bonded systems.
    double mix_beta = 0.3;
    /// @brief When true, write F=−∇E [eV/Å] on launch.
    bool with_forces = true;
    /// @brief When true, write intensive cell stress to `kStress` [eV/Å³].
    ///
    /// Restricted MioSketch: closed-form pair virial of the analytic force
    /// terms. Open-shell / ThreeOb: multi-SCC affine strain-FD of total E.
    /// **Excludes the `d3` dispersion term on both paths** — see
    /// `regressions/semiempirical/dftb2-analytic-stress/REFERENCE.md` Gaps.
    /// Multi-k (`any k_mesh[d] > 1`) does not support stress.
    bool with_stress = false;
    /// @brief Static SK parameter set (member — not a Config bag).
    DftbSkSet sk_set = DftbSkSet::MioSketch;
    /// @brief Net molecular charge [e] (subtracted from valence electron count).
    int charge = 0;
    /// @brief Spin multiplicity 2S+1 (dimensionless). Default 1 = closed shell.
    int multiplicity = 1;
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
    /// @brief Monkhorst–Pack mesh subdivisions per reciprocal axis
    ///        (dimensionless). Default `{1,1,1}` = Γ only (real solver).
    ///
    /// CP2K/DFTB+ convention: \(\kappa_d = (2r - N_d - 1)/(2 N_d)\),
    /// \(r = 1..N_d\). Any \(N_d < 1\) throws. Multi-k requires restricted
    /// MioSketch and energy-only (`with_forces`/`with_stress` false).
    std::array<int, 3> k_mesh{{1, 1, 1}};
    /// @brief Fold \(\pm\mathbf{k}\) pairs (double weight). Default true.
    ///        Exact for real-Hamiltonian closed-shell densities.
    bool k_time_reversal = true;

    DFTB2() = default;

    std::string name() const override { return "dftb2"; }

    void build(const BuildContext& ctx) override;
    void launch(const LaunchContext& ctx) override;

    /// @brief Last total energy E_elec+rep+disp [eV].
    double last_energy() const { return energy_; }

    /// @brief Last D3(BJ) dispersion energy [eV] (0 when `d3` is false).
    double last_e_disp() const { return e_disp_; }

    /// @brief SCC residual max|Δq_out − Δq_in| at exit [e]; below `scc_tol`
    ///        whenever `launch` returned (it throws otherwise).
    double last_scc_residual() const { return scc_residual_; }

    /// @brief SCC iterations the last solve took (dimensionless).
    int last_scc_iterations() const { return scc_iterations_; }

    /// @brief Last intensive Voigt stress [eV/Å³] (zero if stress not requested).
    const std::array<double, 6>& last_stress() const { return last_stress_; }

    /// @brief Per-atom Mulliken charge fluctuation Δq_A [e] after last launch.
    const std::vector<double>& last_mulliken_dq() const { return last_dq_; }

    /// @brief Per-atom Mulliken spin population m_A = Nα−Nβ [e] after last launch.
    const std::vector<double>& last_mulliken_spin() const { return last_spin_; }

  private:
    double energy_ = 0;
    double e_disp_ = 0;
    double scc_residual_ = 0;
    int scc_iterations_ = 0;
    std::array<double, 6> last_stress_{};
    std::vector<double> last_dq_;
    std::vector<double> last_spin_;
};

} // namespace semiempirical
} // namespace cpu
} // namespace ATV
