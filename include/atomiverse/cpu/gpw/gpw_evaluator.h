/// @file cpu/gpw/gpw_evaluator.h
/// @brief GPWEvaluator — GPW DFT energy component for CPU drivers.
#pragma once

#include <array>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "atomiverse/core/role.h"
#include "atomiverse/core/units.h"
#include "atomiverse/cpu/component/gpw.h"
#include "atomiverse/cpu/gpw/cpmd.h"
#include "atomiverse/cpu/gpw/gpw_solver.h"
#include "atomiverse/cpu/gpw/xl_bomd.h"
#include "atomiverse/cpu/scf/aspc_predictor.h"
#include "atomiverse/molrs.hpp"

namespace ATV {
namespace cpu {
namespace gpw {

/// @brief One element's basis + pseudopotential — CP2K's `&KIND` concept.
struct Kind {
    molrs::Element element;
    GaussianBasis basis;
    GthPotential pot;

    Kind(molrs::Element element, GaussianBasis basis, GthPotential pot)
        : element(element), basis(std::move(basis)), pot(std::move(pot)) {}
};

/// @brief GPW total-energy evaluator component.
///
/// Parameters are public members / constructor args — no `GpwEvaluatorConfig`
/// bag.
///
/// **Variable cell:** `launch()` re-targets the solver to the frame's current
/// orthorhombic box on every call (`GpwDft::run` rebuilds all
/// geometry-dependent structure — grids, lattice images, S/T/Vpp — per call),
/// so cell-relaxation drivers may rescale `frame.simbox` between launches.
/// The MO warm start survives a cell change (the AO basis is
/// geometry-independent in size).
class GPWEvaluator : public cpu::GPW {
  public:
    /// @brief Plane-wave cutoff [Ry].
    double cutoff_ry = 300.0;
    /// @brief CP2K-style multigrid levels (`&MGRID NGRIDS`). Default 1 =
    ///        single density grid (CompactBasisGrid). With `rel_cutoff > 0`
    ///        and `ngrids > 1`, product collocation is used (RKS and UKS).
    int ngrids = 1;
    /// @brief CP2K REL_CUTOFF [Ry] for exponent→level mapping. Ignored when
    ///        `ngrids == 1`. Typical twin: 60 Ry with CUTOFF 400 / NGRIDS 4.
    double rel_cutoff = 0.0;
    /// @brief Max SCF iterations.
    int max_iter = 200;
    /// @brief SCF residual tolerance.
    double tol = 1e-7;
    /// @brief Density mixing parameter.
    double mix = 0.3;
    /// @brief XC functional kind.
    XcKind xc = XcKind::PADE;
    /// @brief Non-periodic (cluster) boundary conditions.
    bool nonperiodic = false;
    /// @brief Compute analytic forces (one SCF + one gradient pass).
    bool compute_forces = false;
    /// @brief Compute cell stress (default: analytic fixed-\(P\) path).
    ///        Writes Frame `kStress` [eV/Å³] Voigt; also `last_stress_ha_bohr3()`.
    bool compute_stress = false;
    /// @brief When true, use multi-SCF strain-FD (`gpw_stress_fd`) instead of
    ///        the product analytic path (`gpw_stress`). Default false.
    bool stress_use_fd = false;
    /// @brief Dimensionless strain step for analytic fixed-\(P\) FD and for
    ///        multi-SCF FD when `stress_use_fd` (default 2e-3).
    double stress_strain_h = 2e-3;
    /// @brief Per-axis stress mask (x,y,z): a masked axis skips samples and
    ///        reports σ exactly 0 — pair with the matching
    ///        `AnisotropicCellRelax` `relax_*` mask in variable-cell runs.
    std::array<bool, 3> stress_axes = {true, true, true};
    /// @brief FD displacement [Bohr] for forces (legacy path).
    double fd_step = 1e-3;
    /// @brief Grimme DFT-D3 (zero damping, PBE params). Exclusive vs other disp.
    bool d3 = false;
    /// @brief Grimme DFT-D3(BJ) (PBE params, CP2K CN). Exclusive vs other disp.
    bool d3bj = false;
    /// @brief Grimme DFT-D4 (PBE BJ, two-body). Exclusive vs other disp.
    bool d4 = false;
    /// @brief Tkatchenko–Scheffler (PBE free-atom TS). Exclusive vs other disp.
    bool ts = false;
    /// @brief Free-atom XDM (Becke–Johnson, QE NC-PBE a1/a2). Exclusive vs other disp.
    bool xdm = false;
    /// @brief MBD@rsSCS (PBE β=0.83). Exclusive vs other pairwise / nonlocal disp.
    ///        Default isolated; set `mbd_k_mesh` non-zero for periodic Ewald.
    bool mbd = false;
    /// @brief MBD k-mesh (0,0,0 = free-atom isolated; else periodic Ewald).
    int mbd_k_mesh[3] = {0, 0, 0};
    /// @brief Nonlocal vdW-DF (DRSLL) SC potential. Exclusive vs d3/d3bj/d4/ts/mbd/vdw_df2/rvv10.
    bool vdw_df = false;
    /// @brief Nonlocal vdW-DF2 (LMKLL) SC potential. Exclusive vs d3/d3bj/d4/ts/mbd/vdw_df/rvv10.
    bool vdw_df2 = false;
    /// @brief Nonlocal rVV10 SC potential. Exclusive vs d3/d3bj/d4/ts/mbd/vdw_df/vdw_df2.
    bool rvv10 = false;
    /// @brief Hirshfeld \(V/V_{\mathrm{free}}\) scaling for TS/MBD (requires ts or mbd).
    bool hirshfeld = false;
    /// @brief Orbital Transformation (OT) SCF minimizer instead of
    ///        diagonalization + DIIS. Same gapped ground state; preferred for
    ///        larger AO bases (avoids dense \(O(n_b^3)\) eigensolve each step).
    ///        Units-free switch; energy still [Hartree] on `last_result()`.
    bool use_ot = false;
    /// UKS spin occupations. When both zero, closed-shell RKS from valence/2.
    /// O₂ GTH-q6 triplet: n_alpha=7, n_beta=5.
    int n_alpha = 0;
    int n_beta = 0;
    /// @brief Monkhorst–Pack k mesh (CP2K `&KPOINTS SCHEME MONKHORST-PACK`).
    ///        `{1,1,1}` keeps the Γ-point path. Non-Γ engages multi-k RKS/UKS
    ///        SCF and (with `xl_bomd`) multi-k XL on RKS or UKS.
    int k_mesh[3] = {1, 1, 1};
    /// @brief Fold ±k pairs by time reversal (default true; set false for
    ///        an explicit unreduced mesh, e.g. multi-k XL regression).
    bool k_time_reversal = true;
    /// @brief Extended-Lagrangian BOMD (dissipative XL-BOMD, Niklasson
    ///        JCP 130, 214109 (2009)): propagate an auxiliary density as the
    ///        SCF initial guess across MD steps. Γ RKS (one total \(P\)),
    ///        Γ UKS (independent \(P_\alpha,P_\beta\)), multi-k RKS (one
    ///        complex \(P_k\) per k), or multi-k UKS (dual \(P_k^\alpha,
    ///        P_k^\beta\) per k via `XlBomdCplx`). Diagonalization path
    ///        only — `build()` throws with `use_ot`, `cpmd`, or
    ///        `second_gen_cp`. Units-free switch.
    bool xl_bomd = false;
    /// @brief XL-BOMD dissipation order K (Table I rows K = 3..9).
    int xl_bomd_k = 5;
    /// @brief Car–Parrinello MD: propagate occupied MO coefficients with a
    ///        fictitious mass and electronic RATTLE (not BOMD, not XL-BOMD).
    ///        After one cold SCF seed, each force step is a fixed-density KS
    ///        eval + electronic VV/RATTLE + analytic forces at that density.
    ///        Closed-shell RKS Γ or collinear Γ-UKS (dual α/β propagators) —
    ///        `build()` throws with `use_ot`, `xl_bomd`, or `second_gen_cp`.
    bool cpmd = false;
    /// @brief Fictitious electron mass \(\mu\) [a.u. electron mass]. Typical
    ///        CPMD values \(10^2\)–\(10^3\). Default 400.
    double cpmd_mu = 400.0;
    /// @brief Electronic (and nuclear) timestep [fs]. **Must match** the
    ///        `VelocityVerlet` dt used on the outer MD loop.
    double cpmd_dt = 0.1;
    /// @brief Nosé–Hoover electronic thermostat on fictitious \(\dot C\)
    ///        (requires `cpmd = true`). Regulates \(T_e\) toward
    ///        `cpmd_ekin_target_ha` with coupling time `cpmd_tau_e`.
    ///        Default false (NVE electrons / pure CP).
    bool cpmd_thermostat = false;
    /// @brief Target fictitious electronic kinetic energy \(K_0\) [Hartree].
    ///        Must be > 0 when `cpmd_thermostat` is true.
    double cpmd_ekin_target_ha = 1e-5;
    /// @brief Electronic thermostat coupling time \(\tau_e\) [fs].
    double cpmd_tau_e = 50.0;
    /// @brief Second-generation Car–Parrinello MD (Kühne et al., PRL 98,
    ///        066401 (2007)): Kolafa ASPC density predictor + a **partial**
    ///        SCF corrector (not full BOMD, not fictitious-mass CPMD, not
    ///        Niklasson XL-BOMD).
    ///        - **Γ RKS:** one real total \(P\) (`ASPCPredictor`)
    ///        - **Γ UKS:** dual real \(P_\alpha, P_\beta\)
    ///        - **multi-k RKS:** one complex ASPC per sampled \(k\) on \(P_k\)
    ///        - **multi-k UKS:** dual complex ASPC per \(k\) (\(P_k^\alpha,
    ///          P_k^\beta\))
    ///        **OT corrector** is product at Γ when `use_ot=true` (RKS:
    ///        `second-gen-cpmd-ot`; UKS: `second-gen-cpmd-uks-ot`). Multi-k +
    ///        `use_ot` throws (solver has no \(P_k\) OT warm-start).
    ///        Mutually exclusive with `cpmd` / `xl_bomd`. Requires
    ///        `compute_forces=true`.
    bool second_gen_cp = false;
    /// @brief ASPC extrapolation order \(K\) (Kolafa 2004). Typical 1–4.
    int second_gen_cp_k = 3;
    /// @brief Corrector SCF iteration budget (≪ full BOMD `max_iter`).
    ///        Production steps after ASPC is warm use this cap.
    int second_gen_cp_max_iter = 5;
    /// @brief Dudarev DFT+U (PLUS_U); empty disables. U−J in eV per entry.
    std::vector<PlusU> plus_u;
    /// @brief Per-atom initial spin seed for UKS AFM/FM (0 or n_atoms).
    std::vector<double> atom_spin_guess;
    /// @brief One Kind per distinct element in the system.
    std::vector<Kind> kinds;

    /// @param kinds One Kind per distinct element.
    /// @param cutoff_ry Plane-wave cutoff [Ry].
    /// @param max_iter Max SCF iterations.
    /// @param tol SCF residual tolerance.
    /// @param mix Density mixing (DIIS path; ignored when `use_ot`).
    /// @param xc XC functional.
    /// @param nonperiodic Cluster BC when true.
    /// @param compute_forces Enable analytic forces.
    /// @param fd_step FD step [Bohr].
    /// @param d3 Enable DFT-D3 zero damping.
    /// @param use_ot Enable OT minimizer (product path; no env required).
    explicit GPWEvaluator(std::vector<Kind> kinds, double cutoff_ry = 300.0,
                          int max_iter = 200, double tol = 1e-7,
                          double mix = 0.3, XcKind xc = XcKind::PADE,
                          bool nonperiodic = false, bool compute_forces = false,
                          double fd_step = 1e-3, bool d3 = false,
                          bool use_ot = false)
        : cutoff_ry(cutoff_ry), max_iter(max_iter), tol(tol), mix(mix), xc(xc),
          nonperiodic(nonperiodic), compute_forces(compute_forces),
          fd_step(fd_step), d3(d3), use_ot(use_ot), kinds(std::move(kinds)) {}

    std::string name() const override { return "gpw_evaluator"; }

    void build(const cpu::BuildContext& ctx) override;
    void launch(const cpu::LaunchContext& ctx) override;

    double last_energy_hartree() const { return last_.e_total; }
    /// @brief ⟨S²⟩ from the last UKS run (0 for closed-shell).
    double last_s2() const { return last_.s2; }
    const GpwResult& last_result() const { return last_; }
    /// @brief Last stress Voigt [Hartree/Bohr³] when `compute_stress` ran.
    const std::array<double, 6>& last_stress_ha_bohr3() const {
        return last_stress_ha_;
    }
    /// @brief Electronic kinetic energy \(T_e\) [eV] from the last CPMD
    ///        advance (0 when `cpmd` is off or before the first electronic step).
    double last_electronic_ke_ev() const { return last_te_ev_; }
    /// @brief Nosé–Hoover electronic thermostat energy \(E_{\mathrm{therm}}\)
    ///        [eV] from the last CPMD advance (0 when thermostat is off).
    double last_electronic_thermostat_ev() const { return last_etherm_ev_; }
    /// @brief Car–Parrinello conserved piece \(E_{\mathrm{KS}} + T_e\) [eV]
    ///        (add nuclear KE from `frame.meta.get(kKinetic)` for full NVE
    ///        \(E_{\mathrm{cons}}\); when `cpmd_thermostat`, also add
    ///        `last_electronic_thermostat_ev()` for \(E_{\mathrm{ext}}\)).
    double last_cp_energy_ev() const {
        return last_.e_total * ATV::units::kHartreeInEv + last_te_ev_;
    }

  private:
    GpwResult last_{};
    std::unique_ptr<GpwDft> solver_;
    std::array<double, 6> last_stress_ha_{};
    /// XL-BOMD auxiliary-density propagator (engaged when `xl_bomd` on Γ).
    /// RKS: total \(P\). UKS: \(P_\alpha\). Empty when multi-k XL is active.
    std::optional<XlBomd> xl_;
    /// XL-BOMD β-spin propagator (Γ-UKS only; empty for RKS / multi-k).
    std::optional<XlBomd> xl_beta_;
    /// Multi-k XL α (or RKS total): one complex propagator per sampled k
    /// (empty on Γ).
    std::vector<XlBomdCplx> xl_k_;
    /// Multi-k UKS XL β: one complex propagator per k (empty for RKS / Γ).
    std::vector<XlBomdCplx> xl_k_beta_;
    /// Scratch for multi-k XL `guess_P_k` (lifetime covers `run`).
    std::vector<linalg::CplxMatrix> xl_guess_p_k_;
    /// Scratch for multi-k UKS XL `guess_P_k_beta` (lifetime covers `run`).
    std::vector<linalg::CplxMatrix> xl_guess_p_k_beta_;
    /// CPMD electronic propagator (engaged when `cpmd`). RKS: total occ.
    /// UKS: α channel (occ=1).
    std::optional<Cpmd> cpmd_;
    /// CPMD β-spin propagator (UKS only; empty for RKS).
    std::optional<Cpmd> cpmd_beta_;
    double last_te_ev_ = 0.0;     ///< last electronic KE [eV]
    double last_etherm_ev_ = 0.0; ///< last electronic NH energy [eV]
    /// Second-gen CP ASPC density history (engaged when `second_gen_cp` on Γ).
    /// RKS: total \(P\). UKS: \(P_\alpha\). Empty when multi-k SGCP is active.
    std::optional<scf::ASPCPredictor> sgcp_aspc_;
    /// Second-gen CP β-spin ASPC (Γ-UKS only; empty for RKS / multi-k).
    std::optional<scf::ASPCPredictor> sgcp_aspc_beta_;
    /// Scratch for ASPC → column-major `guess_P` (lifetime covers `run`).
    linalg::RealMatrix sgcp_guess_p_;
    /// Scratch for ASPC β → column-major `guess_P_beta` (UKS).
    linalg::RealMatrix sgcp_guess_p_beta_;
    /// Multi-k SGCP α (or RKS total): one complex ASPC per sampled k (empty on Γ).
    std::vector<scf::ASPCPredictorCplx> sgcp_aspc_k_;
    /// Multi-k SGCP β: one complex ASPC per k (empty for RKS / Γ).
    std::vector<scf::ASPCPredictorCplx> sgcp_aspc_k_beta_;
    /// Scratch for multi-k SGCP `guess_P_k` (lifetime covers `run`).
    std::vector<linalg::CplxMatrix> sgcp_guess_p_k_;
    /// Scratch for multi-k UKS SGCP `guess_P_k_beta`.
    std::vector<linalg::CplxMatrix> sgcp_guess_p_k_beta_;
};

} // namespace gpw
} // namespace cpu
} // namespace ATV
