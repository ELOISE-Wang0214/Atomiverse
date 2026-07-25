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
    /// @brief Grimme DFT-D3 (zero damping, PBE params). Exclusive vs d3bj/d4/ts/mbd.
    bool d3 = false;
    /// @brief Grimme DFT-D3(BJ) (PBE params, CP2K CN). Exclusive vs d3/d4/ts/mbd.
    bool d3bj = false;
    /// @brief Grimme DFT-D4 (PBE BJ, two-body). Exclusive vs d3/d3bj/ts/mbd.
    bool d4 = false;
    /// @brief Tkatchenko–Scheffler (PBE free-atom TS). Exclusive vs d3/d3bj/d4/mbd.
    bool ts = false;
    /// @brief MBD@rsSCS free-atom (PBE β=0.83). Exclusive vs d3/d3bj/d4/ts/vdw_df/vdw_df2/rvv10.
    bool mbd = false;
    /// @brief Nonlocal vdW-DF (DRSLL) post-SCF. Exclusive vs d3/d3bj/d4/ts/mbd/vdw_df2/rvv10.
    bool vdw_df = false;
    /// @brief Nonlocal vdW-DF2 (LMKLL) post-SCF. Exclusive vs d3/d3bj/d4/ts/mbd/vdw_df/rvv10.
    bool vdw_df2 = false;
    /// @brief Nonlocal rVV10 post-SCF. Exclusive vs d3/d3bj/d4/ts/mbd/vdw_df/vdw_df2.
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
    /// @brief Extended-Lagrangian BOMD (dissipative XL-BOMD, Niklasson
    ///        JCP 130, 214109 (2009)): propagate an auxiliary density as the
    ///        SCF initial guess across MD steps. Closed-shell RKS
    ///        diagonalization path only — `build()` throws when combined
    ///        with UKS occupations or `use_ot`. Units-free switch.
    bool xl_bomd = false;
    /// @brief XL-BOMD dissipation order K (Table I rows K = 3..9).
    int xl_bomd_k = 5;
    /// @brief Car–Parrinello MD: propagate occupied MO coefficients with a
    ///        fictitious mass (not BOMD, not XL-BOMD). After one cold SCF seed,
    ///        each force step is a fixed-density KS eval + electronic Verlet +
    ///        analytic forces at that density. Closed-shell RKS Γ only —
    ///        `build()` throws with UKS, `use_ot`, or `xl_bomd`.
    bool cpmd = false;
    /// @brief Fictitious electron mass \(\mu\) [a.u. electron mass]. Typical
    ///        CPMD values \(10^2\)–\(10^3\). Default 400.
    double cpmd_mu = 400.0;
    /// @brief Electronic (and nuclear) timestep [fs]. **Must match** the
    ///        `VelocityVerlet` dt used on the outer MD loop.
    double cpmd_dt = 0.1;
    /// @brief Second-generation Car–Parrinello MD (Kühne et al., PRL 98,
    ///        066401 (2007)): Kolafa ASPC density predictor + a **partial**
    ///        SCF corrector (not full BOMD, not fictitious-mass CPMD, not
    ///        Niklasson XL-BOMD). Closed-shell RKS Γ diagonalization only —
    ///        `build()` throws with UKS, `use_ot`, `cpmd`, or `xl_bomd`.
    ///        Requires `compute_forces=true`.
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
    /// @brief Car–Parrinello conserved piece \(E_{\mathrm{KS}} + T_e\) [eV]
    ///        (add nuclear KE from `frame.meta.get(kKinetic)` for full
    ///        \(E_{\mathrm{cons}}\)).
    double last_cp_energy_ev() const {
        return last_.e_total * ATV::units::kHartreeInEv + last_te_ev_;
    }

  private:
    GpwResult last_{};
    std::unique_ptr<GpwDft> solver_;
    std::array<double, 6> last_stress_ha_{};
    /// XL-BOMD auxiliary-density propagator (engaged when `xl_bomd`).
    std::optional<XlBomd> xl_;
    /// CPMD electronic propagator (engaged when `cpmd`).
    std::optional<Cpmd> cpmd_;
    double last_te_ev_ = 0.0; ///< last electronic KE [eV]
    /// Second-gen CP ASPC density history (engaged when `second_gen_cp`).
    std::optional<scf::ASPCPredictor> sgcp_aspc_;
    /// Scratch for ASPC → column-major `guess_P` (lifetime covers `run`).
    linalg::RealMatrix sgcp_guess_p_;
};

} // namespace gpw
} // namespace cpu
} // namespace ATV
