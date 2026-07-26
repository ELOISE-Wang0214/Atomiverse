/// @file cpu/gpw/cpmd.h
/// @brief Cpmd — Car–Parrinello fictitious electronic DOF propagator (RATTLE +
///        optional Nosé–Hoover electronic thermostat).
///
/// Propagates occupied MO coefficients \(C\) (\(n_b \times n_{\mathrm{occ}}\))
/// with velocity-Verlet and a fictitious electron mass \(\mu\)
/// (Car & Parrinello, Phys. Rev. Lett. 55, 2471 (1985)). The electronic force
/// is the **S-tangent** projection of the closed-shell band-energy gradient
/// (\(P = 2 C C^\top\)):
///
///   \(\mu\,\ddot C = -(I - S C C^\top)\,(4 H C)\)
///
/// so \(\ddot C \to 0\) on the BO surface (\(H C = S C \varepsilon\)).
/// For a collinear UKS spin channel the occupation is 1 and the force is
/// \(\mu\,\ddot C_\sigma = -(I - S C_\sigma C_\sigma^\top)\,(2 H_\sigma C_\sigma)\).
///
/// Holonomic orthonormality \(C^\top S C = I\) is enforced by **electronic
/// RATTLE** (Andersen, J. Comput. Phys. 52, 24 (1983)): a position SHAKE stage
/// with Lagrange multipliers \(C \leftarrow C_{\mathrm{raw}} + C_{\mathrm{ref}}\Lambda\),
/// then a velocity stage that kills \(\mathrm{sym}(C^\top S \dot C)\). This
/// supersedes the Löwdin-only foothold of `cpmd-water`.
///
/// Optional **Nosé–Hoover electronic thermostat** (Blöchl & Parrinello;
/// Sprik; Marx & Hutter) regulates the fictitious kinetic energy
/// \(T_e = \tfrac12\mu\lVert\dot C\rVert_F^2\) toward a target \(K_0\) with
/// coupling time \(\tau_e\). Symmetric NH half-steps wrap the mechanical
/// VV+RATTLE. Extended energy contribution
/// \(E_{\mathrm{therm}} = \tfrac12 Q\xi^2 + 2 K_0\eta\) with
/// \(Q = 2 K_0\tau_e^2\).
///
/// Unit-agnostic dense linear algebra on AO/MO matrices; consumed by
/// `cpu::gpw::GPWEvaluator` (`cpmd = true`) on closed-shell RKS or collinear
/// Γ-UKS (one `Cpmd` per spin channel). \(\mu\) is in atomic mass units
/// (electron-mass units); timesteps enter in fs and are converted to a.u.
/// time internally for the EOM.
#pragma once

#include "atomiverse/cpu/linalg/matrix.h"

namespace ATV {
namespace cpu {
namespace gpw {

/// @brief Car–Parrinello occupied-orbital propagator with electronic RATTLE
///        and optional Nosé–Hoover electronic thermostat.
///
/// Lifecycle: construct with \(\mu > 0\), optionally enable the thermostat
/// via public members, `seed` with S-orthonormal occupied coefficients
/// (typically from a cold SCF), then once per MD force step call
/// `advance(H, S, dt_fs)` after the KS matrix at the current density has
/// been built. Electronic kinetic energy uses the **velocity-RATTLE**
/// \(\dot C\) (after NH rescaling when thermostatted).
class Cpmd {
  public:
    /// @brief Construct with the fictitious electronic mass.
    /// @param mu_au Fictitious mass \(\mu\) in atomic units (electron mass
    ///        units). Typical CPMD values are \(10^2\)–\(10^3\). Default 400.
    /// @throws std::invalid_argument when mu_au ≤ 0.
    explicit Cpmd(double mu_au = 400.0);

    /// @brief Fictitious mass \(\mu\) [a.u.].
    double mu() const { return mu_; }

    /// @brief Whether occupied coefficients have been seeded.
    bool ready() const { return ready_; }

    /// @brief Enable Nosé–Hoover control of fictitious \(T_e\) (default false).
    bool thermostat = false;
    /// @brief Target electronic kinetic energy \(K_0\) [Hartree]. Must be
    ///        finite and > 0 when `thermostat` is true.
    double ekin_target_ha = 1e-5;
    /// @brief Electronic thermostat coupling time \(\tau_e\) [fs]. Must be
    ///        finite and > 0 when `thermostat` is true.
    double tau_e_fs = 50.0;

    /// @brief Current occupied MO coefficients \(C\) (\(n_b \times n_{\mathrm{occ}}\)).
    /// @throws std::logic_error before `seed`.
    const linalg::RealMatrix& C() const;

    /// @brief Current electronic velocity \(\dot C\) (\(n_b \times n_{\mathrm{occ}}\)).
    /// @throws std::logic_error before `seed`.
    const linalg::RealMatrix& Cdot() const;

    /// @brief Electronic kinetic energy \(T_e = \frac12\mu\sum_i\lVert\dot c_i\rVert^2\)
    ///        [Hartree] from the last advance (0 after seed, before any step).
    double kinetic_ha() const { return ke_ha_; }

    /// @brief Nosé–Hoover thermostat energy \(E_{\mathrm{therm}}\) [Hartree]
    ///        (\(0\) when thermostat is off).
    double thermostat_energy_ha() const { return e_therm_ha_; }

    /// @brief NH friction variable \(\xi\) [a.u. time\(^{-1}\)].
    double xi() const { return xi_; }

    /// @brief NH thermostat coordinate \(\eta\) (dimensionless log-scale).
    double eta() const { return eta_; }

    /// @brief Position-RATTLE iterations used by the last `advance` (0 at rest).
    int last_position_iters() const { return last_pos_iters_; }

    /// @brief \(\lVert C^\top S C - I\rVert_\infty\) after the last position stage
    ///        (0 after seed).
    double last_constraint_residual() const { return last_sigma_inf_; }

    /// @brief Seed the propagator at rest: \(C = C_{\mathrm{occ}}\), \(\dot C = 0\),
    ///        \(T_e = 0\). Resets NH state.
    /// @param C_occ Occupied coefficients [n_b × n_occ], S-orthonormal at the
    ///        seed geometry (\(C^\top S C = I\)).
    /// @throws std::invalid_argument when `C_occ` is empty.
    void seed(const linalg::RealMatrix& C_occ);

    /// @brief Seed with a prescribed electronic velocity (restart / tests).
    /// @param C_occ Occupied coefficients [n_b × n_occ].
    /// @param Cdot_occ Electronic velocity, same shape as `C_occ`.
    /// @throws std::invalid_argument on empty / shape mismatch.
    void seed(const linalg::RealMatrix& C_occ, const linalg::RealMatrix& Cdot_occ);

    /// @brief Re-orthonormalize \(C\) under a new metric \(S\) (nuclei moved)
    ///        and project \(\dot C\) onto the velocity constraint manifold.
    /// @param S Overlap at the current geometry, \(n_b \times n_b\).
    /// @throws std::logic_error before `seed`.
    void align_to_overlap(const linalg::RealMatrix& S);

    /// @brief One electronic velocity-Verlet step with RATTLE (pos + vel),
    ///        optionally wrapped by Nosé–Hoover half-steps.
    ///
    /// Uses the S-tangent force
    /// \(F = -(I-SC C^\top)(2\cdot\mathrm{occ}\,H C)\) [Hartree] and
    /// \(\Delta t\) converted from fs to a.u. time. Closed-shell RKS uses
    /// default `occ = 2` (factor 4); each UKS spin channel uses `occ = 1`.
    ///
    /// @param H KS matrix at the current density [Hartree], \(n_b \times n_b\).
    /// @param S Overlap matrix at the current geometry, \(n_b \times n_b\).
    /// @param dt_fs Timestep [fs] — must match the nuclear integrator dt.
    /// @param occ Orbital occupation (2 closed-shell, 1 per UKS spin).
    /// @throws std::logic_error before `seed`.
    /// @throws std::invalid_argument on non-positive dt/occ/thermostat params
    ///         or shape mismatch.
    /// @throws std::runtime_error if position RATTLE fails to converge.
    void advance(const linalg::RealMatrix& H, const linalg::RealMatrix& S,
                 double dt_fs, double occ = 2.0);

    /// @brief Drop electronic history (geometry discontinuity / restart).
    void reset();

  private:
    /// Validate thermostat members when engaged; throw if invalid.
    void check_thermostat_params() const;
    /// One NH half-step of length `dt_half_au` [a.u. time] using current Te.
    void nh_half_step(double dt_half_au);
    /// Recompute \(T_e\) and \(E_{\mathrm{therm}}\) from state.
    void update_energies();

    double mu_ = 400.0;          ///< fictitious mass [a.u.].
    linalg::RealMatrix C_;       ///< occupied coefficients
    linalg::RealMatrix Cdot_;    ///< electronic velocity
    double ke_ha_ = 0.0;         ///< last T_e [Hartree]
    double e_therm_ha_ = 0.0;    ///< last E_therm [Hartree]
    double xi_ = 0.0;            ///< NH friction
    double eta_ = 0.0;           ///< NH coordinate
    int last_pos_iters_ = 0;     ///< last position-RATTLE iteration count
    double last_sigma_inf_ = 0.0;///< last ‖σ‖_∞
    bool ready_ = false;
};

} // namespace gpw
} // namespace cpu
} // namespace ATV
