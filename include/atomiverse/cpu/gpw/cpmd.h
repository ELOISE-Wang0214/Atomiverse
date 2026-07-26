/// @file cpu/gpw/cpmd.h
/// @brief Cpmd / CpmdCplx — Car–Parrinello fictitious electronic DOF
///        propagators (RATTLE + optional Nosé–Hoover electronic thermostat).
///
/// Propagates occupied MO coefficients with velocity-Verlet and a fictitious
/// electron mass \(\mu\) (Car & Parrinello, Phys. Rev. Lett. 55, 2471 (1985)).
///
/// **Real-valued Γ (`Cpmd`):** \(C\) is \(n_b\times n_{\mathrm{occ}}\) real;
/// \(P = 2 C C^\top\) (RKS) or \(P_\sigma = C_\sigma C_\sigma^\top\) (UKS);
/// force \(\mu\,\ddot C = -(I - S C C^\top)\,(2\cdot\mathrm{occ}\,H C)\).
/// Holonomic orthonormality \(C^\top S C = I\) via electronic RATTLE
/// (Andersen, J. Comput. Phys. 52, 24 (1983)).
///
/// **Complex multi-k (`CpmdCplx`):** one propagator per sampled Bloch \(k\);
/// \(C_k\) is complex; constraint \(C_k^\dagger S_k C_k = I\); velocity
/// RATTLE kills \(\mathrm{Herm}(C_k^\dagger S_k \dot C_k)\). The force uses
/// the complex S-tangent gradient (`ot_gradient_k`). Multi-k \(T_e\) is
/// k-weight-summed by the evaluator: \(\sum_k w_k T_e^{(k)}\).
///
/// Optional **Nosé–Hoover electronic thermostat** (Blöchl & Parrinello;
/// Sprik; Marx & Hutter) is product on real `Cpmd` only for this pin
/// (`CpmdCplx` keeps the members for API symmetry but multi-k NH is not
/// regression-pinned).
///
/// Consumed by `cpu::gpw::GPWEvaluator` (`cpmd = true`): Γ RKS/UKS via
/// `Cpmd`; multi-k RKS via a vector of `CpmdCplx`.
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

/// @brief Car–Parrinello propagator on complex Bloch occupied coefficients.
///
/// Hermitian generalization of `Cpmd` for multi-k RKS (one instance per
/// sampled \(k\)). Lifecycle mirrors `Cpmd`: `seed` after a cold multi-k
/// SCF, then `advance(H_k, S_k, dt_fs)` each electronic step. Thermostat
/// members exist for API symmetry; product multi-k CPMD uses NVE electrons
/// (`thermostat = false`).
class CpmdCplx {
  public:
    /// @brief Construct with the fictitious electronic mass.
    /// @param mu_au Fictitious mass \(\mu\) [a.u. electron mass]. Default 400.
    /// @throws std::invalid_argument when mu_au ≤ 0.
    explicit CpmdCplx(double mu_au = 400.0);

    double mu() const { return mu_; }
    bool ready() const { return ready_; }

    bool thermostat = false;
    double ekin_target_ha = 1e-5;
    double tau_e_fs = 50.0;

    /// @brief Occupied Bloch coefficients \(C_k\) (\(n_b \times n_{\mathrm{occ}}\)).
    /// @throws std::logic_error before `seed`.
    const linalg::CplxMatrix& C() const;

    /// @brief Electronic velocity \(\dot C_k\).
    /// @throws std::logic_error before `seed`.
    const linalg::CplxMatrix& Cdot() const;

    double kinetic_ha() const { return ke_ha_; }
    double thermostat_energy_ha() const { return e_therm_ha_; }
    double xi() const { return xi_; }
    double eta() const { return eta_; }
    int last_position_iters() const { return last_pos_iters_; }
    double last_constraint_residual() const { return last_sigma_inf_; }

    /// @brief Seed at rest: \(C = C_{\mathrm{occ}}\), \(\dot C = 0\).
    /// @param C_occ Occupied Bloch block [nb × n_occ], \(C^\dagger S C = I\).
    /// @throws std::invalid_argument when empty.
    void seed(const linalg::CplxMatrix& C_occ);

    /// @brief Seed with a prescribed electronic velocity.
    void seed(const linalg::CplxMatrix& C_occ, const linalg::CplxMatrix& Cdot_occ);

    /// @brief Löwdin-transport \(C\) under a new \(S_k\) and velocity-project.
    /// @param S Hermitian overlap at the current geometry.
    void align_to_overlap(const linalg::CplxMatrix& S);

    /// @brief One electronic VV + RATTLE step on complex Hermitian \(H,S\).
    /// @param H KS matrix at this \(k\) [Hartree], nb × nb Hermitian.
    /// @param S Overlap at this \(k\), nb × nb Hermitian.
    /// @param dt_fs Timestep [fs] (match nuclear VV).
    /// @param occ Orbital occupation (2 closed-shell RKS).
    void advance(const linalg::CplxMatrix& H, const linalg::CplxMatrix& S,
                 double dt_fs, double occ = 2.0);

    void reset();

  private:
    void check_thermostat_params() const;
    void nh_half_step(double dt_half_au);
    void update_energies();

    double mu_ = 400.0;
    linalg::CplxMatrix C_;
    linalg::CplxMatrix Cdot_;
    double ke_ha_ = 0.0;
    double e_therm_ha_ = 0.0;
    double xi_ = 0.0;
    double eta_ = 0.0;
    int last_pos_iters_ = 0;
    double last_sigma_inf_ = 0.0;
    bool ready_ = false;
};

} // namespace gpw
} // namespace cpu
} // namespace ATV
