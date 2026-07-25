/// @file cpu/gpw/cpmd.h
/// @brief Cpmd — Car–Parrinello fictitious electronic DOF propagator.
///
/// Propagates occupied MO coefficients \(C\) (\(n_b \times n_{\mathrm{occ}}\))
/// with a position-Verlet step and a fictitious electron mass \(\mu\)
/// (Car & Parrinello, Phys. Rev. Lett. 55, 2471 (1985)). The electronic force
/// is the **S-tangent** projection of the closed-shell band-energy gradient
/// (\(P = 2 C C^\top\)):
///
///   \(\mu\,\ddot C = -(I - S C C^\top)\,(4 H C)\)
///
/// so \(\ddot C \to 0\) on the BO surface (\(H C = S C \varepsilon\)). A
/// Löwdin re-orthonormalization \(C \leftarrow C (C^\top S C)^{-1/2}\) cleans
/// residual metric drift (SHAKE-free foothold; full RATTLE is a follow-on).
///
/// Unit-agnostic dense linear algebra on AO/MO matrices; consumed by
/// `cpu::gpw::GPWEvaluator` (`cpmd = true`) on the closed-shell RKS Γ path.
/// \(\mu\) is in atomic mass units (electron-mass units); timesteps enter in
/// fs and are converted to a.u. time internally for the EOM.
#pragma once

#include "atomiverse/cpu/linalg/matrix.h"

namespace ATV {
namespace cpu {
namespace gpw {

/// @brief Car–Parrinello occupied-orbital propagator.
///
/// Lifecycle: construct with \(\mu > 0\), `seed` with S-orthonormal occupied
/// coefficients (typically from a cold SCF), then once per MD force step
/// call `advance(H, S, dt_fs)` after the KS matrix at the current density
/// has been built. Electronic kinetic energy is estimated from the central
/// difference \(\dot C \approx (C_{n+1} - C_{n-1}) / (2\Delta t)\).
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

    /// @brief Current occupied MO coefficients \(C\) (\(n_b \times n_{\mathrm{occ}}\)).
    /// @throws std::logic_error before `seed`.
    const linalg::RealMatrix& C() const;

    /// @brief Electronic kinetic energy \(T_e = \frac12\mu\sum_i\lVert\dot c_i\rVert^2\)
    ///        [Hartree] from the last advance (0 after seed, before any step).
    double kinetic_ha() const { return ke_ha_; }

    /// @brief Seed the propagator at rest: \(C_{-1} = C_0 = C_{\mathrm{occ}}\),
    ///        \(T_e = 0\).
    /// @param C_occ Occupied coefficients [n_b × n_occ], S-orthonormal at the
    ///        seed geometry (\(C^\top S C = I\)).
    /// @throws std::invalid_argument when `C_occ` is empty.
    void seed(const linalg::RealMatrix& C_occ);

    /// @brief Re-orthonormalize \(C\) and \(C_{\mathrm{prev}}\) under a new
    ///        metric \(S\) (call when nuclei have moved before the KS eval).
    /// @param S Overlap at the current geometry, \(n_b \times n_b\).
    /// @throws std::logic_error before `seed`.
    void align_to_overlap(const linalg::RealMatrix& S);

    /// @brief One electronic position-Verlet step with Löwdin reorth.
    ///
    /// Uses the S-tangent force \(F = -(I-SC C^\top)(4 H C)\) [Hartree] and
    /// \(\Delta t\) converted from fs to a.u. time:
    /// \(C_{\mathrm{raw}} = 2C - C_{\mathrm{prev}} + (\Delta t_{\mathrm{au}}^2/\mu) F\),
    /// then \(C \leftarrow\) Löwdin\((C_{\mathrm{raw}}, S)\).
    ///
    /// @param H KS matrix at the current density [Hartree], \(n_b \times n_b\).
    /// @param S Overlap matrix at the current geometry, \(n_b \times n_b\).
    /// @param dt_fs Timestep [fs] — must match the nuclear integrator dt.
    /// @throws std::logic_error before `seed`.
    /// @throws std::invalid_argument on non-positive dt or shape mismatch.
    void advance(const linalg::RealMatrix& H, const linalg::RealMatrix& S,
                 double dt_fs);

    /// @brief Drop electronic history (geometry discontinuity / restart).
    void reset();

  private:
    double mu_ = 400.0;           ///< fictitious mass [a.u.].
    linalg::RealMatrix C_;        ///< C_n
    linalg::RealMatrix C_prev_;   ///< C_{n-1}
    double ke_ha_ = 0.0;          ///< last T_e [Hartree]
    bool ready_ = false;
};

} // namespace gpw
} // namespace cpu
} // namespace ATV
