/// @file cpu/scf/tdhf.h
/// @brief Excited-state result types of the TDHF/RPA and TDA (CIS) solvers
///        (public, result-only header — mp2.h precedent).
///
/// User path: `TDA<RHF>{mf}.run(n_states, channel)` after `mf.run`
/// (`pyscf.tdscf.rhf.TDA(mf).kernel`) and `RHF::tdhf(n_states, channel)`.
/// TDA is not a member of `RHF`. The eigenproblem machinery is private
/// (`internal::TDHFPack` over the link-08 `internal::OrbitalHessian` A/B
/// kernel) and is not part of this header.
///
/// References: Dreuw & Head-Gordon, Chem. Rev. 105, 4009 (2005),
/// DOI 10.1021/cr0505627 (TDHF/RPA response formalism); Hirata &
/// Head-Gordon, Chem. Phys. Lett. 314, 291 (1999),
/// DOI 10.1016/S0009-2614(99)01149-5 (TDA); Szabo & Ostlund §7.
#pragma once

#include <vector>

#include <xtensor/containers/xtensor.hpp>

namespace ATV {
namespace cpu {
namespace scf {

/// @brief Spin symmetry of the excitation out of a closed-shell reference
///        (dimensionless selector): `singlet` (spin-conserving, dipole-
///        allowed) or `triplet` (spin-flip manifold, dipole-forbidden from
///        a singlet ground state).
enum class ExcitationChannel { singlet, triplet };

/// @brief Excitation energies and length-gauge transition properties of a
///        TDA or TDHF/RPA solve on a converged closed-shell reference.
struct ExcitedStatesResult {
    /// @brief Excitation energies ω_0n [eV], ascending.
    std::vector<double> energy_ev;

    /// @brief Length-gauge transition dipole moments μ_0n [a.u., e·Bohr],
    ///        shape [n_states, 3]. Singlet channel:
    ///        μ_0n = √2 Σ_ia (X+Y)_ia ⟨φ_i|r|φ_a⟩ (overall sign is the
    ///        arbitrary eigenvector phase); triplet channel: exactly zero
    ///        (spin-forbidden from the singlet ground state).
    xt::xtensor<double, 2> transition_dipole_au;

    /// @brief Length-gauge oscillator strengths (dimensionless):
    ///        f_0n = (2/3)·ω_0n[Ha]·|μ_0n|²[a.u.]. Zero for the triplet
    ///        channel.
    std::vector<double> oscillator_strength;

    /// @brief Which spin channel was solved.
    ExcitationChannel channel = ExcitationChannel::singlet;

    /// @brief true = Tamm-Dancoff approximation (B = 0; CIS on an HF
    ///        reference), false = full TDHF/RPA (dimensionless flag).
    bool tda = false;

    /// @brief true when full MO-product ALDA \(f^{\mathrm{xc}}\) was applied
    ///        to the Casida A matrix (RKS LSDA TDA path).
    bool alda = false;
};

/// @brief Tamm–Dancoff (CIS on HF) excitations of a converged closed-shell
///        reference.
///
/// Mirrors `pyscf.tdscf.rhf.TDA(mf)`: `run` is `kernel`. Call `mf.run`
/// first. Not a method of `RHF` or `MeanFieldMethod`. RKS / UKS TDA stay
/// on those classes (KS linear response is a different kernel).
/// @tparam Method Concrete mean-field type. Instantiated for `RHF`.
template <typename Method>
class TDA {
  public:
    /// @brief Bind a mean-field reference. The reference must outlive `*this`.
    /// @param mf Mean-field method that has already been run (not copied).
    explicit TDA(const Method& mf);

    /// @brief Lowest `n_states` roots of the Hermitian A X = ωX problem.
    /// @param n_states Number of lowest excited states (1 ≤ n_states ≤
    ///                 n_occ·n_virt; dimensionless count).
    /// @param channel  Spin channel. Singlet roots carry length-gauge
    ///                 transition dipoles [a.u.] and oscillator strengths;
    ///                 triplet roots are spin-forbidden (μ = 0, f = 0).
    /// @return Excitation energies [eV], transition dipoles [a.u.],
    ///         oscillator strengths, `tda = true`.
    /// @throws std::runtime_error when called before `mf.run`.
    /// @throws std::invalid_argument when n_states is outside
    ///         [1, n_occ·n_virt] or the system has no virtual orbitals.
    ExcitedStatesResult
    run(int n_states,
        ExcitationChannel channel = ExcitationChannel::singlet) const;

  private:
    const Method* mf_ = nullptr;
};

} // namespace scf
} // namespace cpu
} // namespace ATV
