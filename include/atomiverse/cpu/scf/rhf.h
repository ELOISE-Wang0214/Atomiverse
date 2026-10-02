/// @file cpu/scf/rhf.h
/// @brief Restricted Hartree–Fock (RHF) concrete mean-field method.
#pragma once

#include <string>
#include <vector>

#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/scf/adc.h"
#include "atomiverse/cpu/scf/bse.h"
#include "atomiverse/cpu/scf/eom_cc2.h"
#include "atomiverse/cpu/scf/eom_ea_ccsd.h"
#include "atomiverse/cpu/scf/eom_ip_ccsd.h"
#include "atomiverse/cpu/scf/gf2.h"
#include "atomiverse/cpu/scf/mean_field_method.h"
#include "atomiverse/cpu/scf/rpa.h"

namespace ATV {
namespace cpu {
namespace scf {

/// @brief Restricted Hartree–Fock — a mean-field electronic-structure method
///        (a single Slater determinant in which every electron moves in the
///        averaged field of all the others) whose spin-up (α) and spin-down (β)
///        electrons share one common set of spatial orbitals: the closed-shell,
///        fully paired case.
class RHF : public MeanFieldMethod<RHF> {
  public:
    /// @brief Construct with the RHF defaults (STO-3G; restricted HF; MINAO;
    ///        no density damping; energy tolerance 1e-9 Ha).
    ///
    /// Prefer the basis ctor: `RHF{BasisSet::STO3G}` — **do not** write
    /// `rhf.basis = …` after default construction.
    RHF();

    /// @brief Construct with an explicit basis (primary public form).
    /// @param basis_in Basis set (e.g. `BasisSet::STO3G`).
    explicit RHF(BasisSet basis_in);

    /// @brief Construct with basis and molecular charge.
    /// @param basis_in Basis set.
    /// @param charge_in Molecular charge [e] (electrons = ΣZ − charge).
    RHF(BasisSet basis_in, int charge_in);

    /// @brief Construct with basis, charge, SCF tolerance, and iteration cap.
    /// @param basis_in Basis set.
    /// @param charge_in Molecular charge [e].
    /// @param tol_in Energy / residual tolerance [Ha] (see base `tol`).
    /// @param max_iter_in Maximum SCF iterations [dimensionless].
    RHF(BasisSet basis_in, int charge_in, double tol_in, int max_iter_in = 100);

    std::string name() const { return "rhf"; }

    /// @brief Γ-point periodic RHF. Default false keeps the molecular
    ///        engine bit-identical. When true, Frame.box must be
    ///        orthorhombic; J, K, S, T, V and E_nn share one truncated
    ///        lattice-image window (no Ewald). A cell larger than the
    ///        Gaussian pair reach has only T = 0 and recovers molecular RHF.
    bool periodic = false;

    /// @brief Strict ADC(2) electronic excitation energies of the most recent
    ///        `run()` — Schirmer, Phys. Rev. A **26**, 2395 (1982); Trofimov &
    ///        Schirmer, J. Phys. B **28**, 2299 (1995). Closed-shell singlet
    ///        manifold; matches PySCF `adc.ADC` with `method='adc(2)'` and
    ///        `method_type='ee'`. Builds the second-order singles intermediate
    ///        M_ss, first-order singles–doubles coupling, and zeroth-order
    ///        doubles diagonal; solves the non-symmetric ADC eigenvalue
    ///        problem with a matvec Davidson iteration (no ADC(2)-x or
    ///        ADC(3) terms).
    /// @param n_states Number of lowest ADC(2) roots (1 ≤ n_states ≤
    ///                 n_occ·n_virt; dimensionless count).
    /// @return Excitation energies [eV] and the MP2 correlation energy
    ///         [Hartree] of the RHF reference.
    /// @throws std::runtime_error when called before `run()`.
    /// @throws std::invalid_argument when n_states is outside
    ///         [1, n_occ·n_virt] or the system has no virtual orbitals.
    ADC2Result adc2(int n_states = 1) const;

    /// @brief Extended ADC(2)-x electronic excitation energies of the most
    ///        recent `run()` (closed-shell RHF). Same singles intermediate and
    ///        singles–doubles coupling as strict ADC(2), but the doubles block
    ///        includes first-order two-electron terms (vvvv, oovv, ovvo, oooo)
    ///        — PySCF `adc.ADC` / `method='adc(2)-x'` / `method_type='ee'`.
    /// @param n_states Number of lowest ADC(2)-x roots (1 ≤ n_states ≤
    ///                 n_occ·n_virt; dimensionless count).
    /// @return Excitation energies [eV] and the MP2 correlation energy
    ///         [Hartree] of the RHF reference.
    /// @throws std::runtime_error when called before `run()`.
    /// @throws std::invalid_argument when n_states is outside
    ///         [1, n_occ·n_virt] or the system has no virtual orbitals.
    ADC2Result adc2x(int n_states = 1) const;

    /// @brief Strict ADC(2) ionization potentials of the most recent
    ///        `run()` — Schirmer, Phys. Rev. A **26**, 2395 (1982).
    ///        Closed-shell 1h / 2h1p intermediate-state representation;
    ///        matches PySCF `adc.ADC` with `method='adc(2)'` and
    ///        `method_type='ip'`. Second-order M_ij (occupied–occupied),
    ///        first-order 1h–2h1p coupling via (ia|jk), zeroth-order 2h1p
    ///        diagonal D_ajk = ε_j + ε_k − ε_a; eigenvalues of the
    ///        sign-flipped ISR matrix are the ionization potentials.
    /// @param n_states Number of lowest IP roots (1 ≤ n_states ≤ n_occ).
    /// @return Ionization potentials [eV] and the MP2 correlation energy
    ///         [Hartree] of the RHF reference.
    /// @throws std::runtime_error when called before `run()`.
    /// @throws std::invalid_argument when n_states is outside [1, n_occ]
    ///         or the system has no virtual orbitals.
    IPADC2Result ip_adc2(int n_states = 1) const;

    /// @brief Strict EA-ADC(2) electron-attachment energies of the most recent
    ///        `run()` — Schirmer, Phys. Rev. A **26**, 2395 (1982); matches
    ///        PySCF `adc.ADC` with `method='adc(2)'` and `method_type='ea'`.
    ///        Singles = virtual (1p) space with second-order \(M_{ab}\);
    ///        doubles = 2p–1h with zeroth-order diagonal; first-order
    ///        singles–doubles coupling via \((ia|bc)\). No ADC(2)-x / ADC(3).
    /// @param n_states Number of lowest EA roots (1 ≤ n_states ≤
    ///                 n_virt + n_occ·n_virt²; dimensionless count).
    /// @return Attachment energies ω [eV] and the MP2 correlation energy
    ///         [Hartree] of the N-electron RHF reference.
    /// @throws std::runtime_error when called before `run()`.
    /// @throws std::invalid_argument when n_states is out of range or the
    ///         system has no virtual orbitals.
    EAADC2Result ea_adc2(int n_states = 1) const;

    /// @brief EOM-EA-CCSD electron-attachment energies — drives a fresh RHF
    ///        (materializes ERI), spin-orbital Gauss–Stanton CCSD, then the
    ///        non-symmetric EOM-EA Jacobian (Nooijen & Bartlett, J. Chem.
    ///        Phys. **102**, 3629 (1995)). Matches PySCF
    ///        `cc.CCSD(...).eaccsd()` / `eom_gccsd.EOMEA`.
    /// @param frame Geometry; drives SCF + CCSD + EOM.
    /// @param n_states Number of lowest unique attachment roots (≥ 1).
    /// @return Attachment ω [eV] plus CCSD energetics [Hartree].
    /// @throws std::invalid_argument on open-shell / missing virtuals /
    ///         non-converged SCF or CCSD / n_states out of range.
    EOMEACCSDResult eom_ea_ccsd(Frame& frame, int n_states = 1);

    /// @brief EOM-IP-CCSD ionization potentials — drives a fresh RHF
    ///        (materializes ERI), spin-orbital Gauss–Stanton CCSD, then the
    ///        non-symmetric EOM-IP Jacobian (Tu, Wang & Li, J. Chem. Phys.
    ///        **136**, 174102 (2012)). Matches PySCF
    ///        `cc.CCSD(...).ipccsd()` / `eom_gccsd.EOMIP`.
    /// @param frame Geometry; drives SCF + CCSD + EOM.
    /// @param n_states Number of lowest unique ionization roots (≥ 1).
    /// @return Ionization ω [eV] plus CCSD energetics [Hartree].
    /// @throws std::invalid_argument on open-shell / missing occupieds /
    ///         non-converged SCF or CCSD / n_states out of range.
    EOMIPCCSDResult eom_ip_ccsd(Frame& frame, int n_states = 1);

    /// @brief Direct RPA (dRPA) correlation energy of the most recent `run()`
    ///        (closed-shell RHF reference). Coulomb-only Casida A/B
    ///        (α_hf = 0); correlation via the Furche plasmon formula
    ///        \(E_c = \frac12(\sum_s\Omega_s - \mathrm{Tr} A)\). Exact
    ///        (non-DF) twin of PySCF `gw.rpa.RPA` / `DirectRPA` in the
    ///        complete-auxiliary limit; stock PySCF uses density fitting.
    /// @return Correlation and total energy [Hartree], plus the dRPA
    ///         poles \(\Omega_s\) [Hartree].
    /// @throws std::runtime_error when called before `run()`, on a
    ///         non-restricted reference, or when dRPA is ill-conditioned.
    /// @throws std::invalid_argument when n_occ = 0 or n_virt = 0.
    DirectRPAResult direct_rpa() const;

    /// @brief Singlet BSE@G0W0 excitation energies of the most recent
    ///        `run()` (closed-shell RHF). Spectral G0W0 poles
    ///        (`G0W0<RHF>::run`) enter the ph energy differences and the static
    ///        direct-RPA screened interaction \(W\) (exact 4-centre form).
    ///        Full BSE (default) solves the Casida A/B RPA problem; with
    ///        `tda = true` only A is diagonalized (BSE-TDA).
    ///        Working equations match PySCF `pyscf.gw.bse.BSE` (non-DF ERIs
    ///        + frequency-exact G0W0 poles — see `reg_bse-h2-sto3g`).
    /// @param n_states Number of lowest singlet roots in `[1, n_occ·n_virt]`.
    /// @param tda If true, Tamm–Dancoff (B = 0); else full BSE.
    /// @return Excitation energies [eV] plus the underlying G0W0 result.
    /// @throws std::runtime_error when called before `run()`, on a
    ///         non-restricted reference, or when G0W0 / BSE is ill-conditioned.
    /// @throws std::invalid_argument when `n_states` is out of range or the
    ///         system has no virtual orbitals.
    BSEResult bse(int n_states = 1, bool tda = false) const;

    /// @brief One-shot second-order Green's function (G0F2 / GF2)
    ///        quasiparticle energies of the most recent `run()` (closed-shell
    ///        RHF). Builds the 2h1p+2p1h auxiliary self-energy from the HF
    ///        Green's function and solves the Dyson equation once with
    ///        F = diag(ε_HF). Matches PySCF `agf2.RAGF2` one-shot poles
    ///        (SE from `init_gf` → `se.get_greens_function(F_HF)`). Not
    ///        self-consistent AGF2 (separate roadmap item).
    /// @return QP energies [Hartree], HOMO–LUMO gap / IP / EA [eV],
    ///         Galitskii–Migdal e_corr [Hartree].
    /// @throws std::runtime_error when called before `run()` or on a
    ///         non-restricted reference.
    /// @throws std::invalid_argument when n_occ = 0 or n_virt = 0.
    GF2Result gf2() const;

    /// @brief Singlet EOM-CC2 excitation energies — drives a fresh RHF
    ///        (materializes ERI), spin-orbital Christiansen–Koch–Jørgensen
    ///        CC2 ground state, then the same non-symmetric EOM-EE Jacobian
    ///        evaluated at the CC2 amplitudes (matches PySCF
    ///        `cc.RCCSD` with `cc2=True` followed by `eomee_ccsd_singlet`).
    /// @param frame Geometry; drives SCF + CC2 + EOM.
    /// @param n_states Number of lowest singlet roots (≥ 1).
    /// @return Singlet ω [eV] plus CC2 energetics [Hartree].
    /// @throws std::invalid_argument on open-shell / missing virtuals /
    ///         non-converged SCF or CC2 / n_states out of range.
    EOMCC2Result eom_cc2(Frame& frame, int n_states = 1);

  protected:
    /// @brief Closed-shell occupation: n_alpha = n_beta = n_electrons / 2
    ///        (every occupied spatial orbital holds one α and one β electron).
    void resolve_occupation(int n_electrons) override;
};

} // namespace scf
} // namespace cpu
} // namespace ATV
