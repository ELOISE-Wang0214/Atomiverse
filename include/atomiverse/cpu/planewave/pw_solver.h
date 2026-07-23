/// @file cpu/pw/pw_solver.h
/// @brief cpu::pw::PW — concrete plane-wave DFT method (CPU).
///
/// Parameters are public members (RHF pattern). Frame is the only I/O.
/// The self-consistency loop is owned by this method object (not a
/// public `PWConfig`/`PWSolver` bag). Workflow entry: `PWEvaluator` or
/// direct `build`/`run`.
///
/// @code
///   cpu::pw::PW pw;
///   pw.ecut_ry = 30.0;
///   pw.k_mesh = 4;
///   pw.build(frame);
///   const auto& result = pw.run(frame); // [Ry]/Ha on result
/// @endcode
#pragma once

#include <array>
#include <memory>
#include <string>
#include <vector>

#include "atomiverse/core/precision.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/scf/xc_functional.h"

namespace ATV {
namespace cpu {
namespace pw {

namespace internal {
class PoissonSolver;
class XCBuilder;
class VeffBuilder;
class DavidsonSolver;
class DensityBuilder;
class ConvergenceCheck;
} // namespace internal

/// @brief Metallic occupation smearing kind for plane-wave SCF.
///
/// Default `none` keeps fixed integer occupations (insulators / prior
/// regressions). When non-`none`, set `degauss > 0` [Ry] and enough
/// `extra_bands` for states above the Fermi level.
enum class Smearing : int {
    none = 0,              ///< Fixed closed-shell occupations.
    fermi_dirac = 1,       ///< Fermi–Dirac.
    methfessel_paxton = 2, ///< Methfessel–Paxton of order `mp_order`.
    gaussian = 3,          ///< Gaussian (cold/Marzari–Vanderbilt not here).
};

/// @brief Result of a plane-wave DFT total energy calculation.
struct PWResult {
    double energy = 0;      ///< Total free energy [Ry] (QE-comparable; \(E-TS\)).
    double energy_ha = 0;   ///< Total free energy [Hartree].
    bool converged = false; ///< Whether SCF converged within max_iter.
    int iterations = 0;     ///< Actual SCF iterations taken.
    /// @brief Collinear magnetization N↑ − N↓ [electrons] after SCF.
    ///        Zero when `spin_polarized == false`. Fixed occupations set this
    ///        exactly to `multiplicity - 1` for non-metallic fillings.
    double magnetization = 0;
    double e_fermi = 0;     ///< Fermi level [Ry] (0 if fixed occupations).
    double n_elec_occ = 0;  ///< \(\sum w f\) after last iteration.
    double e_entropy = 0;   ///< \(-TS\) contribution included in energy [Ry].
};

/// @brief Concrete plane-wave DFT method — members + `run(frame)`.
///
/// Construction → set members → `build(reference)` → `run(frame)`.
/// `last_result()` is valid after each `run()`. Geometry-dependent
/// structures rebuild inside each `run()` (AIMD / optimize safe).
class PW {
  public:
    /// @brief UPF v1 pseudopotential path (empty = toy V_loc).
    std::string upf_path;
    /// @brief Wavefunction cutoff [Ry].
    double ecut_ry = 30.0;
    /// @brief Monkhorst-Pack mesh size per dimension (dimensionless).
    int k_mesh = 4;
    /// @brief SCF tolerance on the density-change L2 norm (dimensionless).
    double tol = 1e-6;
    /// @brief Maximum SCF iterations (dimensionless).
    int max_iter = 50;
    /// @brief Density mixing parameter (QE mixing_beta, dimensionless).
    double mix_beta = 0.7;
    /// @brief Pulay mixing history window (dimensionless).
    int pulay_depth = 5;
    /// @brief Exchange–correlation functional (default true GGA PBE).
    ///
    /// `XCType::PBE` uses spectral density gradients (true GGA). Pure LDA is
    /// `XCType::LSDA`. Parameters are members — no public `*Config` bag.
    cpu::scf::XCType xc = cpu::scf::XCType::PBE;
    /// @brief Collinear spin-polarized KS (nspin=2). Default false = restricted.
    bool spin_polarized = false;
    /// @brief Spin multiplicity 2S+1 (used only when `spin_polarized`).
    ///        O₂ triplet ground state uses 3 → N↑−N↓ = 2.
    int multiplicity = 1;
    /// @brief Metallic smearing kind (default: fixed occupations).
    Smearing smearing = Smearing::none;
    /// @brief Smearing width [Ry] (QE `degauss`). Required > 0 when smearing is on.
    double degauss = 0.02;
    /// @brief Methfessel–Paxton order \(N\) (used when `smearing == methfessel_paxton`).
    int mp_order = 1;
    /// @brief Extra bands beyond \(\lceil n_e/2\rceil\) when smearing is on.
    int extra_bands = 4;

    /// @brief Default construct (helpers allocated on first build/run).
    PW();

    ~PW();
    PW(PW&&) noexcept;
    PW& operator=(PW&&) noexcept;
    PW(const PW&) = delete;
    PW& operator=(const PW&) = delete;

    /// @brief One-time setup: parses the UPF (if any) and validates Z.
    /// @param reference Frame used only to validate atom types.
    void build(const cpu::Frame& reference);

    /// @brief Plane-wave SCF at the frame's current geometry.
    /// @param frame Atom positions [Å], types (Z), and simulation cell.
    /// @return Owned `PWResult` (valid until next `run()`).
    PWResult& run(const cpu::Frame& frame);

    /// @brief Most recent SCF result.
    const PWResult& last_result() const { return result_; }

  private:
    bool built_ = false;
    std::unique_ptr<internal::PoissonSolver> poisson_;
    std::unique_ptr<internal::XCBuilder> xc_;
    std::unique_ptr<internal::VeffBuilder> veff_;
    std::unique_ptr<internal::DensityBuilder> density_;
    std::unique_ptr<internal::DavidsonSolver> davidson_;
    int davidson_bands_ = 0;
    PWResult result_;
};

/// @brief Finite-difference ionic forces [eV/Å] (oracle for analytic).
std::vector<std::array<double, 3>>
pw_forces_fd(const cpu::Frame& frame, const PW& method,
             double h_ang = 1e-3);

/// @brief Analytic ionic forces [eV/Å] (local HF + Ewald + KB nonlocal).
std::vector<std::array<double, 3>>
pw_forces_analytic(const cpu::Frame& frame, const PW& method);

/// @brief Cell stress σ [eV/Å³], Voigt order xx,yy,zz,xy,xz,yz.
std::array<double, 6> pw_stress(const cpu::Frame& frame,
                                const PW& method);

} // namespace pw
} // namespace cpu
} // namespace ATV