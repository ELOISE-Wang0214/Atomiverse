/// @file cpu/gpw/gpw_evaluator.h
/// @brief GPWEvaluator — GPW DFT energy component for CPU drivers.
#pragma once

#include <array>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "atomiverse/cpu/component/gpw.h"
#include "atomiverse/cpu/gpw/gpw_solver.h"
#include "atomiverse/core/role.h"
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
    /// @brief Compute cell stress via orthorhombic strain-FD of SCF energy.
    ///        Writes Frame `kStress` [eV/Å³] Voigt; also `last_stress_ha_bohr3()`.
    bool compute_stress = false;
    /// @brief Dimensionless strain step for `compute_stress` (default 2e-3).
    double stress_strain_h = 2e-3;
    /// @brief FD displacement [Bohr] for forces (legacy path).
    double fd_step = 1e-3;
    /// @brief Grimme DFT-D3 (zero damping, PBE params). Mutually exclusive with d3bj.
    bool d3 = false;
    /// @brief Grimme DFT-D3(BJ) (PBE params, CP2K CN). Mutually exclusive with d3.
    bool d3bj = false;
    /// @brief Orbital Transformation (OT) SCF minimizer instead of
    ///        diagonalization + DIIS. Same gapped ground state; preferred for
    ///        larger AO bases (avoids dense \(O(n_b^3)\) eigensolve each step).
    ///        Units-free switch; energy still [Hartree] on `last_result()`.
    bool use_ot = false;
    /// UKS spin occupations. When both zero, closed-shell RKS from valence/2.
    /// O₂ GTH-q6 triplet: n_alpha=7, n_beta=5.
    int n_alpha = 0;
    int n_beta = 0;
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
                          int max_iter = 200, double tol = 1e-7, double mix = 0.3,
                          XcKind xc = XcKind::PADE, bool nonperiodic = false,
                          bool compute_forces = false, double fd_step = 1e-3,
                          bool d3 = false, bool use_ot = false)
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

  private:
    GpwResult last_{};
    std::unique_ptr<GpwDft> solver_;
    std::array<double, 3> built_cell_{};
    std::array<double, 6> last_stress_ha_{};
};

} // namespace gpw
} // namespace cpu
} // namespace ATV