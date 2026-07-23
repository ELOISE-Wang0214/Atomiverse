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
    /// @brief FD displacement [Bohr] for forces (legacy path).
    double fd_step = 1e-3;
    /// @brief Grimme DFT-D3 (zero damping, PBE params).
    bool d3 = false;
    /// @brief Orbital Transformation (OT) SCF minimizer instead of
    ///        diagonalization + DIIS. Same gapped ground state; preferred for
    ///        larger AO bases (avoids dense \(O(n_b^3)\) eigensolve each step).
    ///        Units-free switch; energy still [Hartree] on `last_result()`.
    bool use_ot = false;
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
    /// @param d3 Enable DFT-D3.
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
    const GpwResult& last_result() const { return last_; }

  private:
    GpwResult last_{};
    std::unique_ptr<GpwDft> solver_;
    std::array<double, 3> built_cell_{};
};

} // namespace gpw
} // namespace cpu
} // namespace ATV
