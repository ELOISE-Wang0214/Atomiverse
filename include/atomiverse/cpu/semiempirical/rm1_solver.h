/// @file cpu/semiempirical/rm1_solver.h
/// @brief RM1 concrete method — public members + Frame build/run.
#pragma once

#include <memory>

#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/semiempirical/internal/nddo_solver.h"

namespace ATV {
namespace cpu {
namespace semiempirical {

class RM1Evaluator;

using RM1Result = NDDOResult;

/// @brief RM1 concrete semiempirical method (members + run).
class RM1 {
  public:
    int charge = 0;
    int multiplicity = 1;
    int max_iter = 200;
    double density_tol = 1.0e-8;
    double energy_tol = 1.0e-10;
    /// @brief Finite-difference displacement \(h\) [Å] of `forces()`;
    ///        truncation error \(O(h^{2})\).
    double force_step = 1.0e-4;

    /// @brief Default parameters (charge 0, singlet, max_iter 200, …).
    RM1();
    /// @brief Construct with SCF parameters assigned immediately.
    /// @param charge Total molecular charge [e].
    /// @param multiplicity Spin multiplicity 2S+1 (dimensionless).
    /// @param max_iter SCF iteration cap (dimensionless).
    /// @param density_tol Density residual tolerance (dimensionless).
    /// @param energy_tol Energy residual tolerance [eV].
    RM1(int charge, int multiplicity, int max_iter = 200,
        double density_tol = 1.0e-8, double energy_tol = 1.0e-10);
    ~RM1();
    RM1(RM1&&) noexcept;
    RM1& operator=(RM1&&) noexcept;
    RM1(const RM1&) = delete;
    RM1& operator=(const RM1&) = delete;

    void build(const cpu::Frame& reference);
    RM1Result& run(const cpu::Frame& frame);
    /// @brief `run(const Frame&)`, then writes the total energy [eV] to
    ///        `frame.meta[kEnergy]`.
    /// @param frame Geometry [Å]; energy out.
    /// @return `last_result()`.
    /// @throws std::runtime_error if the SCF does not converge in `max_iter`.
    RM1Result& run(cpu::Frame& frame);

    /// @brief SCF + nuclear forces \(F = -\partial E/\partial R\) [eV/Å] by
    ///        central finite difference of the energy at the converged
    ///        density (step `force_step`).
    ///
    /// See `NDDOSolver::compute_forces` for why holding the density fixed is
    /// exact here and what the difference step costs in accuracy. Writes
    /// `kAtom_fx/fy/fz` [eV/Å] (overwritten, not accumulated) and
    /// `frame.meta[kEnergy]` [eV].
    /// @param frame Geometry [Å]; forces and energy out.
    /// @return `last_result()` with `forces_ev_a` filled (3 · nlocal).
    /// @throws std::runtime_error if the SCF does not converge.
    RM1Result& forces(cpu::Frame& frame);

    const RM1Result& last_result() const;
    bool built() const noexcept;
    bool has_density() const noexcept;
    const NDDOSystem& system() const;
    const NDDOScfState& state() const;

  private:
    friend class RM1Evaluator;

    /// @brief Fill `forces_ev_a` at the converged density without touching
    ///        the Frame (the evaluator accumulates them itself).
    RM1Result& fd_forces_(cpu::Frame& frame);

    void sync_engine_();
    std::unique_ptr<NDDOSolver> engine_;
};

using RM1Solver = RM1;

} // namespace semiempirical
} // namespace cpu
} // namespace ATV
