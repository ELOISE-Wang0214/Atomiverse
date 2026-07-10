/// @file cpu/semiempirical/nddo_solver.h
/// @brief Low-level MNDO/NDDO-family semiempirical solver scaffold.
///
/// This is the semiempirical analogue of `cpu::scf::SCFSolver` and
/// `cpu::pw::PWSolver`: an inner numerical solver owned by an Evaluator, not an
/// outer `Driver<>` subclass.
#pragma once

#include <string>
#include <vector>

#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/semiempirical/nddo_hcore.h"
#include "atomiverse/cpu/semiempirical/nddo_system.h"

namespace ATV {
namespace cpu {
namespace semiempirical {

enum class SemiEmpiricalMethod {
    AM1,
};

struct NDDOConfig {
    SemiEmpiricalMethod method = SemiEmpiricalMethod::AM1;
    int charge = 0;
    int multiplicity = 1;
    int max_iter = 200;
    double density_tol = 1.0e-8;
    double energy_tol = 1.0e-10;
};

struct NDDOResult {
    SemiEmpiricalMethod method = SemiEmpiricalMethod::AM1;
    bool converged = false;
    int iterations = 0;
    double electronic_energy_hartree = 0.0;
    double heat_of_formation_kcal_mol = 0.0;
    std::vector<double> mulliken_charges;
    std::string reference;
};

struct NDDOScfState {
    std::vector<double> initial_diagonal_density;
    std::vector<double> initial_packed_density;
    std::vector<double> initial_spin_density;
    std::vector<double> packed_hcore;
    std::vector<double> initial_one_center_fock;
    double nuclear_repulsion_ev = 0.0;
};

/// @brief Inner solver for MNDO-form NDDO methods.
///
/// `build()` validates and sizes method data for the current atom set.
/// `run()` executes the NDDO reference-shaped AM1 RHF loop for currently ported s/p
/// systems and throws a typed runtime error for unsupported cases.
class NDDOSolver {
  public:
    explicit NDDOSolver(NDDOConfig config);

    void build(const cpu::Frame& reference);
    NDDOResult& run(const cpu::Frame& frame);

    const NDDOResult& last_result() const { return result_; }
    const NDDOConfig& config() const { return config_; }
    const NDDOSystem& system() const { return system_; }
    const NDDOScfState& state() const { return state_; }
    bool built() const noexcept { return built_; }

  protected:
    void validate_supported_atoms_(const cpu::Frame& frame) const;

  private:
    NDDOConfig config_;
    bool built_ = false;
    int n_atoms_ = 0;
    NDDOSystem system_;
    NDDOScfState state_;
    NDDOResult result_;
};

} // namespace semiempirical
} // namespace cpu
} // namespace ATV
