/// @file src/cpu/semiempirical/nddo_solver.cpp
/// @brief NDDO semiempirical inner solver scaffold.

#include "atomiverse/cpu/semiempirical/nddo_solver.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/cpu/semiempirical/am1_parameters.h"
#include "atomiverse/cpu/semiempirical/nddo_density.h"
#include "atomiverse/cpu/semiempirical/nddo_fock.h"
#include "atomiverse/cpu/semiempirical/nddo_hcore.h"
#include "atomiverse/cpu/semiempirical/nddo_matrix.h"
#include "atomiverse/cpu/semiempirical/nddo_system.h"

namespace ATV {
namespace cpu {
namespace semiempirical {

namespace {

std::vector<double>
packed_density_from_diagonal(const std::vector<double>& diagonal) {
    std::vector<double> packed(diagonal.size() * (diagonal.size() + 1u) / 2u,
                               0.0);
    for (std::size_t i = 0; i < diagonal.size(); ++i) {
        packed[packed_index(static_cast<int>(i), static_cast<int>(i))] =
            diagonal[i];
    }
    return packed;
}

std::vector<double> scaled_density(const std::vector<double>& density,
                                   double scale) {
    auto out = density;
    for (auto& value : out) {
        value *= scale;
    }
    return out;
}

double max_abs_delta(const std::vector<double>& a,
                     const std::vector<double>& b) {
    if (a.size() != b.size()) {
        throw std::runtime_error("NDDO density delta size mismatch");
    }
    double out = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        out = std::max(out, std::abs(a[i] - b[i]));
    }
    return out;
}

} // namespace

NDDOSolver::NDDOSolver(NDDOConfig config) : config_(config) {
    result_.method = config_.method;
    result_.reference = "Internal AM1/NDDO reference data";
}

void NDDOSolver::build(const cpu::Frame& reference) {
    validate_supported_atoms_(reference);
    system_ = build_am1_system(reference, config_.charge, config_.multiplicity);
    state_.initial_diagonal_density = initial_diagonal_density_guess(system_);
    state_.initial_packed_density =
        packed_density_from_diagonal(state_.initial_diagonal_density);
    state_.initial_spin_density =
        scaled_density(state_.initial_packed_density, 0.5);
    const auto hcore = NDDOHCoreBuilder{}.build(system_, reference);
    state_.packed_hcore = hcore.packed_hcore;
    state_.nuclear_repulsion_ev = hcore.nuclear_repulsion_ev;
    state_.initial_one_center_fock = NDDOFockBuilder{}.build(
        system_, reference, state_.packed_hcore, state_.initial_packed_density,
        state_.initial_spin_density);

    n_atoms_ = reference.atoms.size();
    result_.converged = false;
    result_.iterations = 0;
    result_.electronic_energy_hartree = 0.0;
    result_.heat_of_formation_kcal_mol = 0.0;
    result_.mulliken_charges.assign(static_cast<std::size_t>(n_atoms_), 0.0);
    built_ = true;
}

NDDOResult& NDDOSolver::run(const cpu::Frame& frame) {
    if (!built_) {
        throw std::runtime_error("NDDOSolver::run called before build()");
    }
    if (frame.atoms.size() != n_atoms_) {
        throw std::runtime_error(
            "NDDOSolver::run requires a fixed atom count after build()");
    }
    validate_supported_atoms_(frame);

    if (system_.n_alpha != system_.n_beta) {
        throw std::runtime_error(
            "NDDOSolver AM1 SCF currently supports closed-shell RHF only");
    }

    const auto hcore = NDDOHCoreBuilder{}.build(system_, frame);
    std::vector<double> density = state_.initial_packed_density;
    std::vector<double> spin_density = scaled_density(density, 0.5);
    std::vector<double> fock = NDDOFockBuilder{}.build(
        system_, frame, hcore.packed_hcore, density, spin_density);
    double previous_energy_ev = std::numeric_limits<double>::infinity();

    for (int iter = 1; iter <= config_.max_iter; ++iter) {
        const auto orbitals =
            NDDODiagonalizer{}.diagonalize_packed(fock, system_.n_orbitals);

        NDDODensityInput density_input;
        density_input.mdim = system_.n_orbitals;
        density_input.norbs = system_.n_orbitals;
        density_input.n_occupied = system_.n_beta;
        density_input.n_fractional = system_.n_beta;
        density_input.full_occupancy = 2.0;
        density_input.fractional_occupancy = 0.0;
        const auto next_density = NDDODensityBuilder{}.build_packed(
            orbitals.coefficients_row_major, density_input);
        const double delta = max_abs_delta(density, next_density);

        density = next_density;
        spin_density = scaled_density(density, 0.5);
        fock = NDDOFockBuilder{}.build(system_, frame, hcore.packed_hcore,
                                       density, spin_density);
        const double electronic_ev = electronic_energy_ev(
            system_.n_orbitals, density, hcore.packed_hcore, fock);
        const double total_ev = electronic_ev + hcore.nuclear_repulsion_ev;
        const double energy_delta = std::abs(total_ev - previous_energy_ev);

        result_.iterations = iter;
        result_.electronic_energy_hartree = total_ev / 27.211386245988;
        result_.mulliken_charges =
            atomic_charges_from_packed_density(system_, density);
        if (delta <= config_.density_tol
            || energy_delta <= config_.energy_tol) {
            result_.converged = true;
            result_.heat_of_formation_kcal_mol = 0.0;
            return result_;
        }
        previous_energy_ev = total_ev;
    }

    result_.converged = false;
    throw std::runtime_error("NDDOSolver AM1 SCF did not converge");
}

void NDDOSolver::validate_supported_atoms_(const cpu::Frame& frame) const {
    if (config_.method != SemiEmpiricalMethod::AM1) {
        throw std::runtime_error(
            "NDDOSolver: unsupported semiempirical method");
    }
    const auto& atoms = frame.atoms;
    const int* __restrict sp_type = frame.span(kAtom_type);
    for (int i = 0; i < atoms.size(); ++i) {
        const int z = sp_type[i];
        if (!am1_has_parameter(z)) {
            throw std::runtime_error(
                "AM1 parameters unavailable for atomic number "
                + std::to_string(z));
        }
    }
}

} // namespace semiempirical
} // namespace cpu
} // namespace ATV
