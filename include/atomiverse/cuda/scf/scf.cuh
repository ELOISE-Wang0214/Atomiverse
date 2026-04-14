/// @file cuda/scf/scf.cuh
/// @brief scf::SCF -- self-consistent field driver.
///
/// Runs an iterative SCF loop through the component lifecycle:
///   build_hamiltonian -> diagonalize -> build_density -> convergence_check -> end_of_step
///
/// Each "step" is one SCF iteration. The user passes max iterations to run().
///
/// Auto-adds ZeroDensity at build_hamiltonian (order -1000) to initialize
/// the density matrix to zero for the first iteration.
///
/// Usage:
/// @code{.cpp}
/// scf::SCF solver;
/// solver.add(scf::Hamiltonian{n_occ, alpha});
/// solver.add(scf::EigenSolver{});
/// solver.add(scf::DensityBuilder{});
/// solver.add(scf::ConvergenceCheck{1e-8});
/// solver.on(scf::Phase::end_of_step).add(Thermo{1});
/// solver.build();
/// solver.run(frame, 100);  // max 100 SCF iterations
/// @endcode
///
/// @see scf::Phase for the SCF phase sequence
/// @see Driver for the inherited fluent API
#pragma once

#include "atomiverse/facade/driver.h"
#include "atomiverse/facade/stage.h"
#include "atomiverse/cuda/scf/hamiltonian.cuh"
#include "atomiverse/cuda/scf/eigen_solver.cuh"
#include "atomiverse/cuda/scf/density_builder.cuh"
#include "atomiverse/cuda/scf/convergence_check.cuh"

namespace ATV {
namespace scf {

/// @brief Zero the density matrix P via cudaMemsetAsync.
///
/// Infrastructure component auto-added by scf::SCF at build_hamiltonian (order -1000).
/// Ensures P starts at zero for the first SCF iteration.
class ZeroDensity : public ComponentBase {
public:
    /// @brief Component name.
    /// @return "zero_density".
    std::string name() const override { return "zero_density"; }

    /// @brief Zero P via cudaMemsetAsync. Only acts on step 0.
    /// @param ctx Launch context.
    void launch(const LaunchContext& ctx) override;
};

/// @brief SCF specification with iterative dense-matrix phases.
///
/// Auto-adds ZeroDensity at build_hamiltonian (order -1000).
class Spec : public ATV::DriverSpec<Phase> {
    using Base = ATV::DriverSpec<Phase>;
public:
    /// @brief Construct SCF spec with auto-added infrastructure.
    Spec() {
        Base::on(Phase::build_hamiltonian).add(ZeroDensity{}, -1000);
    }
};

/// @brief Pure runtime CUDA SCF driver executing a compiled plan.
class Driver : public ATV::Driver<Phase> {
    using Base = ATV::Driver<Phase>;
public:
    using Base::Base;
};

/// @brief Transitional legacy CUDA SCF facade.
class SCF : public ATV::detail::LegacyDriverFacade<Spec, Driver> {};

} // namespace scf
} // namespace ATV
