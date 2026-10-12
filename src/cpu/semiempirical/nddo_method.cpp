/// @file src/cpu/semiempirical/nddo_method.cpp
/// @brief NDDO Driver loop: diagonalize → density → fock → check.

#include "atomiverse/cpu/semiempirical/nddo_method.h"

#include <stdexcept>
#include <utility>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/cpu/semiempirical/internal/nddo_solver.h"
#include "atomiverse/cpu/semiempirical/stages/nddo_stages.h"

namespace ATV {
namespace cpu {
namespace semiempirical {

NddoMethod::NddoMethod() = default;

NddoMethod::~NddoMethod() { reset_built_driver_(); }

void NddoMethod::reset_built_driver_() {
    diagonalize_.reset();
    density_.reset();
    fock_.reset();
    check_.reset();
    converged_.reset();
    this->pending_.clear();
    this->owned_components_.clear();
    this->plan_ = {};
    this->criteria_.clear();
    this->criteria_owned_.clear();
    this->backend_.reset();
}

NDDOResult& NddoMethod::drive_(NDDOSolver& solver, const cpu::Frame& frame) {
    reset_built_driver_();
    solver.scf_prepare(frame);
    loop_converged_ = false;

    diagonalize_ = std::make_unique<NddoStep>(
        ATV::nddo::Phase::diagonalize, "nddo_diagonalize",
        [&solver]() { solver.scf_diagonalize(); });
    density_ =
        std::make_unique<NddoStep>(ATV::nddo::Phase::density, "nddo_density",
                                   [&solver]() { solver.scf_density(); });
    fock_ = std::make_unique<NddoStep>(ATV::nddo::Phase::fock, "nddo_fock",
                                       [&solver]() { solver.scf_fock(); });
    check_ = std::make_unique<NddoStep>(
        ATV::nddo::Phase::check, "nddo_check", [&solver, this]() {
            solver.scf_check();
            loop_converged_ = solver.result().converged;
        });
    converged_ = std::make_unique<NddoConverged>(loop_converged_);

    this->on(ATV::nddo::Phase::diagonalize).add_ref(*diagonalize_);
    this->on(ATV::nddo::Phase::density).add_ref(*density_);
    this->on(ATV::nddo::Phase::fock).add_ref(*fock_);
    this->on(ATV::nddo::Phase::check).add_ref(*check_);
    this->add_criterion_ref(*converged_);

    this->set_serial_on_replica();
    cpu::Driver<ATV::nddo::Phase>::build();

    // The stages close over the solver. Destroy them on every exit, including
    // a throw from setup or from a step, before those unique_ptrs die.
    struct ResetGuard {
        NddoMethod& self;
        ~ResetGuard() { self.reset_built_driver_(); }
    } reset_guard{*this};

    struct PassGuard {
        NddoMethod& self;
        bool armed = false;
        ~PassGuard() {
            if (armed)
                self.teardown_pass();
        }
    } pass_guard{*this};

    cpu::Frame loop_frame;
    setup_pass(loop_frame);
    pass_guard.armed = true;

    const int max_iter = solver.config().max_iter;
    for (int iter = 0; iter < max_iter; ++iter) {
        mutable_outer_status().reset();
        step_once(iter);
        if (poll_criteria_(iter) != ATV::Verdict::Continue)
            break;
    }
    pass_guard.armed = false;
    teardown_pass();

    if (!(this->converged() || loop_converged_))
        throw std::runtime_error("NDDOSolver AM1 SCF did not converge");
    return solver.result();
}

void NddoMethod::publish_(NDDOSolver& solver, double force_step,
                          cpu::Frame& frame, bool forces_requested,
                          bool accumulate) {
    auto& result = drive_(solver, frame);
    if (accumulate) {
        frame.meta.get(kEnergy) += result.total_energy_ev;
        frame.meta.get(kEnergyElectronic) += result.total_energy_ev;
    } else {
        frame.meta.get(kEnergy) = result.total_energy_ev;
        frame.meta.get(kEnergyElectronic) = result.total_energy_ev;
    }
    const int n = static_cast<int>(frame.nlocal());
    StorageReal* __restrict sp_charge = frame.span(kAtom_charge);
    for (int i = 0; i < n; ++i) {
        sp_charge[i] = static_cast<StorageReal>(
            result.mulliken_charges[static_cast<std::size_t>(i)]);
    }
    if (forces_requested) {
        auto& fr = solver.compute_forces(frame, force_step);
        AccReal* __restrict fx = frame.span(kAtom_fx);
        AccReal* __restrict fy = frame.span(kAtom_fy);
        AccReal* __restrict fz = frame.span(kAtom_fz);
        for (int i = 0; i < n; ++i) {
            const std::size_t base = static_cast<std::size_t>(3 * i);
            const AccReal fxi = static_cast<AccReal>(fr.forces_ev_a[base + 0]);
            const AccReal fyi = static_cast<AccReal>(fr.forces_ev_a[base + 1]);
            const AccReal fzi = static_cast<AccReal>(fr.forces_ev_a[base + 2]);
            if (accumulate) {
                fx[i] += fxi;
                fy[i] += fyi;
                fz[i] += fzi;
            } else {
                fx[i] = fxi;
                fy[i] = fyi;
                fz[i] = fzi;
            }
        }
    }
}

} // namespace semiempirical
} // namespace cpu
} // namespace ATV
