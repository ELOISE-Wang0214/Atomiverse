/// @file src/cpu/semiempirical/am1_solver.cpp
#include "atomiverse/cpu/semiempirical/am1_solver.h"

#include <stdexcept>

namespace ATV {
namespace cpu {
namespace semiempirical {

AM1::AM1() = default;
AM1::~AM1() = default;
AM1::AM1(AM1&&) noexcept = default;
AM1& AM1::operator=(AM1&&) noexcept = default;

AM1::AM1(const AM1& other)
    : charge(other.charge), multiplicity(other.multiplicity),
      max_iter(other.max_iter), density_tol(other.density_tol),
      energy_tol(other.energy_tol) {}

AM1& AM1::operator=(const AM1& other) {
    if (this != &other) {
        charge = other.charge;
        multiplicity = other.multiplicity;
        max_iter = other.max_iter;
        density_tol = other.density_tol;
        energy_tol = other.energy_tol;
        engine_.reset();
    }
    return *this;
}

void AM1::sync_engine_() {
    NDDOConfig c;
    c.method = SemiEmpiricalMethod::AM1;
    c.charge = charge;
    c.multiplicity = multiplicity;
    c.max_iter = max_iter;
    c.density_tol = density_tol;
    c.energy_tol = energy_tol;
    engine_ = std::make_unique<NDDOSolver>(c);
}

void AM1::build(const cpu::Frame& reference) {
    sync_engine_();
    engine_->build(reference);
}

AM1Result& AM1::run(const cpu::Frame& frame) {
    if (!engine_)
        build(frame);
    return engine_->run(frame);
}

AM1Result& AM1::compute_forces(cpu::Frame& frame, double step_angstrom) {
    if (!engine_)
        throw std::runtime_error("AM1::compute_forces before build/run");
    return engine_->compute_forces(frame, step_angstrom);
}

const AM1Result& AM1::last_result() const {
    if (!engine_) {
        static AM1Result empty;
        return empty;
    }
    return engine_->last_result();
}

bool AM1::built() const noexcept {
    return engine_ && engine_->built();
}

bool AM1::has_density() const noexcept {
    return engine_ && engine_->has_density();
}

const NDDOSystem& AM1::system() const {
    if (!engine_)
        throw std::runtime_error("AM1::system before build");
    return engine_->system();
}

const NDDOScfState& AM1::state() const {
    if (!engine_)
        throw std::runtime_error("AM1::state before build");
    return engine_->state();
}

} // namespace semiempirical
} // namespace cpu
} // namespace ATV
