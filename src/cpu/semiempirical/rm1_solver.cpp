/// @file src/cpu/semiempirical/rm1_solver.cpp
#include "atomiverse/cpu/semiempirical/rm1_solver.h"
#include <stdexcept>
#include <utility>

#include "atomiverse/core/frame_schema.h"
namespace ATV {
namespace cpu {
namespace semiempirical {
RM1::RM1() = default;

RM1::RM1(int charge_, int multiplicity_, int max_iter_, double density_tol_,
         double energy_tol_)
    : charge(charge_), multiplicity(multiplicity_), max_iter(max_iter_),
      density_tol(density_tol_), energy_tol(energy_tol_) {}

RM1::~RM1() = default;
RM1::RM1(RM1&&) noexcept = default;
RM1& RM1::operator=(RM1&&) noexcept = default;
void RM1::sync_engine_() {
    NDDOConfig c;
    c.method = SemiEmpiricalMethod::RM1;
    c.charge = charge;
    c.multiplicity = multiplicity;
    c.max_iter = max_iter;
    c.density_tol = density_tol;
    c.energy_tol = energy_tol;
    engine_ = std::make_unique<NDDOSolver>(c);
}
void RM1::build(const cpu::Frame& reference) {
    sync_engine_();
    engine_->build(reference);
}
RM1Result& RM1::run(const cpu::Frame& frame) {
    if (!engine_)
        build(frame);
    return engine_->run(frame);
}
RM1Result& RM1::run(cpu::Frame& frame) {
    RM1Result& r = run(std::as_const(frame));
    frame.meta.get(kEnergy) = r.total_energy_ev;
    return r;
}
RM1Result& RM1::forces(cpu::Frame& frame) {
    RM1Result& r = run(std::as_const(frame));
    fd_forces_(frame);
    const int n = frame.nlocal();
    AccReal* __restrict fx = frame.span(kAtom_fx);
    AccReal* __restrict fy = frame.span(kAtom_fy);
    AccReal* __restrict fz = frame.span(kAtom_fz);
    for (int i = 0; i < n; ++i) {
        const std::size_t base = static_cast<std::size_t>(3 * i);
        fx[i] = static_cast<AccReal>(r.forces_ev_a[base + 0]);
        fy[i] = static_cast<AccReal>(r.forces_ev_a[base + 1]);
        fz[i] = static_cast<AccReal>(r.forces_ev_a[base + 2]);
    }
    frame.meta.get(kEnergy) = r.total_energy_ev;
    return r;
}
RM1Result& RM1::fd_forces_(cpu::Frame& frame) {
    if (!engine_)
        throw std::runtime_error("RM1::forces before build/run");
    return engine_->compute_forces(frame, force_step);
}
const RM1Result& RM1::last_result() const {
    if (!engine_) {
        static RM1Result empty;
        return empty;
    }
    return engine_->last_result();
}
bool RM1::built() const noexcept { return engine_ && engine_->built(); }
bool RM1::has_density() const noexcept {
    return engine_ && engine_->has_density();
}
const NDDOSystem& RM1::system() const {
    if (!engine_)
        throw std::runtime_error("RM1::system before build");
    return engine_->system();
}
const NDDOScfState& RM1::state() const {
    if (!engine_)
        throw std::runtime_error("RM1::state before build");
    return engine_->state();
}
} // namespace semiempirical
} // namespace cpu
} // namespace ATV
