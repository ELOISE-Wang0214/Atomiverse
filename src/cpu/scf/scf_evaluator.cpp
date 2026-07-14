/// @file cpu/scf/scf_evaluator.cpp
/// @brief SCFEvaluator implementation — unified AIMD force component.
///
/// As of driver-taxonomy phase 4, the SCF cycle runs through a
/// pre-built `cpu::scf::SCFSolver` member: the n_basis-sized workspaces
/// are allocated once in `build()`; geometry-dependent structures (X,
/// shell-pair lists, KS grid) are rebuilt inside each `run()` because
/// atoms move between launches. The ASPC corrector path runs through a
/// second pre-built solver member (`corrector_solver_`) because it needs
/// a different (tighter `tol`, capped `max_iter`) configuration — it is
/// the steady-state AIMD path once the predictor history fills, so it
/// must not allocate per launch().

#include "atomiverse/cpu/scf/scf_evaluator.h"

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/core/log.h"
#include "atomiverse/core/units.h"

#include <cmath>
#include <utility>

namespace {

/// @brief In-place minimum-image unwrap of an isolated molecule.
///
/// Atom 0 is the anchor; every other atom i is mapped to its
/// minimum-image position relative to atom 0 via `cpu::Box::min_image`
/// (orthorhombic and triclinic cells). After this call, the atoms
/// reassemble into a single
/// physically connected fragment regardless of how the MD integrator
/// wrapped each atom independently across the box.
///
/// This is the fix for the wrap × non-minimum-image SCF bug: the SCF
/// integral builder is not minimum-image aware, so an atom that
/// `VelocityVerlet::wrap_orthorhombic` placed on the opposite face of
/// the box from its molecular partners will look ~one box length away
/// in raw coordinates, and the SCF collapses onto an unphysical
/// fragmented limit.
///
/// @pre  The system is a single isolated fragment (one molecule per
///       SCF call). For multi-fragment SCF (e.g. periodic boundary
///       Ewald or multiple non-bonded molecules in one cell), this
///       anchor-on-atom-0 scheme will need generalisation.
/// @return Snapshot of the wrapped positions so the caller can restore
///         them before launch returns (downstream observers expect
///         wrapped positions in `frame.atoms.x/y/z`).
/// Wrapped positions saved for restoration after the SCF unwrap.
///
/// Plain vectors, restored by an ELEMENT-WISE COPY back into the frame's
/// existing columns — never by swapping a whole array into `frame.atoms.x`.
/// Swapping frees the frame's storage mid-force-phase, dangling every raw
/// pointer a component resolved at its launch boundary (LBFGS's line search
/// replays the force phase and holds exactly such pointers — ASAN caught the
/// use-after-free). Only `grow`/`set_capacity` may reallocate a frame column,
/// and only at `pre_force`; see cpu/frame/atom.h.
struct WrappedPositionsSnapshot {
    std::vector<ATV::StorageReal> x, y, z;
};

WrappedPositionsSnapshot unwrap_for_scf(ATV::cpu::Frame& frame) {
    using namespace ATV;
    const int n = static_cast<int>(frame.atoms.size());
    const StorageReal* __restrict px = frame.span(kAtom_x);
    const StorageReal* __restrict py = frame.span(kAtom_y);
    const StorageReal* __restrict pz = frame.span(kAtom_z);
    WrappedPositionsSnapshot snap{
        std::vector<ATV::StorageReal>(px, px + n),
        std::vector<ATV::StorageReal>(py, py + n),
        std::vector<ATV::StorageReal>(pz, pz + n)};

    if (n <= 1)
        return snap;

    // Gas-phase frame (no/degenerate box): nothing to unwrap, and
    // Box::min_image would divide by a zero edge length.
    if (frame.box.length(0) <= 0.0 || frame.box.length(1) <= 0.0
        || frame.box.length(2) <= 0.0)
        return snap;

    StorageReal* __restrict sp_x = frame.span(kAtom_x);
    StorageReal* __restrict sp_y = frame.span(kAtom_y);
    StorageReal* __restrict sp_z = frame.span(kAtom_z);

    const double x0 = sp_x[0];
    const double y0 = sp_y[0];
    const double z0 = sp_z[0];

    for (int i = 1; i < n; ++i) {
        double dx = sp_x[i] - x0;
        double dy = sp_y[i] - y0;
        double dz = sp_z[i] - z0;
        frame.box.min_image(dx, dy, dz);
        sp_x[i] = x0 + dx;
        sp_y[i] = y0 + dy;
        sp_z[i] = z0 + dz;
    }
    return snap;
}

void restore_wrapped(ATV::cpu::Frame& frame, WrappedPositionsSnapshot&& snap) {
    using namespace ATV;
    // Copy back IN PLACE — the frame's column storage keeps its identity, so
    // pointers other components resolved at their launch boundary stay valid.
    StorageReal* __restrict px = frame.span(kAtom_x);
    StorageReal* __restrict py = frame.span(kAtom_y);
    StorageReal* __restrict pz = frame.span(kAtom_z);
    const int n = static_cast<int>(snap.x.size());
    for (int i = 0; i < n; ++i) {
        px[i] = snap.x[static_cast<std::size_t>(i)];
        py[i] = snap.y[static_cast<std::size_t>(i)];
        pz[i] = snap.z[static_cast<std::size_t>(i)];
    }
}

} // namespace

namespace ATV {
namespace cpu {
namespace scf {

SCFEvaluator::SCFEvaluator(SCFConfig config, BasisSet basis, int charge)
    : config_(std::move(config)), basis_(basis), charge_(charge),
      config_shared_(std::make_shared<SCFConfig>(config_)) {}

// ── Builder methods ─────────────────────────────────────────────────
//
// `predictor()`, `on_evaluate()`, and (in the header) `with_gradient()`
// are ref-qualified so the chain preserves the value category of the
// invocation. When called on an rvalue temporary (e.g.
// `sim.add(SCFEvaluator{cfg}.predictor(3))`), the rvalue overload
// returns `SCFEvaluator&&`, which makes `Driver::add(Comp&& comp)`
// deduce `Comp` as `SCFEvaluator` (rvalue) and move the temporary into
// `unique_ptr` ownership rather than storing a non-owning pointer that
// dangles when the temporary dies. The lvalue overload preserves the
// existing `evaluator.predictor(3); driver.add(std::move(evaluator));`
// pattern.

void SCFEvaluator::set_predictor_(int order, int corrector_steps,
                                  double corrector_tol, bool fallback) {
    use_predictor_ = true;
    corrector_steps_ = corrector_steps;
    corrector_tol_ = corrector_tol;
    fallback_full_scf_ = fallback;
    aspc_alpha_ = std::make_unique<ASPCPredictor>(order);
    aspc_beta_ = std::make_unique<ASPCPredictor>(order);
}

SCFEvaluator& SCFEvaluator::predictor(int order, int corrector_steps,
                                      double corrector_tol, bool fallback) & {
    set_predictor_(order, corrector_steps, corrector_tol, fallback);
    return *this;
}

SCFEvaluator&& SCFEvaluator::predictor(int order, int corrector_steps,
                                       double corrector_tol, bool fallback) && {
    set_predictor_(order, corrector_steps, corrector_tol, fallback);
    return std::move(*this);
}

void SCFEvaluator::set_on_evaluate_(std::function<void(const SCFResult&)> cb) {
    on_evaluate_ = std::move(cb);
}

SCFEvaluator&
SCFEvaluator::on_evaluate(std::function<void(const SCFResult&)> cb) & {
    set_on_evaluate_(std::move(cb));
    return *this;
}

SCFEvaluator&&
SCFEvaluator::on_evaluate(std::function<void(const SCFResult&)> cb) && {
    set_on_evaluate_(std::move(cb));
    return std::move(*this);
}

void SCFEvaluator::build(const cpu::BuildContext& ctx) {
    // Build a representative IntegralPack for the current geometry — used
    // both as the working pack for launch #1 and to size the SCFSolver's
    // persistent workspaces (X, shell-pair lists, KS grid). Subsequent
    // launches reuse those workspaces.
    *last_integrals_ = build_integrals(ctx.frame, basis_, charge_);

    solver_ = std::make_unique<SCFSolver>(config_);
    solver_->build(*last_integrals_);

    if (use_predictor_) {
        // The ASPC corrector runs every steady-state AIMD step once the
        // predictor history fills — pre-build its solver here so launch()
        // never constructs one.
        SCFConfig corrector_config = config_;
        corrector_config.tol = corrector_tol_;
        corrector_config.max_iter = corrector_steps_;
        corrector_solver_ = std::make_unique<SCFSolver>(corrector_config);
        corrector_solver_->build(*last_integrals_);
    }

    ATV_INFO("[scf_eval] build charge={} predictor={} compute_force={} "
             "basis={} n_basis={}",
             charge_, use_predictor_, compute_force_, basis_set_name(basis_),
             last_integrals_->n_basis);
}

void SCFEvaluator::launch(const cpu::LaunchContext& ctx) {
    const auto t_total = LogClock::now();
    auto& frame = ctx.frame;
    const int natm = static_cast<int>(frame.atoms.size());
    ATV_INFO("[scf_eval] step={} start natm={} charge={}", ctx.step, natm,
             charge_);

    // Unwrap the molecule before SCF: build_integrals + scf_gradient are
    // not minimum-image aware, so an atom that the MD integrator wrapped
    // to the opposite face of the box would otherwise look ~box-length
    // away from its partners and collapse the SCF. We anchor atom 0 and
    // map every other atom into its minimum-image position, then restore
    // the wrapped positions before the function returns so downstream
    // observers (recorders, RDF, etc.) see the conventional wrapped form.
    auto snap = unwrap_for_scf(frame);

    const auto t_integrals = LogClock::now();
    *last_integrals_ = build_integrals(frame, basis_, charge_);
    auto& integrals = *last_integrals_;
    ATV_DEBUG("[scf_eval] step={} integrals_ms={:.1f} basis={} n_basis={} "
              "electrons={}",
              ctx.step, elapsed_ms(t_integrals), basis_set_name(basis_),
              integrals.n_basis, integrals.n_electrons);

    SCFResult result;

    const auto t_scf = LogClock::now();
    SCFInitialGuess density_guess;
    if (use_predictor_ && aspc_alpha_ && aspc_beta_) {
        xt::xtensor<double, 2> Pa_guess, Pb_guess;
        const bool have_alpha = aspc_alpha_->predict(Pa_guess);
        const bool have_beta = aspc_beta_->predict(Pb_guess);

        if (have_alpha && have_beta) {
            density_guess.P_alpha = std::move(Pa_guess);
            density_guess.P_beta = std::move(Pb_guess);
            density_guess.has_alpha = true;
            density_guess.has_beta = true;

            // ASPC corrector run: tighter tol, capped iterations — through
            // the corrector solver pre-built in build() (no per-step
            // allocation on this steady-state path).
            ATV_INFO("[scf_eval] step={} scf guess source=aspc_predictor",
                     ctx.step);
            result = corrector_solver_->run(integrals, &density_guess);

            if (!result.converged && fallback_full_scf_)
                result = solver_->run(integrals, &density_guess);

            xt::xtensor<double, 2> Pa_corr, Pb_corr;
            aspc_alpha_->correct(result.P_alpha, density_guess.P_alpha,
                                 Pa_corr);
            aspc_beta_->correct(result.P_beta, density_guess.P_beta, Pb_corr);
            aspc_alpha_->record(Pa_corr);
            aspc_beta_->record(Pb_corr);
            has_last_density_ = true;
            last_P_alpha_ = result.P_alpha;
            last_P_beta_ = result.P_beta;
        } else {
            if (has_last_density_) {
                density_guess.P_alpha = last_P_alpha_;
                density_guess.P_beta = last_P_beta_;
                density_guess.has_alpha = true;
                density_guess.has_beta = true;
                ATV_INFO("[scf_eval] step={} scf guess source=previous_density",
                         ctx.step);
                result = solver_->run(integrals, &density_guess);
            } else {
                result = solver_->run(integrals);
            }
            aspc_alpha_->record(result.P_alpha);
            aspc_beta_->record(result.P_beta);
            has_last_density_ = true;
            last_P_alpha_ = result.P_alpha;
            last_P_beta_ = result.P_beta;
        }
    } else if (has_last_density_) {
        density_guess.P_alpha = last_P_alpha_;
        density_guess.P_beta = last_P_beta_;
        density_guess.has_alpha = true;
        density_guess.has_beta = true;
        ATV_INFO("[scf_eval] step={} scf guess source=previous_density",
                 ctx.step);
        result = solver_->run(integrals, &density_guess);
        has_last_density_ = true;
        last_P_alpha_ = result.P_alpha;
        last_P_beta_ = result.P_beta;
    } else {
        result = solver_->run(integrals);
        has_last_density_ = true;
        last_P_alpha_ = result.P_alpha;
        last_P_beta_ = result.P_beta;
    }
    ATV_DEBUG("[scf_eval] step={} scf_ms={:.1f} converged={} iterations={} "
              "energy_ha={:.10f}",
              ctx.step, elapsed_ms(t_scf), result.converged, result.iterations,
              result.energy);

    // Forces — only when gradient is enabled. Energy-only single-point
    // mode (with_gradient(false)) skips this entirely.
    if (compute_force_) {
        const auto t_gradient = LogClock::now();
        auto grad = scf_gradient(frame, integrals, result, config_);
        ATV_DEBUG("[scf_eval] step={} gradient_ms={:.1f}", ctx.step,
                  elapsed_ms(t_gradient));

        double max_force = 0.0;
        AccReal* __restrict sp_fx = frame.span(kAtom_fx);
        AccReal* __restrict sp_fy = frame.span(kAtom_fy);
        AccReal* __restrict sp_fz = frame.span(kAtom_fz);
        for (int A = 0; A < natm; ++A) {
            sp_fx[A] = grad.forces(A, 0);
            sp_fy[A] = grad.forces(A, 1);
            sp_fz[A] = grad.forces(A, 2);
            const double fx = grad.forces(A, 0);
            const double fy = grad.forces(A, 1);
            const double fz = grad.forces(A, 2);
            const double mag = std::sqrt(fx * fx + fy * fy + fz * fz);
            max_force = std::max(max_force, mag);
        }
        frame.scalars.energy = grad.energy;

        ATV_INFO("[scf_eval] step={} done total_ms={:.1f} energy_eV={:.6f} "
                 "max_force_eVA={:.6f} scf_converged={} scf_iterations={}",
                 ctx.step, elapsed_ms(t_total), frame.scalars.energy, max_force,
                 result.converged, result.iterations);
    } else {
        // Energy-only path. Convert result.energy (Ha) to frame.scalars.energy (eV).
        frame.scalars.energy = result.energy * ATV::kHartree2eV;
        ATV_INFO("[scf_eval] step={} done total_ms={:.1f} energy_eV={:.6f} "
                 "(energy-only) scf_converged={} scf_iterations={}",
                 ctx.step, elapsed_ms(t_total), frame.scalars.energy,
                 result.converged, result.iterations);
    }

    if (on_evaluate_)
        on_evaluate_(result);

    *last_result_ = std::move(result);
    *has_last_result_ = true;

    // Restore the wrapped atomic positions so the rest of the MD loop
    // and external observers see the conventional in-box representation.
    restore_wrapped(frame, std::move(snap));
}

} // namespace scf
} // namespace cpu
} // namespace ATV
