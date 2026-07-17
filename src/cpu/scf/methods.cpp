/// @file cpu/scf/methods.cpp
/// @brief Concrete RHF/UHF/RKS/UKS — Frame API; private integral + loop workspace.

#include "atomiverse/cpu/scf/methods.h"

#include <cmath>
#include <utility>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/core/log.h"
#include "atomiverse/core/units.h"
#include "atomiverse/cpu/scf/aspc_predictor.h"
#include "atomiverse/cpu/scf/internal/fci_pack.h"
#include "atomiverse/cpu/scf/internal/gradient.h"
#include "atomiverse/cpu/scf/internal/integrals.h"
#include "atomiverse/cpu/scf/internal/mp2_pack.h"
#include "atomiverse/cpu/scf/internal/scf_solver.h"
#include "atomiverse/cpu/scf/internal/scf_types.h"

namespace ATV {
namespace cpu {
namespace scf {

class ScfLoop : public SCFSolver {
  public:
    using SCFSolver::SCFSolver;
};

namespace {

struct WrappedPositionsSnapshot {
    std::vector<StorageReal> x, y, z;
};

WrappedPositionsSnapshot unwrap_for_scf(Frame& frame) {
    const int n = static_cast<int>(frame.nlocal());
    const StorageReal* __restrict px = frame.span(kAtom_x);
    const StorageReal* __restrict py = frame.span(kAtom_y);
    const StorageReal* __restrict pz = frame.span(kAtom_z);
    WrappedPositionsSnapshot snap{std::vector<StorageReal>(px, px + n),
                                  std::vector<StorageReal>(py, py + n),
                                  std::vector<StorageReal>(pz, pz + n)};
    if (n <= 1)
        return snap;
    if (frame.simbox.length(0) <= 0.0 || frame.simbox.length(1) <= 0.0
        || frame.simbox.length(2) <= 0.0)
        return snap;
    StorageReal* __restrict sp_x = frame.span(kAtom_x);
    StorageReal* __restrict sp_y = frame.span(kAtom_y);
    StorageReal* __restrict sp_z = frame.span(kAtom_z);
    const double x0 = sp_x[0], y0 = sp_y[0], z0 = sp_z[0];
    for (int i = 1; i < n; ++i) {
        double dx = sp_x[i] - x0, dy = sp_y[i] - y0, dz = sp_z[i] - z0;
        frame.simbox.min_image(dx, dy, dz);
        sp_x[i] = x0 + dx;
        sp_y[i] = y0 + dy;
        sp_z[i] = z0 + dz;
    }
    return snap;
}

void restore_wrapped(Frame& frame, WrappedPositionsSnapshot&& snap) {
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

template <typename D>
SCFConfig make_config(const MeanFieldMethod<D>& m, SpinType spin,
                      MeanFieldType method, XCType xc, int n_alpha,
                      int n_beta) {
    SCFConfig c;
    c.spin = spin;
    c.method = method;
    c.xc_type = xc;
    c.n_alpha = n_alpha;
    c.n_beta = n_beta;
    c.tol = m.tol;
    c.conv_tol_grad = m.conv_tol_grad;
    c.max_iter = m.max_iter;
    c.use_diis = m.use_diis;
    c.diis_max_vectors = m.diis_max_vectors;
    c.diis_start_iter = m.diis_start_iter;
    c.damping_factor = m.damping_factor;
    c.init_guess = m.init_guess;
    c.level_shift = m.level_shift;
    c.lindep_threshold = m.lindep_threshold;
    c.grid_level = m.grid_level;
    c.n_radial = m.n_radial;
    c.n_angular = m.n_angular;
    c.force_native_two_electron_deriv = m.force_native_two_electron_deriv;
    c.xc_grad_project_translation = m.xc_grad_project_translation;
    c.hybrid_alpha_override = m.hybrid_alpha_override;
    c.fixed_grid = m.fixed_grid;
    return c;
}

} // namespace

template <typename D>
MeanFieldMethod<D>::MeanFieldMethod() = default;
template <typename D>
MeanFieldMethod<D>::~MeanFieldMethod() = default;
template <typename D>
MeanFieldMethod<D>::MeanFieldMethod(MeanFieldMethod&&) noexcept = default;
template <typename D>
MeanFieldMethod<D>&
MeanFieldMethod<D>::operator=(MeanFieldMethod&&) noexcept = default;

template <typename D>
void MeanFieldMethod<D>::set_predictor_(int order, int corrector_steps,
                                        double corrector_tol, bool fallback) {
    use_predictor_ = true;
    corrector_steps_ = corrector_steps;
    corrector_tol_ = corrector_tol;
    fallback_full_scf_ = fallback;
    aspc_alpha_ = std::make_unique<ASPCPredictor>(order);
    aspc_beta_ = std::make_unique<ASPCPredictor>(order);
}

template <typename D>
void MeanFieldMethod<D>::publish_config_() {
    *config_sp_ =
        make_config(*this, spin_, method_, xc_type_, n_alpha_, n_beta_);
}

template <typename D>
void MeanFieldMethod<D>::sync_loop_config_() {
    publish_config_();
    loop_ = std::make_unique<ScfLoop>(*config_sp_);
    loop_built_ = false;
    if (use_predictor_) {
        SCFConfig cc = *config_sp_;
        cc.tol = corrector_tol_;
        cc.max_iter = corrector_steps_;
        corrector_loop_ = std::make_unique<ScfLoop>(cc);
    }
}

template <typename D>
void MeanFieldMethod<D>::ensure_loop_() {
    if (!loop_)
        sync_loop_config_();
}

template <typename D>
void MeanFieldMethod<D>::rebuild_integrals_(const Frame& frame) {
    *integrals_sp_ = IntegralPack{frame, basis, charge, compute_eri};
    resolve_occupation(integrals_sp_->n_electrons);
    sync_loop_config_();
}

template <typename D>
SCFResult& MeanFieldMethod<D>::run_on_pack_(const SCFInitialGuess* guess) {
    ensure_loop_();
    if (!loop_built_) {
        loop_->build(*integrals_sp_);
        if (corrector_loop_)
            corrector_loop_->build(*integrals_sp_);
        loop_built_ = true;
    }
    *result_sp_ = loop_->run(*integrals_sp_, guess);
    *has_result_sp_ = true;
    return *result_sp_;
}

template <typename D>
SCFResult& MeanFieldMethod<D>::run(const Frame& frame,
                                   const SCFInitialGuess* guess) {
    rebuild_integrals_(frame);
    return run_on_pack_(guess);
}

template <typename D>
SCFResult& MeanFieldMethod<D>::run(Frame& frame, const SCFInitialGuess* guess,
                                   bool write_energy) {
    auto& res = run(static_cast<const Frame&>(frame), guess);
    if (write_energy)
        frame.meta.get(kEnergy) = res.energy * units::kHartreeInEv;
    return res;
}

template <typename D>
GradientResult MeanFieldMethod<D>::forces(Frame& frame) {
    rebuild_integrals_(frame);
    run_on_pack_(nullptr);
    ensure_loop_();
    auto grad =
        scf_gradient(frame, *integrals_sp_, *result_sp_, loop_->config());
    const int natm = static_cast<int>(frame.nlocal());
    AccReal* __restrict sp_fx = frame.span(kAtom_fx);
    AccReal* __restrict sp_fy = frame.span(kAtom_fy);
    AccReal* __restrict sp_fz = frame.span(kAtom_fz);
    for (int A = 0; A < natm; ++A) {
        sp_fx[A] = grad.forces(A, 0);
        sp_fy[A] = grad.forces(A, 1);
        sp_fz[A] = grad.forces(A, 2);
    }
    frame.meta.get(kEnergy) = grad.energy;
    return grad;
}

template <typename D>
MP2Result MeanFieldMethod<D>::mp2(Frame& frame) {
    const bool prev = compute_eri;
    compute_eri = true;
    run(frame);
    compute_eri = prev;
    return mp2_energy(*integrals_sp_, *result_sp_);
}

template <typename D>
FCIResult MeanFieldMethod<D>::fci(Frame& frame) {
    const bool prev = compute_eri;
    compute_eri = true;
    run(frame);
    compute_eri = prev;
    return fci_energy(*integrals_sp_, *result_sp_);
}

template <typename D>
void MeanFieldMethod<D>::build(const cpu::BuildContext& ctx) {
    rebuild_integrals_(ctx.frame);
    ensure_loop_();
    loop_->build(*integrals_sp_);
    if (corrector_loop_)
        corrector_loop_->build(*integrals_sp_);
    loop_built_ = true;
    ATV_INFO("[scf] build method={} charge={} basis={} n_basis={} "
             "predictor={} compute_force={}",
             name(), charge, basis_set_name(basis), integrals_sp_->n_basis,
             use_predictor_, compute_force);
}

template <typename D>
void MeanFieldMethod<D>::launch(const cpu::LaunchContext& ctx) {
    auto& frame = ctx.frame;
    const int natm = static_cast<int>(frame.nlocal());
    auto snap = unwrap_for_scf(frame);

    rebuild_integrals_(frame);
    ensure_loop_();
    if (!loop_built_) {
        loop_->build(*integrals_sp_);
        if (corrector_loop_)
            corrector_loop_->build(*integrals_sp_);
        loop_built_ = true;
    }

    SCFResult result;
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
            result = corrector_loop_->run(*integrals_sp_, &density_guess);
            if (!result.converged && fallback_full_scf_)
                result = loop_->run(*integrals_sp_, &density_guess);
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
                result = loop_->run(*integrals_sp_, &density_guess);
            } else {
                result = loop_->run(*integrals_sp_);
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
        result = loop_->run(*integrals_sp_, &density_guess);
        has_last_density_ = true;
        last_P_alpha_ = result.P_alpha;
        last_P_beta_ = result.P_beta;
    } else {
        result = loop_->run(*integrals_sp_);
        has_last_density_ = true;
        last_P_alpha_ = result.P_alpha;
        last_P_beta_ = result.P_beta;
    }

    if (compute_force) {
        auto grad =
            scf_gradient(frame, *integrals_sp_, result, loop_->config());
        AccReal* __restrict sp_fx = frame.span(kAtom_fx);
        AccReal* __restrict sp_fy = frame.span(kAtom_fy);
        AccReal* __restrict sp_fz = frame.span(kAtom_fz);
        for (int A = 0; A < natm; ++A) {
            sp_fx[A] = grad.forces(A, 0);
            sp_fy[A] = grad.forces(A, 1);
            sp_fz[A] = grad.forces(A, 2);
        }
        frame.meta.get(kEnergy) = grad.energy;
    } else {
        frame.meta.get(kEnergy) = result.energy * units::kHartreeInEv;
    }

    if (on_evaluate_)
        on_evaluate_(result);

    *result_sp_ = std::move(result);
    *has_result_sp_ = true;
    restore_wrapped(frame, std::move(snap));
}

// explicit instantiations
template class MeanFieldMethod<RHF>;
template class MeanFieldMethod<UHF>;
template class MeanFieldMethod<RKS>;
template class MeanFieldMethod<UKS>;

RHF::RHF() {
    spin_ = SpinType::Restricted;
    method_ = MeanFieldType::HF;
    init_guess = InitGuess::MINAO;
    damping_factor = 0.0;
    tol = 1e-9;
}
void RHF::resolve_occupation(int n_electrons) {
    n_alpha_ = n_electrons / 2;
    n_beta_ = n_electrons / 2;
}

UHF::UHF() {
    spin_ = SpinType::Unrestricted;
    method_ = MeanFieldType::HF;
    init_guess = InitGuess::Hcore;
    damping_factor = 0.7;
    diis_max_vectors = 6;
    diis_start_iter = 5;
    tol = 1e-8;
}
UHF::UHF(int na, int nb) : UHF() {
    n_alpha = na;
    n_beta = nb;
}
void UHF::resolve_occupation(int) {
    n_alpha_ = n_alpha;
    n_beta_ = n_beta;
}

RKS::RKS() {
    spin_ = SpinType::Restricted;
    method_ = MeanFieldType::KS;
    xc_type_ = XCType::PBE;
    xc = XCType::PBE;
    init_guess = InitGuess::MINAO;
    damping_factor = 0.0;
    grid_level = GridLevel::MatchPySCF;
    diis_max_vectors = 8;
    diis_start_iter = 1;
    tol = 1e-9;
}
RKS::RKS(XCType xc_in) : RKS() {
    xc = xc_in;
    xc_type_ = xc_in;
}
void RKS::resolve_occupation(int n_electrons) {
    xc_type_ = xc;
    n_alpha_ = n_electrons / 2;
    n_beta_ = n_electrons / 2;
}

UKS::UKS() {
    spin_ = SpinType::Unrestricted;
    method_ = MeanFieldType::KS;
    xc_type_ = XCType::LSDA;
    xc = XCType::LSDA;
    init_guess = InitGuess::Hcore;
    damping_factor = 0.5;
    diis_max_vectors = 6;
    diis_start_iter = 1;
    tol = 1e-8;
}
UKS::UKS(int na, int nb, XCType xc_in) : UKS() {
    n_alpha = na;
    n_beta = nb;
    xc = xc_in;
    xc_type_ = xc_in;
}
void UKS::resolve_occupation(int) {
    xc_type_ = xc;
    n_alpha_ = n_alpha;
    n_beta_ = n_beta;
}

} // namespace scf
} // namespace cpu
} // namespace ATV
