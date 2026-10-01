/// @file cuda/scf/rhf.cu
/// @brief cuda::scf::RHF — concrete restricted HF on the CUDA backend.
///
/// **Energy path (primary):** STO-3G with l ≤ 1 (H…O, mixed Z, e.g. H₂O):
/// integrals + SCF stay on device (`RHFFockDevice` / `launch_*_sp`) —
/// **CUDA-only HF**, no host ERI tensor. The iteration is the Driver plan
/// (PySCF `scf.hf.kernel` / `scf.diis.CDIIS.update`). Other bases fill
/// S / Hcore / ERI once in `assemble_` (setup), never inside a stage.
///
/// **Forces path (Wave 6):** STO-3G s/p (l ≤ 1, H…O, e.g. H₂O) analytic on
/// device (`launch_rhf_forces_sp`). Out-of-coverage (`compute_force` on a
/// higher-L / other basis) throws with a how-to — no silent host fallback.
/// Use `with_gradient(false)` for energy-only, or `cpu::scf::RHF` for
/// other bases.

#ifdef ATV_HAS_HF

#include "atomiverse/cuda/scf/rhf.cuh"

#include <cuda_runtime.h>

#include <array>
#include <stdexcept>
#include <string>
#include <vector>

#include <xtensor/containers/xtensor.hpp>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/core/units.h"
#include "atomiverse/cpu/frame/transfer.h"
#include "atomiverse/cpu/scf/basis_data/accessors.h"
#include "atomiverse/cpu/scf/dense_matrix.h"
#include "atomiverse/cpu/scf/grid.h"
#include "atomiverse/cpu/scf/internal/integrals.h"
#include "atomiverse/cpu/scf/rhf.h"
#include "atomiverse/cpu/scf/scf_enums.h"
#include "atomiverse/cpu/scf/slater_determinant.h"
#include "atomiverse/cpu/scf/xc_functional.h"
#include "atomiverse/cuda/frame/transfer.cuh"
#include "atomiverse/cuda/frame/types.cuh"
#include "atomiverse/cuda/scf/device_integrals.cuh"
#include "atomiverse/cuda/scf/internal/rhf_scf_cycle.cuh"

namespace ATV {
namespace cuda {
namespace scf {
namespace {

/// @brief True when every atom is s/p STO-3G (l ≤ 1) with uniform n_prim.
///
/// Wave 4: mixed Z (H₂O) OK — molecular shell list built from the table.
/// Numbers come from `basis_data::sto3g()`, never retyped here.
bool device_sp_sto3g_ok(cpu::scf::BasisSet basis, const int* h_Z, int natm) {
    if (basis != cpu::scf::BasisSet::STO3G || natm <= 0 || h_Z == nullptr)
        return false;
    int n_prim0 = -1;
    for (int i = 0; i < natm; ++i) {
        const int Z = h_Z[i];
        if (!cpu::scf::basis_data::is_sp_sto3g(Z))
            return false;
        const int np = cpu::scf::basis_data::sp_sto3g_n_prim(Z);
        if (n_prim0 < 0)
            n_prim0 = np;
        else if (np != n_prim0)
            return false;
    }
    return true;
}

/// @brief Reconstruct a host Frame from device positions + nuclear Z.
///
/// Used for force assembly and the transitional host-integral energy path.
cpu::Frame host_frame_from_device(const cuda::Frame& dframe) {
    const int n = dframe.n_atoms();
    if (n <= 0) {
        throw std::invalid_argument(
            "cuda::scf::RHF: device frame has no atoms");
    }
    std::vector<real_t> hx(static_cast<size_t>(n));
    std::vector<real_t> hy(static_cast<size_t>(n));
    std::vector<real_t> hz(static_cast<size_t>(n));
    std::vector<int> htype(static_cast<size_t>(n));
    cudaMemcpy(hx.data(), dframe.atoms.x,
               sizeof(real_t) * static_cast<size_t>(n), cudaMemcpyDeviceToHost);
    cudaMemcpy(hy.data(), dframe.atoms.y,
               sizeof(real_t) * static_cast<size_t>(n), cudaMemcpyDeviceToHost);
    cudaMemcpy(hz.data(), dframe.atoms.z,
               sizeof(real_t) * static_cast<size_t>(n), cudaMemcpyDeviceToHost);
    cudaMemcpy(htype.data(), dframe.atoms.type,
               sizeof(int) * static_cast<size_t>(n), cudaMemcpyDeviceToHost);

    cpu::Frame host(n);
    auto x = host.span(kAtom_x);
    auto y = host.span(kAtom_y);
    auto z = host.span(kAtom_z);
    auto type = host.span(kAtom_type);
    for (int i = 0; i < n; ++i) {
        x[i] = static_cast<double>(hx[static_cast<size_t>(i)]);
        y[i] = static_cast<double>(hy[static_cast<size_t>(i)]);
        z[i] = static_cast<double>(hz[static_cast<size_t>(i)]);
        type[i] = htype[static_cast<size_t>(i)];
    }
    // Isolated-molecule box placeholder (integrals do not use the cell).
    host.box = cpu::Box::cubic(100.0);
    return host;
}

/// @brief Transitional host fill of S/Hcore/ERI → device SCF (p shells).
std::unique_ptr<internal::RhfScfCycle>
make_host_integral_solver(const cpu::Frame& frame, cpu::scf::BasisSet basis,
                          int charge, double tol, int max_iter,
                          double level_shift) {
    xt::xtensor<double, 2> S;
    xt::xtensor<double, 2> Hcore;
    xt::xtensor<double, 4> ERI;
    double V_nuc = 0.0;
    int n_basis = 0;
    int n_electrons = 0;
    cpu::scf::internal::IntegralWorkspace::fill_rhf_matrices(
        frame, basis, charge, /*compute_eri=*/true, S, Hcore, ERI, V_nuc,
        n_basis, n_electrons);
    return std::make_unique<internal::RhfScfCycle>(
        std::move(S), std::move(Hcore), std::move(ERI), V_nuc, n_electrons, tol,
        max_iter, level_shift);
}

/// @brief CUDA-only RHF: molecular s/p STO-3G → device 1e/ERI + SCF (Wave 4).
///
/// Builds a full-molecule shell list from `basis_data::sp_sto3g_shells(Z)`
/// per atom (mixed Z OK). n_AO = Σ n_cart(L) with s→1, p→3.
/// @param level_shift Virtual level shift λ [Hartree]; wired on device Fock.
std::unique_ptr<internal::RhfScfCycle>
make_device_sp_sto3g_solver(const int* h_Z, int natm, int charge, double tol,
                            int max_iter, double level_shift) {
    int n_electrons = -charge;
    for (int i = 0; i < natm; ++i)
        n_electrons += h_Z[i];
    if (n_electrons < 0 || n_electrons % 2 != 0) {
        throw std::invalid_argument(
            "cuda::scf::RHF: closed-shell CUDA path requires even non-negative "
            "electron count (n_electrons="
            + std::to_string(n_electrons) + ", charge=" + std::to_string(charge)
            + ")");
    }
    const int n_occ = n_electrons / 2;

    // Count shells / AOs, then fill molecular SoA from the basis table.
    int n_shells = 0;
    int n_ao = 0;
    int n_prim = -1;
    for (int a = 0; a < natm; ++a) {
        const auto& shells = cpu::scf::basis_data::sp_sto3g_shells(h_Z[a]);
        if (n_prim < 0)
            n_prim = cpu::scf::basis_data::sp_sto3g_n_prim(h_Z[a]);
        for (const auto& sh : shells) {
            ++n_shells;
            n_ao += cpu::scf::basis_data::n_cart_l(sh.l);
        }
    }

    std::vector<double> alpha_mol(static_cast<std::size_t>(n_shells)
                                  * static_cast<std::size_t>(n_prim));
    std::vector<double> coeff_mol(alpha_mol.size());
    std::vector<int> L_mol(static_cast<std::size_t>(n_shells));
    std::vector<int> atom_of(static_cast<std::size_t>(n_shells));

    int s = 0;
    for (int a = 0; a < natm; ++a) {
        const auto& shells = cpu::scf::basis_data::sp_sto3g_shells(h_Z[a]);
        for (const auto& sh : shells) {
            L_mol[static_cast<std::size_t>(s)] = sh.l;
            atom_of[static_cast<std::size_t>(s)] = a;
            for (int p = 0; p < n_prim; ++p) {
                const std::size_t k = static_cast<std::size_t>(s * n_prim + p);
                alpha_mol[k] = sh.alpha[static_cast<std::size_t>(p)];
                coeff_mol[k] = sh.coeff[static_cast<std::size_t>(p)];
            }
            ++s;
        }
    }

    return std::make_unique<internal::RhfScfCycle>(
        n_shells, n_ao, n_prim, alpha_mol.data(), coeff_mol.data(),
        L_mol.data(), atom_of.data(), n_occ, tol, max_iter, level_shift);
}

} // namespace

struct EngineSlot {
    std::unique_ptr<internal::RhfScfCycle> cycle;
    cuda::Context* ctx = nullptr;
    cuda::Frame* frame = nullptr;
};

class RhfStage : public cuda::Component {
  public:
    explicit RhfStage(EngineSlot* slot) : slot_(slot) {}

  protected:
    EngineSlot* slot_ = nullptr;
};

class RhfFockStage : public RhfStage {
  public:
    using RhfStage::RhfStage;
    std::string name() const override { return "scf_fock_stage"; }
    void launch(const cuda::LaunchContext& ctx) override {
        if (!ctx.at(ATV::scf::Phase::fock) || !slot_ || !slot_->cycle)
            return;
        slot_->cycle->step_fock(ctx.runtime.stream, ctx.frame, ctx.step);
    }
};

class RhfMixStage : public RhfStage {
  public:
    using RhfStage::RhfStage;
    std::string name() const override { return "scf_mix_stage"; }
    void launch(const cuda::LaunchContext& ctx) override {
        if (!ctx.at(ATV::scf::Phase::mix) || !slot_ || !slot_->cycle)
            return;
        slot_->cycle->step_mix(ctx.runtime.stream, ctx.frame, ctx.step);
    }
};

class RhfDiagStage : public RhfStage {
  public:
    using RhfStage::RhfStage;
    std::string name() const override { return "scf_diagonalize_stage"; }
    void launch(const cuda::LaunchContext& ctx) override {
        if (!ctx.at(ATV::scf::Phase::diagonalize) || !slot_ || !slot_->cycle)
            return;
        slot_->cycle->step_diagonalize(ctx.runtime.stream, ctx.frame, ctx.step);
    }
};

class RhfDensityStage : public RhfStage {
  public:
    using RhfStage::RhfStage;
    std::string name() const override { return "scf_density_stage"; }
    void launch(const cuda::LaunchContext& ctx) override {
        if (!ctx.at(ATV::scf::Phase::density) || !slot_ || !slot_->cycle)
            return;
        slot_->cycle->step_density(ctx.runtime.stream, ctx.frame, ctx.step);
    }
};

class RhfCheckStage : public RhfStage {
  public:
    using RhfStage::RhfStage;
    std::string name() const override { return "scf_check_stage"; }
    void launch(const cuda::LaunchContext& ctx) override {
        if (!ctx.at(ATV::scf::Phase::check) || !slot_ || !slot_->cycle)
            return;
        slot_->cycle->step_check(ctx.runtime.stream, ctx.frame, ctx.step);
    }
};

class RhfConverged : public cuda::Criterion {
  public:
    explicit RhfConverged(EngineSlot* slot) : slot_(slot) {}
    const char* name() const override { return "scf_converged"; }
    ATV::Verdict after_iteration(const cuda::LoopState& state) override {
        if (!slot_ || !slot_->cycle)
            return ATV::Verdict::Continue;
        // PySCF kernel: stop when |ΔE| and the orbital gradient are inside
        // tolerance. This path judges the density residual the check stage
        // wrote (device flag or host `workspace.converged`).
        if (slot_->cycle->poll_converged(state.runtime.stream, state.iteration))
            return ATV::Verdict::Converged;
        return ATV::Verdict::Continue;
    }

  private:
    EngineSlot* slot_ = nullptr;
};

struct RhfLoopState {
    EngineSlot slot;
    RhfFockStage fock;
    RhfMixStage mix;
    RhfDiagStage diag;
    RhfDensityStage density;
    RhfCheckStage check;

    RhfLoopState()
        : fock(&slot), mix(&slot), diag(&slot), density(&slot), check(&slot) {}
};

RHF::RHF(cpu::scf::BasisSet basis_in, int charge_in, double tol_in,
         int max_iter_in, double level_shift_in, bool compute_force_in)
    : basis(basis_in), charge(charge_in), tol(tol_in), max_iter(max_iter_in),
      level_shift(level_shift_in), compute_force(compute_force_in),
      loop_(std::make_unique<RhfLoopState>()) {
    wire_stages_();
}

RHF::~RHF() { release_engine_(); }

RHF::RHF(RHF&& other) noexcept
    : cuda::Driver<ATV::scf::Phase>(std::move(other)), basis(other.basis),
      charge(other.charge), tol(other.tol), max_iter(other.max_iter),
      level_shift(other.level_shift), compute_force(other.compute_force),
      loop_(std::move(other.loop_)), energy_ha_(other.energy_ha_),
      iterations_(other.iterations_),
      last_force_device_(other.last_force_device_),
      last_device_sp_(other.last_device_sp_),
      host_Z_(std::move(other.host_Z_)) {
    other.energy_ha_ = 0.0;
    other.iterations_ = 0;
    other.compute_force = false;
    other.last_force_device_ = false;
    other.last_device_sp_ = false;
}

void RHF::wire_stages_() {
    static_assert(static_cast<int>(ATV::scf::Phase::fock) == 0);
    static_assert(static_cast<int>(ATV::scf::Phase::mix) == 1);
    static_assert(static_cast<int>(ATV::scf::Phase::diagonalize) == 2);
    static_assert(static_cast<int>(ATV::scf::Phase::density) == 3);
    static_assert(static_cast<int>(ATV::scf::Phase::check) == 4);

    this->on(ATV::scf::Phase::fock).add_ref(loop_->fock);
    this->on(ATV::scf::Phase::mix).add_ref(loop_->mix);
    this->on(ATV::scf::Phase::diagonalize).add_ref(loop_->diag);
    this->on(ATV::scf::Phase::density).add_ref(loop_->density);
    this->on(ATV::scf::Phase::check).add_ref(loop_->check);
    this->add_criterion(RhfConverged{&loop_->slot});
    this->build();
}

void RHF::cache_host_types(const cpu::Frame& frame) {
    const int n = static_cast<int>(frame.nlocal());
    if (n <= 0) {
        host_Z_.clear();
        return;
    }
    host_Z_.resize(static_cast<std::size_t>(n));
    const auto type = frame.span(kAtom_type);
    for (int i = 0; i < n; ++i)
        host_Z_[static_cast<std::size_t>(i)] = type[i];
}

void RHF::release_engine_() {
    if (!loop_ || !loop_->slot.cycle)
        return;
    cuda::Context* ctx = loop_->slot.ctx;
    cuda::Frame* frame = loop_->slot.frame;
    if (ctx && frame)
        loop_->slot.cycle->destroy(ctx->stream, *frame);
    loop_->slot.cycle.reset();
    loop_->slot.ctx = nullptr;
    loop_->slot.frame = nullptr;
}

void RHF::assemble_(cuda::Context& ctx, cuda::Frame& frame) {
    const int natm = frame.n_atoms();
    if (natm <= 0) {
        throw std::invalid_argument(
            "cuda::scf::RHF: device frame has no atoms");
    }
    if (static_cast<int>(host_Z_.size()) != natm) {
        host_Z_.resize(static_cast<std::size_t>(natm));
        auto* stream = static_cast<cudaStream_t>(ctx.stream);
        cudaMemcpyAsync(host_Z_.data(), frame.atoms.type,
                        sizeof(int) * static_cast<std::size_t>(natm),
                        cudaMemcpyDeviceToHost, stream);
        cudaStreamSynchronize(stream);
    }
    const int* h_Z = host_Z_.data();

    last_force_device_ = false;
    if (compute_force && !device_sp_sto3g_ok(basis, h_Z, natm)) {
        throw std::invalid_argument(
            "cuda::scf::RHF::forces: device analytic gradients cover "
            "STO-3G l<=1 (H-O) only. Use cpu::scf::RHF{basis}.forces(frame) "
            "for other bases or higher angular momentum.");
    }

    release_engine_();
    loop_->slot.ctx = &ctx;
    loop_->slot.frame = &frame;
    if (device_sp_sto3g_ok(basis, h_Z, natm)) {
        loop_->slot.cycle = make_device_sp_sto3g_solver(h_Z, natm, charge, tol,
                                                        max_iter, level_shift);
        last_device_sp_ = true;
    } else {
        cpu::Frame host = host_frame_from_device(frame);
        loop_->slot.cycle = make_host_integral_solver(host, basis, charge, tol,
                                                      max_iter, level_shift);
        last_device_sp_ = false;
    }
    loop_->slot.cycle->build(ctx.stream, frame);
}

void RHF::assemble_forces_(cuda::Context& ctx, cuda::Frame& frame) {
    const auto& ws = loop_->slot.cycle->workspace();
    auto* stream = static_cast<cudaStream_t>(ctx.stream);
    const bool device_forces = ws.device_sp_shells && ws.n > 0 && ws.P && ws.C
                               && ws.eig && ws.d_shell_alpha && ws.d_shell_atom
                               && ws.d_shell_L && ws.d_ao_shell && ws.d_ao_cart
                               && ws.d_atom_pos && frame.atoms.fx
                               && frame.atoms.fy && frame.atoms.fz;
    if (!device_forces) {
        throw std::runtime_error(
            "cuda::scf::RHF: compute_force is set but the device "
            "s/p force workspace is incomplete (P/C/eig/shell maps).");
    }
    prec_t* d_W = ws.d_W;
    const bool own_W = (d_W == nullptr);
    const int nn = ws.n * ws.n;
    if (own_W)
        cudaMalloc(&d_W, static_cast<std::size_t>(nn) * sizeof(prec_t));
    launch_build_W_rhf(d_W, ws.C, ws.eig, ws.n, ws.n_occ, stream);
    launch_rhf_forces_sp(frame.atoms.fx, frame.atoms.fy, frame.atoms.fz, ws.P,
                         d_W, ws.d_shell_alpha, ws.d_shell_coeff,
                         ws.d_shell_center, ws.d_shell_atom, ws.d_shell_L,
                         ws.d_ao_shell, ws.d_ao_cart, ws.d_atom_Z,
                         ws.d_atom_pos, ws.n, ws.n_prim, ws.natm, stream);
    last_force_device_ = true;
    if (own_W) {
        cudaStreamSynchronize(stream);
        cudaFree(d_W);
    }
}

void RHF::pull_energy_(cuda::Context& ctx, cuda::Frame& frame) {
    prec_t e_ev = 0;
    if (frame.scalars) {
        auto* stream = static_cast<cudaStream_t>(ctx.stream);
        cudaMemcpyAsync(&e_ev, &frame.scalars->energy, sizeof(prec_t),
                        cudaMemcpyDeviceToHost, stream);
        cudaStreamSynchronize(stream);
    }
    energy_ha_ = static_cast<double>(e_ev) / units::kHartreeInEv;
}

void RHF::step_stages_(cuda::Context& ctx, cuda::Frame& frame, int step) {
    struct Guard {
        bool& flag;
        explicit Guard(bool& f) : flag(f) { flag = true; }
        ~Guard() { flag = false; }
    } guard(this->launching_);
    for (const auto& plan_slot : this->plan().slots) {
        if (plan_slot.ops.empty())
            continue;
        for (const auto& op : plan_slot.ops) {
            if (!op.component)
                continue;
            cuda::LaunchContext lc{ctx, frame, plan_slot.phase, step};
            op.component->launch(lc);
        }
    }
}

int RHF::drive_loop_(cuda::Context& ctx, cuda::Frame& frame, int budget) {
    struct Guard {
        RHF* self;
        ~Guard() { self->release_engine_(); }
    } guard{this};

    assemble_(ctx, frame);
    const bool device = loop_->slot.cycle->device_integrals();
    loop_->slot.cycle->begin(ctx.stream);

    const cuda::LoopServices saved = ctx.loop;
    ctx.loop.status = &this->status_;
    ctx.loop.metrics = &this->metrics_;
    this->status_.reset();
    this->metrics_ = LoopMetrics{};

    const int cap = budget > 0 ? budget : max_iter;
    int iters = 0;
    for (int i = 0; i < cap; ++i) {
        this->mutable_outer_status().reset();
        step_stages_(ctx, frame, i);
        ++iters;
        if (this->poll_handlers_(i, frame, ctx) != ATV::Verdict::Continue)
            break;
    }
    loop_->slot.cycle->finish(ctx.stream);
    if (compute_force)
        assemble_forces_(ctx, frame);
    pull_energy_(ctx, frame);
    this->iterations_ = iters;
    last_device_sp_ = device;
    ctx.loop = saved;
    return iters;
}

void RHF::build_on(cuda::Context&, cuda::Frame&) {
    if (this->plan().slots.empty())
        this->build();
}

int RHF::run_nested_erased(cuda::Context& ctx, cuda::Frame& frame, int budget) {
    return drive_loop_(ctx, frame, budget);
}

void RHF::destroy_on(cuda::Context&, cuda::Frame&) { release_engine_(); }

void RHF::run(cpu::Frame& frame) {
    cache_host_types(frame);
    this->setup_pass(frame);
    try {
        drive_loop_(this->backend_->context(), this->backend_->frame(),
                    max_iter);
    } catch (...) {
        release_engine_();
        this->teardown_pass(frame);
        throw;
    }
    this->teardown_pass(frame);
}

cpu::scf::GradientResult RHF::forces(cpu::Frame& frame) {
    const bool prev = compute_force;
    compute_force = true;
    try {
        run(frame);
    } catch (...) {
        compute_force = prev;
        throw;
    }
    compute_force = prev;

    const int natm = static_cast<int>(frame.nlocal());
    xt::xtensor<double, 2> f = xt::zeros<double>(
        std::array<std::size_t, 2>{static_cast<std::size_t>(natm), 3});
    auto fx = frame.span(kAtom_fx);
    auto fy = frame.span(kAtom_fy);
    auto fz = frame.span(kAtom_fz);
    for (int A = 0; A < natm; ++A) {
        f(A, 0) = static_cast<double>(fx[A]);
        f(A, 1) = static_cast<double>(fy[A]);
        f(A, 2) = static_cast<double>(fz[A]);
    }
    const double e_ev = static_cast<double>(frame.meta.get(kEnergy));
    return cpu::scf::GradientResult{std::move(f), e_ev};
}

const SCFWorkspace& RHF::workspace() const {
    if (!loop_ || !loop_->slot.cycle) {
        throw std::runtime_error(
            "cuda::scf::RHF::workspace: no live SCF workspace. Call run() or "
            "forces(); the workspace is released when that call returns.");
    }
    return loop_->slot.cycle->workspace();
}

} // namespace scf
} // namespace cuda
} // namespace ATV

#endif // ATV_HAS_HF
