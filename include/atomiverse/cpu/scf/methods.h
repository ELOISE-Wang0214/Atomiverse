/// @file cpu/scf/methods.h
/// @brief Concrete Gaussian mean-field methods: RHF, UHF, RKS, UKS.
///
/// There is no bag named `SCFSolver` and no user-facing `SCFConfig` /
/// `IntegralPack` step. Each method is a concrete class; parameters are
/// **members**; the only I/O currency is `Frame`.
///
/// SCF self-consistency is a loop owned by the method (Driver-shaped).
/// Integrals and AO workspaces are private.
///
/// @code
///   cpu::scf::RHF rhf;
///   rhf.basis = BasisSet::STO3G;
///   rhf.tol = 1e-10;
///   rhf.run(frame);                    // energy / orbitals on last_result()
///   // forces: rhf.with_gradient(true) then use as HF component, or
///   //         GradientResult g = rhf.forces(frame);
/// @endcode
#pragma once

#include <functional>
#include <memory>
#include <string>
#include <utility>
#include <xtensor/containers/xtensor.hpp>

#include "atomiverse/cpu/component/hf.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/scf/basis_set.h"
#include "atomiverse/cpu/scf/fci.h"
#include "atomiverse/cpu/scf/grid.h"
#include "atomiverse/cpu/scf/mp2.h"
#include "atomiverse/cpu/scf/scf_result.h"
#include "atomiverse/cpu/scf/xc_functional.h"

namespace ATV {
namespace cpu {
namespace scf {

// Forward — private loop / workspace (not user types).
class ScfLoop;
struct IntegralPack;
struct SCFConfig;

/// @brief Shared view for observers/recorders after a method launch.
///
/// `integrals` / `config` are opaque private workspace published only so
/// cube writers can read shells. Users never construct those types.
struct SCFResultHandle {
    std::shared_ptr<const SCFResult> result;
    std::shared_ptr<const IntegralPack>
        integrals; // incomplete outside internal
    std::shared_ptr<const bool> ready;
    std::shared_ptr<const SCFConfig> config; // incomplete outside internal
};

/// @brief Shared members and Frame API for concrete mean-field methods.
///
/// CRTP (`Derived` = `RHF`/`UHF`/`RKS`/`UKS`) so chaining
/// (`RHF{}.predictor(3)`) keeps the concrete type for `Driver::add`.
template <typename Derived>
class MeanFieldMethod : public HF {
  public:
    // ── Geometry / basis (inputs from Frame + these members) ──────────
    BasisSet basis = BasisSet::STO3G;
    int charge = 0;

    // ── Convergence (members, not a Config bag) ───────────────────────
    double tol = 1e-9;         ///< |ΔE| threshold [Hartree].
    double conv_tol_grad = -1; ///< DIIS-error tol; <0 → sqrt(tol).
    int max_iter = 100;
    bool use_diis = true;
    int diis_max_vectors = 8;
    int diis_start_iter = 1;
    double damping_factor = 0.0;
    InitGuess init_guess = InitGuess::MINAO;
    double level_shift = 0.0;
    double lindep_threshold = 1e-8;

    // ── KS grid (ignored for pure HF) ─────────────────────────────────
    GridLevel grid_level = GridLevel::SG0;
    int n_radial = 50;
    int n_angular = 110;

    // ── Debug / regression knobs ──────────────────────────────────────
    bool force_native_two_electron_deriv = false;
    bool xc_grad_project_translation = false;
    double hybrid_alpha_override = -1.0;
    const MolecularGrid* fixed_grid = nullptr;

    /// @brief When true, materialize the dense ERI tensor (MP2/FCI/CUDA).
    bool compute_eri = false;

    /// @brief When true (default), component `launch()` writes forces.
    bool compute_force = true;

    MeanFieldMethod();
    ~MeanFieldMethod() override;
    MeanFieldMethod(MeanFieldMethod&&) noexcept;
    MeanFieldMethod& operator=(MeanFieldMethod&&) noexcept;
    MeanFieldMethod(const MeanFieldMethod&) = delete;
    MeanFieldMethod& operator=(const MeanFieldMethod&) = delete;

    /// @brief Run SCF on `frame`. Integrals are private workspace.
    /// @param frame  Geometry (Å) + Z in type column.
    /// @param guess  Optional density guess.
    /// @param write_energy  If true, write eV energy to `frame.meta[kEnergy]`.
    /// @return Reference to owned `last_result()` (Hartree internals).
    SCFResult& run(Frame& frame, const SCFInitialGuess* guess = nullptr,
                   bool write_energy = false);
    /// @brief Const-geometry energy (does not write Frame meta).
    SCFResult& run(const Frame& frame, const SCFInitialGuess* guess = nullptr);

    /// @brief Analytic forces [eV/Å] after (or including) a fresh SCF.
    /// @param frame  Geometry; forces written to fx/fy/fz and energy to meta.
    GradientResult forces(Frame& frame);

    /// @brief Closed-shell MP2 on top of a fresh RHF-style run (materializes ERI).
    /// @pre Derived method is closed-shell HF (use `RHF`).
    MP2Result mp2(Frame& frame);

    /// @brief Tiny full-CI on top of a fresh RHF-style run (materializes ERI).
    FCIResult fci(Frame& frame);

    const SCFResult& last_result() const { return *result_sp_; }
    bool has_result() const { return *has_result_sp_; }

    /// @brief Observer handle (shared_ptrs survive driver ownership moves).
    SCFResultHandle result_handle() const {
        return {result_sp_, integrals_sp_, has_result_sp_, config_sp_};
    }

    /// @brief Component path — size loop from a representative geometry.
    void build(const cpu::BuildContext& ctx) override;
    /// @brief Component path — SCF (+ optional forces) each force step.
    void launch(const cpu::LaunchContext& ctx) override;

    Derived& with_gradient(bool enabled) & noexcept {
        compute_force = enabled;
        return static_cast<Derived&>(*this);
    }
    Derived&& with_gradient(bool enabled) && noexcept {
        compute_force = enabled;
        return static_cast<Derived&&>(*this);
    }

    Derived& on_evaluate(std::function<void(const SCFResult&)> cb) & {
        on_evaluate_ = std::move(cb);
        return static_cast<Derived&>(*this);
    }
    Derived&& on_evaluate(std::function<void(const SCFResult&)> cb) && {
        on_evaluate_ = std::move(cb);
        return static_cast<Derived&&>(*this);
    }

    Derived& predictor(int order = 3, int corrector_steps = 2,
                       double corrector_tol = 1e-4, bool fallback = true) & {
        set_predictor_(order, corrector_steps, corrector_tol, fallback);
        return static_cast<Derived&>(*this);
    }
    Derived&& predictor(int order = 3, int corrector_steps = 2,
                        double corrector_tol = 1e-4, bool fallback = true) && {
        set_predictor_(order, corrector_steps, corrector_tol, fallback);
        return static_cast<Derived&&>(*this);
    }

  protected:
    void set_predictor_(int order, int corrector_steps, double corrector_tol,
                        bool fallback);
    /// @brief Fixed by the concrete subclass (RHF/UHF/RKS/UKS).
    SpinType spin_ = SpinType::Restricted;
    MeanFieldType method_ = MeanFieldType::HF;
    XCType xc_type_ = XCType::LSDA;
    int n_alpha_ = 0;
    int n_beta_ = 0;

    /// @brief Set n_alpha/n_beta from electron count (and subclass rules).
    virtual void resolve_occupation(int n_electrons) = 0;

    /// @brief Snapshot members into the private loop config.
    void sync_loop_config_();

  private:
    void ensure_loop_();
    void rebuild_integrals_(const Frame& frame);
    SCFResult& run_on_pack_(const SCFInitialGuess* guess);
    void publish_config_();

    std::unique_ptr<ScfLoop> loop_;
    std::unique_ptr<ScfLoop> corrector_loop_;
    bool loop_built_ = false;

    // Shared workspace (also published via result_handle for observers).
    std::shared_ptr<IntegralPack> integrals_sp_{
        std::make_shared<IntegralPack>()};
    std::shared_ptr<SCFResult> result_sp_{std::make_shared<SCFResult>()};
    std::shared_ptr<bool> has_result_sp_{std::make_shared<bool>(false)};
    std::shared_ptr<SCFConfig> config_sp_{std::make_shared<SCFConfig>()};

    // ASPC (component path)
    bool use_predictor_ = false;
    int corrector_steps_ = 2;
    double corrector_tol_ = 1e-4;
    bool fallback_full_scf_ = true;
    std::unique_ptr<class ASPCPredictor> aspc_alpha_;
    std::unique_ptr<class ASPCPredictor> aspc_beta_;
    bool has_last_density_ = false;
    xt::xtensor<double, 2> last_P_alpha_;
    xt::xtensor<double, 2> last_P_beta_;
    std::function<void(const SCFResult&)> on_evaluate_;
};

/// @brief Restricted Hartree–Fock.
class RHF : public MeanFieldMethod<RHF> {
  public:
    RHF();
    std::string name() const override { return "rhf"; }

  protected:
    void resolve_occupation(int n_electrons) override;
};

/// @brief Unrestricted Hartree–Fock.
class UHF : public MeanFieldMethod<UHF> {
  public:
    int n_alpha = 0; ///< Alpha electrons (must be set before run).
    int n_beta = 0;  ///< Beta electrons.

    UHF();
    explicit UHF(int n_alpha, int n_beta);
    std::string name() const override { return "uhf"; }

  protected:
    void resolve_occupation(int n_electrons) override;
};

/// @brief Restricted Kohn–Sham (default PBE / MatchPySCF grid).
class RKS : public MeanFieldMethod<RKS> {
  public:
    XCType xc = XCType::PBE;

    RKS();
    explicit RKS(XCType xc);
    std::string name() const override { return "rks"; }

  protected:
    void resolve_occupation(int n_electrons) override;
};

/// @brief Unrestricted Kohn–Sham.
class UKS : public MeanFieldMethod<UKS> {
  public:
    int n_alpha = 0;
    int n_beta = 0;
    XCType xc = XCType::LSDA;

    UKS();
    UKS(int n_alpha, int n_beta, XCType xc = XCType::LSDA);
    std::string name() const override { return "uks"; }

  protected:
    void resolve_occupation(int n_electrons) override;
};

// ── Compatibility aliases during migration (component name) ───────────
// Historical `cpu::scf::SCF` was a method-agnostic bag. Prefer RHF/UHF/RKS/UKS.
using SCF = RKS; // production default: closed-shell KS/PBE

} // namespace scf
} // namespace cpu
} // namespace ATV
