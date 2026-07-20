/// @file cpu/scf/mean_field_method.h
/// @brief Shared CRTP (curiously recurring template pattern) base for the
///        concrete Gaussian mean-field methods: `MeanFieldMethod<D>` (base of
///        RHF / UHF / RKS / UKS).
///
/// There is no bag named `SCFSolver` and no user-facing `SCFConfig` /
/// `IntegralPack` step. Each method is a concrete class; parameters are
/// **members**; the only I/O currency is `Frame`.
///
/// SCF self-consistency is a **Driver loop**: `MeanFieldMethod<D>` derives
/// `cpu::CPUDriver<ATV::scf::Phase>` and drives the five sub-spec-02 stage
/// components (Fock / Mix / Diagonalize / Density / ConvergenceMeasure) through
/// `fock → mix → diagonalize → density → check`, budget = `max_iter`, with the
/// default `ScfConverged` criterion registered via `add_criterion`. All entries
/// (both `run` overloads, `forces`, the component `build`/`launch` path,
/// `mp2`/`fci`) funnel through one private member `drive_scf_`. Integrals and AO
/// workspaces are private. The concrete methods live in `rhf.h` / `uhf.h` /
/// `rks.h` / `uks.h`; `methods.h` is their umbrella.
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
#include <utility>
#include <vector>
#include <xtensor/containers/xtensor.hpp>

#include "atomiverse/cpu/component/hf.h"
#include "atomiverse/cpu/driver/cpu_driver.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/scf/basis_set.h"
#include "atomiverse/cpu/scf/dense_matrix.h"
#include "atomiverse/cpu/scf/fci.h"
#include "atomiverse/cpu/scf/grid.h"
#include "atomiverse/cpu/scf/mp2.h"
#include "atomiverse/cpu/scf/scf_result.h"
#include "atomiverse/cpu/scf/stages/density_stage.h"
#include "atomiverse/cpu/scf/stages/diagonalize_stage.h"
#include "atomiverse/cpu/scf/stages/fock_stage.h"
#include "atomiverse/cpu/scf/stages/mix_stage.h"
#include "atomiverse/cpu/scf/stages/scf_converged.h"
#include "atomiverse/cpu/scf/xc_functional.h"

namespace ATV {
namespace cpu {
namespace scf {

/// @brief Shared view for observers/recorders after a method launch.
///
/// `integrals` is opaque private workspace published only so cube writers can
/// read shells; users never construct that type. Spin treatment is read from
/// the result (`result->spin`) — the handle no longer snapshots a `SCFConfig`.
struct SCFResultHandle {
    std::shared_ptr<const SCFResult> result;
    std::shared_ptr<const IntegralPack> integrals;
    std::shared_ptr<const bool> ready;
};

/// @brief Shared members and Frame API for concrete mean-field methods.
///
/// CRTP (`Derived` = `RHF`/`UHF`/`RKS`/`UKS`) so chaining
/// (`RHF{}.predictor(3)`) keeps the concrete type for `Driver::add`.
///
/// Derives BOTH `cpu::CPUDriver<ATV::scf::Phase>` (the self-consistency loop
/// node) and `HF` (still a Component, so an outer `MD`/`Evaluate` can bind it at
/// its force/energy phase and call `launch`). `run` (driver) and `launch`
/// (component) COEXIST this sub-spec — a deliberate, temporary break of the
/// "never coexist" rule; both funnel into `drive_scf_`. Sub-spec 04 deletes
/// `build`/`launch` and narrows `HF` to a pure tag.
template <typename Derived>
class MeanFieldMethod : public cpu::CPUDriver<ATV::scf::Phase>, public HF {
  public:
    /// Un-hide the loop node's no-argument Plan `build()` (the Component base's
    /// virtual `build(const BuildContext&)` below would otherwise hide it).
    using cpu::CPUDriver<ATV::scf::Phase>::build;

    // ── Geometry / basis (inputs from Frame + these members) ──────────
    BasisSet basis = BasisSet::STO3G;
    int charge = 0;

    // ── Convergence (members, not a Config bag) ───────────────────────
    // "DIIS" below = direct inversion of the iterative subspace: the standard
    // SCF convergence accelerator that extrapolates the Fock matrix from a short
    // history of the previous iterations' error vectors.
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

    /// @brief Run SCF on `frame` (driver loop). Integrals are private workspace.
    /// @param frame  Geometry (Å) + Z in type column.
    /// @param guess  Optional density guess.
    /// @param write_energy  If true, write eV energy to `frame.meta[kEnergy]`.
    /// @return Reference to owned `last_result()` (Hartree internals).
    SCFResult& run(Frame& frame, const SCFInitialGuess* guess = nullptr,
                   bool write_energy = false);
    /// @brief Const-geometry energy (does not write Frame meta or move atoms).
    /// @param frame  Geometry; left untouched.
    /// @param guess  Optional density guess.
    /// @return Reference to owned `last_result()` (Hartree internals).
    SCFResult& run(const Frame& frame, const SCFInitialGuess* guess = nullptr);

    /// @brief Analytic forces [eV/Å] after (or including) a fresh SCF.
    /// @param frame  Geometry; forces written to fx/fy/fz and energy to meta.
    /// @return Forces + energy in eV / eV·Å⁻¹.
    GradientResult forces(Frame& frame);

    /// @brief Closed-shell MP2 on top of a fresh RHF-style run (materializes ERI).
    /// @param frame  Geometry; drives a fresh SCF then adds the correction.
    /// @return MP2 correlation + total energy.
    /// @pre Derived method is closed-shell HF (use `RHF`).
    MP2Result mp2(Frame& frame);

    /// @brief Tiny full-CI on top of a fresh RHF-style run (materializes ERI).
    /// @param frame  Geometry; drives a fresh SCF then adds the correction.
    /// @return Full-CI correlation + total energy.
    FCIResult fci(Frame& frame);

    /// @brief Most recent converged (or last) SCF result [Hartree internals].
    const SCFResult& last_result() const { return *result_sp_; }
    /// @brief Whether a result has been produced by any run/launch.
    bool has_result() const { return *has_result_sp_; }

    /// @brief Observer handle (shared_ptrs survive driver ownership moves).
    /// @return Shared view of result / integrals for cube writers.
    SCFResultHandle result_handle() const {
        return {result_sp_, integrals_sp_, has_result_sp_};
    }

    /// @brief Component path — one-time setup (debt; 04 deletes build/launch).
    /// @param ctx Build context (unused; drive_scf_ owns per-run setup).
    void build(const cpu::BuildContext& ctx) override;
    /// @brief Component path — SCF (+ optional forces) each force step.
    /// @param ctx Launch context carrying the frame; funnels into drive_scf_.
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

    /// @brief Stages are fixed in `ensure_driver_built_`; nothing is
    ///        user-suppliable, so there is nothing to validate.
    void validate_components() const override {}

    /// @brief Fixed by the concrete subclass (RHF/UHF/RKS/UKS).
    SpinType spin_ = SpinType::Restricted;
    MeanFieldType method_ = MeanFieldType::HF;
    XCType xc_type_ = XCType::LSDA;
    int n_alpha_ = 0;
    int n_beta_ = 0;

    /// @brief Set n_alpha/n_beta from electron count (and subclass rules).
    virtual void resolve_occupation(int n_electrons) = 0;

  private:
    // ── Per-run SCF iteration state — plainly-named private members ────────────
    //
    // Not a state bag: each field is wired into a stage component by ctor
    // reference in `ensure_driver_built_` (const& reads, & writes), per the
    // driver-world ruling against `*Workspace` state bags. xtensor identity is
    // stable, so the per-run resize inside `drive_scf_` never invalidates a
    // bound reference; a move of the method drops the wiring
    // (`reset_built_driver_`) so the next `drive_scf_` re-binds against the
    // moved-to object's members — no dangling reference.
    RealMatrix S_;                ///< Overlap matrix (dimensionless).
    RealMatrix Hcore_;            ///< Core Hamiltonian [Hartree].
    RealMatrix X_;                ///< Orthogonalization transform.
    RealMatrix F_alpha_;          ///< α Fock matrix [Hartree].
    RealMatrix F_beta_;           ///< β Fock matrix [Hartree].
    RealMatrix J_;                ///< Coulomb matrix [Hartree].
    RealMatrix err_alpha_;        ///< α DIIS error channel.
    RealMatrix err_beta_;         ///< β DIIS error channel.
    RealMatrix C_alpha_;          ///< α MO coefficients.
    RealMatrix C_beta_;           ///< β MO coefficients.
    RealVector eig_alpha_;        ///< α orbital energies [Hartree].
    RealVector eig_beta_;         ///< β orbital energies [Hartree].
    RealMatrix P_alpha_;          ///< α density (old-in / new-out).
    RealMatrix P_beta_;           ///< β density (old-in / new-out).
    double e_total_ = 0;          ///< Total SCF energy [Hartree].
    double e_one_ = 0;            ///< One-electron energy [Hartree].
    double e_coulomb_ = 0;        ///< Coulomb (J) energy [Hartree].
    double e_xc_grid_ = 0;        ///< Grid XC energy [Hartree].
    double e_exact_exchange_ = 0; ///< Hybrid exact-exchange energy [Hartree].
    double int_n_alpha_ = 0;      ///< Grid-integrated α electron count.
    double int_n_beta_ = 0;       ///< Grid-integrated β electron count.
    double v_nuc_ = 0;            ///< Nuclear-repulsion energy [Hartree].
    int n_basis_ = 0;             ///< Basis dimension.
    int n_dropped_ = 0;           ///< Overlap eigen-directions dropped.
    int iter_ = 0;                ///< Current SCF iteration index.
    int diis_pairs_ = 0;          ///< DIIS history-size mirror.
    bool diis_active_ = false;    ///< DIIS extrapolating (≥ 2 vectors).
    bool restricted_ = false;     ///< spin_ == Restricted (per-run cache).
    bool is_hf_ = false;          ///< method_ == HF (per-run cache).

    /// @brief Options selecting which of the shared entry behaviors apply.
    ///        Private, per-call — routing flags, not a config bag.
    struct DriveOpts {
        bool write_energy = false; ///< Write eV energy to `frame.meta`.
        bool forces = false;       ///< Compute + write analytic forces.
        bool unwrap = true;        ///< Unwrap positions (mutating paths only).
        int max_iter = -1;         ///< Budget override; <0 ⇒ member `max_iter`.
        double tol_override = -1;  ///< |ΔE| tol override; <0 ⇒ member `tol`.
    };

    /// @brief RAII guard: unwrap positions for the integral build on entry,
    ///        restore the wrapped storage on exit (absorbs the old free
    ///        `unwrap_for_scf`/`restore_wrapped`; no free helper functions).
    class PositionGuard {
      public:
        PositionGuard(Frame& frame, bool active);
        ~PositionGuard();
        PositionGuard(const PositionGuard&) = delete;
        PositionGuard& operator=(const PositionGuard&) = delete;

      private:
        Frame& frame_;
        bool active_ = false;
        std::vector<StorageReal> x_, y_, z_;
    };

    /// @brief THE funnel: prologue (per-run assembly) → CPUDriver loop →
    ///        epilogue (de-shift / spin / result fill / energy / forces).
    /// @param frame Geometry; mutated only on the mutating entries.
    /// @param guess Optional density guess.
    /// @param opts  Which shared behaviors apply on this entry.
    /// @return Reference to the owned `SCFResult`.
    SCFResult& drive_scf_(Frame& frame, const SCFInitialGuess* guess,
                          const DriveOpts& opts);

    /// @brief Lazily wire the five stages + `ScfConverged` and compile the Plan.
    ///        Deferred to first `drive_scf_` so member config (grid, level
    ///        shift, xc) is finalized before the stages bind references to the
    ///        method's members. Idempotent, and self-healing after a move: if
    ///        the object was relocated since it wired (detected via
    ///        `wired_self_`), it drops the stale wiring and rebinds to `this`.
    void ensure_driver_built_();

    /// @brief Drop the wired stages + their Plan / criterion registration and
    ///        clear `driver_built_`, so the next `drive_scf_` rebuilds. Invoked
    ///        when a move leaves the stages bound to a stale address.
    void reset_built_driver_();

    // Shared workspace (also published via result_handle for observers).
    std::shared_ptr<IntegralPack> integrals_sp_{
        std::make_shared<IntegralPack>()};
    std::shared_ptr<SCFResult> result_sp_{std::make_shared<SCFResult>()};
    std::shared_ptr<bool> has_result_sp_{std::make_shared<bool>(false)};

    // Stage components (heap-owned; borrowed into the Plan by add_ref). Null
    // until ensure_driver_built_ wires them to this method's members.
    std::unique_ptr<FockStage> fock_stage_;
    std::unique_ptr<MixStage> mix_stage_;
    std::unique_ptr<DiagonalizeStage> diag_stage_;
    std::unique_ptr<DensityStage> density_stage_;
    std::unique_ptr<ScfConverged> converged_stage_;
    bool driver_built_ = false;
    /// Address the stages were wired to; a mismatch on `this` flags a move.
    const void* wired_self_ = nullptr;

    // Most recent analytic gradient (forces() return + component force write).
    GradientResult last_grad_;

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

} // namespace scf
} // namespace cpu
} // namespace ATV
