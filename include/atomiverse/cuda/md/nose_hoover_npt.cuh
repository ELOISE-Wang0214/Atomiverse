/// @file cuda/md/nose_hoover_npt.cuh
/// @brief CUDA Martyna-Tobias-Klein NPT integrator (LAMMPS-equivalent
///        `fix npt`).
///
/// CUDA sibling of `cpu::md::NoseHooverNPT`. Same algorithm (Tuckerman
/// 2006 / Shinoda 2004 MTK) and Trotter ordering (P1–P11), same
/// LAMMPS-aligned `mtk_term1` / `mtk_term2`. The barostat is a single
/// chain over the active `ε̇_α` set with per-axis `couple_axis[]`
/// guards; non-coupled axes have `ε̇_α ≡ 0`. Pressure-coupling convention
/// matches LAMMPS `fix npt`: positive ⇔ compressive; per-axis driving
/// pressure `P_α = (σ_αα + 2·KE_α)/V` averaged over coupled axes when
/// `pdim==3` (the legacy isotropic case is bit-identical at
/// equipartition).
///
/// **Replaces (do not add alongside)** `cuda::md::VelocityVerlet` +
/// `cuda::md::NoseHoover` + `cuda::md::Barostat`. Rewrites
/// `frame.scalars->kinetic` and `frame.scalars->conserved` at
/// `end_of_step` so downstream observers (`Thermo`, recorders) see
/// post-NPT values.
///
/// Single-GPU scope (`msu-pe-validation` fixed-scope #1). Pressure-sign
/// is LAMMPS-standard (same as `cuda::VirialStress` /
/// `frame.scalars->stress`) — NOT the inverted convention used by
/// `cuda::md::Barostat` (Berendsen). Trotter ordering across the three
/// MD phases mirrors `cpu::md::NoseHooverNPT::launch`:
///
/// - `initial_integrate` (P1–P6): nhc_baro_half → nhc_atoms_half →
///   update_eps_dot → scale_v_by_baro → kick_v_by_force → propagate_x.
/// - `final_integrate` (P7): kick_v_by_force.
/// - `end_of_step` (P8–P11, after ReduceEnergy): recompute KE_global →
///   scale_v_by_baro → update_eps_dot → nhc_atoms_half (reversed
///   atoms-then-baro) → nhc_baro_half → write kinetic + conserved.
///
/// @see cpu::md::NoseHooverNPT for the CPU reference + full algorithm
///      documentation block.
/// @see LAMMPS `fix_nh.cpp` lines 2541-2569 (`mtk_term1` / `mtk_term2`).
#pragma once

#include <array>
#include <cstdint>
#include <string>
#include <vector>

#include "atomiverse/cuda/component.cuh"
#include "atomiverse/facade/component_role.h"
#include "atomiverse/facade/stage.h"

namespace ATV {
namespace cuda {
namespace md {

/// @brief CUDA MTK NPT integrator. Same API surface and Trotter ordering
///        as `cpu::md::NoseHooverNPT`.
class NoseHooverNPT : public ComponentBase {
  public:
    /// @brief Configuration aggregate (mirror of `cpu::md::NoseHooverNPT::Config`).
    struct Config {
        double T_target;   ///< Target temperature [K].
        double tau_T;      ///< Thermostat coupling time [fs] (LAMMPS `Tdamp`).
        double P_target;   ///< Target pressure [eV/Å³] (LAMMPS sign, +).
        double tau_P;      ///< Barostat coupling time [fs] (LAMMPS `Pdamp`).
        double dt;         ///< Timestep [fs].
        double drag = 0.0; ///< Empirical NHC drag (LAMMPS `drag`).
        int mtchain = 3;   ///< Atom-NHC chain length.
        int mpchain = 3;   ///< Baro-NHC chain length.
        std::array<bool, 3> couple_axis = {true, true, true};
        double T_stop = -1.0; ///< Ramp end T [K]; <0 folds to T_target.
    };

    explicit NoseHooverNPT(Config cfg);

    /// @brief Set the temperature-ramp run span `[nfirst, nlast]`.
    /// @return `*this`, chainable.
    NoseHooverNPT& set_run_span(std::int64_t nfirst, std::int64_t nlast) {
        nfirst_ = nfirst;
        nlast_ = nlast;
        return *this;
    }

    /// @return Live thermostat setpoint [K] — for observers.
    double T_target() const noexcept { return T_live_; }

    std::string name() const override { return "nose_hoover_npt"; }
    ComponentRole role() const override { return ComponentRole::Integrator; }
    /// @brief Reads — consumes `Stress` to drive ε̇ via `update_eps_dot_`.
    ///        Triggers `cuda::md::MD::install_infra_` to auto-inject
    ///        `cuda::VirialStress` if none is registered. The Integrator
    ///        role alone could not express this dependency.
    Capability reads() const override { return Capability::Stress; }

    std::vector<DefaultBinding> default_bindings() const override {
        return {{static_cast<int>(ATV::md::Phase::initial_integrate), 0},
                {static_cast<int>(ATV::md::Phase::final_integrate), 0},
                {static_cast<int>(ATV::md::Phase::end_of_step), -40}};
    }

    void build(const BuildContext& ctx) override;
    void launch(const LaunchContext& ctx) override;
    void destroy(const DestroyContext& ctx) override;

  private:
    // Configuration
    double T_target_, T_stop_;
    double tau_T_, P_target_, tau_P_, dt_, drag_;
    int mtchain_, mpchain_, pdim_ = 3;
    std::array<bool, 3> couple_axis_;
    std::int64_t nfirst_ = 0, nlast_ = 0;

    // Derived in build()
    int n_dof_ = 0;
    int n_atoms_ = 0;
    double kT_ = 0.0;
    double T_live_ = 0.0;
    double tdrag_factor_ = 1.0;
    double pdrag_factor_ = 1.0;
    std::vector<double> Q_t_, Q_b_;
    double W_ = 0.0;

    // NHC state (host-side; chains are short, M ≤ a few)
    std::vector<double> eta_t_, eta_t_dot_;
    std::vector<double> eta_b_, eta_b_dot_;
    std::array<double, 3> eps_dot_ = {0.0, 0.0, 0.0};

    // Host-cached box diagonal (read once at build, updated after each
    // propagate_x_and_box). Used by the host helpers to compute volume
    // without a per-step D2H of the box.
    double box_lx_ = 0.0, box_ly_ = 0.0, box_lz_ = 0.0;

    // Cached global kinetic energy carried between phases.
    double KE_global_ = 0.0;

    // Device scratch (cudaMalloc'd in build()).
    void* d_ke_buf_ = nullptr;  // 4 doubles: [0]=global, [1..3]=per-axis

    // Helpers — see cpu::md::NoseHooverNPT for the canonical algorithm.
    double compute_global_KE_(const LaunchContext& ctx);
    std::array<double, 3> compute_axis_KE_(const LaunchContext& ctx);
    void compute_temp_target_(int step);
    void nhc_atoms_half_(const LaunchContext& ctx, double half_dt);
    void nhc_baro_half_(double half_dt);
    void update_eps_dot_(const LaunchContext& ctx, double half_dt);
    void scale_v_by_baro_(const LaunchContext& ctx, double half_dt);
    void kick_v_by_force_(const LaunchContext& ctx, double half_dt);
    void propagate_x_and_box_(const LaunchContext& ctx);
    double extended_hamiltonian_(const LaunchContext& ctx) const;
};

} // namespace md
} // namespace cuda
} // namespace ATV
