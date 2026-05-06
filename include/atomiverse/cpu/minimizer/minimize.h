/// @file cpu/minimizer/minimize.h
/// @brief cpu::minimizer::Minimize — concrete energy-minimizer driver.
///
/// `Minimize` is the user-facing CPU outer-driver under the
/// `OptimizeDriver` category. It iterates an optimizer component
/// (LBFGS / CG / future FIRE) until the registered convergence-check
/// component reports converged, or until `max_iter` is reached.
///
/// Exposes the same LAMMPS-shaped categorical setter API as
/// `cpu::md::MD` (minus the time-evolution setters): `set_pair`,
/// `set_bond`, `set_angle`, `set_dihedral`, `set_neighbor`,
/// `set_optimizer`, `set_observer`, `set_record`, `set_thermo`.
///
/// `add()` and `on(phase).add()` from the underlying `Driver<>`
/// remain public for user-defined hooks.
#pragma once

#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>

#include "atomiverse/core/real.h"
#include "atomiverse/cpu/component.h"
#include "atomiverse/cpu/pair_evaluator.h"
#include "atomiverse/cpu/minimizer/init.h"
#include "atomiverse/cpu/optimize_driver.h"
#include "atomiverse/cpu/observer/thermo_component.h"
#include "atomiverse/cpu/snapshot_publisher.h"
#include "atomiverse/facade/component_role.h"
#include "atomiverse/facade/stage.h"

namespace ATV {
namespace cpu {
namespace minimizer {

/// @brief Energy-minimization driver — convergence-aware Optimize loop.
///
/// Auto-adds `ZeroForces` and `ZeroEnergy` at force phase. Wires the
/// shared `cpu::SnapshotPublisher` into `end_of_step @ +10` so async
/// observers see snapshots.
///
/// Validation rules inherited from `OptimizeDriver`: requires ≥ 1
/// `Optimizer` component (e.g. `LBFGS`, `CG`) and ≥ 1 `Evaluator`
/// with `Capability::Energy` or `Capability::Force` (e.g.
/// `SCFEvaluator`).
class Minimize : public OptimizeDriver {
public:
    Minimize() {
        this->add(ZeroForces{});
        this->add(ZeroEnergy{});
        this->on(ATV::minimizer::Phase::end_of_step)
             .add(cpu::SnapshotPublisher{}, 10);
    }

    // ─── Categorical setters (LAMMPS-shaped) ──────────────────────────

    /// @brief Configure the pair-cutoff force evaluator (singular).
    /// @tparam P Subclass of `cpu::md::PairEvaluator`.
    template <typename P>
    Minimize& set_pair(P&& p) {
        using D = std::decay_t<P>;
        static_assert(std::is_base_of_v<ATV::cpu::PairEvaluator, D>,
                      "set_pair requires a cpu::PairEvaluator subclass");
        this->set_singular_at_(std::forward<P>(p), pair_set_, "set_pair",
                               ATV::minimizer::Phase::force, 0);
        return *this;
    }

    template <typename B>
    Minimize& set_bond(B&& b) {
        this->set_singular_at_(std::forward<B>(b), bond_set_, "set_bond",
                               ATV::minimizer::Phase::force, 0);
        return *this;
    }
    template <typename A>
    Minimize& set_angle(A&& a) {
        this->set_singular_at_(std::forward<A>(a), angle_set_, "set_angle",
                               ATV::minimizer::Phase::force, 0);
        return *this;
    }
    template <typename D>
    Minimize& set_dihedral(D&& d) {
        this->set_singular_at_(std::forward<D>(d), dihedral_set_, "set_dihedral",
                               ATV::minimizer::Phase::force, 0);
        return *this;
    }

    /// @brief Configure the neighbor list. Bound at `force @ -10` so
    /// it runs before force evaluators.
    template <typename N>
    Minimize& set_neighbor(N&& n_comp) {
        this->set_singular_at_(std::forward<N>(n_comp), neighbor_set_, "set_neighbor",
                               ATV::minimizer::Phase::force, -10);
        return *this;
    }

    /// @brief Configure the optimization algorithm (LAMMPS `min_style`).
    template <typename O>
    Minimize& set_optimizer(O&& o) {
        this->set_singular_default_(std::forward<O>(o), optimizer_set_,
                                    "set_optimizer");
        return *this;
    }

    /// @brief Add a derived-quantity observer (plural).
    template <typename Obs>
    Minimize& set_observer(Obs&& obs) {
        this->add(std::forward<Obs>(obs));
        return *this;
    }

    /// @brief Add a recorder / I/O sink (plural).
    template <typename R>
    Minimize& set_record(R&& r) {
        this->add(std::forward<R>(r));
        return *this;
    }

    /// @brief Configure the thermo print cadence (singular). Binds
    /// `cpu::Thermo{every}` at `end_of_step @ +20`.
    Minimize& set_thermo(int every) {
        if (thermo_every_) {
            if (*thermo_every_ == every) return *this;
            throw std::invalid_argument(
                "set_thermo already configured with " +
                std::to_string(*thermo_every_) +
                ", rejected re-call with " + std::to_string(every));
        }
        thermo_every_ = every;
        this->on(ATV::minimizer::Phase::end_of_step)
             .add(ATV::cpu::Thermo{every}, +20);
        return *this;
    }

protected:
    /// @brief Missing-NL throw only. (No Stress auto-injection on
    /// `Minimize` — variable-cell minimization is not yet supported.
    /// MPI halo on `Minimize` waits for a real multi-rank optimization
    /// workflow.)
    void install_infra() override {
        const auto s = this->scan_pending_infra_();
        this->throw_if_missing_neighbor_list_(
            s, "Minimize",
            "Add `opt.set_neighbor(cpu::md::NeighborList<>{cutoff, skin})`.");
    }

private:
    std::optional<std::string> pair_set_;
    std::optional<std::string> bond_set_;
    std::optional<std::string> angle_set_;
    std::optional<std::string> dihedral_set_;
    std::optional<std::string> neighbor_set_;
    std::optional<std::string> optimizer_set_;
    std::optional<int>         thermo_every_;
};

} // namespace minimizer
} // namespace cpu
} // namespace ATV
