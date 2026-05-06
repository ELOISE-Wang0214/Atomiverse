/// @file cpu/evaluate.h
/// @brief cpu::Evaluate — concrete one-pass evaluator driver.
///
/// `cpu::Evaluate` is the user-facing CPU outer-driver under the
/// `EvaluateDriver` category. It runs the registered evaluator
/// components exactly once on the current `cpu::Frame` and returns —
/// no time advance, no minimization loop.
///
/// Single-point energy / force / property workflows go here. Exposes
/// the same LAMMPS-shaped categorical setter API as `cpu::md::MD`
/// (minus the time-evolution and optimizer setters): `set_pair`,
/// `set_bond`, `set_angle`, `set_dihedral`, `set_neighbor`,
/// `set_observer`, `set_record`. `add()` and `on(phase).add()` from
/// the underlying `Driver<>` remain public for user-defined hooks
/// (e.g. SCF evaluators registered via `add(SCFEvaluator{...})`).
///
/// Auto-adds `cpu::md::ZeroForces` and `cpu::md::ZeroEnergy` at
/// `evaluate::Phase::prepare` (orders -1000 / -999) so frame fields
/// start clean before evaluators write to them. The MD-flavoured Zero*
/// components are reused via explicit `on(phase).add()` binding —
/// their launch bodies are phase-agnostic (they just zero the frame
/// fields).
///
/// @code{.cpp}
///   // SCF energy + forces (default)
///   cpu::Evaluate eval;
///   eval.add(cpu::scf::SCFEvaluator{cpu::scf::default_config(integrals)});
///   eval.build();
///   eval.run(frame);
///
///   // Force-field force evaluation (LAMMPS-shaped)
///   cpu::Evaluate eval2;
///   eval2.set_pair(cpu::md::LJ{eps, sigma, cutoff});
///   eval2.set_neighbor(cpu::md::NeighborList<>{cutoff, skin});
///   eval2.build();
///   eval2.run(frame);
/// @endcode
#pragma once

#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>

#include "atomiverse/cpu/evaluate_driver.h"
#include "atomiverse/cpu/md/init.h"
#include "atomiverse/cpu/pair_evaluator.h"
#include "atomiverse/facade/stage.h"

namespace ATV {
namespace cpu {

/// @brief One-pass evaluator driver.
class Evaluate : public EvaluateDriver {
public:
    Evaluate() {
        // Override the MD-flavoured default_bindings (which target
        // md::Phase::pre_force) by binding explicitly to evaluate::Phase.
        this->on(ATV::evaluate::Phase::prepare).add(cpu::md::ZeroForces{}, -1000);
        this->on(ATV::evaluate::Phase::prepare).add(cpu::md::ZeroEnergy{},  -999);
    }

    // ─── Categorical setters (LAMMPS-shaped) ──────────────────────────

    /// @brief Configure the pair-cutoff force evaluator (singular).
    /// @tparam P Subclass of `cpu::md::PairEvaluator`.
    template <typename P>
    Evaluate& set_pair(P&& p) {
        using D = std::decay_t<P>;
        static_assert(std::is_base_of_v<ATV::cpu::PairEvaluator, D>,
                      "set_pair requires a cpu::PairEvaluator subclass");
        this->set_singular_at_(std::forward<P>(p), pair_set_, "set_pair",
                               ATV::evaluate::Phase::evaluate, 0);
        return *this;
    }

    template <typename B>
    Evaluate& set_bond(B&& b) {
        this->set_singular_at_(std::forward<B>(b), bond_set_, "set_bond",
                               ATV::evaluate::Phase::evaluate, 0);
        return *this;
    }
    template <typename A>
    Evaluate& set_angle(A&& a) {
        this->set_singular_at_(std::forward<A>(a), angle_set_, "set_angle",
                               ATV::evaluate::Phase::evaluate, 0);
        return *this;
    }
    template <typename D>
    Evaluate& set_dihedral(D&& d) {
        this->set_singular_at_(std::forward<D>(d), dihedral_set_, "set_dihedral",
                               ATV::evaluate::Phase::evaluate, 0);
        return *this;
    }

    /// @brief NeighborList runs earlier in `evaluate::Phase::evaluate`
    /// to populate the pair list before pair / bonded forces.
    template <typename N>
    Evaluate& set_neighbor(N&& n_comp) {
        this->set_singular_at_(std::forward<N>(n_comp), neighbor_set_, "set_neighbor",
                               ATV::evaluate::Phase::evaluate, -10);
        return *this;
    }

    template <typename Obs>
    Evaluate& set_observer(Obs&& obs) {
        this->add(std::forward<Obs>(obs));
        return *this;
    }

    template <typename R>
    Evaluate& set_record(R&& r) {
        this->add(std::forward<R>(r));
        return *this;
    }

protected:
    /// @brief Missing-NL throw only. No Stress auto-injection (no
    /// integrator). No MPI halo (single-rank only today).
    void install_infra() override {
        const auto s = this->scan_pending_infra_();
        this->throw_if_missing_neighbor_list_(
            s, "Evaluate",
            "Add `eval.set_neighbor(cpu::md::NeighborList<>{cutoff, skin})`.");
    }

private:
    std::optional<std::string> pair_set_;
    std::optional<std::string> bond_set_;
    std::optional<std::string> angle_set_;
    std::optional<std::string> dihedral_set_;
    std::optional<std::string> neighbor_set_;
};

} // namespace cpu
} // namespace ATV
