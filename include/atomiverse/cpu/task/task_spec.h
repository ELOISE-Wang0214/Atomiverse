/// @file cpu/task/task_spec.h
/// @brief TaskSpec -- the versioned, in-memory description of one dispatchable
///        task (mpi-launcher-embedded-python-02).
///
/// A TaskSpec is the byte-broadcastable payload of the launcher's task
/// protocol: it enumerates EXACTLY what the Python bindings can express --
/// nine MD components (mirroring `bind_md.cpp`'s bindable constructors) plus
/// the AM1 single-point evaluator (mirroring `cpu::semiempirical::AM1Config`)
/// -- as plain parameter aggregates. Components outside this enumeration are
/// unrepresentable by construction (the variant IS the allow-list).
///
/// Internal wire protocol only: never written to disk, never user-visible.
/// Frame data never enters a TaskSpec (the Frame travels separately via
/// `cpu::comm::bcast_frame`).
///
/// Two validating constructors carry the task-kind invariants:
///  - MD: `TaskSpec(entries, n_steps)` -- `n_steps >= 1`, no evaluator-only
///    entry (AM1) in the list.
///  - SinglePoint: `TaskSpec(AM1EvaluatorParams)` -- exactly one entry, no
///    explicit binding (the materializer pins `evaluate::Phase::evaluate`).
///
/// @see cpu/task/task_codec.h        pack/unpack to a byte buffer
/// @see cpu/task/task_materializer.h TaskSpec -> runnable driver
#pragma once

#include <cstddef>
#include <cstdint>
#include <optional>
#include <variant>
#include <vector>

#include "atomiverse/facade/stage.h"

namespace ATV {
namespace cpu {
namespace task {

/// @brief Which driver pipeline a TaskSpec materializes into.
enum class TaskKind : std::uint8_t {
    MD = 0,          ///< `cpu::md::MD` run over a step budget.
    SinglePoint = 1, ///< one `cpu::Evaluate` pass (v1: AM1 only).
};

// ── Component parameter mirrors (one per bindable constructor) ────────────

/// @brief Mirrors `cpu::PairLJCut(epsilon, sigma, cutoff)`.
struct PairLJCutParams {
    double epsilon; ///< LJ well depth [eV].
    double sigma;   ///< LJ diameter [Angstrom].
    double cutoff;  ///< Pair cutoff [Angstrom].
};
/// @brief Field-wise equality.
/// @return True iff every field compares equal.
inline bool operator==(const PairLJCutParams& a, const PairLJCutParams& b) {
    return a.epsilon == b.epsilon && a.sigma == b.sigma && a.cutoff == b.cutoff;
}

/// @brief Mirrors `cpu::md::NeighborList<>(cutoff, skin = -1)`.
struct NeighborListParams {
    double cutoff; ///< Neighbor cutoff [Angstrom].
    double skin;   ///< Skin distance [Angstrom]; -1 = component default.
};
/// @brief Field-wise equality.
/// @return True iff every field compares equal.
inline bool operator==(const NeighborListParams& a,
                       const NeighborListParams& b) {
    return a.cutoff == b.cutoff && a.skin == b.skin;
}

/// @brief Mirrors `cpu::BondHarmonic(K, r0)`.
struct BondHarmonicParams {
    std::vector<double> K;  ///< Force constants per bond type [eV/Angstrom^2].
    std::vector<double> r0; ///< Equilibrium distances per bond type [Angstrom].
};
/// @brief Field-wise equality.
/// @return True iff every field compares equal.
inline bool operator==(const BondHarmonicParams& a,
                       const BondHarmonicParams& b) {
    return a.K == b.K && a.r0 == b.r0;
}

/// @brief Mirrors `cpu::AngleHarmonic(K, theta0)`.
struct AngleHarmonicParams {
    std::vector<double> K;      ///< Force constants per angle type [eV/rad^2].
    std::vector<double> theta0; ///< Equilibrium angles per angle type [rad].
};
/// @brief Field-wise equality.
/// @return True iff every field compares equal.
inline bool operator==(const AngleHarmonicParams& a,
                       const AngleHarmonicParams& b) {
    return a.K == b.K && a.theta0 == b.theta0;
}

/// @brief Mirrors `cpu::DihedralMultiHarmonic(Cn)`.
struct DihedralMultiHarmonicParams {
    std::vector<std::vector<double>> Cn; ///< Coefficients per dihedral type [eV].
};
/// @brief Field-wise equality.
/// @return True iff every field compares equal.
inline bool operator==(const DihedralMultiHarmonicParams& a,
                       const DihedralMultiHarmonicParams& b) {
    return a.Cn == b.Cn;
}

/// @brief Mirrors `cpu::md::VelocityVerlet(dt)`.
struct VelocityVerletParams {
    double dt; ///< Timestep [fs].
};
/// @brief Field-wise equality.
/// @return True iff every field compares equal.
inline bool operator==(const VelocityVerletParams& a,
                       const VelocityVerletParams& b) {
    return a.dt == b.dt;
}

/// @brief Mirrors `cpu::md::InitVelocity(temperature, seed)`.
struct InitVelocityParams {
    double temperature; ///< Target temperature [K].
    unsigned seed;      ///< RNG seed (dimensionless).
};
/// @brief Field-wise equality.
/// @return True iff every field compares equal.
inline bool operator==(const InitVelocityParams& a,
                       const InitVelocityParams& b) {
    return a.temperature == b.temperature && a.seed == b.seed;
}

/// @brief Mirrors `cpu::md::NoseHoover(T_target, tau_T, dt)`.
struct NoseHooverParams {
    double T_start; ///< Target temperature [K].
    double tau_T;   ///< Coupling time (LAMMPS `Tdamp`) [fs].
    double dt;      ///< Timestep [fs].
};
/// @brief Field-wise equality.
/// @return True iff every field compares equal.
inline bool operator==(const NoseHooverParams& a, const NoseHooverParams& b) {
    return a.T_start == b.T_start && a.tau_T == b.tau_T && a.dt == b.dt;
}

/// @brief Mirrors `cpu::Thermo(freq)`.
struct ThermoParams {
    int freq; ///< Output every N steps (dimensionless).
};
/// @brief Field-wise equality.
/// @return True iff every field compares equal.
inline bool operator==(const ThermoParams& a, const ThermoParams& b) {
    return a.freq == b.freq;
}

/// @brief Mirrors `cpu::semiempirical::AM1Config` (= `NDDOConfig`) minus
///        `method`, which is pinned to AM1 by this variant alternative
///        itself (and by `AM1Solver::normalize_`).
struct AM1EvaluatorParams {
    int charge = 0;             ///< Total molecular charge [e].
    int multiplicity = 1;       ///< Spin multiplicity (dimensionless).
    int max_iter = 200;         ///< SCF iteration cap (dimensionless).
    double density_tol = 1e-8;  ///< Density-norm convergence (dimensionless).
    double energy_tol = 1e-10;  ///< Energy convergence [Hartree].
};
/// @brief Field-wise equality.
/// @return True iff every field compares equal.
inline bool operator==(const AM1EvaluatorParams& a,
                       const AM1EvaluatorParams& b) {
    return a.charge == b.charge && a.multiplicity == b.multiplicity &&
           a.max_iter == b.max_iter && a.density_tol == b.density_tol &&
           a.energy_tol == b.energy_tol;
}

/// @brief The closed component allow-list. AM1 is deliberately the LAST
///        alternative (index 9) so the nine MD component_ids stay stable
///        on the wire.
using ComponentParams =
    std::variant<PairLJCutParams, NeighborListParams, BondHarmonicParams,
                 AngleHarmonicParams, DihedralMultiHarmonicParams,
                 VelocityVerletParams, InitVelocityParams, NoseHooverParams,
                 ThermoParams, AM1EvaluatorParams>;

/// @brief Human-readable component name for error messages, by variant index.
/// @param index `ComponentParams::index()` value.
/// @return Static name string; "unknown" past the enumeration.
inline const char* component_name(std::size_t index) {
    static const char* const kNames[] = {
        "PairLJCut",      "NeighborList",   "BondHarmonic",
        "AngleHarmonic",  "DihedralMultiHarmonic", "VelocityVerlet",
        "InitVelocity",   "NoseHoover",     "Thermo",
        "AM1Evaluator"};
    static_assert(sizeof(kNames) / sizeof(kNames[0]) ==
                      std::variant_size_v<ComponentParams>,
                  "kNames must track ComponentParams");
    return index < std::variant_size_v<ComponentParams> ? kNames[index]
                                                        : "unknown";
}

/// @brief Explicit `md.on(phase).add(comp, order)` registration record.
///        The phase stays the TYPED enum in memory (stage-typing rule);
///        it crosses the wire as `ATV::slot_of(phase)`.
struct ExplicitBinding {
    ATV::md::Phase phase; ///< Target MD phase.
    int order;            ///< In-phase ordering key (dimensionless).
};
/// @brief Field-wise equality.
/// @return True iff every field compares equal.
inline bool operator==(const ExplicitBinding& a, const ExplicitBinding& b) {
    return a.phase == b.phase && a.order == b.order;
}

/// @brief One registered component: its parameters plus HOW it was
///        registered (`nullopt` = `md.add(comp)` default-binding path).
struct ComponentEntry {
    ComponentParams params;                ///< The component's ctor args.
    std::optional<ExplicitBinding> binding; ///< Explicit phase/order, if any.
};
/// @brief Field-wise equality.
/// @return True iff every field compares equal.
inline bool operator==(const ComponentEntry& a, const ComponentEntry& b) {
    return a.params == b.params && a.binding == b.binding;
}

/// @brief The complete, validated description of one dispatchable task.
class TaskSpec {
  public:
    /// @brief Wire-format version; bumped when the bindable set grows.
    static constexpr std::uint32_t kWireVersion = 1;

    /// @brief MD task: the component list plus a step budget.
    /// @param entries Registered components (default- or explicit-binding).
    /// @param n_steps Step budget (steps); must be >= 1 and fit an int.
    /// @throws std::invalid_argument if `n_steps < 1`, if `n_steps`
    ///         exceeds the int step-budget limit (checked narrowing), or
    ///         if any entry carries an evaluator-only component
    ///         (AM1Evaluator).
    TaskSpec(std::vector<ComponentEntry> entries, std::int64_t n_steps);

    /// @brief AM1 single-point task: one evaluator parameter set, no
    ///        binding (the materializer pins `evaluate::Phase::evaluate`),
    ///        no step budget.
    /// @param params AM1 evaluator configuration.
    explicit TaskSpec(AM1EvaluatorParams params);

    /// @brief The task kind this spec materializes into.
    /// @return TaskKind::MD or TaskKind::SinglePoint.
    TaskKind kind() const noexcept { return kind_; }
    /// @brief The registered component entries (SinglePoint: exactly one).
    /// @return Entry list; reference valid for the spec's lifetime.
    const std::vector<ComponentEntry>& entries() const noexcept {
        return entries_;
    }
    /// @brief Step budget [steps]; fixed to 1 for SinglePoint.
    /// @return The validated step budget.
    std::int64_t n_steps() const noexcept { return n_steps_; }

  private:
    TaskKind kind_;
    std::vector<ComponentEntry> entries_;
    std::int64_t n_steps_;
};

} // namespace task
} // namespace cpu
} // namespace ATV
