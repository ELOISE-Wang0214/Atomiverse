/// @file cpu/path/pimd.h
/// @brief PIMD skeleton — ring-polymer MD over IdentityReplication.
///
/// Budget-terminating preset (like Remd): expands one Frame into P identical
/// beads, runs P owned `cpu::md::MD` trees with ring-polymer springs, collapses
/// to the **centroid**. Second real consumer of IdentityReplication after REMD.
///
/// Free-particle path: program builder may install nothing (springs only).
/// Optional external isotropic harmonic trap k_ext (eV/Å²) about the origin.
///
/// Out of scope (v1 skeleton): RPMD estimators, PILE thermostat, ring-polymer
/// contraction, multi-rank bead parallel, Trotter-factor tables.
///
/// Spring constant (Tuckerman / Ceriotti form):
///   k = (m / β² ℏ²) · P    [converted to eV/Å² via mass_v2_to_energy]
/// with ℏ = h/(2π), β = 1/(k_B T).
#pragma once

#include <cstdint>
#include <memory>
#include <vector>

#include "atomiverse/cpu/dynamics/md.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/path/ring_polymer_spring.h"
#include "atomiverse/cpu/replica/replica_expansion.h"
#include "atomiverse/cpu/replica/replica_set.h"

namespace ATV {
namespace cpu {
namespace path {

/// @brief Configuration for one PIMD ring polymer.
struct PimdConfig {
    int n_beads = 4;             ///< Trotter number P (≥ 2).
    double temperature = 300.0;  ///< Physical temperature [K] (> 0).
    double dt = 0.5;             ///< Nested MD timestep [fs] (> 0).
    int steps_per_block = 10;    ///< MD steps between outer bookkeeping.
    std::uint64_t seed = 1;      ///< InitVelocity seeds.
    /// @brief Optional external harmonic trap k |r|²/2 [eV/Å²]; 0 = free particle.
    double k_external = 0.0;
};

/// @brief Ring-polymer MD preset over IdentityReplication (bead indices).
///
/// Not a new Driver category — owns N `MD` programs + `ReplicaSet`, like Remd.
class Pimd {
  public:
    explicit Pimd(PimdConfig cfg);

    /// @brief Install physical force program (called once per bead at setup).
    /// @param md Nested MD for this bead (integrator already set).
    /// @param bead_index Bead index in [0, P).
    /// @param user Opaque user pointer.
    using ProgramBuilder = void (*)(cpu::md::MD& md, int bead_index, void* user);
    void set_program_builder(ProgramBuilder builder, void* user = nullptr);

    /// @brief Expand, build P MD programs, run `n_blocks` blocks, collapse centroid.
    /// @param frame Input geometry (copied to every bead); overwritten with centroid.
    /// @param n_blocks Number of outer blocks (≥ 1).
    void run(cpu::Frame& frame, int n_blocks);

    /// @brief Number of beads P.
    int n_beads() const { return cfg_.n_beads; }

    /// @brief Spring constant k used for the ring [eV/Å²].
    double spring_constant() const { return k_spring_; }

    /// @brief Mean squared bead–centroid displacement after last run [Å²]
    ///        (diagnostic; free particle should grow then stabilize vs springs).
    double mean_bead_spread() const { return mean_spread_; }

  private:
    void ensure_built_(cpu::Frame& frame);
    void propagate_();
    void collapse_centroid_(cpu::Frame& out);
    static double compute_k_spring_(double temperature_K, int n_beads,
                                    double mass_amu);

    PimdConfig cfg_;
    ProgramBuilder builder_ = nullptr;
    void* builder_user_ = nullptr;
    std::unique_ptr<IdentityReplication> expansion_;
    std::unique_ptr<ReplicaSet> band_;
    // springs_ must outlive programs_ (add_ref borrow) — declared first so
    // it is destroyed last.
    std::vector<std::unique_ptr<RingPolymerSpring>> springs_;
    std::vector<std::unique_ptr<cpu::md::MD>> programs_;
    double k_spring_ = 0;
    double mean_spread_ = 0;
    bool built_ = false;
};

} // namespace path
} // namespace cpu
} // namespace ATV
