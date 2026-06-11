/// @file cpu/potentials/coul_long.h
/// @brief Real-space screened Coulomb (LAMMPS `pair_style coul/long`) — CPU.
///
/// The short-range half of an Ewald/PME split. Evaluates
///
///   E_real = kCoulomb · Σ_{i<j} q_i q_j · erfc(α r_ij) / r_ij
///
/// over the neighbor list, where α is the Ewald splitting parameter shared
/// with the reciprocal-space `KSpaceEwald` / `KSpacePME` component. The
/// complementary `erf` part of the full 1/r interaction is supplied by the
/// reciprocal sum; the two must use the **same** α and the neighbor-list
/// cutoff must be ≥ this component's cutoff.
///
/// This is a cross-driver force component (added via `.on(force).add()`),
/// not a `set_pair` style, so it coexists with a Lennard-Jones (or other)
/// `Pair` evaluator in the same force phase. It reads the shared
/// `runtime.pair_list` built by `NeighborList`.
#pragma once

#include <vector>

#include "atomiverse/core/real.h"
#include "atomiverse/cpu/component/component.h"
#include "atomiverse/facade/component_role.h"

namespace ATV {
namespace cpu {

/// @brief Real-space screened Coulomb force/energy (Ewald/PME real part).
class CoulLong : public cpu::ComponentBase {
  public:
    /// @param alpha  Ewald splitting parameter α [1/Å]. Must match the
    ///               reciprocal-space component.
    /// @param cutoff Real-space cutoff [Å]. The neighbor list must be built
    ///               with a cutoff ≥ this value.
    CoulLong(double alpha, double cutoff);

    std::string name() const override { return "coul/long"; }
    ComponentRole role() const override { return ComponentRole::Evaluator; }
    Capability capabilities() const override {
        return Capability::Energy | Capability::Force | Capability::Stress;
    }
    bool needs_neighbor_list() const override { return true; }

    /// @brief Real-space cutoff [Å] (mirrors `PairEvaluator::pair_cutoff`).
    Real pair_cutoff() const { return cutoff_; }

    void build(const BuildContext& ctx) override;
    void launch(const LaunchContext& ctx) override;
    void destroy(const DestroyContext& ctx) override;

  private:
    double alpha_;
    double cutoff_, cutoff_sq_;

    // Per-thread force reduction buffers (Newton's 3rd law), mirroring
    // PairLJCut. Index range covers owned + ghost atoms under MPI. Ghost
    // charges arrive via the build_halo border message (x,y,z,tag,q), so
    // f.atoms.charge[j] is valid for ghost j — no tag→charge lookup table.
    int n_threads_ = 1;
    std::vector<std::vector<double>> fx_thr_, fy_thr_, fz_thr_;
};

} // namespace cpu
} // namespace ATV
