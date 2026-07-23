/// @file cpu/potentials/pair_born.h
/// @brief Born–Mayer–Huggins pair: U = A exp((σ−r)/ρ) − C/r⁶ + D/r⁸ (CPU).
///
/// Matches LAMMPS `pair_style born` with energy shift at cutoff
/// (`pair_modify shift yes`): \(U(r_c)=0\). Units: energy [eV], length [Å].
#pragma once

#include <string>

#include "atomiverse/cpu/potentials/pair.h"

namespace ATV {
namespace cpu {

/// @brief Born–Mayer–Huggins pair potential (single-type, shifted at cutoff).
///
/// \[
/// U(r)=A\,e^{(\sigma-r)/\rho}-\frac{C}{r^{6}}+\frac{D}{r^{8}}-U_{\mathrm{shift}},
/// \quad r < r_c
/// \]
/// with
/// \(U_{\mathrm{shift}}=A e^{(\sigma-r_c)/\rho}-C/r_c^{6}+D/r_c^{8}\).
///
/// Force: \(F=-\nabla U\) with
/// \(\mathrm{d}U/\mathrm{d}r=-(A/\rho)\,e^{(\sigma-r)/\rho}+6C/r^{7}-8D/r^{9}\)
/// (unshifted; shift is constant).
///
/// O(N²) all-pairs walk (NeighborConsumer withdrawn) — same foothold shape as
/// `PairBuck` / `PairYukawa`. Production multi-type + NL is deferred to
/// Fumi–Tosi / salt productization.
class PairBorn : public Pair {
  public:
    /// @param A      Born–Mayer prefactor [eV] (must be ≥ 0).
    /// @param rho    Softness length [Å] (must be > 0).
    /// @param sigma  Ionic-pair radius parameter [Å].
    /// @param C      Dispersion coefficient [eV·Å⁶] (must be ≥ 0).
    /// @param D      \(r^{-8}\) coefficient [eV·Å⁸].
    /// @param cutoff Pair cutoff [Å] (must be > 0).
    PairBorn(double A, double rho, double sigma, double C, double D,
             double cutoff);

    std::string name() const override { return "born"; }

    /// @brief O(N²) path: withdraw neighbor-list gate.
    const NeighborConsumer* as_neighbor_consumer() const final {
        return nullptr;
    }

  protected:
    void compute(const LaunchContext& ctx) override;

  private:
    double A_ = 0;
    double rho_ = 0;
    double sigma_ = 0;
    double C_ = 0;
    double D_ = 0;
    double cutoff_sq_ = 0;
    double v_shift_ = 0; ///< U(r_c) so shifted potential is zero at cutoff
};

} // namespace cpu
} // namespace ATV
