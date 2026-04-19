/// @file cpu/md/lj.h
/// @brief CPU Lennard-Jones 12-6 (shifted) force component.
#pragma once

#include <string>
#include <vector>

#include "atomiverse/cpu/component.h"

namespace ATV {
namespace cpu {
namespace md {

/// @brief Lennard-Jones 12-6 with shifted potential.
///
/// V(r) = 4*eps * [(sig/r)^12 - (sig/r)^6] - V_shift
/// F(r) = 24*eps/r * [2*(sig/r)^12 - (sig/r)^6] * r_hat
///
/// Optional special-bonds exclusions (LAMMPS convention): pairs that are
/// 1-2 (bonded), 1-3 (across an angle vertex) or 1-4 (across a dihedral)
/// are skipped during the pair loop. Enable via `set_special_bonds()`. Lists
/// are derived from `frame.bonds`/`frame.angles`/`frame.dihedrals` at
/// `build()` time and consulted with a per-atom binary search at
/// `launch()` time. Default: no exclusions (matches the original LJ
/// behaviour for plain non-bonded fluids).
///
/// Thread-parallel with per-thread force reduction buffers (Newton's 3rd law).
class LJ : public cpu::ComponentBase {
public:
    /// @param epsilon Well depth [eV].
    /// @param sigma Size parameter [Angstrom].
    /// @param cutoff Cutoff distance [Angstrom].
    LJ(double epsilon, double sigma, double cutoff);

    /// @brief Enable/disable LAMMPS-style "special bonds" exclusions.
    ///
    /// `s12 == 0` skips 1-2 (bonded) pairs entirely. `s12 == 1` keeps full
    /// LJ on them. Other scalings are not yet implemented and trigger an
    /// assertion -- production force fields with 1-4 scaling (e.g. OPLS-AA
    /// `s14 = 0.5`) are tracked under msu-pe-validation Phase A.5.
    /// @param s12 Scale for 1-2 (bond) pairs. 0 or 1.
    /// @param s13 Scale for 1-3 (angle) pairs. 0 or 1.
    /// @param s14 Scale for 1-4 (dihedral) pairs. 0 or 1.
    /// @return `*this` for fluent setup.
    LJ& set_special_bonds(double s12, double s13, double s14);

    std::string name() const override { return "lj"; }
    void build(const cpu::BuildContext& ctx) override;
    void launch(const cpu::LaunchContext& ctx) override;
    void destroy(const cpu::DestroyContext& ctx) override;

private:
    double epsilon_, sigma_, cutoff_;
    double sigma6_, cutoff_sq_, v_shift_;

    // Special-bonds: when any of these is < 1, the corresponding bonded
    // pair list is built at build() time and consulted at launch() time.
    double s12_ = 1.0, s13_ = 1.0, s14_ = 1.0;

    // Per-atom sorted exclusion lists (only populated when at least one of
    // s12_/s13_/s14_ is < 1 -- otherwise empty and the launch loop skips
    // the lookup entirely).
    std::vector<std::vector<int>> exclusions_;

    // Per-thread force reduction buffers
    int n_threads_ = 1;
    std::vector<std::vector<double>> fx_thr_, fy_thr_, fz_thr_;
};

} // namespace md
} // namespace cpu
} // namespace ATV
