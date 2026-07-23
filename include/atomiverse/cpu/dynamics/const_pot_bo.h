/// @file cpu/dynamics/const_pot_bo.h
/// @brief Constant-potential Born–Oppenheimer electrode charge solve (E1).
///
/// Gaussian elastance matrix \(A\), dense inverse, optional electroneutrality
/// (LAMMPS ELECTRODE `symm on` / MetalWalls BO). Component on `cpu::md::MD`
/// at `pre_force`; writes electrode charges to `kAtom_charge`. Not a Driver
/// category and not part of `KSpace*`.
#pragma once

#include <cstddef>
#include <string>
#include <utility>
#include <vector>

#include "atomiverse/core/binding.h"
#include "atomiverse/core/stage.h"
#include "atomiverse/cpu/component/component.h"

namespace ATV {
namespace cpu {
namespace md {

/// @brief One electrode group: all atoms of a given type at fixed \(\Psi\) [V].
struct ElectrodeGroup {
    int type = 0;     ///< Atom type (matches `kAtom_type`).
    double psi = 0.0; ///< Applied electrostatic potential [V].
};

/// @brief Constant-potential BO charge solver (dense matrix inverse).
///
/// Parameters are constructor arguments and public members — no Config bag.
///
/// Parallel-plate (E1):
/// @code
/// cpu::md::ConstPotBO conp(/*eta [1/Å]=*/2.0);
/// conp.add_electrode(/*type=*/1, /*Psi [V]=*/0.0);
/// conp.add_electrode(/*type=*/2, /*Psi [V]=*/1.0);
/// conp.electroneutrality = true;
/// conp.xy_image_range = 200;
/// sim.on(md::Phase::pre_force).add(std::move(conp));
/// @endcode
///
/// Finite-field + single electrode (E3, Dufils-style): one wall, free total
/// charge, external field on the BO RHS \((AQ)_i=\Psi_i-b_i+\mathbf{E}\cdot\mathbf{r}_i\).
/// Compose with `FixEfield` for forces:
/// @code
/// cpu::md::ConstPotBO conp(2.0);
/// conp.add_electrode(1, 0.0);
/// conp.Ez = 0.1;                 // [V/Å]
/// conp.electroneutrality = false;
/// sim.add(std::move(conp));
/// sim.add(cpu::md::FixEfield{0.0, 0.0, 0.1});
/// @endcode
class ConstPotBO : public cpu::Component {
  public:
    /// @param eta Gaussian charge width \(\eta\) [1/Å] (shared by all
    ///            electrode sites in v1).
    explicit ConstPotBO(double eta);

    /// @brief Register one electrode group (type + applied potential).
    /// @param type Atom type of electrode sites.
    /// @param psi Applied potential [V].
    ConstPotBO& add_electrode(int type, double psi) &;
    ConstPotBO add_electrode(int type, double psi) &&;

    std::string name() const override { return "const_pot_bo"; }

    std::vector<DefaultBinding> default_bindings() const override {
        // Before pair/kspace force so they see updated charges; after Zero*.
        return {ATV::make_binding(ATV::md::Phase::pre_force, -50)};
    }

    void build(const BuildContext& ctx) override;
    void launch(const LaunchContext& ctx) override;
    void destroy(const DestroyContext& ctx) override;

    /// @brief Total charge on electrode group `ie` [e] after last launch.
    /// @param ie Electrode index in registration order.
    double total_charge(std::size_t ie) const;

    /// @brief \(\Psi_{\mathrm{last}}-\Psi_{\mathrm{first}}\) [V] (two-electrode convenience).
    double delta_psi() const;

    /// @brief \(\lvert\sum_i Q_i\rvert\) [e] over all electrode sites.
    double neutrality_residual() const;

    /// @brief \(\lVert AQ-(\Psi-b+\mathbf{E}\cdot\mathbf{r})\rVert_\infty\) after last
    ///        solve (includes Lagrange shift when electroneutrality is on).
    double solve_residual() const { return residual_; }

    /// @brief Number of electrode atoms in the last solve.
    int n_electrode() const { return static_cast<int>(elec_atoms_.size()); }

    // ── Public parameters (mutable after construction) ───────────────
    double eta = 1.979;           ///< Gaussian width [1/Å].
    bool electroneutrality = true; ///< Project \(\sum Q_i = 0\).
    /// Half-range of 2D \(xy\) lattice images for planar PBC (0 = primary
    /// cell only). Matches LAMMPS ELECTRODE madelung / plate_cap sum when
    /// large (e.g. 200+). Cost \(\propto N_{\mathrm{elec}}^2 \times R^2\).
    int xy_image_range = 0;
    /// Rebuild \(A\) every launch (movable electrodes). v1 walls are frozen;
    /// leave false for geometry-static electrodes.
    bool rebuild_matrix_every_step = false;
    /// Finite-field components [V/Å] on the BO RHS (Dufils single-electrode /
    /// 3D PBC). Zero leaves the E1 parallel-plate path unchanged. Compose
    /// with `FixEfield` for electrolyte forces — this member only updates \(Q\).
    double Ex = 0.0;
    double Ey = 0.0;
    double Ez = 0.0;

  private:
    std::vector<ElectrodeGroup> groups_;

    // Maps: electrode local index → frame atom row; group id per electrode atom.
    std::vector<int> elec_atoms_;
    std::vector<int> elec_group_;
    std::vector<double> psi_vec_;   ///< \(\Psi_i\) per electrode atom [V].
    std::vector<double> Q_;         ///< Solved charges [e].
    std::vector<double> b_;         ///< Electrolyte potential [V].
    std::vector<double> A_;         ///< Row-major elastance [V/e].
    std::vector<double> A_work_;    ///< Factorization workspace.
    std::vector<double> total_Q_;   ///< Per-group totals [e].
    bool matrix_ready_ = false;
    double residual_ = 0.0;

    void collect_electrodes_(const cpu::Frame& frame);
    void build_matrix_(const cpu::Frame& frame);
    void build_rhs_(const cpu::Frame& frame);
    void solve_();
    void write_charges_(cpu::Frame& frame) const;
    double pair_A_(double r, double eta_i, double eta_j) const;
    double pair_b_(double r, double q, double eta_i) const;
};

} // namespace md
} // namespace cpu
} // namespace ATV
