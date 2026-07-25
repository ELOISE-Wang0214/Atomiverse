/// @file cpu/dynamics/const_pot_bo.h
/// @brief Constant-potential Born–Oppenheimer electrode charge solve (E1/E2).
///
/// Gaussian elastance matrix \(A\), dense inverse **or** CG, optional
/// electroneutrality (LAMMPS ELECTRODE `symm on` / MetalWalls BO), and
/// optional Kolafa ASPC charge predictor for CG warm-start (MetalWalls E2).
/// Component on `cpu::md::MD` at `pre_force`; writes electrode charges to
/// `kAtom_charge`. Not a Driver category and not part of `KSpace*`.
#pragma once

#include <cstddef>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "atomiverse/core/binding.h"
#include "atomiverse/core/stage.h"
#include "atomiverse/cpu/component/component.h"
#include "atomiverse/cpu/electrodes/charge_aspc.h"

namespace ATV {
namespace cpu {
namespace md {

/// @brief One electrode group: all atoms of a given type at fixed \(\Psi\) [V].
struct ElectrodeGroup {
    int type = 0;     ///< Atom type (matches `kAtom_type`).
    double psi = 0.0; ///< Applied electrostatic potential [V].
};

/// @brief Constant-potential BO charge solver (dense inverse or CG + ASPC).
///
/// Parameters are constructor arguments and public members — no Config bag.
/// Optional Thomas–Fermi metallicity (map E7 / Scalfi et al. 2020): set
/// `l_tf` [Å] and `site_density` [1/Å³] to harden the elastance diagonal.
///
/// Parallel-plate (E1):
/// @code
/// cpu::md::ConstPotBO conp(/*eta [1/Å]=*/2.0);
/// conp.add_electrode(/*type=*/1, /*Psi [V]=*/0.0);
/// conp.add_electrode(/*type=*/2, /*Psi [V]=*/1.0);
/// conp.electroneutrality = true;
/// conp.use_cg = true;            // iterative path
/// conp.predictor_order = 2;      // ASPC warm-start (0 = cold)
/// conp.xy_image_range = 200;
/// conp.l_tf = 2.0;            // [Å]; 0 = perfect metal
/// conp.site_density = 0.0593; // [1/Å³] Au FCC
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

    /// @brief Integral capacitance \(|Q_{\mathrm{last}}| / |\Delta\Psi|\) [e/V]
    ///        for a two-electrode capacitor (uses first/last electrode groups).
    double capacitance() const;

    /// @brief \(\lvert\sum_i Q_i\rvert\) [e] over all electrode sites.
    double neutrality_residual() const;

    /// @brief \(\lVert AQ-(\Psi-b+\mathbf{E}\cdot\mathbf{r})\rVert_\infty\) after last
    ///        solve (includes Lagrange shift when electroneutrality is on).
    double solve_residual() const { return residual_; }

    /// @brief Number of electrode atoms in the last solve.
    int n_electrode() const { return static_cast<int>(elec_atoms_.size()); }

    /// @brief CG iterations used on the last launch (0 if dense Cholesky).
    /// @return Iteration count of the primary \(A Q = \mathrm{rhs}\) CG.
    int last_iters() const { return last_iters_; }

    /// @brief Whether the last CG launch used an ASPC warm-start guess.
    /// @return True if `predict()` succeeded before the last solve.
    bool last_used_predictor() const { return last_used_predictor_; }

    /// @brief Target \(\sum Q\) used on the last neutrality projection [e]
    ///        (`0` for electrode-only; \(-\sum q_{\mathrm{elyt}}\) for global).
    double neutrality_target() const { return neutrality_target_; }

    /// @brief Frobenius condition number proxy \(\kappa_2(A)\approx
    ///        \lambda_{\max}/\lambda_{\min}\) from power iteration (empty-cell
    ///        elastance debug / MetalWalls matrix-condition tool).
    double matrix_condition() const;

    /// @brief Copy of last electrode charges [e] (length `n_electrode()`).
    const std::vector<double>& charges() const { return Q_; }

    // ── Public parameters (mutable after construction) ───────────────
    double eta = 1.979; ///< Gaussian width [1/Å].
    bool electroneutrality =
        true; ///< Project electrode net charge (see global_*).
    /// When true with `electroneutrality`, enforce global neutrality
    /// \(\sum_i Q_i = -\sum_j q_j^{\mathrm{electrolyte}}\) (MetalWalls /
    /// Scalfi). When false (default), enforce electrode-only \(\sum Q_i=0\).
    bool global_electroneutrality = false;
    /// When true, after the charge solve accumulate Hellmann–Feynman-style
    /// forces on electrode atoms from electrolyte point charges (fixed-\(Q\)
    /// slice of \(-\nabla U\)). Required for movable walls / piston coupling.
    bool write_electrode_forces = false;
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

    /// @brief Use conjugate-gradient on \(A Q = \Psi-b\) instead of dense
    ///        Cholesky (MetalWalls iterative path). Required for ASPC gain.
    bool use_cg = false;
    /// @brief CG residual tolerance on \(\lVert A Q - \mathrm{rhs}\rVert_2
    ///        / \lVert\mathrm{rhs}\rVert_2\).
    double cg_tol = 1e-10;
    /// @brief Hard cap on CG iterations; throws if exceeded.
    int max_cg_iters = 500;
    /// @brief Kolafa ASPC order \(K\) for charge warm-start. `0` disables
    ///        prediction (cold start every step). Typical product value: 2.
    int predictor_order = 0;
    /// Thomas–Fermi screening length \(l_{\mathrm{TF}}\) [Å]. Zero recovers
    /// perfect-metal BO (MetalWalls / Scalfi). Requires `site_density` > 0
    /// when positive.
    double l_tf = 0.0;
    /// Electrode site number density \(d\) [1/Å³] for the TF hardness
    /// \(A_{ii}\mathrel{+}= 4\pi k_e l_{\mathrm{TF}}^{2} d\). Au FCC with
    /// \(a_0=4.07\) Å gives \(d = 4/a_0^3 \approx 0.0593\).
    double site_density = 0.0;

  private:
    std::vector<ElectrodeGroup> groups_;

    // Maps: electrode local index → frame atom row; group id per electrode atom.
    std::vector<int> elec_atoms_;
    std::vector<int> elec_group_;
    std::vector<double> psi_vec_; ///< \(\Psi_i\) per electrode atom [V].
    std::vector<double> Q_;       ///< Solved charges [e].
    std::vector<double> b_;       ///< Electrolyte potential [V].
    std::vector<double> A_;       ///< Row-major elastance [V/e].
    std::vector<double> A_work_;  ///< Factorization workspace.
    std::vector<double> total_Q_; ///< Per-group totals [e].
    bool matrix_ready_ = false;
    double residual_ = 0.0;
    int last_iters_ = 0;
    bool last_used_predictor_ = false;

    std::unique_ptr<electrodes::ChargeASPC> aspc_;
    /// Unconstrained \(A Q = \Psi-b\) solution recorded for ASPC (before
    /// electroneutrality projection). Empty when ASPC is off.
    std::vector<double> aspc_record_;

    void collect_electrodes_(const cpu::Frame& frame);
    void build_matrix_(const cpu::Frame& frame);
    void build_rhs_(const cpu::Frame& frame);
    void solve_();
    void solve_dense_();
    void solve_cg_();
    /// @brief CG solve \(A x = b\) (row-major SPD \(A\)), starting from \(x\).
    /// @return Iteration count.
    int cg_solve_(const std::vector<double>& A, int n,
                  const std::vector<double>& b, std::vector<double>& x) const;
    void write_charges_(cpu::Frame& frame) const;
    void write_electrode_forces_(cpu::Frame& frame) const;
    double electrolyte_charge_(const cpu::Frame& frame) const;
    double pair_A_(double r, double eta_i, double eta_j) const;
    double pair_b_(double r, double q, double eta_i) const;
    void ensure_aspc_();
    void record_aspc_();

    double neutrality_target_ = 0.0;
};

} // namespace md
} // namespace cpu
} // namespace ATV
