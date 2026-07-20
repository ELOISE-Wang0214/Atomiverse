/// @file cpu/scf/stages/fock_stage.cpp
/// @brief FockStage — wraps internal::HamiltonianBuilder at scf::Phase::fock.

#include "atomiverse/cpu/scf/stages/fock_stage.h"

#include "atomiverse/facade/stage.h"

namespace ATV {
namespace cpu {
namespace scf {

FockStage::FockStage(bool restricted, bool is_hf, XCType xc_type,
                     GridLevel grid_level, int n_radial, int n_angular,
                     double hybrid_alpha_override,
                     const MolecularGrid* fixed_grid,
                     const IntegralPack& integrals, const RealMatrix& P_alpha,
                     const RealMatrix& P_beta, RealMatrix& F_alpha,
                     RealMatrix& F_beta, RealMatrix& J, double& e_xc_grid,
                     double& e_exact_exchange, double& int_n_alpha,
                     double& int_n_beta)
    : restricted_(restricted), is_hf_(is_hf), xc_type_(xc_type),
      grid_level_(grid_level), n_radial_(n_radial), n_angular_(n_angular),
      hybrid_alpha_override_(hybrid_alpha_override), fixed_grid_(fixed_grid),
      integrals_(integrals), P_alpha_(P_alpha), P_beta_(P_beta),
      F_alpha_(F_alpha), F_beta_(F_beta), J_(J), e_xc_grid_(e_xc_grid),
      e_exact_exchange_(e_exact_exchange), int_n_alpha_(int_n_alpha),
      int_n_beta_(int_n_beta) {}

std::vector<DefaultBinding> FockStage::default_bindings() const {
    return {ATV::make_binding(ATV::scf::Phase::fock, 0)};
}

void FockStage::build(const cpu::BuildContext& ctx) {
    (void)ctx;
    // One-time allocation of shell-pair lists + (KS) the integration grid,
    // via HamiltonianBuilder's config-free entry.
    hamiltonian_.prepare(integrals_, is_hf_, grid_level_, n_radial_, n_angular_,
                         fixed_grid_);
}

void FockStage::launch(const cpu::LaunchContext& ctx) {
    if (!ctx.at(ATV::scf::Phase::fock))
        return;

    // Same core, same order, no extra arithmetic: the discrete step() forwards
    // to the shared HamiltonianBuilder core, so the outputs are bit-identical to
    // the shipping solver path.
    const internal::FockResult result =
        hamiltonian_.step(integrals_, restricted_, is_hf_, xc_type_,
                          hybrid_alpha_override_, P_alpha_, P_beta_);

    F_alpha_ = result.F_alpha;
    F_beta_ = result.F_beta;
    J_ = result.J;
    e_xc_grid_ = result.E_xc_grid;
    e_exact_exchange_ = result.E_exact_exchange;
    int_n_alpha_ = result.int_n_alpha;
    int_n_beta_ = result.int_n_beta;
}

} // namespace scf
} // namespace cpu
} // namespace ATV
