/// @file cpu/scf/uhf.cpp
/// @brief UHF constructor defaults, occupation resolution, and analytic
///        nuclear Hessian (unrestricted CPHF). The shared driver loop lives
///        in mean_field_method.cpp (explicitly instantiated there).

#include "atomiverse/cpu/scf/uhf.h"

#include <stdexcept>
#include <string>

#include <xtensor/containers/xtensor.hpp>

#include "atomiverse/core/units.h"
#include "atomiverse/cpu/scf/dense_matrix.h"
#include "atomiverse/cpu/scf/internal/integral_derivatives.h"
#include "atomiverse/cpu/scf/internal/integrals.h"
#include "atomiverse/cpu/scf/internal/ucphf_pack.h"
#include "atomiverse/cpu/scf/slater_determinant.h"

namespace ATV {
namespace cpu {
namespace scf {

UHF::UHF() {
    spin_ = SpinType::Unrestricted;
    method_ = MeanFieldType::HF;
    init_guess = InitGuess::Hcore;
    damping_factor = 0.7;
    diis_max_vectors = 6;
    diis_start_iter = 5;
    tol = 1e-8;
}
UHF::UHF(int na, int nb) : UHF() {
    n_alpha = na;
    n_beta = nb;
}
void UHF::resolve_occupation(int) {
    n_alpha_ = n_alpha;
    n_beta_ = n_beta;
}

vib::HessianResult UHF::hessian(Frame& frame, const SCFInitialGuess* guess) {
    if (ri != RI::None)
        throw std::runtime_error(
            "UHF::hessian: the analytic Hessian is not wired for density "
            "fitting (RI) — derivative three-center integrals are a "
            "follow-up chain link; use ri = RI::None (no silent fallback)");

    // Fresh SCF at this geometry (writes the eV energy to frame.meta).
    run(frame, guess);
    if (!converged())
        throw std::runtime_error(
            "UHF::hessian: SCF did not converge — no reference "
            "wavefunction to differentiate (" + stop_reason() + ")");

    const SlaterDeterminant& det = determinant();
    const IntegralPack& pack = integral_workspace_();
    if (pack.spherical())
        throw std::runtime_error(
            "UHF::hessian: analytic second derivatives are not implemented "
            "for spherical-harmonic basis sets (def2 family) yet — the "
            "derivative c2s transform is a follow-up chain link. Use a "
            "Cartesian-convention basis or the finite-difference "
            "vib::Hessian driver.");
    if (pack.has_ecp)
        throw std::runtime_error(
            "UHF::hessian: analytic second derivatives are not implemented "
            "for ECP packs — d²U/dR² of the semi-local potential is a "
            "follow-up chain link. Use the finite-difference vib::Hessian "
            "driver.");

    const xt::xtensor<double, 4> eri =
        (pack.ERI.size() > 0) ? pack.ERI : pack.eri_working();

    using Spin = SlaterDeterminant::Spin;
    const RealMatrix Pa = det.density(Spin::alpha);
    const RealMatrix Pb = det.density(Spin::beta);

    // Skeleton first derivatives: dS plus per-spin derivative Focks
    // F^{σ(x)} = dT + dV + J^{(x)}[P] − K^{(x)}[P^σ].
    const IntegralDerivPack dpack{pack};
    xt::xtensor<double, 4> dFa = dpack.dT + dpack.dV;
    xt::xtensor<double, 4> dFb = dpack.dT + dpack.dV;
    accumulate_two_electron_fock_derivative_uhf(dFa, dFb, pack, Pa, Pb);

    // Coupled α/β CPHF solve (3N right-hand sides) + Hessian assembly.
    internal::UCphfPack cphf{det, pack, eri, dpack.dS, dFa, dFb};
    cphf.solve();
    const xt::xtensor<double, 2> H_ha = cphf.hessian();

    const int natm = static_cast<int>(frame.nlocal());
    const double to_ev_ang2 =
        units::kHartreeInEv / (units::kBohrInAngstrom * units::kBohrInAngstrom);
    vib::HessianResult out;
    out.n_atoms = natm;
    out.hessian = H_ha * to_ev_ang2;
    return out;
}

} // namespace scf
} // namespace cpu
} // namespace ATV
