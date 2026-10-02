/// @file cpu/scf/rhf.cpp
/// @brief RHF constructor defaults + occupation resolution + Seeger-Pople
///        stability analysis. The shared driver loop lives in
///        mean_field_method.cpp (explicitly instantiated there).

#include "atomiverse/cpu/scf/rhf.h"

#include <cmath>
#include <stdexcept>
#include <vector>

#include <xtensor-blas/xlinalg.hpp>

#include "atomiverse/cpu/scf/internal/integrals.h"
#include "atomiverse/cpu/scf/internal/one_electron.h"
#include "atomiverse/cpu/scf/internal/orbital_hessian.h"
#include "atomiverse/cpu/scf/internal/stability.h"

namespace ATV {
namespace cpu {
namespace scf {

namespace {
/// Numerical-zero guard on the stability eigenvalues [Hartree]: a lowest
/// eigenvalue above −kStabilityZeroTol counts as stable (guards round-off
/// on genuinely zero modes; physical instabilities are orders larger).
constexpr double kStabilityZeroTol = 1e-8;

/// Rotation amplitude [rad] of the leading component of the external
/// instability mode — matches the `apply_homo_lumo_mixing` default (20°)
/// so the rotated guess bridges the existing broken-symmetry machinery.
constexpr double kRotationTheta = 0.3490658503988659;
} // namespace

RHF::RHF() {
    spin_ = SpinType::Restricted;
    method_ = MeanFieldType::HF;
    init_guess = InitGuess::MINAO;
    damping_factor = 0.0;
    tol = 1e-9;
}

RHF::RHF(BasisSet basis_in) : RHF() { basis = basis_in; }

RHF::RHF(BasisSet basis_in, int charge_in) : RHF(basis_in) {
    charge = charge_in;
}

RHF::RHF(BasisSet basis_in, int charge_in, double tol_in, int max_iter_in)
    : RHF(basis_in, charge_in) {
    tol = tol_in;
    max_iter = max_iter_in;
}
void RHF::resolve_occupation(int n_electrons) {
    n_alpha_ = n_electrons / 2;
    n_beta_ = n_electrons / 2;
}

StabilityResult RHF::stability() const {
    if (!has_result())
        throw std::runtime_error(
            "RHF::stability: no SCF result — call run(frame) first");

    const SlaterDeterminant& det = determinant();
    const internal::IntegralView pack = integral_view_();

    // Dense working-dimension ERI: reuse the materialized tensor when the
    // run kept one (compute_eri), otherwise rebuild it from the private
    // shell workspace (Cartesian eri_dense + spherical sandwich).
    const xt::xtensor<double, 4> eri =
        (pack.ERI.size() > 0) ? pack.ERI : pack.eri_working();

    const internal::OrbitalHessian hessian{det, eri};
    const internal::StabilityAnalysis analysis{hessian};

    StabilityResult out;
    out.internal_eigenvalue = analysis.lowest(internal::Channel::singlet);
    out.external_eigenvalue = analysis.lowest(internal::Channel::triplet);
    out.internal_stable = out.internal_eigenvalue >= -kStabilityZeroTol;
    out.external_stable = out.external_eigenvalue >= -kStabilityZeroTol;

    const RealMatrix& C = det.coefficients(SlaterDeterminant::Spin::alpha);
    out.rotated_C = C;
    if (out.external_stable)
        return out;

    // Rotate the coefficients along the external instability mode: build the
    // antisymmetric occ-virt generator K from the mode amplitudes X_ai
    // (scaled so the leading amplitude rotates by kRotationTheta), then
    // C_rot = C · exp(K). Same occupied/virtual split as OrbitalHessian
    // (integer occupations, f = 1 ⇔ occupied).
    const RealVector mode = analysis.lowest_vector(internal::Channel::triplet);
    const RealVector& f = det.occupations(SlaterDeterminant::Spin::alpha);
    const int n_mo = static_cast<int>(C.shape()[1]);
    std::vector<int> occ, virt;
    for (int p = 0; p < n_mo; ++p)
        (f(p) > 0.0 ? occ : virt).push_back(p);
    const int no = static_cast<int>(occ.size());
    const int nv = static_cast<int>(virt.size());

    double x_max = 0.0;
    for (std::size_t k = 0; k < mode.size(); ++k)
        x_max = std::max(x_max, std::abs(mode(k)));
    const double scale = kRotationTheta / x_max;

    RealMatrix K = xt::zeros<double>(
        {static_cast<std::size_t>(n_mo), static_cast<std::size_t>(n_mo)});
    for (int a = 0; a < nv; ++a)
        for (int i = 0; i < no; ++i) {
            const double k_ai = scale * mode(a * no + i);
            K(virt[static_cast<std::size_t>(a)],
              occ[static_cast<std::size_t>(i)]) = k_ai;
            K(occ[static_cast<std::size_t>(i)],
              virt[static_cast<std::size_t>(a)]) = -k_ai;
        }

    // exp(K) by plain Taylor series: ‖K‖ ≤ kRotationTheta < 1, so the series
    // converges to machine precision in a handful of terms, and the exact
    // exponential of an antisymmetric K is orthogonal.
    RealMatrix U = xt::eye<double>(static_cast<std::size_t>(n_mo));
    RealMatrix term = U;
    for (int k = 1; k <= 40; ++k) {
        term = xt::linalg::dot(term, K) / static_cast<double>(k);
        U = U + term;
        double t_max = 0.0;
        for (int r = 0; r < n_mo; ++r)
            for (int c = 0; c < n_mo; ++c)
                t_max = std::max(t_max, std::abs(term(r, c)));
        if (t_max < 1e-16)
            break;
    }
    out.rotated_C = xt::linalg::dot(C, U);
    return out;
}

} // namespace scf
} // namespace cpu
} // namespace ATV
