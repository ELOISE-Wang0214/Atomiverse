/// @file cpu/scf/rks.cpp
/// @brief RKS constructor defaults, occupation resolution, and nuclear
///        Hessian via central finite difference of analytic RKS forces.
///        The shared SCF driver lives in mean_field_method.cpp.

#include "atomiverse/cpu/scf/rks.h"

#include <stdexcept>
#include <string>
#include <vector>

#include <xtensor/containers/xtensor.hpp>

#include "atomiverse/core/units.h"
#include "atomiverse/cpu/scf/gradient_result.h"

namespace ATV {
namespace cpu {
namespace scf {

RKS::RKS() {
    spin_ = SpinType::Restricted;
    method_ = MeanFieldType::KS;
    xc_type_ = XCType::PBE;
    xc = XCType::PBE;
    init_guess = InitGuess::MINAO;
    damping_factor = 0.0;
    grid_level = GridLevel::MatchPySCF;
    diis_max_vectors = 8;
    diis_start_iter = 1;
    tol = 1e-9;
}
RKS::RKS(XCType xc_in) : RKS() {
    xc = xc_in;
    xc_type_ = xc_in;
}
void RKS::resolve_occupation(int n_electrons) {
    xc_type_ = xc;
    n_alpha_ = n_electrons / 2;
    n_beta_ = n_electrons / 2;
}

vib::HessianResult RKS::hessian(Frame& frame) {
    if (ri != RI::None)
        throw std::runtime_error(
            "RKS::hessian: the analytic Hessian is not wired for density "
            "fitting (RI) — use ri = RI::None (no silent fallback)");

    // Reference SCF + energy (writes eV energy to frame.meta).
    run(frame);
    if (!converged())
        throw std::runtime_error(
            "RKS::hessian: SCF did not converge — no reference "
            "wavefunction to differentiate (" + stop_reason() + ")");

    const int natm = static_cast<int>(frame.nlocal());
    const int n3 = 3 * natm;
    // Central FD of analytic RKS forces. Match the offline PySCF pin step
    // (h = 1e-4 Bohr); convert forces eV/Å → Ha/Bohr for the FD, then
    // scale the Hessian to eV/Å².
    //   g = −F;  H_ha = ∂g/∂R [Ha/Bohr²];  H_ev = H_ha * scale.
    // Full CPKS with analytic XC second derivatives is a follow-up (the
    // ALDA/response pieces are exercised by TDDFT; the nuclear XC skeleton
    // second-derivative grid response is the remaining hard cut). Fresh
    // RKS per displacement avoids DIIS/MO state bleed.
    const double h_bohr = 1e-4;
    const double h_ang = h_bohr * units::kBohrInAngstrom;
    const double f_to_ha_bohr =
        units::kBohrInAngstrom / units::kHartreeInEv;
    const double hess_to_ev_ang2 =
        units::kHartreeInEv
        / (units::kBohrInAngstrom * units::kBohrInAngstrom);

    auto x = frame.span(kAtom_x);
    auto y = frame.span(kAtom_y);
    auto z = frame.span(kAtom_z);
    auto coord = [&](int A, int d) -> StorageReal& {
        if (d == 0)
            return x[A];
        if (d == 1)
            return y[A];
        return z[A];
    };

    std::vector<double> r0(static_cast<std::size_t>(n3));
    for (int A = 0; A < natm; ++A)
        for (int d = 0; d < 3; ++d)
            r0[static_cast<std::size_t>(3 * A + d)] =
                static_cast<double>(coord(A, d));

    xt::xtensor<double, 2> H = xt::zeros<double>({n3, n3});

    auto make_scratch = [&]() {
        RKS s(xc);
        s.basis = basis;
        s.tol = tol;
        s.max_iter = max_iter;
        s.init_guess = init_guess;
        s.damping_factor = damping_factor;
        s.diis_max_vectors = diis_max_vectors;
        s.diis_start_iter = diis_start_iter;
        s.grid_level = grid_level;
        s.n_radial = n_radial;
        s.n_angular = n_angular;
        s.conv_tol_grad = conv_tol_grad;
        return s;
    };

    for (int A = 0; A < natm; ++A) {
        for (int d = 0; d < 3; ++d) {
            const int col = 3 * A + d;
            for (int B = 0; B < natm; ++B)
                for (int e = 0; e < 3; ++e)
                    coord(B, e) = static_cast<StorageReal>(
                        r0[static_cast<std::size_t>(3 * B + e)]);
            coord(A, d) = static_cast<StorageReal>(
                r0[static_cast<std::size_t>(col)] + h_ang);
            RKS plus = make_scratch();
            const GradientResult gp = plus.forces(frame);
            if (!plus.converged())
                throw std::runtime_error(
                    "RKS::hessian: SCF failed on +h displacement of atom "
                    + std::to_string(A) + " axis " + std::to_string(d));

            coord(A, d) = static_cast<StorageReal>(
                r0[static_cast<std::size_t>(col)] - h_ang);
            RKS minus = make_scratch();
            const GradientResult gm = minus.forces(frame);
            if (!minus.converged())
                throw std::runtime_error(
                    "RKS::hessian: SCF failed on −h displacement of atom "
                    + std::to_string(A) + " axis " + std::to_string(d));

            // g = −F; H = ∂g/∂R = −∂F/∂R.
            for (int B = 0; B < natm; ++B)
                for (int e = 0; e < 3; ++e) {
                    const int row = 3 * B + e;
                    const double Fp =
                        static_cast<double>(gp.forces(B, e)) * f_to_ha_bohr;
                    const double Fm =
                        static_cast<double>(gm.forces(B, e)) * f_to_ha_bohr;
                    const double dF_dR = (Fp - Fm) / (2.0 * h_bohr);
                    H(row, col) = (-dF_dR) * hess_to_ev_ang2;
                }
        }
    }

    // Restore geometry and re-run reference SCF.
    for (int A = 0; A < natm; ++A)
        for (int d = 0; d < 3; ++d)
            coord(A, d) = static_cast<StorageReal>(
                r0[static_cast<std::size_t>(3 * A + d)]);
    run(frame);

    // Symmetrize (FD noise).
    for (int i = 0; i < n3; ++i)
        for (int j = 0; j < i; ++j) {
            const double s = 0.5 * (H(i, j) + H(j, i));
            H(i, j) = s;
            H(j, i) = s;
        }

    vib::HessianResult out;
    out.n_atoms = natm;
    out.hessian = std::move(H);
    return out;
}

} // namespace scf
} // namespace cpu
} // namespace ATV
