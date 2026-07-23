/// @file cpu/scf/rohf.cpp
/// @brief ROHF constructor defaults, occupation resolution, and nuclear
///        Hessian via central finite difference of analytic ROHF forces.
///        The shared driver loop lives in mean_field_method.cpp; the
///        Roothaan effective-Fock coupling lives in stages/rohf_couple_stage.h.

#include "atomiverse/cpu/scf/rohf.h"

#include <stdexcept>
#include <string>
#include <vector>

#include <xtensor/containers/xtensor.hpp>

#include "atomiverse/core/units.h"
#include "atomiverse/cpu/scf/gradient_result.h"

namespace ATV {
namespace cpu {
namespace scf {

ROHF::ROHF() {
    spin_ = SpinType::RestrictedOpen;
    method_ = MeanFieldType::HF;
    init_guess = InitGuess::Hcore;
    damping_factor = 0.7;
    diis_max_vectors = 6;
    diis_start_iter = 5;
    tol = 1e-8;
}
ROHF::ROHF(int na, int nb) : ROHF() {
    n_alpha = na;
    n_beta = nb;
}
void ROHF::resolve_occupation(int) {
    if (n_alpha < n_beta)
        throw std::invalid_argument(
            "ROHF::run: n_alpha (" + std::to_string(n_alpha)
            + ") < n_beta (" + std::to_string(n_beta)
            + ") — ROHF defines the open shell as excess alpha electrons; "
              "swap the counts (no fallback).");
    n_alpha_ = n_alpha;
    n_beta_ = n_beta;
}

vib::HessianResult ROHF::hessian(Frame& frame) {
    if (ri != RI::None)
        throw std::runtime_error(
            "ROHF::hessian: the analytic Hessian is not wired for density "
            "fitting (RI) — use ri = RI::None (no silent fallback)");

    // Reference SCF + energy (writes eV energy to frame.meta).
    run(frame);
    if (!converged())
        throw std::runtime_error(
            "ROHF::hessian: SCF did not converge — no reference "
            "wavefunction to differentiate (" + stop_reason() + ")");

    const int natm = static_cast<int>(frame.nlocal());
    const int n3 = 3 * natm;
    // Central FD of analytic ROHF forces. Match the PySCF pin step
    // (h = 1e-4 Bohr); convert forces eV/Å → Ha/Bohr for the FD, then
    // scale the Hessian to eV/Å².
    //   g = −F;  H_ha = ∂g/∂R [Ha/Bohr²];  H_ev = H_ha * scale.
    // Full three-space ROHF CPHF is a follow-up (PySCF also lacks
    // hessian.ROHF). Fresh ROHF per displacement avoids DIIS/MO state bleed.
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
        ROHF s(n_alpha, n_beta);
        s.basis = basis;
        s.tol = tol;
        s.max_iter = max_iter;
        s.init_guess = init_guess;
        s.damping_factor = damping_factor;
        s.diis_max_vectors = diis_max_vectors;
        s.diis_start_iter = diis_start_iter;
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
            ROHF plus = make_scratch();
            const GradientResult gp = plus.forces(frame);
            if (!plus.converged())
                throw std::runtime_error(
                    "ROHF::hessian: SCF failed on +h displacement of atom "
                    + std::to_string(A) + " axis " + std::to_string(d));

            coord(A, d) = static_cast<StorageReal>(
                r0[static_cast<std::size_t>(col)] - h_ang);
            ROHF minus = make_scratch();
            const GradientResult gm = minus.forces(frame);
            if (!minus.converged())
                throw std::runtime_error(
                    "ROHF::hessian: SCF failed on −h displacement of atom "
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

    for (int A = 0; A < natm; ++A)
        for (int d = 0; d < 3; ++d)
            coord(A, d) = static_cast<StorageReal>(
                r0[static_cast<std::size_t>(3 * A + d)]);
    run(frame);

    for (int i = 0; i < n3; ++i)
        for (int j = 0; j < i; ++j) {
            const double s = 0.5 * (H(i, j) + H(j, i));
            H(i, j) = s;
            H(j, i) = s;
        }

    vib::HessianResult out;
    out.n_atoms = natm;
    out.hessian = H;
    return out;
}

} // namespace scf
} // namespace cpu
} // namespace ATV
