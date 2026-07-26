/// @file regressions/gpw/gpw-fd-phonon-si/reg_gpw-fd-phonon-si.cpp
/// @brief END-TO-END: finite-displacement vibrational spectrum from GPW
///        analytic forces, Si₂ (UKS triplet, isolated) vs offline CP2K.
///
/// Rebuilt 2026-07-26 on corrected physics — see REFERENCE.md "History". The
/// three defects that made the first attempt produce garbage were in the SETUP,
/// not in the force code:
///   1. Si₂'s ground state is \(^3\Sigma_g^-\); running it closed-shell (RKS)
///      lets OT report `converged` on a saddle of the closed-shell functional,
///      where the frozen-density force is not the energy gradient.
///   2. An atom at (0,0,0) with Martyna–Tuckerman cluster boundaries clips the
///      collocated density, breaking translation invariance of the energy.
///   3. An isolated-cluster spectrum was compared against a `FULLY_PERIODIC T`
///      crystal pin. An isolated diatomic has ONE stretch plus five near-zero
///      modes (3 translations + 2 rotations), not three optical modes.
///
/// This case therefore runs the dimer as a UKS triplet, centred in the cell, and
/// compares against a CP2K `VIBRATIONAL_ANALYSIS` run on the identical system.
///
/// Public API only. No live CP2K. See REFERENCE.md.

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <vector>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/core/stage.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/gpw/gaussian_basis.h"
#include "atomiverse/cpu/gpw/gpw_evaluator.h"
#include "atomiverse/cpu/gpw/gth_potential.h"
#include "atomiverse/cpu/vibration/phonon.h"

using namespace ATV;

namespace {

// ── Offline CP2K 2026.1 pins (see REFERENCE.md) ────────────────────────────
// Si2 UKS triplet, isolated (PERIODIC NONE + PSOLVER MT), r = 2.246 A, 10 A
// cube, PADE, CUTOFF 400 / REL_CUTOFF 40 / NGRIDS 4, EPS_SCF 1e-7, DX = 1e-3.
// NOTE: this file is NOT registered in regressions/CMakeLists.txt yet. The pins
// below are real and were taken this wave, but Atomiverse's SCF does not
// converge at these settings for the open-shell dimer (linear/DIIS/OT vs CP2K's
// Broyden), and 12 displacements at ~38 s per SCF would exceed the 60 s e2e tier
// anyway. Re-enabling is a one-line CMake change once the mixer lands — the full
// convergence/cost matrix is in REFERENCE.md.
constexpr double kOmegaCp2k = 463.489404; // cm^-1, the single stretch
constexpr double kOmegaTol = 25.0;         // cm^-1 (see REFERENCE tolerances)
constexpr double kZeroModeTol = 40.0;      // cm^-1 for translations/rotations
constexpr double kForceXCp2k = 1.99658355e-02; // Ha/Bohr, atom 0, +x
constexpr double kForceTol = 5.0e-3;           // Ha/Bohr, same-method twin

constexpr int kNAtoms = 2;
constexpr int kNModes = 3 * kNAtoms; // 6: one stretch + 5 near-zero
constexpr int kNZeroModes = 5;
constexpr double kL = 10.0;      // A, cubic cell
constexpr double kBond = 2.246;  // A, Si-Si
constexpr double kWallBudget = 60.0;

} // namespace

int main() {
    using clock = std::chrono::steady_clock;

    // Centred on the cell: the MT cluster boundary condition is not translation
    // invariant, and GpwDft::build rejects off-centre nonperiodic systems.
    const double c = 0.5 * kL;
    const double half = 0.5 * kBond;
    const double xs[kNAtoms] = {c - half, c + half};

    cpu::Frame frame(kNAtoms);
    for (int i = 0; i < kNAtoms; ++i) {
        frame.span(kAtom_type)[i] = 14;
        frame.span(kAtom_mass)[i] = 28.0855;
        frame.span(kAtom_tag)[i] = static_cast<TagInt>(i + 1);
        frame.span(kAtom_x)[i] = static_cast<StorageReal>(xs[i]);
        frame.span(kAtom_y)[i] = static_cast<StorageReal>(c);
        frame.span(kAtom_z)[i] = static_cast<StorageReal>(c);
        frame.span(kAtom_fx)[i] = frame.span(kAtom_fy)[i] =
            frame.span(kAtom_fz)[i] = 0;
    }
    frame.simbox = cpu::Box::cubic(kL);

    const cpu::gpw::Kind Si{
        molrs::Element::Si,
        cpu::gpw::GaussianBasis(molrs::Element::Si, "DZVP-MOLOPT-GTH"),
        cpu::gpw::GthPotential(molrs::Element::Si, "GTH-PBE-q4")};

    // UKS triplet: 8 valence e- as n_alpha = 5, n_beta = 3 (CP2K LSD +
    // MULTIPLICITY 3). This is Si2's actual ground state; forcing the
    // closed-shell solution is what broke the first attempt.
    cpu::gpw::GPWEvaluator eval({Si});
    eval.cutoff_ry = 400.0;
    eval.max_iter = 400;
    eval.tol = 1e-7;
    eval.mix = 0.4;
    eval.xc = cpu::gpw::XcKind::PADE;
    eval.nonperiodic = true;
    eval.n_alpha = 5;
    eval.n_beta = 3;
    eval.compute_forces = true;

    // No NeighborList: GPW is not a NeighborConsumer.
    cpu::vib::Phonon phonon(/*na=*/1, /*nb=*/1, /*nc=*/1);
    phonon.on(ATV::vib::Phase::force).add_ref(eval);
    phonon.enable_timing(false);
    phonon.build();

    const auto t0 = clock::now();
    const cpu::vib::PhononResult r = phonon.run(frame);
    const auto t1 = clock::now();
    const double wall = std::chrono::duration<double>(t1 - t0).count();

    if (r.n_atoms != kNAtoms) {
        std::printf("FAIL: n_atoms got=%d want=%d\n", r.n_atoms, kNAtoms);
        return 1;
    }
    if (static_cast<int>(r.frequencies.size()) != kNModes) {
        std::printf("FAIL: n_modes got=%zu want=%d\n",
                    static_cast<std::size_t>(r.frequencies.size()), kNModes);
        return 1;
    }

    std::vector<double> w(kNModes);
    for (int k = 0; k < kNModes; ++k) {
        w[static_cast<std::size_t>(k)] = r.frequencies(k);
        if (!std::isfinite(w[static_cast<std::size_t>(k)])) {
            std::printf("FAIL: non-finite omega[%d]=%g\n", k,
                        w[static_cast<std::size_t>(k)]);
            return 1;
        }
        std::printf("  omega[%d] = % .4f cm^-1\n", k,
                    w[static_cast<std::size_t>(k)]);
    }

    // Five modes must be near zero (3 translations + 2 rotations of a linear
    // molecule); the remaining one is the stretch.
    std::vector<double> by_abs = w;
    std::sort(by_abs.begin(), by_abs.end(),
              [](double a, double b) { return std::abs(a) < std::abs(b); });
    int n_zero = 0;
    for (int k = 0; k < kNZeroModes; ++k)
        if (std::abs(by_abs[static_cast<std::size_t>(k)]) < kZeroModeTol)
            ++n_zero;
    if (n_zero != kNZeroModes) {
        std::printf("FAIL: got %d of %d near-zero modes (|omega| < %.0f cm^-1);"
                    " translations/rotations are not clean\n",
                    n_zero, kNZeroModes, kZeroModeTol);
        return 1;
    }

    const double omega = by_abs[static_cast<std::size_t>(kNModes - 1)];
    const double d_omega = std::abs(omega - kOmegaCp2k);
    std::printf("  stretch: atv=%.4f  cp2k=%.4f  |d|=%.3f cm^-1\n", omega,
                kOmegaCp2k, d_omega);
    if (d_omega > kOmegaTol) {
        std::printf("FAIL: |omega_atv - omega_cp2k| = %.3f > %.1f cm^-1\n",
                    d_omega, kOmegaTol);
        return 1;
    }

    // Equilibrium-geometry force cross-check against the same CP2K twin: the
    // Hessian is only as good as the forces it differentiates.
    const auto fx = frame.span(kAtom_fx);
    const double f0_ev_a = static_cast<double>(fx[0]);
    const double f0_ha_bohr =
        f0_ev_a / ATV::units::kHartreePerBohrInEvPerAngstrom;
    const double d_f = std::abs(std::abs(f0_ha_bohr) - kForceXCp2k);
    std::printf("  F_x(atom 0): atv=% .6e  cp2k=% .6e  |d|=%.3e Ha/Bohr\n",
                f0_ha_bohr, kForceXCp2k, d_f);
    if (d_f > kForceTol) {
        std::printf("FAIL: |F_x| mismatch %.3e > %.1e Ha/Bohr\n", d_f,
                    kForceTol);
        return 1;
    }

    if (wall > kWallBudget) {
        std::printf("FAIL: wall=%.2f s exceeds %.0f s e2e budget\n", wall,
                    kWallBudget);
        return 1;
    }

    std::printf("PASS: gpw-fd-phonon-si  modes=%d  near-zero=%d  "
                "stretch=%.4f cm^-1 (cp2k %.4f, |d|=%.3f)  "
                "|dF_x|=%.3e Ha/Bohr  wall=%.2f s\n",
                kNModes, n_zero, omega, kOmegaCp2k, d_omega, d_f, wall);
    return 0;
}
