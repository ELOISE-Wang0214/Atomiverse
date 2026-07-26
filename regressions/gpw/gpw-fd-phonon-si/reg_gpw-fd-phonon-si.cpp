/// @file regressions/gpw/gpw-fd-phonon-si/reg_gpw-fd-phonon-si.cpp
/// @brief END-TO-END: Γ-point FD phonons with GPW analytic forces (Si₂).
///
/// See REFERENCE.md. Public API only. Requires ATV_BUILD_GPW (+ HF + PW) and
/// vibration. Exit non-zero with got/want/tol. No live CP2K.

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

// Offline CP2K 2026.1 VIBRATIONAL_ANALYSIS optical pins (cm^-1),
// nonperiodic MT twin, FULLY_PERIODIC T → 3N−3 modes. See REFERENCE.md.
constexpr double kOmega_cp2k_opt[3] = {-263.48092846, -194.72181159,
                                       383.78258256};
constexpr double kOptTol = 50.0;      // cm^-1
constexpr double kAcousticTol = 80.0; // cm^-1
constexpr int kNAtoms = 2;
constexpr int kNModes = 3 * kNAtoms; // 6
constexpr double kA0 = 5.43;         // Å

} // namespace

int main() {
    cpu::Frame frame(kNAtoms);
    const double pos[2][3] = {{0.0, 0.0, 0.0},
                              {kA0 / 4.0, kA0 / 4.0, kA0 / 4.0}};
    for (int i = 0; i < kNAtoms; ++i) {
        frame.span(kAtom_type)[i] = 14;
        frame.span(kAtom_mass)[i] = 28.0855;
        frame.span(kAtom_tag)[i] = static_cast<TagInt>(i + 1);
        frame.span(kAtom_x)[i] = static_cast<StorageReal>(pos[i][0]);
        frame.span(kAtom_y)[i] = static_cast<StorageReal>(pos[i][1]);
        frame.span(kAtom_z)[i] = static_cast<StorageReal>(pos[i][2]);
        frame.span(kAtom_fx)[i] = frame.span(kAtom_fy)[i] =
            frame.span(kAtom_fz)[i] = 0;
    }
    frame.simbox = cpu::Box::cubic(kA0);

    const cpu::gpw::Kind Si{
        molrs::Element::Si,
        cpu::gpw::GaussianBasis(molrs::Element::Si, "DZVP-MOLOPT-GTH"),
        cpu::gpw::GthPotential(molrs::Element::Si, "GTH-PBE-q4")};

    // nonperiodic=true: MT/cluster Poisson. Periodic Γ MOLOPT image sums
    // exceed the ≤60 s e2e budget (same product constraint as reg_gpw-pbesol-si).
    cpu::gpw::GPWEvaluator eval(
        /*kinds*/ {Si}, /*cutoff_ry*/ 150.0, /*max_iter*/ 80, /*tol*/ 1e-7,
        /*mix*/ 0.4, /*xc*/ cpu::gpw::XcKind::PADE, /*nonperiodic*/ true,
        /*compute_forces*/ true);

    // No NeighborList: GPW is not a NeighborConsumer.
    cpu::vib::Phonon phonon(/*na=*/1, /*nb=*/1, /*nc=*/1);
    phonon.on(ATV::vib::Phase::force).add_ref(eval);
    phonon.enable_timing(false);
    phonon.build();

    const auto t0 = std::chrono::steady_clock::now();
    const cpu::vib::PhononResult r = phonon.run(frame);
    const auto t1 = std::chrono::steady_clock::now();
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
    int n_ac = 0;
    for (int k = 0; k < kNModes; ++k) {
        w[static_cast<std::size_t>(k)] = r.frequencies(k);
        if (!std::isfinite(w[static_cast<std::size_t>(k)])) {
            std::printf("FAIL: non-finite ω[%d]=%g\n", k,
                        w[static_cast<std::size_t>(k)]);
            return 1;
        }
        if (std::abs(w[static_cast<std::size_t>(k)]) < kAcousticTol)
            ++n_ac;
        std::printf("  ω[%d]=% .6f cm^-1\n", k, w[static_cast<std::size_t>(k)]);
    }
    if (n_ac < 3) {
        std::printf("FAIL: acoustic residual count got=%d want>=3 "
                    "(|ω|<%.0f cm^-1; Γ ASR)\n",
                    n_ac, kAcousticTol);
        return 1;
    }

    // Optical pins: drop 3 near-zero (by |ω|), sort remaining signed ascending.
    std::vector<double> by_abs = w;
    std::sort(by_abs.begin(), by_abs.end(),
              [](double a, double b) { return std::abs(a) < std::abs(b); });
    std::vector<double> opt(by_abs.begin() + 3, by_abs.end());
    std::sort(opt.begin(), opt.end());

    double max_dopt = 0.0;
    for (int k = 0; k < 3; ++k) {
        const double d = std::abs(opt[static_cast<std::size_t>(k)]
                                  - kOmega_cp2k_opt[k]);
        max_dopt = std::max(max_dopt, d);
        std::printf("  opt[%d] atv=% .4f  cp2k=% .4f  |d|=%.3f cm^-1\n", k,
                    opt[static_cast<std::size_t>(k)], kOmega_cp2k_opt[k], d);
    }
    if (max_dopt > kOptTol) {
        std::printf("FAIL: max|ω_opt-cp2k|=%.3f > tol %.1f cm^-1\n", max_dopt,
                    kOptTol);
        return 1;
    }

    if (wall > 60.0) {
        std::printf("FAIL: wall=%.2f s exceeds 60 s e2e budget\n", wall);
        return 1;
    }

    std::printf("PASS: GPW FD phonon Si2 (MT)  a=%.3f A  N=%d  modes=%d  "
                "n_acoustic(|ω|<%.0f)=%d  max|dω_opt|=%.3f cm^-1  wall=%.2f s\n",
                kA0, kNAtoms, kNModes, kAcousticTol, n_ac, max_dopt, wall);
    return 0;
}
