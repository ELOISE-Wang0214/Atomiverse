/// @file tests/cpu/gpw/test_cpu_gpw_evaluator.cpp
/// @brief GPWEvaluator plugs into the Evaluate driver and writes frame energy.
///
/// The `GpwEvaluatorSolver` suite (arch-refactor-04-engines Design §2) pins the
/// evaluator behaviour across the switch to an owned GPWSolver: the energy is
/// bit-invariant vs a direct gpw_scf solve of the identically-assembled system,
/// and launching after the frame box changed from the build-time cell throws
/// std::runtime_error (GPW does not support variable cells). These cases use a
/// fully inline basis/pseudopotential — no CP2K data files required.

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <gtest/gtest.h>

#include "atomiverse/core/units.h"
#include "atomiverse/cpu/component/component.h"
#include "atomiverse/cpu/driver/single_point.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/gpw/gpw_evaluator.h"
#include "atomiverse/cpu/gpw/gpw_solver.h"
#include "atomiverse/cpu/gpw/internal/cp2k_basis.h"
#include "atomiverse/cpu/gpw/internal/cp2k_io.h"
#include "atomiverse/cpu/gpw/basis_data/gpw_data.h"
#include "atomiverse/facade/component_role.h"

using namespace ATV::cpu::gpw;

namespace {
constexpr double kBohr2Ang = 1.0 / 1.8897259885789;
constexpr double kAng2Bohr =
    1.8897259885789; // same constant as gpw_evaluator.cpp

// Inline H element: one L2-normalized s primitive (alpha = 0.5) + GTH-PADE-q1
// local pseudopotential — the minimal converging H2 system of
// test_cpu_gpw_solver.cpp, expressed as the GaussianBasis/GthPotential structs
// the evaluator consumes. No file parsing.
GpwElement inline_hydrogen_element() {
    GaussianBasis basis;
    basis.element = "H";
    basis.name = "INLINE-S";
    GaussianBasisSet set;
    set.n = 1;
    set.lmin = 0;
    set.lmax = 0;
    set.exps = {0.5};
    set.nshell = {1};
    set.ncoef = 1;
    set.coef = {1.0};
    basis.sets = {set};

    GthPotential pot;
    pot.element = "H";
    pot.name = "INLINE-q1";
    pot.zion = 1;
    pot.rloc = 0.2;
    pot.cloc = {-4.17890044, 0.72446331};

    GpwElement e;
    e.z = 1;
    e.basis = basis;
    e.pot = pot;
    return e;
}

// H2 at 1.4 Bohr in a cubic 8 Bohr box, expressed in the frame's Angstrom.
ATV::cpu::Frame inline_h2_frame() {
    const double box = 8.0; // [Bohr]
    ATV::cpu::Frame frame(2);
    frame.box = ATV::cpu::Box::orthorhombic(box * kBohr2Ang, box * kBohr2Ang,
                                            box * kBohr2Ang);
    frame.atoms.x(0) = (box / 2 - 0.7) * kBohr2Ang;
    frame.atoms.y(0) = (box / 2) * kBohr2Ang;
    frame.atoms.z(0) = (box / 2) * kBohr2Ang;
    frame.atoms.x(1) = (box / 2 + 0.7) * kBohr2Ang;
    frame.atoms.y(1) = (box / 2) * kBohr2Ang;
    frame.atoms.z(1) = (box / 2) * kBohr2Ang;
    frame.atoms.type(0) = 1; // H
    frame.atoms.type(1) = 1;
    return frame;
}

GpwEvaluatorConfig inline_h2_evaluator_config() {
    GpwEvaluatorConfig cfg;
    cfg.cutoff_ry = 120.0;
    cfg.max_iter = 100;
    cfg.tol = 1e-7;
    cfg.mix = 0.5;
    cfg.xc = XcKind::PADE;
    cfg.nonperiodic = false;
    cfg.elements = {inline_hydrogen_element()};
    return cfg;
}

void pin_single_thread() {
#ifdef _OPENMP
    omp_set_num_threads(1); // bit-stable FP sums; see test_cpu_gpw_solver.cpp
#endif
}
} // namespace

TEST(GpwEvaluator, H2SinglePointThroughDriver) {
    // H2 at 1 Bohr separation in a 16 Bohr cubic box (large enough that the
    // diffuse DZV-GTH-PADE images are negligible — see the gpw-periodic-lattice-sum
    // note; the tight 8x4x4 box over-images this very diffuse PADE basis). This
    // test exercises the evaluator → Evaluate-driver → frame.energy plumbing;
    // periodic CP2K parity for the production basis is in the MOLOPT zno/nacl tests.
    ATV::cpu::Frame frame(2);
    frame.box = ATV::cpu::Box::orthorhombic(16.0 * kBohr2Ang, 16.0 * kBohr2Ang,
                                            16.0 * kBohr2Ang);
    frame.atoms.x(0) = 3.5 * kBohr2Ang;
    frame.atoms.y(0) = 2.0 * kBohr2Ang;
    frame.atoms.z(0) = 2.0 * kBohr2Ang;
    frame.atoms.x(1) = 4.5 * kBohr2Ang;
    frame.atoms.y(1) = 2.0 * kBohr2Ang;
    frame.atoms.z(1) = 2.0 * kBohr2Ang;
    frame.atoms.type(0) = 1; // H
    frame.atoms.type(1) = 1;

    GpwElement H;
    H.z = 1;
    H.basis = builtin_gth_basis("H", "DZV-GTH-PADE");
    H.pot = builtin_gth_potential("H", "GTH-PADE-q1");

    GpwEvaluatorConfig cfg;
    cfg.cutoff_ry = 140.0;
    cfg.xc = XcKind::PADE;
    cfg.nonperiodic = false;
    cfg.mix = 0.3;
    cfg.elements = {H};

    GPWEvaluator eval(cfg);
    ATV::cpu::single_point(frame, eval);

    // The evaluator computed n_occ=1 from 2 valence electrons, ran the SCF, and
    // wrote eV onto the frame. Compare (loosely) to converged isolated CP2K H2
    // (-1.10676 Ha): the DZV-GTH-PADE/PADE-LDA path sits ~30 mHa above CP2K
    // (a pre-existing basis/grid offset, documented in gpw-periodic-lattice-sum);
    // the gate brackets that offset so the test validates plumbing + rough physics.
    const double e_ev = frame.scalars.energy;
    const double ref_ev = -1.10676359396999 * 27.211386245988;
    EXPECT_TRUE(eval.last_result().converged);
    EXPECT_NEAR(e_ev, ref_ev,
                1.0) // ~37 mHa: brackets the pre-existing PADE-path offset
        << "frame.energy=" << e_ev << " eV  ref=" << ref_ev << " eV";
    // Frame energy must equal the solver's Hartree result times the conversion.
    EXPECT_NEAR(e_ev, eval.last_energy_hartree() * 27.211386245988, 1e-9);
}

// ── Evaluator-owns-GPWSolver contract (arch-refactor-04-engines Design §2) ──

// Numerical anchor (must stay green across the refactor): the evaluator's
// energy is BIT-IDENTICAL to a direct gpw_scf solve of the same system,
// assembled exactly the way gpw_evaluator.cpp assembles it (Angstrom → Bohr,
// n_occ = Σ zion / 2, cold start). Launch keeps calling the same inner solve —
// switching to an owned GPWSolver must not move a single bit.
TEST(GpwEvaluatorSolver, InlineH2EnergyBitMatchesDirectSolve) {
    pin_single_thread();
    ATV::cpu::Frame frame = inline_h2_frame();
    const GpwEvaluatorConfig ecfg = inline_h2_evaluator_config();

    GPWEvaluator eval(ecfg);
    ATV::cpu::single_point(frame, eval);
    ASSERT_TRUE(eval.last_result().converged);

    // Reference: mirror the evaluator's own atom/config assembly bit-for-bit.
    const GpwElement H = inline_hydrogen_element();
    std::vector<GpwAtom> atoms;
    for (int i = 0; i < 2; ++i) {
        const double pos[3] = {frame.atoms.x(i) * kAng2Bohr,
                               frame.atoms.y(i) * kAng2Bohr,
                               frame.atoms.z(i) * kAng2Bohr};
        atoms.push_back(make_gpw_atom(pos, H.basis, H.pot, 1));
    }
    GpwConfig cfg;
    cfg.cell[0] = frame.box.length(0) * kAng2Bohr;
    cfg.cell[1] = frame.box.length(1) * kAng2Bohr;
    cfg.cell[2] = frame.box.length(2) * kAng2Bohr;
    cfg.cutoff_ry = ecfg.cutoff_ry;
    cfg.n_occ = 1; // Σ zion = 2 valence electrons, closed shell
    cfg.max_iter = ecfg.max_iter;
    cfg.tol = ecfg.tol;
    cfg.mix = ecfg.mix;
    cfg.xc = ecfg.xc;
    cfg.nonperiodic = ecfg.nonperiodic;
    const GpwResult ref = gpw_scf(atoms, cfg);
    ASSERT_TRUE(ref.converged);

    EXPECT_DOUBLE_EQ(eval.last_energy_hartree(), ref.e_total);
    EXPECT_DOUBLE_EQ(frame.scalars.energy, ref.e_total * ATV::kHartree2eV);
}

// Variable-cell guard: GPW does not support a cell that changes after build()
// (the build-time cell fixes the grids). Launching with a frame box that
// differs from the build-time cell must throw std::runtime_error instead of
// silently solving in the new box (previously undefined behaviour).
TEST(GpwEvaluatorSolver, LaunchAfterCellChangeThrowsRuntimeError) {
    ATV::cpu::Frame frame = inline_h2_frame();

    GPWEvaluator eval(inline_h2_evaluator_config());
    ATV::cpu::Context runtime;
    eval.build(ATV::cpu::BuildContext{runtime, frame});

    // NPT-style cell change between build and launch.
    const double a = frame.box.length(0);
    const double b = frame.box.length(1);
    const double c = frame.box.length(2);
    frame.box = ATV::cpu::Box::orthorhombic(1.10 * a, b, c);

    ATV::cpu::LaunchContext ctx{runtime, frame, 0, 0};
    EXPECT_THROW(eval.launch(ctx), std::runtime_error);
}
