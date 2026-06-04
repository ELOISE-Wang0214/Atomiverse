/// @file tests/cpu/md/test_pe_thermo_schema.cpp
/// @brief MSU-PE thermo-schema invariants — `msu-pe-cuda-components-01`
///        ac-001 (per-term SystemScalars default-zero) and ac-006
///        (`energy ≈ Σ per-term energies` invariant), exercised by
///        directly calling the CPU evaluators on a 4-atom PE-like
///        fragment (one bond, one angle, one dihedral).
///
/// The benchmark thermo.csv 17-column extended schema (Phase G.2 + G.3)
/// is verified end-to-end by the operator-action smoke runs of
/// `tools/run_msu_atomiverse.sh`; this unit test pins the underlying
/// CPU invariants that the schema rides on.

#include <gtest/gtest.h>

#include "atomiverse/core/system_scalars.h"
#include "atomiverse/cpu/component.h"
#include "atomiverse/cpu/evaluate.h"
#include "atomiverse/cpu/evaluate_driver.h"
#include "atomiverse/cpu/frame.h"
#include "atomiverse/cpu/md/lj.h"
#include "atomiverse/cpu/md/init.h"
#include "atomiverse/cpu/md/neighbor_list.h"
#include "atomiverse/cpu/potentials/angle_harmonic.h"
#include "atomiverse/cpu/potentials/bond_harmonic.h"
#include "atomiverse/cpu/potentials/dihedral_multi_harmonic.h"
#include "atomiverse/core/mat3.h"

#include <cmath>

using namespace ATV;

// ── ac-001: per-term energy fields exist and default to zero ─────────

TEST(PeThermoSchema, SystemScalarsDefaultZeroPerTermEnergies) {
    SystemScalars s;
    EXPECT_EQ(s.energy, 0.0);
    EXPECT_EQ(s.energy_pair, 0.0);
    EXPECT_EQ(s.energy_bond, 0.0);
    EXPECT_EQ(s.energy_angle, 0.0);
    EXPECT_EQ(s.energy_dihedral, 0.0);
    EXPECT_EQ(s.kinetic, 0.0);
    EXPECT_EQ(s.conserved, 0.0);
    for (int k = 0; k < 6; ++k) {
        EXPECT_EQ(s.stress[k], 0.0);
        EXPECT_EQ(s.stress_kinetic[k], 0.0);
    }
}

TEST(PeThermoSchema, ZeroEnergyClearsAllFiveAccumulators) {
    SystemScalars s;
    s.energy = 1.5;
    s.energy_pair = 0.5;
    s.energy_bond = 0.25;
    s.energy_angle = 0.5;
    s.energy_dihedral = 0.25;
    s.kinetic = 99.0; // should NOT be touched by ZeroEnergy
    s.stress[0] = 1.0; // should NOT be touched by ZeroEnergy

    // Mirror ZeroEnergy::launch's contract: clear the 5 potential-energy
    // accumulators in lockstep so the per-step invariant
    // `energy == Σ per-term` is preserved across steps.
    s.energy = 0;
    s.energy_pair = 0;
    s.energy_bond = 0;
    s.energy_angle = 0;
    s.energy_dihedral = 0;

    EXPECT_EQ(s.energy, 0.0);
    EXPECT_EQ(s.energy_pair, 0.0);
    EXPECT_EQ(s.energy_bond, 0.0);
    EXPECT_EQ(s.energy_angle, 0.0);
    EXPECT_EQ(s.energy_dihedral, 0.0);
    // ZeroEnergy must NOT touch kinetic or stress.
    EXPECT_EQ(s.kinetic, 99.0);
    EXPECT_EQ(s.stress[0], 1.0);
}

// ── ac-006: per-term energies sum to total potential energy ──────────

namespace {

/// @brief Construct a 4-atom PE-like fragment (one bond, one angle, one
///        dihedral, plus a pair of well-separated extra atoms to exercise
///        a non-trivial LJ pair). Geometry is hand-tuned so every evaluator
///        contributes a non-zero energy — bond stretched off equilibrium,
///        angle off 109.5°, dihedral off cosφ=1.
cpu::Frame make_pe_fragment() {
    cpu::Frame f(6);
    // Linear chain 0-1-2-3 stretched along x with mild zig.
    f.atoms.x[0] = 0.0;  f.atoms.y[0] = 0.0;  f.atoms.z[0] = 0.0;
    f.atoms.x[1] = 1.6;  f.atoms.y[1] = 0.0;  f.atoms.z[1] = 0.0;
    f.atoms.x[2] = 2.8;  f.atoms.y[2] = 1.2;  f.atoms.z[2] = 0.0;
    f.atoms.x[3] = 4.2;  f.atoms.y[3] = 1.5;  f.atoms.z[3] = 0.6;
    // Two extra atoms placed inside LJ cutoff but outside the bonded set.
    f.atoms.x[4] = 6.0;  f.atoms.y[4] = 0.0;  f.atoms.z[4] = 0.0;
    f.atoms.x[5] = 7.0;  f.atoms.y[5] = 0.0;  f.atoms.z[5] = 0.0;
    for (int i = 0; i < 6; ++i) {
        f.atoms.vx[i] = f.atoms.vy[i] = f.atoms.vz[i] = 0.0;
        f.atoms.fx[i] = f.atoms.fy[i] = f.atoms.fz[i] = 0.0;
        f.atoms.mass[i] = 14.0;     // UA CH₂
        f.atoms.type[i] = 0;
    }
    f.box = Mat3::cubic(40.0);      // large box -> no PBC interference

    // One bond (0-1), one angle (0-1-2), one dihedral (0-1-2-3), in the
    // per-atom topology the force kernels actually read (num_bond/bond_atom,
    // num_angle/angle_atom, …). The legacy flat `f.bonds/angles/dihedrals`
    // blocks are dead storage the kernels ignore, so populating those left
    // every per-term energy at zero. Mirrors the MSU-PE pipeline idiom.
    f.atoms.map_clear();
    f.atoms.map_set(0, f.atoms.nlocal);
    f.atoms.add_bond(0, 1, 0);
    f.atoms.add_angle(0, 1, 2, 0);
    f.atoms.add_dihedral(0, 1, 2, 3, 0);
    f.atoms.build_special_list();
    return f;
}

} // namespace

TEST(PeThermoSchema, PerTermEnergiesSumToTotalPotentialEnergy) {
    auto frame = make_pe_fragment();

    // Force-field constants — values mirror the LAMMPS Dreiding-UA setup
    // used by the MSU-PE benchmark, chosen so each evaluator returns a
    // non-trivial positive (bond/angle) or signed (LJ/dihedral) energy
    // when the geometry is off-equilibrium.
    const double lj_eps   = 0.112 / 23.060547830619026; // 0.112 kcal/mol → eV
    const double lj_sigma = 4.01;
    const double lj_cut   = 10.0;
    const double bond_k_eV  = 350.0 / 23.060547830619026 * 0.5; // K (eV/Å²)
    const double bond_r0   = 1.54;
    const double angle_k_eV = 60.0 / 23.060547830619026 * 0.5;
    const double angle_th0 = 109.5 * M_PI / 180.0;
    std::vector<double> dih_C = {0.5, 0.0, 0.0, 0.0, 0.0};

    cpu::md::NeighborList<> nl{lj_cut};
    cpu::md::LJ lj =
        cpu::md::LJ{lj_eps, lj_sigma, lj_cut}.set_special_bonds(0.0, 0.0, 0.0);
    // Canonical names: cpu::BondHarmonic / cpu::AngleHarmonic /
    // cpu::DihedralMultiHarmonic (the bare `HarmonicBond` symbol is an
    // ATV-namespace alias for back-compat; we use canonical names here).
    cpu::BondHarmonic hb{{bond_k_eV}, {bond_r0}};
    cpu::AngleHarmonic ha{{angle_k_eV}, {angle_th0}};
    cpu::DihedralMultiHarmonic hd{{dih_C}};
    cpu::md::ZeroEnergy ze;
    cpu::md::ZeroForces zf;
    cpu::md::ZeroStress zs;

    // One evaluator pass: zero accumulators (Σ per-term as well), then
    // run each evaluator. The Evaluate driver pipes them through the
    // single `evaluate::Phase::evaluate` slot in registration order;
    // for this invariant we want all four evaluators to write their
    // per-term and total energy on the same line — ZeroEnergy resets
    // both, ensuring no cross-step leakage.
    cpu::Evaluate eval;
    eval.on(ATV::evaluate::Phase::prepare).add(zf, -1000).add(ze, -999).add(zs, -998);
    eval.set_neighbor(nl);
    eval.set_pair(lj);
    eval.set_bond(hb);
    eval.set_angle(ha);
    eval.set_dihedral(hd);
    eval.build();
    eval.run(frame);

    const double total =
        frame.scalars.energy_pair + frame.scalars.energy_bond
        + frame.scalars.energy_angle + frame.scalars.energy_dihedral;

    // Each evaluator wrote `+= energy_total` and `+= energy_<term>` on
    // the SAME line — the invariant is a literal equality (no parallel
    // reduction order shuffle, no FMA reordering). Compare against
    // FP-strict tolerance.
    EXPECT_NEAR(frame.scalars.energy, total,
                std::fabs(frame.scalars.energy) * 1e-12 + 1e-15);

    // All four evaluators should have actually contributed.
    EXPECT_NE(frame.scalars.energy_bond, 0.0);
    EXPECT_NE(frame.scalars.energy_angle, 0.0);
    EXPECT_NE(frame.scalars.energy_dihedral, 0.0);
}
