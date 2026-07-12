/// @file tests/cpu/md/test_fix_rattle.cpp
/// @brief Unit: RATTLE holonomic constraints — the invariants, the DOF accounting,
///        and the role/validation gate.
///
/// A constraint is not judged by a number it reproduces; it is judged by an
/// INVARIANT it never breaks. So: the bond lengths must be constant to machine
/// precision for the whole trajectory, and the relative velocity along every bond
/// must vanish after the velocity stage.

#include <cmath>
#include <vector>

#include <gtest/gtest.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "atomiverse/core/mat3.h"
#include "atomiverse/cpu/component/context.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/md/fix_rattle.h"
#include "atomiverse/cpu/md/md.h"
#include "atomiverse/cpu/md/velocity_verlet.h"
#include "atomiverse/cpu/neighbor/neighbor_list.h"
#include "atomiverse/cpu/potentials/angle_harmonic.h"
#include "atomiverse/cpu/potentials/bond_harmonic.h"
#include "atomiverse/cpu/potentials/pair_lj_cut.h"

using namespace ATV;

namespace {

/// Pin the thread count. These are 24-atom systems, and on a 48-core node the
/// OpenMP fork/join overhead alone costs ~4 s per test while the physics costs
/// ~2 ms (measured: 48 threads 3974 ms, 4 threads 5 ms, 1 thread 2 ms). A unit
/// test must not be at the mercy of the machine's core count.
///
/// Noted as an API gap: the engine takes its thread count from
/// `omp_get_max_threads()` and exposes no public setter, so a user with a small
/// system has no in-code way to ask for fewer threads.
struct PinThreads {
    PinThreads() {
#ifdef _OPENMP
        omp_set_num_threads(2);
#endif
    }
};
const PinThreads kPinThreads;

// TIP3P geometry: r(OH) = 0.9572 A, theta(HOH) = 104.52 deg.
constexpr double kROH = 0.9572;
constexpr double kThetaDeg = 104.52;
// The H-H distance that encodes the angle: d = 2 r sin(theta/2).
const double kDHH = 2.0 * kROH * std::sin(0.5 * kThetaDeg * M_PI / 180.0);

/// A box of rigid-geometry waters on a lattice — several call sites below.
cpu::Frame water_box(int n_side, double spacing) {
    const int n_mol = n_side * n_side * n_side;
    const int n = 3 * n_mol;
    cpu::Frame f(n);

    const double half = 0.5 * kThetaDeg * M_PI / 180.0;
    int idx = 0, mol = 0;
    for (int ix = 0; ix < n_side; ++ix)
        for (int iy = 0; iy < n_side; ++iy)
            for (int iz = 0; iz < n_side; ++iz) {
                const double ox = (ix + 0.5) * spacing;
                const double oy = (iy + 0.5) * spacing;
                const double oz = (iz + 0.5) * spacing;
                // Rotate each molecule a little so they are not all aligned.
                const double phi = 0.7 * double(mol);

                const int o = idx, h1 = idx + 1, h2 = idx + 2;
                f.atoms.x[o] = ox; f.atoms.y[o] = oy; f.atoms.z[o] = oz;
                f.atoms.mass[o] = 15.9994; f.atoms.type[o] = 8;

                f.atoms.x[h1] = ox + kROH * std::sin(half) * std::cos(phi);
                f.atoms.y[h1] = oy + kROH * std::sin(half) * std::sin(phi);
                f.atoms.z[h1] = oz + kROH * std::cos(half);
                f.atoms.mass[h1] = 1.008; f.atoms.type[h1] = 1;

                f.atoms.x[h2] = ox - kROH * std::sin(half) * std::cos(phi);
                f.atoms.y[h2] = oy - kROH * std::sin(half) * std::sin(phi);
                f.atoms.z[h2] = oz + kROH * std::cos(half);
                f.atoms.mass[h2] = 1.008; f.atoms.type[h2] = 1;

                idx += 3;
                ++mol;
            }
    f.box = Mat3::cubic(n_side * spacing);
    f.atoms.map_set(0, n);

    // Topology through the loader API (it owns the LAMMPS ownership rules: the
    // smaller-tag atom owns a bond; the centre atom owns an angle).
    for (int m = 0; m < n_mol; ++m) {
        const int o = 3 * m, h1 = o + 1, h2 = o + 2;
        f.atoms.add_bond(f.atoms.tag(o), f.atoms.tag(h1), 0);
        f.atoms.add_bond(f.atoms.tag(o), f.atoms.tag(h2), 0);
        f.atoms.add_angle(f.atoms.tag(h1), f.atoms.tag(o), f.atoms.tag(h2), 0);
    }
    // Intramolecular pairs MUST be excluded from the pair style. The O and its H's
    // sit 0.96 A apart — well inside sigma — so an unexcluded LJ would blast the
    // molecule apart and no constraint solver could hold it. (This is exactly what
    // my first draft of this test did: the position stage failed to converge, and
    // the bug was in the test, not in RATTLE.)
    f.atoms.build_special_list();
    return f;
}

/// Worst relative deviation of any constrained distance from its target.
double worst_geometry_error(const cpu::Frame& f, int n_mol) {
    double worst = 0.0;
    for (int m = 0; m < n_mol; ++m) {
        const int o = 3 * m, h1 = o + 1, h2 = o + 2;
        const std::pair<std::pair<int, int>, double> checks[3] = {
            {{o, h1}, kROH}, {{o, h2}, kROH}, {{h1, h2}, kDHH}};
        for (const auto& c : checks) {
            double dx = f.atoms.x[c.first.first] - f.atoms.x[c.first.second];
            double dy = f.atoms.y[c.first.first] - f.atoms.y[c.first.second];
            double dz = f.atoms.z[c.first.first] - f.atoms.z[c.first.second];
            f.box.min_image(dx, dy, dz);
            const double r = std::sqrt(dx * dx + dy * dy + dz * dz);
            worst = std::max(worst, std::fabs(r - c.second) / c.second);
        }
    }
    return worst;
}

} // namespace

/// The DOF accounting: 3 constraints per rigid water, so 9 - 3 = 6 DOF per molecule.
TEST(FixRattle, ReportsThreeConstraintsPerWater) {
    cpu::Frame f = water_box(2, 6.0); // 8 molecules
    cpu::md::FixRattle r{1.0};
    r.constrain_bond(0, kROH);
    r.constrain_angle(0, kDHH);
    const cpu::Context rt{}; // single-rank
    EXPECT_EQ(r.n_constraints(rt, f), 3 * 8);
}

/// The gate the ComponentRole::Constraint entry exists to keep shut: a Plan with a
/// constraint but NO integrator must not run. Before Constraint existed, tagging the
/// component as an Integrator would have satisfied DynamicsDriver's validation and
/// the system would have "run" with nothing propagating the equations of motion.
TEST(FixRattle, ConstraintDoesNotSatisfyTheIntegratorRequirement) {
    cpu::Frame f = water_box(2, 6.0);
    cpu::md::MD sim;
    sim.set_pair(cpu::PairLJCut{0.1, 3.0, 8.0});
    sim.set_neighbor(cpu::md::NeighborList<>{8.0, 0.3});
    cpu::md::FixRattle rattle{1.0};
    rattle.constrain_bond(0, kROH);
    rattle.constrain_angle(0, kDHH);
    sim.add(rattle);
    // No integrator was registered.
    EXPECT_THROW(sim.build(), std::runtime_error);
}

/// An unresolvable partner THROWS. A bonded force style silently drops such a tuple
/// (a small missing energy term); a dropped CONSTRAINT is a molecule that is simply
/// not rigid — a physics corruption, not an approximation.
TEST(FixRattle, ThrowsOnUnresolvableConstraint) {
    cpu::Frame f = water_box(1, 6.0);
    f.atoms.bond_atom(0, 0) = 999999; // a tag that is not in the frame
    f.atoms.map_set(0, f.atoms.size());

    cpu::md::MD sim;
    sim.set_pair(cpu::PairLJCut{0.0067, 3.15, 6.0}.set_special_bonds(0.0, 0.0, 0.0));
    sim.set_neighbor(cpu::md::NeighborList<>{6.0, 0.5});
    sim.set_integrator(cpu::md::VelocityVerlet{1.0});
    cpu::md::FixRattle rattle{1.0};
    rattle.constrain_bond(0, kROH);
    sim.add(rattle);
    sim.enable_timing(false);
    sim.build();
    // The tag has no mass in the (global) table: no rank owns such an atom.
    EXPECT_THROW(sim.run(f, 1), std::runtime_error);
}

/// THE INVARIANT: the geometry stays rigid for the whole trajectory.
TEST(FixRattle, GeometryStaysRigidOverATrajectory) {
    const int n_side = 2, n_mol = n_side * n_side * n_side; // 8 waters — unit
    cpu::Frame f = water_box(n_side, 5.0);          // tests stay seconds-fast
    const double dt = 2.0; // fs — larger than a flexible water could survive

    cpu::md::MD sim;
    sim.set_pair(cpu::PairLJCut{0.0067, 3.15, 6.0}.set_special_bonds(0.0, 0.0, 0.0));
    sim.set_neighbor(cpu::md::NeighborList<>{6.0, 0.5});
    sim.set_integrator(cpu::md::VelocityVerlet{dt});
    cpu::md::FixRattle rattle{dt};
    rattle.constrain_bond(0, kROH);
    rattle.constrain_angle(0, kDHH);
    sim.add(rattle);
    sim.enable_timing(false);
    sim.build();

    // Kick the molecules so the constraints have something to fight.
    for (int i = 0; i < f.atoms.size(); ++i) {
        f.atoms.vx[i] = 0.004 * std::sin(double(3 * i + 1));
        f.atoms.vy[i] = 0.004 * std::cos(double(5 * i + 2));
        f.atoms.vz[i] = 0.004 * std::sin(double(7 * i + 3));
    }

    EXPECT_LT(worst_geometry_error(f, n_mol), 1e-12) << "start geometry is not rigid";

    for (int chunk = 0; chunk < 5; ++chunk) {
        sim.run(f, 20);
        EXPECT_LT(worst_geometry_error(f, n_mol), 1e-8)
            << "constrained geometry drifted after " << 20 * (chunk + 1) << " steps";
    }
}

/// THE OTHER INVARIANT: after the velocity stage, the relative velocity along every
/// constrained bond vanishes — d/dt |r_ij|^2 = 0. Skip this stage and the lengths
/// stay right while the molecule silently absorbs energy along its own bonds.
TEST(FixRattle, RelativeVelocityAlongBondsVanishes) {
    const int n_side = 2, n_mol = 8;
    cpu::Frame f = water_box(n_side, 6.0);
    const double dt = 2.0;

    cpu::md::MD sim;
    sim.set_pair(cpu::PairLJCut{0.0067, 3.15, 6.0}.set_special_bonds(0.0, 0.0, 0.0));
    sim.set_neighbor(cpu::md::NeighborList<>{6.0, 0.5});
    sim.set_integrator(cpu::md::VelocityVerlet{dt});
    cpu::md::FixRattle rattle{dt};
    rattle.constrain_bond(0, kROH);
    rattle.constrain_angle(0, kDHH);
    sim.add(rattle);
    sim.enable_timing(false);
    sim.build();

    for (int i = 0; i < f.atoms.size(); ++i) {
        f.atoms.vx[i] = 0.003 * std::sin(double(2 * i + 1));
        f.atoms.vy[i] = 0.003 * std::cos(double(4 * i));
        f.atoms.vz[i] = 0.003 * std::sin(double(6 * i + 5));
    }
    sim.run(f, 25);

    for (int m = 0; m < n_mol; ++m) {
        const int o = 3 * m, h1 = o + 1, h2 = o + 2;
        const std::pair<int, int> pairs[3] = {{o, h1}, {o, h2}, {h1, h2}};
        for (const auto& p : pairs) {
            double dx = f.atoms.x[p.first] - f.atoms.x[p.second];
            double dy = f.atoms.y[p.first] - f.atoms.y[p.second];
            double dz = f.atoms.z[p.first] - f.atoms.z[p.second];
            f.box.min_image(dx, dy, dz);
            const double dvx = f.atoms.vx[p.first] - f.atoms.vx[p.second];
            const double dvy = f.atoms.vy[p.first] - f.atoms.vy[p.second];
            const double dvz = f.atoms.vz[p.first] - f.atoms.vz[p.second];
            const double rv = dx * dvx + dy * dvy + dz * dvz;
            EXPECT_LT(std::fabs(rv), 1e-8)
                << "molecule " << m << ": relative velocity along a constrained bond";
        }
    }
}
