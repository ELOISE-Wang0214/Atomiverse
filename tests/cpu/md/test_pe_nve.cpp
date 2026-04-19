/// @file tests/cpu/md/test_pe_nve.cpp
/// @brief Phase A.4 (msu-pe-validation) -- single-chain PE NVE energy conservation.
///
/// Runs an all-trans C20H42 chain with a deterministic, tiny initial-velocity
/// kick for 1000 Velocity-Verlet steps at Δt = 0.5 fs. The three bonded
/// potentials (HarmonicBond, HarmonicAngle, MultiHarmonicDihedral) are on;
/// LJ is intentionally off so this test isolates the integrator × bonded
/// pipeline. With dihedral forces fixed (see phase-A.2 commit), NVE drift
/// should stay below 1e-4 relative to the initial total energy over 500 fs.

#include <gtest/gtest.h>

#include <cmath>
#include <vector>

#include "atomiverse/cpu/bonded/harmonic_angle.h"
#include "atomiverse/cpu/bonded/harmonic_bond.h"
#include "atomiverse/cpu/bonded/multi_harmonic_dihedral.h"
#include "atomiverse/cpu/frame.h"
#include "atomiverse/cpu/md/md.h"
#include "atomiverse/cpu/md/velocity_verlet.h"
#include "atomiverse/core/mat3.h"
#include "atomiverse/core/units.h"

using namespace ATV;

namespace {

// ── OPLS-AA UA parameters (same as test_pe_single_chain.cpp) ──────────

inline constexpr double kKcalPerMolToEv = 1.0 / 23.060547830619026;
inline constexpr double kBondK_kcal  = 268.0;
inline constexpr double kBondR0      = 1.529;
inline constexpr double kAngleK_kcal = 58.35;
inline constexpr double kAngleTheta0 = 112.7 * kPI / 180.0;
inline constexpr double kV1_kcal = 1.411;
inline constexpr double kV2_kcal = -0.271;
inline constexpr double kV3_kcal = 3.145;

inline double bond_K_eV()  { return kBondK_kcal  * kKcalPerMolToEv; }
inline double angle_K_eV() { return kAngleK_kcal * kKcalPerMolToEv; }
inline std::vector<double> dihedral_Cn_eV() {
    const double v1 = kV1_kcal * kKcalPerMolToEv;
    const double v2 = kV2_kcal * kKcalPerMolToEv;
    const double v3 = kV3_kcal * kKcalPerMolToEv;
    return {
        v1 / 2.0 + v2 + v3 / 2.0,
        v1 / 2.0 - 1.5 * v3,
        -v2,
         2.0 * v3,
         0.0,
    };
}

cpu::Frame make_all_trans_alkane(int n_carbons, double box_len) {
    cpu::Frame frame(n_carbons);

    const double half_apex = (kPI - kAngleTheta0) / 2.0;
    const double dx_step = kBondR0 * std::cos(half_apex);
    const double dy_step = kBondR0 * std::sin(half_apex);

    const double x0 = (box_len - (n_carbons - 1) * dx_step) / 2.0;
    const double y0 = box_len / 2.0;
    const double z0 = box_len / 2.0;

    for (int i = 0; i < n_carbons; ++i) {
        frame.atoms.x[i] = x0 + i * dx_step;
        frame.atoms.y[i] = y0 + (i % 2 == 0 ? 0.0 : dy_step);
        frame.atoms.z[i] = z0;
        frame.atoms.mass[i] = 14.027;  // CH2 united atom (amu)
        frame.atoms.type[i] = 0;
    }

    const int n_bonds = n_carbons - 1;
    frame.bonds.atomi = xt::zeros<int>({n_bonds});
    frame.bonds.atomj = xt::zeros<int>({n_bonds});
    frame.bonds.type  = xt::zeros<int>({n_bonds});
    for (int k = 0; k < n_bonds; ++k) {
        frame.bonds.atomi[k] = k;
        frame.bonds.atomj[k] = k + 1;
        frame.bonds.type[k]  = 0;
    }

    const int n_angles = n_carbons - 2;
    frame.angles.atomi = xt::zeros<int>({n_angles});
    frame.angles.atomj = xt::zeros<int>({n_angles});
    frame.angles.atomk = xt::zeros<int>({n_angles});
    frame.angles.type  = xt::zeros<int>({n_angles});
    for (int k = 0; k < n_angles; ++k) {
        frame.angles.atomi[k] = k;
        frame.angles.atomj[k] = k + 1;
        frame.angles.atomk[k] = k + 2;
        frame.angles.type[k]  = 0;
    }

    const int n_dihedrals = n_carbons - 3;
    frame.dihedrals.atomi = xt::zeros<int>({n_dihedrals});
    frame.dihedrals.atomj = xt::zeros<int>({n_dihedrals});
    frame.dihedrals.atomk = xt::zeros<int>({n_dihedrals});
    frame.dihedrals.atoml = xt::zeros<int>({n_dihedrals});
    frame.dihedrals.type  = xt::zeros<int>({n_dihedrals});
    for (int k = 0; k < n_dihedrals; ++k) {
        frame.dihedrals.atomi[k] = k;
        frame.dihedrals.atomj[k] = k + 1;
        frame.dihedrals.atomk[k] = k + 2;
        frame.dihedrals.atoml[k] = k + 3;
        frame.dihedrals.type[k]  = 0;
    }

    frame.box = Mat3::cubic(box_len);
    return frame;
}

} // namespace

// ── NVE conservation over 1000 Velocity-Verlet steps ────────────────

TEST(PeNve, EnergyConserves1000Steps) {
    auto frame = make_all_trans_alkane(20, 50.0);

    // Give a tiny deterministic transverse kick so the chain samples a
    // non-trivial bond-angle-dihedral trajectory without drifting far from
    // the equilibrium geometry. Amplitude chosen so the excursion is well
    // inside the harmonic regime of every bonded potential.
    for (int i = 0; i < frame.atoms.size(); ++i) {
        const double sign = (i % 2 == 0) ? +1.0 : -1.0;
        frame.atoms.vz[i] = sign * 1e-3;   // Å / fs
    }

    cpu::md::MD sim;
    sim.add(cpu::md::VelocityVerlet{0.5});   // Δt = 0.5 fs
    sim.on(ATV::md::Phase::force)
        .add(cpu::bonded::HarmonicBond{{bond_K_eV()}, {kBondR0}})
        .add(cpu::bonded::HarmonicAngle{{angle_K_eV()}, {kAngleTheta0}})
        .add(cpu::bonded::MultiHarmonicDihedral{{dihedral_Cn_eV()}});
    sim.build();

    sim.run(frame, 1);
    const double e0 = frame.energy + frame.kinetic;
    ASSERT_GT(std::abs(e0), 0.0) << "initial energy must be nonzero";

    sim.run(frame, 999);
    const double e1 = frame.energy + frame.kinetic;

    const double drift = std::abs(e1 - e0) / std::abs(e0);
    EXPECT_LT(drift, 1e-4)
        << "NVE drift " << drift
        << " exceeds 1e-4 over 1000 Velocity-Verlet steps (E0=" << e0
        << ", E1=" << e1 << ")";
}
