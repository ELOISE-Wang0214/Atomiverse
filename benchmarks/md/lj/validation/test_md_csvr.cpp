/// @file tests/cpu/md/test_cpu_md_csvr.cpp
/// @brief CSVR (Bussi) thermostat — canonical sampling (ac-008) and the
///        drift-free conserved quantity (ac-009).
///
/// Drives a 256-atom FCC argon LJ gas with VelocityVerlet + FixCSVR and
/// checks (a) the kinetic-energy distribution is the canonical one — mean
/// T at the setpoint AND relative variance Var(K)/<K>^2 ≈ 2/N_f (the
/// signature that distinguishes CSVR from Berendsen, which collapses the
/// fluctuations); (b) the published H~ = E_pot + E_kin + H~_therm
/// (frame.meta.get(ATV::kConserved)) shows no secular drift.

#include <gtest/gtest.h>

#include <cmath>
#include <vector>

#include "atomiverse/core/mat3.h"
#include "atomiverse/core/units.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/md/fix_csvr.h"
#include "atomiverse/cpu/md/init.h"
#include "atomiverse/cpu/md/md.h"
#include "atomiverse/cpu/md/velocity_verlet.h"
#include "atomiverse/cpu/neighbor/neighbor_list.h"
#include "atomiverse/cpu/potentials/pair_lj_cut.h"

using namespace ATV;

namespace {

cpu::Frame make_fcc_argon(int n_cells, double a) {
    const int n = n_cells * n_cells * n_cells * 4;
    cpu::Frame frame(n);
    const double basis[4][3] = {
        {0.0, 0.0, 0.0}, {0.5, 0.5, 0.0}, {0.5, 0.0, 0.5}, {0.0, 0.5, 0.5}};
    int idx = 0;
    for (int ix = 0; ix < n_cells; ++ix)
        for (int iy = 0; iy < n_cells; ++iy)
            for (int iz = 0; iz < n_cells; ++iz)
                for (int b = 0; b < 4; ++b) {
                    frame.span(ATV::kAtom_x)[idx] = (ix + basis[b][0]) * a;
                    frame.span(ATV::kAtom_y)[idx] = (iy + basis[b][1]) * a;
                    frame.span(ATV::kAtom_z)[idx] = (iz + basis[b][2]) * a;
                    frame.span(ATV::kAtom_mass)[idx] = 39.948;
                    frame.span(ATV::kAtom_type)[idx] = 18;
                    ++idx;
                }
    frame.simbox = ATV::cpu::Box::cubic(n_cells * a);
    return frame;
}

// Kinetic energy [eV] from the live velocities (independent of any observer).
double kinetic_ev(const cpu::Frame& f) {
    double s = 0.0;
    for (int i = 0; i < f.nlocal(); ++i)
        s += f.span(ATV::kAtom_mass)[i]
             * (f.span(ATV::kAtom_vx)[i] * f.span(ATV::kAtom_vx)[i]
                + f.span(ATV::kAtom_vy)[i] * f.span(ATV::kAtom_vy)[i]
                + f.span(ATV::kAtom_vz)[i] * f.span(ATV::kAtom_vz)[i]);
    return 0.5 * units::mass_v2_to_energy * s;
}

} // namespace

TEST(FixCSVR, CanonicalSamplingAndConservedQuantity) {
    const double T = 300.0;   // K
    const double dt = 5.0;    // fs
    const double tau = 100.0; // fs (CP2K TIMECON-style coupling)
    const double eps = 0.0104, sigma = 3.4, cutoff = 8.5;

    auto frame = make_fcc_argon(4, 5.26); // 256 atoms
    const int natoms = frame.nlocal();
    const int n_dof = 3 * natoms - 3; // COM removed

    cpu::md::MD sim;
    sim.add(cpu::md::VelocityVerlet{dt});
    sim.add(cpu::md::InitVelocity{T, /*seed=*/2024});
    sim.add_potential(cpu::PairLJCut{eps, sigma, cutoff});
    sim.set_neighbor(cpu::md::NeighborList<>{cutoff, 1.0});
    sim.set_thermostat(cpu::md::FixCSVR{
        T, tau, dt, /*seed=*/7}); // plugs in via set_thermostat
    sim.build();

    sim.run(frame, 3000); // equilibrate

    const int n_samples = 1500;
    const int interval = 15; // ~tau/dt decorrelation
    std::vector<double> ke, cons;
    ke.reserve(n_samples);
    cons.reserve(n_samples);
    for (int s = 0; s < n_samples; ++s) {
        sim.run(frame, interval);
        ke.push_back(kinetic_ev(frame));
        cons.push_back(frame.meta.get(ATV::kConserved));
    }

    // ── ac-008: canonical kinetic-energy distribution ──────────────────
    double mean_ke = 0.0;
    for (double k : ke)
        mean_ke += k;
    mean_ke /= n_samples;
    double var_ke = 0.0;
    for (double k : ke)
        var_ke += (k - mean_ke) * (k - mean_ke);
    var_ke /= n_samples;

    const double mean_T = 2.0 * mean_ke / (n_dof * units::kB_ev_per_k);
    const double rel_var = var_ke / (mean_ke * mean_ke);
    const double canonical_rel_var = 2.0 / n_dof; // exact for the canonical K

    std::printf("[CSVR] N=%d N_f=%d  mean_T=%.2f K (target %.0f)  "
                "rel_var(K)=%.4e  canonical 2/N_f=%.4e  ratio=%.3f\n",
                natoms, n_dof, mean_T, T, rel_var, canonical_rel_var,
                rel_var / canonical_rel_var);

    // Mean temperature at the setpoint (CSVR is unbiased).
    EXPECT_NEAR(mean_T, T, 0.06 * T) << "mean temperature off setpoint";
    // Canonical fluctuation signature: relative variance ≈ 2/N_f. Berendsen
    // would give << this; a runaway would give >>. Allow a generous
    // statistical band but tight enough to be a real canonical check.
    EXPECT_GT(rel_var, 0.5 * canonical_rel_var)
        << "kinetic fluctuations too small — not canonical (Berendsen-like)";
    EXPECT_LT(rel_var, 1.8 * canonical_rel_var)
        << "kinetic fluctuations too large";

    // ── ac-009: conserved quantity H~ is drift-free ────────────────────
    // Compare the mean of the first vs the second half of the trajectory;
    // a correctly book-kept H~ has no secular drift. Scale the tolerance by
    // the kinetic energy (the natural energy unit of the run).
    double first = 0.0, second = 0.0;
    const int half = n_samples / 2;
    for (int s = 0; s < half; ++s)
        first += cons[s];
    for (int s = half; s < n_samples; ++s)
        second += cons[s];
    first /= half;
    second /= (n_samples - half);
    const double drift = std::abs(second - first);

    std::printf("[CSVR] H~ first-half=%.6f second-half=%.6f  drift=%.3e eV  "
                "(mean_ke=%.4f eV)\n",
                first, second, drift, mean_ke);

    // H~ should be conserved to the velocity-Verlet accuracy: any secular
    // drift must be a tiny fraction of the kinetic energy scale.
    EXPECT_LT(drift, 0.01 * mean_ke)
        << "CSVR conserved quantity drifts — bookkeeping (Eq. 15) is wrong";
}
