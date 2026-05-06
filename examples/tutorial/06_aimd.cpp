/// @file examples/tutorial/06_aimd.cpp
/// @brief Tutorial example for ab initio molecular dynamics on H2.
///
/// AIMD = Born-Oppenheimer molecular dynamics. The classical MD loop runs
/// unchanged; the only difference is that the per-step force evaluator is
/// `cpu::scf::SCFEvaluator` (a self-consistent SCF + analytic gradient)
/// instead of a classical pair potential. This is the same Driver/Component
/// composition that drives geometry optimization in
/// `examples/tutorial/05_geom_opt.cpp`, only with the integration loop
/// (`cpu::md::MD`) instead of the optimization loop (`cpu::minimizer::
/// Minimize`).
///
/// Demonstrates:
///   1. Construct an H2 frame and integrals at STO-3G.
///   2. Compose `cpu::md::MD` with InitVelocity + VelocityVerlet +
///      NoseHoover (NVT ensemble).
///   3. Bind `SCFEvaluator{config}.predictor(3)` at md::Phase::force —
///      ASPC predictor warm-starts the density across MD steps.
///   4. Run 200 steps and report initial / final energy + temperature.
///
/// Build:
///   cmake -B build -DATV_BUILD_CUDA=OFF -DATV_BUILD_RHF=ON \
///         -DATV_BUILD_EXAMPLES=ON
///   cmake --build build --target 06_aimd
///
/// Run:
///   ./build/examples/06_aimd
///
/// Expected: H2 vibrates around its equilibrium bond length (~0.74 Å at
/// RHF/STO-3G), kinetic temperature drifts toward 100 K under the
/// Nose-Hoover thermostat, total potential energy stays near
/// -31.7 eV (-1.165 Ha).
///
/// We pick T = 100 K rather than 300 K because RHF cannot describe
/// homolytic H2 dissociation correctly — a hot RHF trajectory will
/// drift onto the unphysical RHF ionic curve and the molecule appears
/// to fly apart. UKS-PBE or CASSCF would handle this correctly; this
/// tutorial keeps things RHF for clarity, so we keep the kinetic
/// energy comfortably below the (incorrect RHF) barrier.

#include <cmath>

#include "atomiverse/core/log.h"
#include "atomiverse/core/units.h"
#include "atomiverse/cpu/frame.h"
#include "atomiverse/cpu/md/md.h"
#include "atomiverse/cpu/md/init.h"
#include "atomiverse/cpu/md/velocity_verlet.h"
#include "atomiverse/cpu/md/nose_hoover.h"
#include "atomiverse/cpu/scf/integrals.h"
#include "atomiverse/cpu/scf/scf_evaluator.h"
#include "atomiverse/cpu/scf/scf_solve.h"

using namespace ATV;

namespace {

/// @brief H2 starting geometry. Placement inside the simulation box is
///        irrelevant because `SCFEvaluator::launch` reassembles the
///        molecule via minimum-image relative to atom 0 before each
///        integral build, so atoms split across a periodic boundary
///        are recomposed transparently.
cpu::Frame make_h2() {
    cpu::Frame f(2);
    f.atoms.x = {0.0, 0.74};   // 0.74 Å — RHF/STO-3G equilibrium
    f.atoms.y = {0.0, 0.0};
    f.atoms.z = {0.0, 0.0};
    f.atoms.mass = {1.008, 1.008};
    f.atoms.type = {1, 1};
    f.box = Mat3::cubic(100.0);  // SCF treats the molecule as isolated
    return f;
}

/// @brief Approximate kinetic temperature from frame.kinetic.
///        T = 2·KE / (n_dof · kB), with COM removed (n_dof = 3·N - 3).
double temperature_K(const cpu::Frame& frame) {
    const int n_dof = 3 * frame.atoms.size() - 3;
    if (n_dof <= 0) return 0.0;
    return 2.0 * static_cast<double>(frame.kinetic) / (n_dof * kBoltzmann);
}

} // namespace

int main() {
    ATV::log_init(spdlog::level::info);
    ATV_INFO("=== Atomiverse: AIMD on H2 (NVT, Nose-Hoover) ===");
    ATV_INFO("    Molecule:  H2");
    ATV_INFO("    Method:    RHF / STO-3G");
    ATV_INFO("    Ensemble:  NVT @ 100 K (Nose-Hoover, tau = 200 fs)");
    ATV_INFO("    Predictor: ASPC order 3 (density warm-start)");

    auto frame = make_h2();

    auto integrals = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto config = cpu::scf::rhf_config(integrals);

    const double dt        = 0.5;     // fs
    const double T_target  = 100.0;   // K — well below RHF/STO-3G unphysical barrier
    const double tau       = 200.0;   // fs (loose coupling — H2 has only 3 DOFs)
    const int    n_steps   = 100;     // 50 fs total — captures ~7 vib periods (~7.6 fs)
    // Note: instantaneous T fluctuates strongly for a 2-atom system because
    // there are only 3 vibrational/rotational DOFs after COM removal. The
    // ensemble average over many trajectories converges to T_target;
    // a single short trajectory does not.

    // Build the simulation. `SCFEvaluator{config}.predictor(3)` is a
    // chained rvalue builder: ref-qualified overloads in SCFEvaluator
    // make the temporary itself become the rvalue argument to add(),
    // which moves it into driver-owned storage.
    cpu::md::MD sim;
    sim.add(cpu::md::InitVelocity{T_target, /*seed=*/42u});
    sim.add(cpu::md::VelocityVerlet{dt});
    sim.add(cpu::md::NoseHoover{T_target, tau, dt});
    sim.on(md::Phase::force).add(
        cpu::scf::SCFEvaluator{config}.predictor(/*order=*/3));
    sim.build();

    // Take one step to populate frame.kinetic and frame.energy from the
    // integrator's first force evaluation.
    sim.run(frame, 1);
    const double E_init  = static_cast<double>(frame.energy);
    const double KE_init = static_cast<double>(frame.kinetic);
    const double T_init  = temperature_K(frame);
    ATV_INFO(
        "[aimd] step=1   PE_eV={:>10.4f}  KE_eV={:>8.4f}  T_K={:>7.2f}",
        E_init, KE_init, T_init);

    sim.run(frame, n_steps - 1);

    const double E_final  = static_cast<double>(frame.energy);
    const double KE_final = static_cast<double>(frame.kinetic);
    const double T_final  = temperature_K(frame);
    ATV_INFO(
        "[aimd] step={:<3d} PE_eV={:>10.4f}  KE_eV={:>8.4f}  T_K={:>7.2f}",
        n_steps, E_final, KE_final, T_final);

    // Bond length via minimum image — atoms may have wrapped through
    // the box independently during the trajectory, but we want the
    // physical bond.
    const double lx = frame.box(0, 0);
    auto min_image = [](double d, double L) {
        return (L > 0.0) ? (d - L * std::round(d / L)) : d;
    };
    const double dx = min_image(static_cast<double>(frame.atoms.x(1) - frame.atoms.x(0)), lx);
    const double dy = min_image(static_cast<double>(frame.atoms.y(1) - frame.atoms.y(0)), lx);
    const double dz = min_image(static_cast<double>(frame.atoms.z(1) - frame.atoms.z(0)), lx);
    const double bond = std::sqrt(dx * dx + dy * dy + dz * dz);
    ATV_INFO("[aimd] final H-H bond = {:.4f} Å (RHF/STO-3G eq ≈ 0.74 Å)", bond);

    ATV_INFO("AIMD run complete. {} SCF + gradient evaluations performed.",
             n_steps);
    return 0;
}
