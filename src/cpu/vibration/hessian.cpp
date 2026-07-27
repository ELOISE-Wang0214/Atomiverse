/// @file src/cpu/vib/hessian.cpp
/// @brief cpu::vib::Hessian — constructor wiring and run() implementation.

#include "atomiverse/cpu/vibration/hessian.h"

#include "atomiverse/cpu/component/zero.h"
#include "atomiverse/cpu/vibration/frequency_report.h"
#include "atomiverse/cpu/vibration/hessian_diagonalizer.h"
#include "atomiverse/cpu/vibration/mass_weighting.h"
#include "atomiverse/cpu/vibration/normal_mode_writer.h"
#include "atomiverse/cpu/vibration/thermochem.h"

namespace ATV {
namespace cpu {
namespace vib {

namespace {

/// @brief Bind the driver-owned accumulator zeroing for the force phase.
///
/// Zeroing the force/energy accumulators is the hosting driver's job, never a
/// force source's: a source that zeroed would erase a co-bound source. Without
/// it an accumulating (`+=`) classical Potential integrates its forces across
/// the `6N` displacement passes and the central difference returns `k/2`
/// instead of the curvature `k`. Orders `-1000`/`-999` mirror `md::MD` and
/// `vib::Phonon`, so every driver's force slot has the same head.
///
/// Called from **every** construction path (`wire_full_chain` and the
/// `WithoutAccumulator` ctor) exactly once, so the two wiring sites cannot
/// drift apart.
///
/// @param driver Driver to register the two Zero components into.
template <class Driver>
void wire_zeroing(Driver& driver) {
    driver.on(ATV::vib::Phase::force)
        .add(cpu::ZeroForces{}, -1000)
        .add(cpu::ZeroEnergy{}, -999);
}

/// @brief Register the full vibrational-analysis component chain.
///
/// Auto-registers the four numerical-Hessian structural components, the
/// two frequency components, and the three thermochemistry-link
/// components, all bound to the driver-owned shared workspace. On
/// `vib::Phase::diagonalize` the ordering is Symmetrize (-1) ->
/// MassWeighting (0) -> HessianDiagonalizer (+1); on `vib::Phase::thermochem`
/// the Thermochem Observer runs; on `vib::Phase::end_of_step` the
/// FrequencyReport (0) and NormalModeWriter (+1) Recorders run.
///
/// @param driver Driver to register the components into.
/// @param ws     Driver-owned shared workspace.
template <class Driver>
void wire_full_chain(Driver& driver, HessianWorkspace* ws, double temperature,
                     double pressure, int multiplicity, int symmetry_number,
                     const std::string& frequency_tsv,
                     const std::string& mode_xyz_prefix) {
    wire_zeroing(driver);
    driver.add(ReferenceGeometry{ws});
    driver.add(DisplaceStep{ws});
    driver.add(HessianAccumulator{ws});
    driver.add(Symmetrize{ws});
    driver.add(MassWeighting{ws});
    driver.add(HessianDiagonalizer{ws});
    driver.add(
        Thermochem{ws, temperature, pressure, multiplicity, symmetry_number});
    driver.add(FrequencyReport{ws, frequency_tsv});
    driver.add(NormalModeWriter{ws, mode_xyz_prefix});
}

} // namespace

Hessian::Hessian() {
    // Default thermochemistry: standard conditions, closed-shell singlet.
    wire_full_chain(*this, workspace_.get(), 298.15, 101325.0, 1, 1,
                    "frequencies.tsv", "mode");
}

Hessian::Hessian(double temperature, double pressure, int multiplicity,
                 int symmetry_number) {
    wire_full_chain(*this, workspace_.get(), temperature, pressure,
                    multiplicity, symmetry_number, "frequencies.tsv", "mode");
}

Hessian::Hessian(double temperature, double pressure, int multiplicity,
                 int symmetry_number, std::string frequency_tsv,
                 std::string mode_xyz_prefix) {
    wire_full_chain(*this, workspace_.get(), temperature, pressure,
                    multiplicity, symmetry_number, frequency_tsv,
                    mode_xyz_prefix);
}

Hessian::Hessian(WithoutAccumulator /*tag*/) {
    // Test-only configuration: omit HessianAccumulator so the build-time
    // validation throws for the missing Capability::Hessian component.
    wire_zeroing(*this);
    add(ReferenceGeometry{workspace_.get()});
    add(DisplaceStep{workspace_.get()});
    add(Symmetrize{workspace_.get()});
}

HessianResult Hessian::run(cpu::Frame& frame) {
    // Size the shared workspace from the frame before the component loop
    // runs — no allocation happens inside any component launch().
    const int n_atoms = frame.nlocal();
    workspace_->allocate(n_atoms);

    // The 6N-displacement loop: two displacements per Cartesian coordinate. A fixed
    // budget, not a policy — every displacement must complete for the Hessian to be
    // defined, so there is nothing for a criterion to decide and none is polled.
    const int n_disp = 6 * n_atoms;
    setup_pass(frame);
    for (int d = 0; d < n_disp; ++d)
        step_once(d);
    // The finalize pass: Symmetrize, MassWeighting and the diagonalizer run here;
    // DisplaceStep restores the equilibrium geometry.
    step_once(n_disp);
    teardown_pass();

    // Publish the assembled, symmetrized Hessian, the vibrational
    // frequencies (MassWeighting + HessianDiagonalizer), and the harmonic
    // thermochemistry (Thermochem).
    last_result_.hessian = workspace_->hessian;
    last_result_.n_atoms = n_atoms;
    last_result_.units = "eV/Angstrom^2";
    last_result_.frequencies = workspace_->frequencies;
    last_result_.normal_modes = workspace_->normal_modes;
    last_result_.n_imaginary = workspace_->n_imaginary;
    last_result_.is_linear = workspace_->is_linear;
    last_result_.zpe = workspace_->zpe;
    last_result_.u_tot = workspace_->u_tot;
    last_result_.h_tot = workspace_->h_tot;
    last_result_.s_tot = workspace_->s_tot;
    last_result_.s_elec = workspace_->s_elec;
    last_result_.s_trans = workspace_->s_trans;
    last_result_.s_rot = workspace_->s_rot;
    last_result_.s_vib = workspace_->s_vib;
    last_result_.g_tot = workspace_->g_tot;
    last_result_.g_thermal = workspace_->g_thermal;
    return last_result_;
}

} // namespace vib
} // namespace cpu
} // namespace ATV
