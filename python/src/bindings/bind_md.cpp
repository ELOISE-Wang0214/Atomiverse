/// @file python/bind_md.cpp
/// @brief MD driver, components, potentials, and observer bindings.

#include "cpu_frame_wrapper.h"

#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include "atomiverse/cpu/md/md.h"
#include "atomiverse/cpu/md/velocity_verlet.h"
#include "atomiverse/cpu/md/init.h"
#include "atomiverse/cpu/md/nose_hoover.h"
#include "atomiverse/cpu/md/lj.h"
#include "atomiverse/cpu/md/neighbor_list.h"
#include "atomiverse/cpu/bonded/harmonic_bond.h"
#include "atomiverse/cpu/bonded/harmonic_angle.h"
#include "atomiverse/cpu/bonded/multi_harmonic_dihedral.h"
#include "atomiverse/cpu/observer/thermo_component.h"
#include "atomiverse/facade/stage.h"

void bind_md(nb::module_& m) {
    // ── Phase enum ──
    nb::enum_<md::Phase>(m, "Phase", "MD phase sequence (LAMMPS-style).")
        .value("initial_integrate", md::Phase::initial_integrate)
        .value("pre_force",         md::Phase::pre_force)
        .value("force",             md::Phase::force)
        .value("post_force",        md::Phase::post_force)
        .value("final_integrate",   md::Phase::final_integrate)
        .value("end_of_step",       md::Phase::end_of_step)
    ;

    // ── MD Driver (CPU backend) ──
    using CpuMD = cpu::md::MD;
    using PhaseBuilder = PhaseBuilder<md::Phase, cpu::device_tag>;

    nb::class_<CpuMD>(m, "MD", "Molecular dynamics driver (CPU backend).")
        .def(nb::init<>())
        // add() with default bindings — one overload per component type
        .def("add_velocity_verlet", [](CpuMD& self, cpu::md::VelocityVerlet comp) {
            self.add(std::move(comp));
        }, nb::arg("comp"))
        .def("add_init_velocity", [](CpuMD& self, cpu::md::InitVelocity comp) {
            self.add(std::move(comp));
        }, nb::arg("comp"))
        .def("add_nose_hoover", [](CpuMD& self, cpu::md::NoseHoover comp) {
            self.add(std::move(comp));
        }, nb::arg("comp"))
        // on(phase) returns a PhaseBuilder for explicit binding
        .def("on", [](CpuMD& self, md::Phase phase) -> PhaseBuilder {
            return self.on(phase);
        }, nb::arg("phase"), nb::rv_policy::move)
        .def("build", &CpuMD::build, "Compile the execution plan.")
        .def("run", [](CpuMD& self, nb::object frame_obj, int n_steps) {
            // Extract CpuFrame from Python object
            auto& cf = nb::cast<CpuFrame&>(frame_obj);
            {
                nb::gil_scoped_release release;
                self.run(cf.frame, n_steps);
            }
        }, nb::arg("frame"), nb::arg("n_steps"),
           "Run simulation. GIL is released during computation.")
        .def("enable_timing", &CpuMD::enable_timing)
        .def("summary", &CpuMD::summary)
    ;

    // ── PhaseBuilder ──
    nb::class_<PhaseBuilder>(m, "PhaseBuilder",
        "Fluent builder for explicit phase binding.")
        .def("add_lj", [](PhaseBuilder& self, cpu::md::LJ comp, int order) -> PhaseBuilder& {
            return self.add(std::move(comp), order);
        }, nb::arg("comp"), nb::arg("order") = 0, nb::rv_policy::reference)
        .def("add_neighbor_list", [](PhaseBuilder& self, cpu::md::NeighborList<> comp, int order) -> PhaseBuilder& {
            return self.add(std::move(comp), order);
        }, nb::arg("comp"), nb::arg("order") = 0, nb::rv_policy::reference)
        .def("add_harmonic_bond", [](PhaseBuilder& self, cpu::bonded::HarmonicBond comp, int order) -> PhaseBuilder& {
            return self.add(std::move(comp), order);
        }, nb::arg("comp"), nb::arg("order") = 0, nb::rv_policy::reference)
        .def("add_harmonic_angle", [](PhaseBuilder& self, cpu::bonded::HarmonicAngle comp, int order) -> PhaseBuilder& {
            return self.add(std::move(comp), order);
        }, nb::arg("comp"), nb::arg("order") = 0, nb::rv_policy::reference)
        .def("add_multi_harmonic_dihedral", [](PhaseBuilder& self, cpu::bonded::MultiHarmonicDihedral comp, int order) -> PhaseBuilder& {
            return self.add(std::move(comp), order);
        }, nb::arg("comp"), nb::arg("order") = 0, nb::rv_policy::reference)
        .def("add_thermo", [](PhaseBuilder& self, cpu::Thermo comp, int order) -> PhaseBuilder& {
            return self.add(std::move(comp), order);
        }, nb::arg("comp"), nb::arg("order") = 0, nb::rv_policy::reference)
    ;

    // ── MD Components ──

    nb::class_<cpu::md::VelocityVerlet>(m, "VelocityVerlet",
        "Symplectic velocity-Verlet integrator.")
        .def(nb::init<Real>(), nb::arg("dt"),
             "Construct with timestep dt [fs].")
    ;

    nb::class_<cpu::md::InitVelocity>(m, "InitVelocity",
        "Maxwell-Boltzmann velocity initialization.")
        .def(nb::init<Real, unsigned>(),
             nb::arg("temperature"), nb::arg("seed"),
             "Construct with target temperature [K] and RNG seed.")
    ;

    nb::class_<cpu::md::NoseHoover>(m, "NoseHoover",
        "Nose-Hoover thermostat.")
        .def(nb::init<Real, Real, Real>(),
             nb::arg("T_target"), nb::arg("tau"), nb::arg("dt"),
             "Construct with target temperature [K], coupling time [fs], and timestep [fs].")
    ;

    // ── Potentials ──

    nb::class_<cpu::md::LJ>(m, "LJ",
        "Lennard-Jones 12-6 pair potential.")
        .def(nb::init<Real, Real, Real>(),
             nb::arg("epsilon"), nb::arg("sigma"), nb::arg("cutoff"),
             "Construct with epsilon [eV], sigma [A], and cutoff [A].")
    ;

    nb::class_<cpu::md::NeighborList<>>(m, "NeighborList",
        "Cell-list based neighbor list component.")
        .def(nb::init<Real, Real>(),
             nb::arg("cutoff"), nb::arg("skin") = Real(-1),
             "Construct with cutoff [A] and optional skin distance [A].")
    ;

    nb::class_<cpu::bonded::HarmonicBond>(m, "HarmonicBond",
        "Harmonic bond potential V = K(r - r0)^2.")
        .def(nb::init<std::vector<Real>, std::vector<Real>>(),
             nb::arg("K"), nb::arg("r0"),
             "Construct with force constants K and equilibrium distances r0 per bond type.")
    ;

    nb::class_<cpu::bonded::HarmonicAngle>(m, "HarmonicAngle",
        "Harmonic angle potential V = K(theta - theta0)^2.")
        .def(nb::init<std::vector<Real>, std::vector<Real>>(),
             nb::arg("K"), nb::arg("theta0"),
             "Construct with force constants K and equilibrium angles theta0 per angle type.")
    ;

    nb::class_<cpu::bonded::MultiHarmonicDihedral>(m, "MultiHarmonicDihedral",
        "Multi-harmonic dihedral potential.")
        .def(nb::init<std::vector<std::vector<Real>>>(),
             nb::arg("Cn"),
             "Construct with coefficients Cn per dihedral type.")
    ;

    // ── Observers ──

    nb::class_<cpu::Thermo>(m, "Thermo",
        "Thermodynamic output at end_of_step.")
        .def(nb::init<int>(), nb::arg("freq"),
             "Construct with output frequency (every N steps).")
    ;
}
